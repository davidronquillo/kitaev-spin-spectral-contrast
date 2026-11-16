"""Selected canonical quotient honeycomb geometry and hopping construction.

No geometry or numerical acquisition is constructed at import.
The background sampler and external labeling-map reader are not included. Links are stored
as s=-i b_A b_B, in lexicographic cell-key order followed by x,y,z. The physical
matter parity is derived from that actual labeled Clifford word.
"""
import hashlib
from collections import deque
from functools import lru_cache, reduce
from operator import xor

import numpy as np

LABELS = ("C3n02", "C3n04", "C3n08", "R24x32", "C3n12", "C3n16")
COLUMNS = {"C3n01": ((2, 2), (-2, 4)), "C3n02": ((4, 4), (-4, 8)),
           "C3n04": ((8, 8), (-8, 16)), "C3n08": ((16, 16), (-16, 32)),
           "R12x16": ((12, 0), (0, 16)), "R24x32": ((24, 0), (0, 32)),
           "C3n12": ((24,24),(-24,48)), "C3n16": ((32,32),(-32,64))}
BOND_NAMES = ("x", "y", "z")
DISPLACEMENTS = ((0, 0), (-1, 1), (-1, 0))
FLIP_THRESHOLD = 66871911161826185187153288884631757060


def parity_permutation(word):
    """Sign of a permutation of distinct integer Clifford generators."""
    if sorted(word) != list(range(len(word))):
        raise ValueError("Clifford word does not contain every generator once")
    seen = 0
    inversions = 0
    for generator in word:
        inversions ^= bin(seen >> (generator + 1)).count("1") & 1
        seen |= 1 << generator
    return -1 if inversions else 1


def clifford_word(word):
    mask = 0
    sign = 1
    for generator in word:
        if bin(mask >> (generator+1)).count("1") % 2:
            sign = -sign
        mask ^= 1 << generator
    return sign,mask


def rank_binary(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)




class Geometry:
    def __init__(self, label):
        self.label = label
        self.columns = COLUMNS[label]
        (a, c), (b, d) = self.columns
        self.cells = a*d-b*c
        self.N = 2*self.cells
        self.adj = ((d, -b), (-c, a))
        reps = {(0, 0): (0, 0)}
        pending = deque([(0, 0)])
        while pending:
            key = pending.popleft()
            x, y = reps[key]
            for dx, dy in ((1, 0), (0, 1)):
                coord = (x+dx, y+dy)
                new = self.key(*coord)
                if new not in reps:
                    reps[new] = coord
                    pending.append(new)
        self.keys = tuple(sorted(reps))
        self.coordinates = tuple(reps[key] for key in self.keys)
        self.index = {key: i for i, key in enumerate(self.keys)}
        if len(self.keys) != self.cells or self.keys[0] != (0, 0):
            raise ArithmeticError("canonical quotient cardinality/source mismatch")
        self.edges = tuple((name, 2*k, self.site(x+dx, y+dy, 1))
                           for k, (x, y) in enumerate(self.coordinates)
                           for name, (dx, dy) in zip(BOND_NAMES, DISPLACEMENTS))
        self.edge_index = {frozenset((a, b)): i for i, (_, a, b) in enumerate(self.edges)}
        self.ends = np.asarray([(a, b) for _, a, b in self.edges], dtype=int)
        self.source_site = 0
        self.source_edge = 2
        self.face_sites = []
        self.face_edges = []
        nnn = []
        for x, y in self.coordinates:
            face = (self.site(x,y,0), self.site(x,y,1), self.site(x+1,y-1,0),
                    self.site(x,y-1,1), self.site(x,y-1,0), self.site(x-1,y,1))
            edges = tuple(self.edge_index[frozenset((face[i],face[(i+1)%6]))]
                          for i in range(6))
            self.face_sites.append(face)
            self.face_edges.append(edges)
            nnn.extend((face[i],face[(i+2)%6],edges[i],edges[(i+1)%6]) for i in range(6))
        self.nnn = np.asarray(nnn, dtype=int)
        self.faces = tuple(sum(1 << i for i in edges) for edges in self.face_edges)
        self.loops = tuple(self.period_loop(*period) for period in self.columns)
        word = [4*i+BOND_NAMES.index(name) for name,a,b in self.edges for i in (a,b)]
        word.extend(4*i+3 for i in range(self.N))
        self.clifford_permutation = parity_permutation(word)
        phase = (1,-1j,-1,1j)[(len(self.edges)+self.N//2)%4]
        if phase not in (-1, 1):
            raise ArithmeticError("nonreal physical parity coefficient")
        self.parity_sign = int(self.clifford_permutation*complex(phase).real)
        self.check_graph()

    def key(self, x, y):
        return tuple((u*x+v*y) % self.cells for u,v in self.adj)

    def site(self, x, y, sub):
        return 2*self.index[self.key(x,y)]+sub

    def period_loop(self, dx, dy):
        x = y = mask = 0
        steps = [1 if dx >= 0 else -1]*abs(dx)+[2 if dy >= 0 else -2]*abs(dy)
        for step in steps:
            start = self.site(x,y,0)
            if abs(step) == 1:
                if step > 0:
                    mid = self.site(x,y,1)
                    x += 1
                else:
                    x -= 1
                    mid = self.site(x,y,1)
            else:
                if step > 0:
                    mid = self.site(x-1,y+1,1)
                    y += 1
                else:
                    y -= 1
                    mid = self.site(x-1,y+1,1)
            finish = self.site(x,y,0)
            mask ^= (1 << self.edge_index[frozenset((start,mid))])
            mask ^= (1 << self.edge_index[frozenset((finish,mid))])
        if (x,y) != (dx,dy) or self.site(x,y,0) != self.source_site:
            raise ArithmeticError("period path lift/closure mismatch")
        return mask

    def check_graph(self):
        if len(self.edge_index) != 3*self.cells or len(self.edges) != 3*self.cells:
            raise ArithmeticError("duplicate nearest bond")
        adjacency = [set() for _ in range(self.N)]
        for _,a,b in self.edges:
            adjacency[a].add(b)
            adjacency[b].add(a)
        seen = {0}
        pending = [0]
        while pending:
            for j in adjacency[pending.pop()]:
                if j not in seen:
                    seen.add(j)
                    pending.append(j)
        if len(seen) != self.N or any(len(v) != 3 for v in adjacency):
            raise ArithmeticError("connectivity/degree check")
        if any(len(set(v)) != 6 for v in self.face_sites):
            raise ArithmeticError("nonhexagonal face")
        if (rank_binary(self.faces) != self.cells-1 or reduce(xor,self.faces) != 0
                or rank_binary(self.faces+self.loops) != self.cells+1):
            raise ArithmeticError("face/Wilson cycle completeness")
        if len({frozenset((int(i),int(j))) for i,j,_,_ in self.nnn}) != 6*self.cells:
            raise ArithmeticError("duplicate oriented next-neighbor pair")
        self.changed_faces = [i for i,f in enumerate(self.faces) if f & (1 << self.source_edge)]
        if len(self.changed_faces) != 2 or sum(self.source_edge in (e,f) for _,_,e,f in self.nnn) != 4:
            raise ArithmeticError("source does not flip two faces/four NNN paths")
        for sites,edges in zip(self.face_sites,self.face_edges):
            missing = [next(name for name in BOND_NAMES if name not in
                            (self.edges[edges[i-1]][0],self.edges[edges[i]][0]))
                       for i in range(6)]
            spin_D = [4*s+j for s,name in zip(sites,missing) for j in (BOND_NAMES.index(name),3)]
            spin_D.extend(4*s+j for s in sites for j in range(4))
            link_word = [4*s+BOND_NAMES.index(self.edges[e][0]) for e in edges for s in self.edges[e][1:]]
            if clifford_word(spin_D) != clifford_word(link_word):
                raise ArithmeticError("signed physical face/Clifford link identity")


    def flipped(self, links):
        result = np.asarray(links, dtype=np.int8).copy()
        result[self.source_edge] *= -1
        return result

    def sparse(self, links):
        from scipy.sparse import coo_matrix
        links=np.asarray(links)
        if links.shape!=(3*self.cells,) or not np.all((links==1)|(links==-1)):
            raise ValueError("invalid complete link assignment")
        a,b=self.ends.T; i,j,e,f=self.nnn.T
        rows=np.r_[a,b,i,j]; cols=np.r_[b,a,j,i]
        data=np.r_[2.*links,-2.*links,.5*links[e]*links[f],-.5*links[e]*links[f]]
        result=coo_matrix((data,(rows,cols)),shape=(self.N,self.N)).tocsr()
        result.sum_duplicates();result.sort_indices()
        return result

    def matrix(self, links):
        return self.sparse(links).toarray()

    def link_record(self, links):
        links = np.asarray(links, dtype=np.int8)
        return dict(sha256=hashlib.sha256(links.tobytes()).hexdigest(),
                    encoding="int8 signs; canonical cell key then x,y,z",
                    negative_count=int(np.count_nonzero(links < 0)),
                    required_parity=int(self.parity_sign*np.prod(links, dtype=np.int64)),
                    fluxes=[int(np.prod(links[list(edges)])) for edges in self.face_edges],
                    wilson=[int(np.prod(links[[i for i in range(len(links)) if (mask>>i)&1]])) for mask in self.loops])




@lru_cache(maxsize=8)
def geometry(label):
    return Geometry(label)


