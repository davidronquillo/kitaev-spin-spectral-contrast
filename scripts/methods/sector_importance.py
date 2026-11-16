"""Reference sector importance estimator for an odd pure fermionic Gaussian chart.

Energy units are J_sigma. Uses no repository data. Run this file to check against
finite exact proposal sums and independent Slater-minor amplitudes. --benchmark
adds a reproducible synthetic throughput observation.
The implementation produces individual bounded, unbiased draws for a fixed
represented chart. Statistical and physical/numerical enclosures are external.
"""
import argparse
import itertools
import math
import time
from functools import lru_cache

import numpy as np


@lru_cache(maxsize=None)
def pf_matchings(vertices):
    """Signed perfect matchings in increasing vertex order."""
    if not vertices:
        return ((1, ()),)
    out = []
    first = vertices[0]
    for position in range(1, len(vertices)):
        second = vertices[position]
        rest = vertices[1:position] + vertices[position + 1:]
        for sign, pairs in pf_matchings(rest):
            parity = 1 if position % 2 else -1
            out.append((parity * sign, ((first, second),) + pairs))
    return tuple(out)


def newton_remainder(tail5_upper, z3_lower, first_omitted):
    """Upper mass at odd q >= first_omitted from T(q>=5) and Z3.

    Ordinary double precision for reference only. Widen bounds in production.
    """
    if first_omitted < 5 or first_omitted % 2 != 1:
        raise ValueError("first omitted sector must be odd and >= 5")
    if z3_lower <= 0:
        return tail5_upper
    ratio = tail5_upper / z3_lower
    result = tail5_upper
    for k in range(2, (first_omitted - 1) // 2):
        result *= min(1.0, 2 * ratio / (k + 1))
    return result


def proposal_table(weights):
    """Positive integer masses summing to 2**52, including sub-resolution atoms.

    Drawing an integer against the integer CDF induces the *recorded* probabilities
    exactly. Keep those probabilities in the importance denominator; they can differ
    slightly from the floating-point singleton/edge weights used to select them.
    """
    resolution = 2**52
    total = math.fsum(weights)
    if total <= 0 or not math.isfinite(total):
        raise ValueError("empty or nonfinite proposal mass")
    masses = np.floor((resolution - len(weights))*weights/total).astype(np.int64) + 1
    first = int(np.argmax(weights))
    masses[first] += resolution - int(masses.sum(dtype=np.int64))
    if np.min(masses) <= 0:
        raise ArithmeticError("integer proposal normalization failed")
    cdf = np.cumsum(masses, dtype=np.int64)
    if cdf[-1] != resolution:
        raise ArithmeticError("integer CDF support")
    return cdf, masses.astype(float) / resolution


class SectorProposal:
    """Draw coherent sector weights using a singleton and independent mode pairs."""

    def __init__(self, x, pair, z1, energies, q):
        if q not in (5, 7):
            raise ValueError("the supported sectors are q=5 and q=7")
        self.q = q
        self.k = (q - 1) // 2
        self.x = np.asarray(x)
        self.pair = np.asarray(pair)
        self.energies = np.asarray(energies)
        self.z1 = float(z1)
        n = len(self.x)
        if n < q or self.pair.shape != (n, n) or len(self.energies) != n:
            raise ValueError("incompatible chart dimension or particle count")
        if not np.allclose(self.pair, -self.pair.T, rtol=1e-12, atol=1e-14):
            raise ValueError("pair must be skew-symmetric in ordinary transpose")
        self.i, self.j = np.triu_indices(n, 1)
        singleton_mass = abs(self.x) ** 2
        edge_mass = abs(self.pair[self.i, self.j]) ** 2
        self.norm_x2 = math.fsum(singleton_mass)
        self.b = math.fsum(edge_mass)
        if min(self.z1, self.norm_x2, self.b) <= 0:
            raise ValueError("singular or empty proposal; break out zero-weight cases")
        self.singleton_cdf, self.singleton_prob = proposal_table(singleton_mass)
        self.edge_cdf, self.edge_prob = proposal_table(edge_mass)
        self.edge_prob_matrix = np.zeros((n,n))
        self.edge_prob_matrix[self.i, self.j] = self.edge_prob
        self.terms = pf_matchings(tuple(range(q + 1)))
        scale1 = float(np.max(singleton_mass / self.singleton_prob))
        scale2 = float(np.max(edge_mass / self.edge_prob))
        self.weight_cap = self.z1 * scale1 * scale2**self.k * len(self.terms) / math.factorial(self.k)

    def subset_quantities(self, indices):
        """For increasing, distinct rows return amplitudes, weights and proposal."""
        indices = np.asarray(indices)
        if indices.ndim != 2 or indices.shape[1] != self.q:
            raise ValueError("wrong mode tuple shape")
        if np.any(np.diff(indices, axis=1) <= 0):
            raise ValueError("expected strictly increasing mode indices")
        # Vertex zero is the auxiliary blocked orbital; vertices 1..q are modes.
        factors = {(0, a + 1): self.x[indices[:, a]] for a in range(self.q)}
        prob_factors = {(0, a + 1): self.singleton_prob[indices[:, a]] for a in range(self.q)}
        for a, b in itertools.combinations(range(self.q), 2):
            factors[(a + 1, b + 1)] = self.pair[indices[:, a], indices[:, b]]
            prob_factors[(a + 1, b + 1)] = self.edge_prob_matrix[indices[:, a], indices[:, b]]
        amplitude = np.zeros(len(indices), complex)
        envelope = np.zeros(len(indices))
        for sign, matching in self.terms:
            terms = factors[matching[0]].copy()
            proposal_term = prob_factors[matching[0]].copy()
            for edge in matching[1:]:
                terms *= factors[edge]
                proposal_term *= prob_factors[edge]
            amplitude += sign * terms
            envelope += proposal_term
        amplitude *= math.sqrt(self.z1)
        probabilities = math.factorial(self.k) * envelope
        return amplitude, abs(amplitude) ** 2, probabilities

    def draw(self, rng, size):
        """Return (importance weights, total energies). Collisions stay as draws."""
        single = np.searchsorted(self.singleton_cdf, rng.integers(0, 2**52, size=size), side="right")
        edges = np.searchsorted(self.edge_cdf, rng.integers(0, 2**52, size=(size, self.k)), side="right")
        indices = np.column_stack((single, self.i[edges], self.j[edges]))
        indices.sort(axis=1)
        distinct = np.all(np.diff(indices, axis=1) > 0, axis=1)
        valid = indices[distinct]
        _, w, p = self.subset_quantities(valid)
        if np.any(p <= 0) or np.any(~np.isfinite(p)):
            raise ArithmeticError("unresolved proposal mass")
        weights = np.zeros(size)
        energies = np.zeros(size)
        weights[distinct] = w / p
        energies[distinct] = np.sum(self.energies[valid], axis=1)
        if np.any(weights > self.weight_cap * (1 + 1e-10)):
            raise ArithmeticError("importance weight exceeds Cauchy cap")
        return weights, energies


def kernels(energy):
    width = 0.05
    center, left, right = [
        width ** 2 / (width ** 2 + (energy - c) ** 2) for c in (2, 1, 3)
    ]
    return np.column_stack((np.ones_like(energy), center, center - (left + right) / 2))


def toy_chart(n, pair_odds, seed):
    """Use an isometry so canonical pairs and blocked mode mix coherently."""
    rng = np.random.default_rng(seed)
    m = 1 + 2 * len(pair_odds)
    u, _ = np.linalg.qr(rng.normal(size=(n, m)) + 1j * rng.normal(size=(n, m)))
    a = u[:, 1::2] * pair_odds
    b = a @ u[:, 2::2].T
    b = b - b.T
    z1 = 1 / np.prod(1 + abs(pair_odds) ** 2)
    return u, b, z1


def check_exact():
    n = 7
    odds = np.array([0.24 + 0.53j, -0.38 + 0.19j, 0.27 - 0.16j])
    u, b, z1 = toy_chart(n, odds, 16601)
    energies = np.array([0.11, 0.18, 0.24, 0.33, 0.46, 0.72, 0.93])
    sector_weights = []
    for k in range(len(odds) + 1):
        z = z1 * sum(np.prod(abs(odds[list(s)]) ** 2)
                     for s in itertools.combinations(range(len(odds)), k))
        sector_weights.append(z)
    np.testing.assert_allclose(sum(sector_weights), 1, atol=1e-14)
    for q in (5, 7):
        proposal = SectorProposal(u[:, 0], b, z1, energies, q)
        subsets = list(itertools.combinations(range(n), q))
        indices = np.array(subsets)
        amplitudes, w, p = proposal.subset_quantities(indices)
        minors = np.zeros(len(subsets), complex)
        for s in itertools.combinations(range(len(odds)), proposal.k):
            columns = [0] + [j for t in s for j in (2 * t + 1, 2 * t + 2)]
            coefficient = math.sqrt(z1) * np.prod(odds[list(s)])
            for i, subset in enumerate(subsets):
                minors[i] += coefficient * np.linalg.det(u[np.ix_(subset, columns)])
        np.testing.assert_allclose(amplitudes, minors, rtol=1e-12, atol=2e-15)
        qsums = dict.fromkeys(subsets, 0.0)
        invalid = 0.0
        edges = list(zip(proposal.i.tolist(), proposal.j.tolist()))
        edge_prob = proposal.edge_prob
        for single in range(n):
            single_prob = proposal.singleton_prob[single]
            for matching in itertools.product(range(len(edges)), repeat=proposal.k):
                tup = tuple(sorted([single] + [a for e in matching for a in edges[e]]))
                prob = single_prob * math.prod(edge_prob[j] for j in matching)
                if len(set(tup)) < q:
                    invalid += prob
                else:
                    qsums[tup] += prob
        p_exact = np.array([qsums[tup] for tup in subsets])
        np.testing.assert_allclose(p, p_exact, atol=2e-15, rtol=1e-12)
        np.testing.assert_allclose(sum(p) + invalid, 1., atol=2e-13, rtol=0.)
        y = (w / p)[:, None] * kernels(np.sum(energies[indices], axis=1))
        np.testing.assert_allclose(p_exact @ y, abs(minors) ** 2 @ kernels(
            np.sum(energies[indices], axis=1)), atol=1e-15)
        assert np.max(w / p) <= proposal.weight_cap * (1 + 1e-13)
        np.testing.assert_allclose(sum(w), sector_weights[proposal.k], atol=1e-15)
        print(f"q={q}: exact matching signs, Slater-minor amplitudes, all ordered "
              f"proposal probabilities, collisions and kernel expectations passed")
    for first in (5, 7, 9):
        actual = sum(sector_weights[(first - 1) // 2:])
        cap = newton_remainder(sum(sector_weights[2:]), sector_weights[1], first)
        assert actual <= cap + 2e-15
    print("Gaussian number mass and Newton remainder controls passed")


def benchmark():
    import os
    n = 900
    odds = np.exp(-np.arange(12)/3.) * (0.45+0.2j)
    u, b, z1 = toy_chart(n, odds, 16602)
    energies = np.linspace(.003, 6.2, n)
    rng = np.random.default_rng(16603)
    for q in (5, 7):
        proposal = SectorProposal(u[:, 0], b, z1, energies, q)
        draws = {5: 200000, 7: 40000}[q]
        total, count = np.zeros(3), 0
        start = time.perf_counter()
        while count < draws:
            weights, e = proposal.draw(rng, min(2048, draws-count))
            total += np.sum(weights[:, None] * kernels(e), axis=0)
            count += len(weights)
        elapsed = time.perf_counter()-start
        print(f"q={q}, n={n}, {draws} proposal draws, {elapsed:.3f}s, "
              f"{draws/elapsed:.0f} draws/s, cap {proposal.weight_cap:.4g}, "
              f"checksum {total[0]:.6g}, pid {os.getpid()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    check_exact()
    if args.benchmark:
        benchmark()
