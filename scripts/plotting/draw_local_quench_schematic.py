#!/usr/bin/env python3
"""Draw a vector-only conceptual figure; no numerical scientific data are used.

The patch is a rotated honeycomb patch. The four affected next-nearest-neighbor
paths are j-i-a, j-i-b, i-j-c, and i-j-d. Bar partition and all marker heights, including the paired-depth examples, are
illustrative layout choices, not evaluated observables. Original quantitative
figure files are never read or modified.

Run: python3 scripts/plotting/draw_local_quench_schematic.py
Requires Matplotlib and NumPy. Outputs PDF, SVG, and a review PNG.
"""
from pathlib import Path
import argparse
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, default=ROOT/'outputs'/'plots')
OUT = parser.parse_args().output_dir
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    'font.family': 'DejaVu Serif', 'font.size': 8.5,
    'mathtext.fontset': 'dejavuserif', 'pdf.fonttype': 42,
    'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.linewidth': 0.7, 'savefig.facecolor': 'white',
})
ink='#20252b'; gray='#66717a'; pale='#eef2f5'
one='#216b83'; higher='#efc17b'; changed='#c96b25'; nnn='#5d78a0'
TOTAL_HEIGHT=4.90
TOP_HEIGHT=3.30
BOTTOM_HEIGHT=TOTAL_HEIGHT-TOP_HEIGHT
TOP_BASE=BOTTOM_HEIGHT/TOTAL_HEIGHT
TOP_SCALE=TOP_HEIGHT/TOTAL_HEIGHT
fig=plt.figure(figsize=(7.4,TOTAL_HEIGHT),facecolor='white')
canvas=fig.add_axes([0,TOP_BASE,1,TOP_SCALE]);canvas.set_xlim(0,1);canvas.set_ylim(0,1);canvas.axis('off')
def label(x,y,s,**kw):
    kw.setdefault('color',ink)
    return canvas.text(x,y,s,transform=canvas.transAxes,**kw)
label(.025,.955,'(a)  Local spin quench',fontsize=10,fontweight='bold',va='top')
label(.345,.955,'(b)  Full response and control',fontsize=10,fontweight='bold',va='top')
label(.755,.955,'(c)  Measured contrast',fontsize=10,fontweight='bold',va='top')
for x in [.324,.735]:
    canvas.plot([x,x],[.075,.91],color='#d8dde1',lw=.7,transform=canvas.transAxes)

# (a) Two hexagons share the vertical selected z link i--j.
patch=fig.add_axes([.025,TOP_BASE+.285*TOP_SCALE,.285,.57*TOP_SCALE])
patch.set_aspect('equal');patch.set_xlim(-2.02,2.02);patch.set_ylim(-1.35,1.60);patch.axis('off')
s3=math.sqrt(3)
pos={'i':(0,.5),'j':(0,-.5),'a':(-s3/2,1),'b':(s3/2,1),
     'c':(-s3/2,-1),'d':(s3/2,-1),'l0':(-s3,.5),'l1':(-s3,-.5),
     'r0':(s3,.5),'r1':(s3,-.5)}
faces=[['i','a','l0','l1','c','j'],['i','j','d','r1','r0','b']]
edges=set()
for face in faces:
    patch.add_patch(Polygon([pos[v] for v in face],facecolor=pale,edgecolor='none',zorder=0))
    for a,b in zip(face,face[1:]+face[:1]):edges.add(tuple(sorted((a,b))))
for a,b in sorted(edges):
    if {a,b}=={'i','j'}:continue
    patch.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color='#7d8790',lw=1.1,zorder=1)
paths=[('j','i','a'),('j','i','b'),('i','j','c'),('i','j','d')]
affected=[]
for a,mid,b in paths:
    affected.append([a,b])
    patch.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color=nnn,lw=1.3,
               linestyle=(0,(3,2)),zorder=2)
patch.plot([0,0],[-.5,.5],color=changed,lw=3,zorder=3)
subA={'i','c','d','l0','r0'}
for n,(x,y) in pos.items():
    patch.plot(x,y,'o',ms=4.0,mfc=ink if n in subA else 'white',mec=ink,mew=.8,zorder=4)
patch.text(-.20,.77,r'$i$',ha='right',va='bottom',fontsize=9)
patch.text(.20,-.78,r'$j$',ha='left',va='top',fontsize=9)
# Identify the orange z bond in the keyed label below the patch, clear of NNN lines.
patch.annotate(r'$S_i^z$',xy=(0,.55),xytext=(0,1.36),ha='center',va='bottom',
               fontsize=10,arrowprops={'arrowstyle':'-|>','lw':.9,'color':ink})
patch.text(-1.07,0,r'$W_L\mapsto -W_L$',ha='center',va='center',fontsize=7.4)
patch.text(1.07,0,r'$W_R\mapsto -W_R$',ha='center',va='center',fontsize=7.4)
label(.168,.858,'Frozen flux background',ha='center',fontsize=8.5)
label(.168,.277,r'$z$ bond: $u_{ij}\mapsto-u_{ij}$',ha='center',fontsize=9.3,color=changed)
canvas.add_line(Line2D([.042,.080],[.226,.226],color=nnn,lw=1.3,linestyle=(0,(3,2)),transform=canvas.transAxes))
label(.088,.226,'Four affected next-nearest',fontsize=8,va='center')
label(.088,.182,'hoppings change sign',fontsize=8,va='center')
label(.168,.080,'Same background and quench\nfor both predictions',ha='center',fontsize=8.3,linespacing=1.25)

# (b) Total unfiltered measures, never filtered ordinates, have mass one.
label(.535,.855,'Eligible odd sources',ha='center',fontsize=9,fontweight='bold')
label(.535,.792,'Total unfiltered intensity = 1',ha='center',fontsize=8.7)
x=.425;width=.285;height=.064;split=.40 # purely illustrative layout fraction
for y in [.680,.550]:
    canvas.add_patch(Rectangle((x,y),width,height,transform=canvas.transAxes,facecolor=one,edgecolor=ink,lw=.6))
canvas.add_patch(Rectangle((x+split*width,.680),(1-split)*width,height,transform=canvas.transAxes,facecolor=higher,edgecolor=ink,lw=.6))
label(.410,.712,'Full',ha='right',va='center',fontsize=8.6)
label(.410,.582,'Control',ha='right',va='center',fontsize=8.6)
label(x+split*width/2,.712,r'$w_1$',ha='center',va='center',color='white',fontsize=9)
label(x+(1+split)*width/2,.712,r'$1-w_1$',ha='center',va='center',fontsize=9)
label(x+width/2,.582,'retained one-particle shape',ha='center',va='center',color='white',fontsize=8.2)
label(.535,.465,r'$S_\gamma=w_1P_\gamma+h_\gamma$',ha='center',fontsize=10)
label(.535,.382,r'$S_\gamma-P_\gamma=h_\gamma-(1-w_1)P_\gamma$',ha='center',fontsize=9.1)
label(.535,.303,r'$h$: not separately unit-normalized',ha='center',fontsize=8.5)
label(.535,.225,r'$P$ retains actual one-particle energies',ha='center',fontsize=8.3)
label(.535,.180,'and relative intensities',ha='center',fontsize=8.3)
label(.535,.105,'Even-source responses are shared.',ha='center',fontsize=8.3)

# (c) Generic discrete markers illustrate the definition, not a spectrum.
probe=fig.add_axes([.785,TOP_BASE+.405*TOP_SCALE,.195,.36*TOP_SCALE])
probe.set_xlim(.72,3.28);probe.set_ylim(.03,1.03)
probe.spines[['top','right','left']].set_visible(False)
probe.spines['bottom'].set_color(gray)
probe.set_yticks([]);probe.set_xticks([1,2,3]);probe.tick_params(axis='x',length=3,pad=2,labelsize=8)
probe.set_xlabel(r'$c/J_\sigma$',fontsize=8.5,labelpad=2)
y1,y2,y3=.71,.29,.89 # illustrative positions only, not physical observations
mean=(y1+y3)/2
probe.plot([1,3],[mean,mean],color=gray,lw=.8,linestyle=(0,(2,2)))
probe.plot([1,1],[y1,mean],color=gray,lw=.7)
probe.plot([3,3],[mean,y3],color=gray,lw=.7)
probe.plot([1,2,3],[y1,y2,y3],linestyle='none',marker='o',ms=4,mfc='white',mec=ink,mew=1)
probe.annotate('',xy=(2.13,mean),xytext=(2.13,y2),arrowprops={'arrowstyle':'<->','lw':.9,'color':changed})
probe.text(2.27,(mean+y2)/2,r'$d(U)$',fontsize=8.6,va='center',color=changed)
probe.text(2,mean+.06,'mean shoulder',ha='center',fontsize=8.1,color=gray)
probe.text(2,y2-.08,'center',ha='center',va='top',fontsize=8.1,color=gray)
label(.865,.827,'Generic filtered response',ha='center',fontsize=8.3)
label(.865,.235,r'$d(U)=\frac{U(1)+U(3)}{2}-U(2)$',ha='center',fontsize=9)
label(.865,.115,r'$D=d(P)-d(S)$',ha='center',fontsize=10)


# (d) Two invented triplets illustrate only the ordering of paired depths.
# These are layout choices, not sampled sources, averages, median spectra, or data.
# Scientific result assignments and schematic qualifications are in the caption.
bottom=fig.add_axes([0,0,1,TOP_BASE]);bottom.set_xlim(0,1);bottom.set_ylim(0,1);bottom.axis('off')
def dlabel(x,y,s,**kw):
    kw.setdefault('color',ink)
    return bottom.text(x,y,s,transform=bottom.transAxes,**kw)
# Place panel (d) in physical inches so removing the result bands cannot shrink
# the plot area, endpoint guides, or type. Panels (a)--(c) retain their 3.30 in.
def dy(inches):
    return inches/BOTTOM_HEIGHT
bottom.plot([.025,.98],[dy(1.58),dy(1.58)],color='#d8dde1',lw=.7)
dlabel(.025,dy(1.51),'(d)  Sign of the paired correction',fontsize=10,fontweight='bold',va='top')
# Legend shape/color is shared by both examples.
bottom.plot(.67,dy(1.46),'o',ms=3.8,mfc='white',mec=one,mew=1,transform=bottom.transAxes)
dlabel(.685,dy(1.46),r'$P$: matched control',fontsize=8.2,va='center')
bottom.plot(.872,dy(1.46),'s',ms=3.8,mfc=changed,mec=changed,mew=.8,transform=bottom.transAxes)
dlabel(.887,dy(1.46),r'$S$: full',fontsize=8.2,va='center')

dlabel(.25,dy(1.18),r'Shallowing: $d(S)<d(P)$',ha='center',fontsize=9.1,fontweight='bold')
dlabel(.75,dy(1.18),r'Deepening: $d(S)>d(P)$',ha='center',fontsize=9.1,fontweight='bold')
bottom.plot([.50,.50],[dy(.32),dy(1.17)],color='#d8dde1',lw=.7)

examples=[(.075,(.85,.40,.85),(.60,.30,.60)),
          (.575,(.65,.40,.65),(.95,.50,.95))]
for left,P,S in examples:
    # 0.7564 in = the unchanged 0.305 x 2.48 in height of Figure Revision 2.
    a=fig.add_axes([left,.31/TOTAL_HEIGHT,.35,.7564/TOTAL_HEIGHT])
    a.set_xlim(.10,3.90);a.set_ylim(.12,1.04)
    a.spines[['top','right','left']].set_visible(False)
    a.spines['bottom'].set_color(gray)
    a.set_yticks([]);a.set_xticks([1,2,3]);a.tick_params(axis='x',length=2.5,pad=2,labelsize=7.7)
    a.set_xlabel(r'$c/J_\sigma$',fontsize=8,labelpad=1)
    # Put depth dimensions in dedicated outer lanes. Mean-shoulder and
    # central-level extensions meet each arrow at its exact endpoint height.
    # Opposite-color guides never cross an arrow or its label.
    for values,color,marker,filled,arrow_x,label_x,align,side in [
          (P,one,'o',False,.72,.52,'right','left'),
          (S,changed,'s',True,3.28,3.48,'left','right')]:
        shoulder=(values[0]+values[2])/2
        a.plot([1,3],[shoulder,shoulder],color=color,lw=.85,linestyle=(0,(2,2)),alpha=.85)
        shoulder_anchor=1 if side=='left' else 3
        a.plot([arrow_x,shoulder_anchor],[shoulder,shoulder],color=color,lw=.65)
        a.plot([arrow_x,2],[values[1],values[1]],color=color,lw=.65)
        a.plot([1,2,3],values,linestyle='none',marker=marker,ms=4.0,
               mfc=color if filled else 'white',mec=color,mew=.9,zorder=3)
        a.annotate('',xy=(arrow_x,shoulder),xytext=(arrow_x,values[1]),
                   arrowprops={'arrowstyle':'<->','lw':.9,'color':color,
                               'shrinkA':0,'shrinkB':0,'mutation_scale':8.5})
        term='P' if marker=='o' else 'S'
        a.text(label_x,(shoulder+values[1])/2,r'$d('+term+r')$',ha=align,va='center',
               fontsize=8.4,color=color)
    # The guides are separate channel-specific mean shoulders, not fitted spectra.

base=OUT/'local_quench_schematic'
fig.savefig(str(base)+'.pdf',metadata={'Title':'Local spin quench and intensity-matched comparison','Author':'David C. Ronquillo','Subject':'Conceptual vector schematic; no numerical data'})
fig.savefig(str(base)+'.svg',metadata={'Title':'Local spin quench and intensity-matched comparison','Description':'Conceptual vector schematic; no numerical data'})
fig.savefig(str(base)+'.png',dpi=200)
plt.close(fig)
print(json.dumps({'pdf':str(base)+'.pdf','svg':str(base)+'.svg','review_png':str(base)+'.png','size_inches':[7.4,TOTAL_HEIGHT],'affected_NNN_paths':paths,'affected_NNN_endpoints':affected,'scientific_data_used':False,'illustrative_depth_triplets':examples}))
