"""Draw five diagnostic views from saved scalar tables; these are not the approved manuscript artwork."""
import csv
import json
from pathlib import Path
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"data"/"diagnostics"
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output-dir",type=Path,default=ROOT/"outputs"/"plots")
FIG=parser.parse_args().output_dir
FIG.mkdir(parents=True,exist_ok=True)
A=json.loads((DATA/"plot_summary.json").read_text())
TEAL="#09697A"
ORANGE="#B54D27"
NAVY="#243C5A"
PURPLE="#785383"
GREY="#849097"
INK="#192D36"
G005=("g005",TEAL,"o",-.12)
G010=("g010",ORANGE,"s",.12)
plt.rcParams.update({"font.family":"DejaVu Sans","mathtext.fontset":"dejavusans",
    "font.size":8.3,"axes.labelsize":8.3,"axes.titlesize":9.0,"axes.titleweight":"semibold",
    "axes.labelcolor":INK,"text.color":INK,"xtick.color":INK,"ytick.color":INK,
    "axes.spines.top":False,"axes.spines.right":False,"axes.edgecolor":"#94a1a7",
    "axes.linewidth":.65,"xtick.major.width":.6,"ytick.major.width":.6,
    "xtick.major.size":3.0,"ytick.major.size":3.0,"pdf.fonttype":42,"savefig.facecolor":"white"})


def rows(name):
    with (DATA/name).open() as f:return list(csv.DictReader(f))


def interval(ax,x,lo,hi,colour,marker=None,filled=True,ms=4):
    ax.errorbar(x,(lo+hi)/2,yerr=(hi-lo)/2,fmt=marker or "none",
        mec=colour,mfc=colour if filled else "white",color=colour,lw=1.2,
        ms=ms,capsize=2.1,capthick=1.05,zorder=4)


def marker(c,m,txt,filled=True):
    return Line2D([],[],ls="none",marker=m,ms=4.7,mec=c,mfc=c if filled else "white",label=txt)


def title(ax,letter,txt):
    ax.set_title(f"({letter})  {txt}",loc="left",pad=6)


def save(f,name):
    f.savefig(FIG/(name+".pdf"),bbox_inches="tight",pad_inches=.045,
              metadata={"Title":name.replace("_"," "),"Creator":"Matplotlib"})
    plt.close(f)


def figure1():
    f,(a,b)=plt.subplots(1,2,figsize=(7.12,2.82),gridspec_kw={"width_ratios":[1.38,1]})
    f.subplots_adjust(left=.087,right=.982,bottom=.248,top=.818,wspace=.44)
    title(a,"a","Odd-source median depth correction")
    a.axhline(0,color=INK,lw=.8)
    a.axvline(5.53,color="#b7c1c5",ls=":",lw=.85)
    geoms=("C3n02","C3n04","C3n08","C3n12","C3n16","R24x32")
    for w,col,m,dx in [G005,G010]:
        l,u=A["exact_N24_medians"][w]["D"];interval(a,dx,1000*l,1000*u,col,"D",ms=3.7)
        for x,geom in zip([1,2,3,4,5,6.17],geoms):
            v=next(v for v in A["medians"] if (v["width"],v["geometry"],v["endpoint"])==(w,geom,"D"))
            interval(a,x+dx,1000*v["lo"],1000*v["hi"],col)
    a.set_xticks([0,1,2,3,4,5,6.17],["24","96","384","1536","3456","6144","1536\nrect."])
    a.tick_params(axis="x",labelsize=7)
    a.set(xlim=(-.4,6.62),xlabel="Sites, $N$",ylabel=r"$D$  ($10^{-3}$ of total intensity)")
    title(b,"b","Signs on the same backgrounds")
    fs=rows("N6144_sign_frequencies.csv")
    b.axvline(85,color=GREY,ls=":",lw=1)
    for w,col,m,dy in [G005,G010]:
        for y,event in zip([2,1,0],["X>0","D<0","joint"]):
            r=next(r for r in fs if (r["width"],r["event"])==(w,event))
            lo,hi=int(r["definite"]),int(r["possible"])
            b.plot([lo,hi],[y-dy,y-dy],lw=2.2,color=col)
            b.plot(lo,y-dy,marker=m,ms=3.6,mfc=col,mec=col)
            b.plot(hi,y-dy,marker="|",ms=6,mec=col)
    b.set(yticks=[2,1,0],yticklabels=[r"$X(2)>0$",r"$D<0$","Both"],xlim=(62,148),ylim=(-.48,2.52),
          xlabel="Count among 170 odd sources",xticks=[70,85,100,115,130,145])
    b.tick_params(axis="y",length=0)
    b.text(.98,.97,"$N=6144$",transform=b.transAxes,ha="right",va="top",fontsize=7.6)
    f.legend(handles=[marker(TEAL,"o",r"$\gamma/J_\sigma=0.05$"),marker(ORANGE,"s",r"$\gamma/J_\sigma=0.10$")],
             loc="upper center",bbox_to_anchor=(.51,1.035),ncol=2,frameon=False,columnspacing=2)
    save(f,"depth_reversal")


def figure2():
    f,ax=plt.subplots(2,3,figsize=(7.15,3.95),sharey="row",sharex=True)
    f.subplots_adjust(left=.085,right=.985,top=.82,bottom=.15,wspace=.15,hspace=.25)
    ss=rows("absolute_model_ordinates.csv")
    types=[("C3n01","exact_odd_mean","24 sites\nOdd ensemble"),
           ("C3n16","fixed_odd_mean","6,144 sites\nOdd cohort"),
           ("C3n01","exact_all_mean","24 sites\nAll backgrounds")]
    for j,(g,stat,label) in enumerate(types):
        for i,w in enumerate(("g005","g010")):
            a=ax[i,j]
            if i==0:a.set_title(label,pad=5)
            a.text(.035,.07,"("+chr(97+i*3+j)+")",transform=a.transAxes,ha="left",va="bottom",fontsize=7.6,color=INK)
            gains = []
            for c in (1,2,3):
                r=next(r for r in ss if (r["geometry"],r["estimand"],r["width"],int(r["c_Jsigma"]))==(g,stat,w,c))
                s_mid = (float(r["S_abs_lo"]) + float(r["S_abs_hi"])) / 2
                p_mid = (float(r["P_abs_lo"]) + float(r["P_abs_hi"])) / 2
                gains.append(s_mid - p_mid)
                for name,dx,col,m,fill in [("P",-.09,NAVY,"s",False),("S",.09,ORANGE,"o",True)]:
                    interval(a,c+dx,1000*float(r[name+"_abs_lo"]),1000*float(r[name+"_abs_hi"]),col,m,fill,4.1)
            depth = 1000 * (gains[1] - (gains[0] + gains[2]) / 2)
            a.text(.99, .075, rf"$10^3\overline{{D}}_{{\rm abs}}\simeq {depth:+.2f}$",
                   transform=a.transAxes, ha="right", va="bottom", fontsize=7.1,
                   bbox=dict(facecolor="white", edgecolor="none", pad=.4, alpha=.92))
            a.set(xticks=[1,2,3],xlim=(.58,3.42),ylim=(0,11.1 if i==0 else 20.6))
            a.grid(axis="y",color="#e5eaec",lw=.5);a.set_axisbelow(True)
            if j==0:a.set_ylabel(r"$10^3 I_{\rm abs}$")
            if j==2:a.text(.97,.955,r"$\gamma/J_\sigma=$"+("0.05" if i==0 else "0.10"),transform=a.transAxes,
                           ha="right",va="top",fontsize=7.1)
    f.text(.54,.013,r"Filter center  $c/J_\sigma$",ha="center")
    f.legend(handles=[marker(NAVY,"s","Intensity-matched",False),marker(ORANGE,"o","Full response")],
             loc="upper center",bbox_to_anchor=(.52,1.038),ncol=2,frameon=False,columnspacing=2)
    save(f,"absolute_spin_intensities")


def figure3():
    f,(a,b)=plt.subplots(1,2,figsize=(7.12,2.95))
    f.subplots_adjust(left=.082,right=.981,bottom=.215,top=.88,wspace=.365)
    title(a,"a","Paired width ratios")
    r=A["width_diagnostic"]["c200"]["c200_rows_for_figure"]
    p=np.array([x["P"] for x in r]);h=np.array([x["h"] for x in r])
    pm,hm=p.mean(1),h.mean(1)
    lines=[[(x,l),(x,u)] for x,(l,u) in zip(pm,h)]+[[(l,y),(u,y)] for y,(l,u) in zip(hm,p)]
    a.add_collection(LineCollection(lines,colors="#648b97",alpha=.42,linewidths=.55))
    a.scatter(pm,hm,s=5,color=TEAL,alpha=.70,lw=0,zorder=3)
    a.plot([1,4],[1,4],":",color=GREY,lw=1);a.axhline(2,ls="--",color=ORANGE,lw=.85)
    a.set(xlim=(.91,4.1),ylim=(.93,4.1),xticks=[1,2,3,4],yticks=[1,2,3,4],xlabel=r"One particle  $R_P$",ylabel=r"Higher particles  $R_h$")
    title(b,"b","Within-window fraction difference")
    ww=sorted(rows("window_fraction_bounds_N6144.csv"),key=lambda r:float(r["f_h_minus_f_P_lo"]))
    xx=np.arange(1,len(ww)+1)
    for tag,color in [("unresolved","#92a3ac"),("higher",TEAL),("one",ORANGE)]:
        ind=[i for i,r in enumerate(ww) if r["forced_order"]==tag]
        if not ind:continue
        lows=np.array([float(ww[i]["f_h_minus_f_P_lo"]) for i in ind])
        highs=np.array([float(ww[i]["f_h_minus_f_P_hi"]) for i in ind])
        b.vlines(xx[ind],lows,highs,color=color,alpha=.53,lw=.72,zorder=1)
        b.plot(xx[ind],lows,ls="none",marker=".",mec=color,mfc=color,ms=1.65,zorder=2)
    b.axhline(0,color=INK,lw=.9,zorder=3)
    b.set(xlim=(0,171),ylim=(-1.04,1.06),xticks=[1,50,100,150,170],yticks=[-1,-.5,0,.5,1],
          xlabel="Background, ordered by lower bound",ylabel=r"$f_h-f_P$")
    b.text(.98,.075,"37 strictly positive\n133 overlap zero",transform=b.transAxes,ha="right",fontsize=7.1,
           bbox=dict(boxstyle="square,pad=.19",facecolor="white",edgecolor="none",alpha=.92))
    save(f,"filter_ratio_and_window")


def figure4():
    f,a=plt.subplots(figsize=(3.42,2.55))
    f.subplots_adjust(left=.18,right=.98,bottom=.215,top=.85)
    ss=rows("physical_050_comparison.csv")
    for name,col,m,dx,fill in [("odd_q1_ablation",NAVY,"s",-.14,False),("full_frozen",ORANGE,"o",0,True),
                              ("physical",TEAL,"^",.14,True)]:
        for r in (r for r in ss if r["quantity"]==name):
            interval(a,int(r["center_Jsigma"])+dx,1000*float(r["lower"]),1000*float(r["upper"]),col,m,fill,4.4)
    a.set(xticks=[1,2,3],xlim=(.55,3.45),ylim=(2.2,8.8),xlabel=r"Filter center  $c/J_\sigma$",ylabel=r"$10^3 I_{\rm abs}$")
    a.grid(axis="y",color="#e1e8eb",lw=.55);a.set_axisbelow(True)
    f.legend(handles=[marker(NAVY,"s","Matched",False),marker(ORANGE,"o","Full"),marker(TEAL,"^","Physical")],
             loc="upper center",bbox_to_anchor=(.53,1.035),ncol=3,frameon=False,columnspacing=.7,handletextpad=.15,fontsize=7.5)
    save(f,"prespecified_model_benchmark")


def figure5():
    f,(a,b)=plt.subplots(1,2,figsize=(7.12,2.7),gridspec_kw={"width_ratios":[1.1,1]})
    f.subplots_adjust(left=.078,right=.982,bottom=.21,top=.81,wspace=.37)
    title(a,"a","Coherent physical response")
    tab=rows("physical_field_curves.csv"); minima=rows("physical_minimum_candidates.csv")
    for field,col,sym in [("050",TEAL,"o"),("055",PURPLE,"D")]:
        t=[r for r in tab if r["field"]==field]; x=np.array([float(r["center_Jsigma"]) for r in t])
        l=np.array([1000*float(r["Y_lo"]) for r in t]);u=np.array([1000*float(r["Y_hi"]) for r in t])
        a.fill_between(x,l,u,color=col,alpha=.28)
        a.plot(x,(l+u)/2,color=col,lw=1.1)
        ca=np.array([float(r["center_Jsigma"]) for r in minima if r["field"]==field and r["layer"]=="raw"])
        yy=[(l+u)[np.argmin(abs(x-xc))]/2 for xc in ca]
        a.scatter(ca,yy,s=11,marker=sym,facecolors="white",edgecolors=col,lw=.7,zorder=4)
    a.set(xlim=(1.5,3.5),xticks=[1.5,2,2.5,3,3.5],xlabel=r"Filter center  $c/J_\sigma$",ylabel=r"$10^3 Y_{\rm abs}$")
    title(b,"b","Admissible fine-grid shifts")
    for val in (0,.20):b.axvline(val,color=GREY,lw=1,ls=":")
    for tag,y,m,col in [("raw",1,"o",INK),("inelastic",0,"D",INK)]:
        x=A["physical_fine_displacements"][tag]
        b.scatter(x,[y]*len(x),s=14,marker=m,facecolors=col if y else "white",edgecolors=col,lw=.8,zorder=3)
    b.set(xlim=(-.008,.208),ylim=(-.62,1.65),yticks=[0,1],yticklabels=["Elastic\nremoved","Raw"],
          xticks=[0,.05,.1,.15,.2],xticklabels=["0",".05",".10",".15",".20"],
          xlabel=r"Displacement  $\Delta c/J_\sigma$")
    b.text(.2,1.46,"Bare Zeeman",ha="right",fontsize=7.15,color="#596a72");b.text(0,1.46,"Stationary",ha="left",fontsize=7.15,color="#596a72")
    b.spines["left"].set_visible(False);b.tick_params(axis="y",length=0)
    f.legend(handles=[Line2D([],[],color=TEAL,lw=1.4,label=r"$b/J=0.50$"),
                      Line2D([],[],color=PURPLE,lw=1.4,label=r"$b/J=0.55$")],frameon=False,ncol=2,
             loc="upper center",bbox_to_anchor=(.28,1.04),fontsize=7.7)
    save(f,"physical_field_displacement")


if __name__=="__main__":
    for plot in (figure1,figure2,figure3,figure4,figure5):plot()
