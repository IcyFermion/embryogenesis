"""Terminal Figure S3: four tracking geometries with near-natural insets.

``render`` is copied verbatim from ``terminal_pareto/figS3_cross_geometry.py``.
That script computes and replay-checks the cross-geometry projections and its
source hash is pinned by the pooled-migration validator, so it stays
byte-identical (a documented compatibility exception). This module draws from
its validated projection table and references.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from publication import style

SOURCES = ["embryo1", "embryo2", "embryo3", "pooled"]
LABELS = {"embryo1":"Embryo 1", "embryo2":"Embryo 2", "embryo3":"Embryo 3", "pooled":"Pooled (all three)"}
COLORS = dict(zip(SOURCES,["#4477AA", "#CC6677", "#228833", "#202020"]))
INSET_OFFSETS = {"travel": [-0.026, 0.036], "state": [-0.09, 0.09]}


def render(frame,references,zoom,insets=False,*,out_dir,publication=False):
    fig,axes = plt.subplots(2,2,figsize=(8.0,7.6))
    if publication:
        fig.subplots_adjust(left=.105,right=.97,bottom=.06,top=.82,hspace=.43,wspace=.29)
    else:
        fig.subplots_adjust(left=.105,right=.97,bottom=.15,top=.73,hspace=.49,wspace=.29)
        fig.text(.105,.975,"Optimize in one geometry; measure in all four",fontsize=14,fontweight="bold",va="top")
        fig.text(.105,.929,"Each panel measures the same four sets of assignments using its own travel distances.",fontsize=9,color="#444444")
    handles = [Line2D([],[],color=COLORS[s],lw=2,label=LABELS[s]) for s in SOURCES]
    fig.legend(handles=handles,title="Assignments optimized using",loc="upper center",
               bbox_to_anchor=(.53,.992 if publication else .89),ncol=4,
               frameon=False,fontsize=8.5,title_fontsize=9,
               columnspacing=1.8,handletextpad=.55,borderaxespad=0)
    style_handles = [
        Line2D([],[],color="#333333",lw=1.8,ls="-",label="Native front"),
        Line2D([],[],color="#333333",lw=1.8,ls=(0,(4,2)),
               label="Transferred assignments"),
        Line2D([],[],color="#222222",marker="X",markersize=7,ls="None",
               label="Natural lineage"),
    ]
    fig.legend(handles=style_handles,loc="upper center",
               bbox_to_anchor=(.53,.917 if publication else .80),ncol=3,
               frameon=False,fontsize=8,handlelength=2.2,
               columnspacing=2.1,handletextpad=.55,borderaxespad=0)
    all_x,all_y = frame.x.to_numpy(),frame.y.to_numpy()
    for letter,ax,target in zip("ABCD",axes.flat,SOURCES):
        inset = ax.inset_axes([.40,.43,.56,.51]) if insets else None
        # Draw the in-place front last; keep every source's color fixed.
        for source in [s for s in SOURCES if s != target]+[target]:
            d = frame[(frame.target==target)&(frame.source==source)].sort_values("step",ascending=False)
            for plot_ax in [ax]+([inset] if inset is not None else []):
                plot_ax.plot(d.x,d.y,color=COLORS[source],lw=1.9 if source==target else 1.15,
                        ls="-" if source==target else (0,(4,2)),alpha=1 if source==target else .85)
        nx,ny = references[target]["natural_display"]
        ax.scatter([nx],[ny],marker="X",s=55,color="#222222",edgecolor="white",lw=.6,zorder=8)
        if inset is not None:
            inset.scatter([nx],[ny],marker="X",s=36,color="#222222",edgecolor="white",lw=.5,zorder=8)
            inset.set_xlim(nx+INSET_OFFSETS["travel"][0],
                           nx+INSET_OFFSETS["travel"][1])
            inset.set_ylim(ny+INSET_OFFSETS["state"][0],
                           ny+INSET_OFFSETS["state"][1])
            inset.set_title("Near natural lineage",fontsize=7.5,pad=5)
            inset.tick_params(labelsize=6.5,pad=2,length=2)
            inset.xaxis.set_major_locator(plt.MaxNLocator(3))
            inset.yaxis.set_major_locator(plt.MaxNLocator(3))
            inset.grid(alpha=.15,lw=.4)
            for spine in inset.spines.values():
                spine.set_color("#888888")
                spine.set_linewidth(.6)
            ax.indicate_inset_zoom(inset,edgecolor="#888888",alpha=.65,lw=.6)
        ax.text(-.17,1.09,letter,transform=ax.transAxes,fontsize=12,fontweight="bold")
        title = f"Travel measured in embryo {target[-1]}" if target != "pooled" else "Travel measured in pooled geometry"
        ax.set_title(title,loc="left",fontsize=9.5,fontweight="semibold",pad=11)
        ax.set_xlabel("Travel cost (endpoint-normalized)",fontsize=8.5)
        ax.set_ylabel("Cell-state cost\n(endpoint-normalized)",fontsize=8.5)
        ax.tick_params(labelsize=8)
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(alpha=.18,lw=.5)
        if zoom:
            # Fixed declared window around nature, shared as relative offsets;
            # each panel retains its own reference endpoint transform.
            ax.set_xlim(nx+INSET_OFFSETS["travel"][0],
                        nx+INSET_OFFSETS["travel"][1])
            ax.set_ylim(ny+INSET_OFFSETS["state"][0],
                        ny+INSET_OFFSETS["state"][1])
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))
            ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        else:
            ax.set_xlim(min(-.015,float(all_x.min())-.015),max(1.025,float(all_x.max())+.015))
            ax.set_ylim(min(-.015,float(all_y.min())-.015),max(1.025,float(all_y.max())+.015))
    if not publication:
        scope = "Full fronts; insets zoom near natural lineage" if insets else ("Zoom around natural lineage" if zoom else "Full sampled fronts")
        fig.text(.105,.094,scope+". All curves within a panel share that geometry's reference endpoints.",fontsize=8)
        fig.text(.105,.066,"Pooled travel = mean of the three travel costs, each divided by its embryo's natural-lineage total.",fontsize=8)
        fig.text(.105,.038,"275 shared edges; fixed protein costs; 301 weighted solutions per source. Lines connect sampled assignments.",fontsize=8)
    stem = "figS3_four_geometries_with_insets" if insets else ("figS3_four_geometries_near_natural" if zoom else "figS3_four_geometries_full")
    for ext in ("pdf","png","svg"):
        fig.savefig(out_dir/f"{stem}.{ext}",dpi=220,bbox_inches="tight",facecolor="white")
    plt.close(fig)


# Typography the S3 script applies before rendering (kept for parity).
S3_RC = {"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"}


def render_publication(frame, references, *, out_dir):
    with style.matplotlib_defaults(**S3_RC):
        render(frame, references, zoom=False, insets=True, out_dir=out_dir, publication=True)
