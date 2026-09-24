"""Project four saved assignment sweeps into four tracking cost geometries.

The default command renders the pooled publication layout in its run directory.
No optimization is performed here.
"""
import argparse
from pathlib import Path
import hashlib
import json
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from terminal_pareto.analysis_context import build_analysis_context, validate_existing_context_manifest
from terminal_pareto.global_analysis import load_global_analysis
from terminal_pareto.front_coordinates import EndpointTransform

BASE = ROOT / "terminal_pareto/output"
TRACK = BASE / "tracking_geometry_sensitivity/assignments_and_costs.npz"
SOURCES = ["embryo1", "embryo2", "embryo3", "pooled"]
LABELS = {"embryo1":"Embryo 1", "embryo2":"Embryo 2", "embryo3":"Embryo 3", "pooled":"Pooled (all three)"}
COLORS = dict(zip(SOURCES,["#4477AA", "#CC6677", "#228833", "#202020"]))
INSET_OFFSETS = {"travel": [-0.026, 0.036], "state": [-0.09, 0.09]}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_projections(context, tracking_path=TRACK):
    """Validate source caches, then score all 4 x 4 x 301 assignments."""
    validate_existing_context_manifest(context)
    pooled = load_global_analysis(context)
    if context.profile != "pooled_tracking_v1" or context.spec.sweep_intervals != 300:
        raise ValueError("Cross-geometry Figure S3 requires the pooled 300-interval run")
    tracking_path = Path(tracking_path)
    historical = json.loads((tracking_path.parent / "provenance.json").read_text())
    if historical.get("n") != len(context.terminal_nodes) or historical.get("intervals") != 300:
        raise ValueError("Historical tracking assignment cache has different sampling settings")
    for source_path, recorded_hash in historical.get("sha256", {}).items():
        if source_path.endswith(("/data/embryo1/tracks.txt", "/data/embryo2/tracks.txt",
                                 "/data/embryo3/tracks.txt", "/data/cell_lineage.json")):
            if digest(source_path) != recorded_hash:
                raise ValueError(f"Historical tracking source changed: {source_path}")
    with np.load(tracking_path, allow_pickle=False) as cache:
        np.testing.assert_array_equal(cache["children"], context.terminal_nodes)
        np.testing.assert_array_equal(cache["natural_parents"], context.terminal_parents)
        np.testing.assert_allclose(cache["expression_matrix"],context.expression_matrix,atol=1e-12,rtol=0)
        alpha = cache["alpha"].copy()
        np.testing.assert_array_equal(alpha, np.linspace(0,1,context.spec.sweep_intervals+1))
        # Convert historical row->column permutations to child->parent-slot.
        historical_assignments = {
            r:cache[r+"_assignments"].copy() for r in SOURCES[:3]}
        assignments = {r:np.argsort(historical_assignments[r],axis=1)
                       for r in SOURCES[:3]}
        historical_matrices = {r:cache[r+"_travel_matrix"].copy() for r in SOURCES[:3]}
        historical_state = cache["expression_matrix"].copy()
        for r in SOURCES[:3]:
            np.testing.assert_allclose(historical_matrices[r],context.component_travel_matrices[r],atol=1e-12,rtol=0)
    assignments["pooled"] = pooled.assignments
    matrices = dict(context.component_travel_matrices, pooled=context.travel_matrix)
    children = np.arange(len(context.terminal_nodes))
    natural_parent_slots = np.asarray(context.terminal_parents)
    expected_parent_counts = dict(zip(*np.unique(natural_parent_slots, return_counts=True)))
    state_costs = {}
    for name, a in assignments.items():
        assert a.shape == (301,275)
        np.testing.assert_array_equal(np.sort(a,axis=1),np.broadcast_to(children,a.shape))
        for assigned in a:
            observed = dict(zip(*np.unique(natural_parent_slots[assigned], return_counts=True)))
            if observed != expected_parent_counts:
                raise ValueError(f"{name} assignment violates biological-parent capacity")
        state_costs[name] = context.expression_matrix[a,children].sum(axis=1)
        if name != "pooled":
            np.testing.assert_allclose(
                state_costs[name], historical_state[
                    children,historical_assignments[name]].sum(axis=1),
                atol=1e-10,rtol=0)
    np.testing.assert_allclose(state_costs["pooled"],pooled.raw_front_state,atol=1e-10,rtol=0)
    records, references = [], {}
    max_projection_error = 0.0
    for target, matrix in matrices.items():
        own_a = assignments[target]
        own_t = matrix[own_a,children].sum(axis=1)
        own_s = state_costs[target]
        rr,cc = linear_sum_assignment(matrix)
        np.testing.assert_allclose(own_t[-1],matrix[rr,cc].sum(),atol=1e-10,rtol=0)
        rr,cc = linear_sum_assignment(context.expression_matrix)
        np.testing.assert_allclose(own_s[0],context.expression_matrix[rr,cc].sum(),atol=1e-10,rtol=0)
        if target == "pooled":
            np.testing.assert_allclose(own_t,pooled.raw_front_travel,atol=1e-12,rtol=0)
        transform = EndpointTransform.from_endpoints(
            reference_analysis_id=f"{context.cache_key}:{target}",
            travel_optimum_assignment_id=f"{target}:sweep:300",
            state_optimum_assignment_id=f"{target}:sweep:0",
            travel_optimum_costs=(own_t[-1],own_s[-1]),
            state_optimum_costs=(own_t[0],own_s[0]))
        natural = (float(np.trace(matrix)),float(np.trace(context.expression_matrix)))
        nx,ny = transform.transform(*natural)
        anchors = transform.transform(np.array([own_t[-1],own_t[0]]),np.array([own_s[-1],own_s[0]]))
        np.testing.assert_allclose(anchors,[[0,1],[1,0]],atol=1e-12,rtol=0)
        references[target] = dict(transform.metadata(), natural_raw=list(natural),
                                   natural_display=[float(nx),float(ny)])
        for source in SOURCES:
            travel = matrix[assignments[source],children].sum(axis=1)
            state = state_costs[source]
            if target != "pooled":
                saved_costs = historical_matrices[target][
                    assignments[source],children].sum(axis=1)
                max_projection_error = max(
                    max_projection_error, float(np.max(np.abs(travel-saved_costs))))
                np.testing.assert_allclose(travel,saved_costs,atol=1e-10,rtol=0)
            else:
                saved_costs = np.mean([
                    historical_matrices[r][assignments[source],children].sum(axis=1)
                    / context.component_natural_totals[r]
                    for r in SOURCES[:3]],axis=0)
                max_projection_error = max(
                    max_projection_error, float(np.max(np.abs(travel-saved_costs))))
                np.testing.assert_allclose(travel,saved_costs,atol=1e-10,rtol=0)
            x,y = transform.transform(travel,state)
            tx,sy = transform.inverse(x,y)
            np.testing.assert_allclose(tx,travel,atol=1e-10,rtol=0)
            np.testing.assert_allclose(sy,state,atol=1e-10,rtol=0)
            records.extend(dict(target=target,source=source,step=k,alpha=alpha[k],
                travel=travel[k],state=state[k],x=x[k],y=y[k],in_place=target==source)
                for k in range(len(alpha)))
    frame = pd.DataFrame(records)
    # Every source assignment is evaluated unchanged in each target geometry.
    assert len(frame) == 4*4*301
    assert (frame.groupby(["source","step"]).state.nunique() == 1).all()
    for source in SOURCES:
        means = sum(frame[(frame.source==source)&(frame.target==r)].travel.to_numpy()
                    / context.component_natural_totals[r] for r in SOURCES[:3])/3
        actual = frame[(frame.source==source)&(frame.target=="pooled")].travel.to_numpy()
        np.testing.assert_allclose(means,actual,atol=1e-12,rtol=0)
    return frame,references,assignments,max_projection_error


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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="migration_candidate_20260920")
    parser.add_argument("--output-root", type=Path, default=BASE)
    parser.add_argument("--tracking-cache", type=Path, default=TRACK)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--publication", action="store_true", default=True,
                        help="Use manuscript-size art and keep the full/local views optional.")
    parser.add_argument("--diagnostics", action="store_true",
                        help="Also export separate full and near-natural views.")
    args = parser.parse_args(argv)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"pdf.fonttype":42,"svg.fonttype":"none"})
    context = build_analysis_context(
        "pooled_tracking_v1", run_id=args.run_id, output_root=args.output_root)
    if args.out is None:
        args.out = context.run_paths.display("endpoint")
    args.out.mkdir(parents=True,exist_ok=True)
    frame,references,assignments,max_projection_error=load_projections(
        context, args.tracking_cache)
    if args.diagnostics or not args.publication:
        render(frame,references,zoom=True,out_dir=args.out,
               publication=args.publication)
        render(frame,references,zoom=False,out_dir=args.out,
               publication=args.publication)
    render(frame,references,zoom=False,insets=True,out_dir=args.out,
           publication=args.publication)
    coordinates = args.out / "cross_geometry_coordinates.csv"
    frame.to_csv(coordinates,index=False)
    inputs=[args.tracking_cache,args.tracking_cache.parent/"provenance.json",
            context.run_paths.analysis/"analysis_manifest.json",
            context.run_paths.analysis/"global_terminal_analysis.npz",
            context.run_paths.analysis/"global_terminal_analysis.json",
            Path(__file__),ROOT/"terminal_pareto/front_coordinates.py"]
    inset_limits = {
        target: {
            "travel": [natural[0]+offset for offset in INSET_OFFSETS["travel"]],
            "state": [natural[1]+offset for offset in INSET_OFFSETS["state"]],
        }
        for target,ref in references.items()
        for natural in [ref["natural_display"]]
    }
    provenance = {
        "status":"pooled publication build",
        "profile": context.profile,
        "run_id": context.run_paths.run_id,
        "context_cache_key":context.cache_key,
        "projection_count": len(frame),
        "max_saved_cost_projection_error": max_projection_error,
        "projection_table_sha256": digest(coordinates),
        "references":references,
        "source_assignments_sha256":{s:hashlib.sha256(a.tobytes()).hexdigest() for s,a in assignments.items()},
        "input_sha256":{str(p.resolve()):digest(p) for p in inputs},
        "checks":["source-derived cohort and matrix equivalence","historical tracking source hashes and 300-interval sweep",
                  "4,816 target/source/sweep costs verified against saved matrices",
                  "301 feasible slot permutations and biological-parent capacities per source",
                  "independent endpoint optimum checks","pooled production cache cost reconstruction",
                  "state costs invariant across evaluation geometries","pooled costs equal mean component ratios",
                  "shared target-native endpoints and inverse transformations"],
        "display":{
            "asset":"figS3_four_geometries_with_insets.pdf",
            "main_axes":"full sampled fronts and all transferred points",
            "inset_offsets_from_natural":INSET_OFFSETS,
            "inset_limits_by_target":inset_limits,
            "reference_policy":"each target native front anchors all four sources and natural lineage",
            "line_policy":"sweep order retained; transferred points not re-filtered",
            "panel_order":SOURCES,
        },
    }
    (args.out/"provenance.json").write_text(json.dumps(
        provenance,indent=2,sort_keys=True)+"\n")
    if not args.publication:
        subprocess.run(["pdfunite",str(args.out/"figS3_four_geometries_near_natural.pdf"),
                        str(args.out/"figS3_four_geometries_full.pdf"),
                        str(args.out/"four_geometry_review.pdf")],check=True)
    print(args.out/"figS3_four_geometries_with_insets.pdf")


if __name__=="__main__":
    main()
