"""Matched terminal-only comparison pilot; no accepted publication outputs touched.

Run from the repository root in the dev environment. Raw 3D uses the existing
loader convention; C. briggsae axial calibration is unresolved. XY is a separate
sensitivity, not a corrected 3D reconstruction. All outputs go to output/pilot_cross_dataset.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist, cosine

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from terminal_pareto import data_loader as dl
from terminal_pareto import pareto_engine as pe
from terminal_pareto import plot_style as ps
from terminal_pareto.fig5_table1_ce_canonical_metrics import compute_front_landmarks
from terminal_pareto.subtree_analysis import exact_cousin_stats, build_type_map

OUT = ROOT / "terminal_pareto/output/pilot_cross_dataset"
STYLES = {
    "ce_protein": ("C. elegans protein", "#0072B2", "-"),
    "ce_rna": ("C. elegans RNA", "#E69F00", "-"),
    "cb_af16_p1": ("C. briggsae RNA · AF16 p1", "#009E73", "-"),
    "cb_af16_p5": ("C. briggsae RNA · AF16 p5", "#006B50", "--"),
}
GEOMETRIES = {"raw3d": "Existing 3D convention — C. briggsae axial scale unverified",
              "xy": "XY-only sensitivity — axial coordinates omitted"}


def load_inputs(out):
    lineage = dl.load_json(str(ROOT / "data/cell_lineage.json"))
    prot, ce, cb = dl.load_protein_expression(), dl.load_ce_rna(), dl.load_cb_rna()
    psels, rsels = dl.load_prot_sel(), dl.load_rna_sel()
    xyz, valid = dl.load_elegans_tracking(255)
    configs = {
        "ce_protein": dict(xyz=xyz, valid=valid, exp=prot, features=psels,
                           path=dl.CE_TRACKS_PATH, cutoff=255),
        "ce_rna": dict(xyz=xyz, valid=valid, exp=ce, features=rsels,
                       path=dl.CE_TRACKS_PATH, cutoff=255),
    }
    metadata = pd.read_csv(ROOT / "data/c_briggsae/yiming/metadata.csv").set_index("ID")
    for key, rep in [("cb_af16_p1", "210519ZZY0874p1"),
                     ("cb_af16_p5", "210519ZZY0874p5")]:
        _, path, cutoff = next(row for row in dl.CB_REPLICATES_NEW if row[0].startswith(rep))
        assert metadata.loc[rep, "Type"] == "AF16"
        xyz, valid = dl.load_briggsae_tracking(cutoff, path)
        configs[key] = dict(xyz=xyz, valid=valid, exp=cb, features=rsels,
                            path=path, cutoff=cutoff)
    for cfg in configs.values():
        exp = cfg["exp"]
        assert exp.index.is_unique and exp.columns.is_unique
        assert len(cfg["features"]) == len(set(cfg["features"])) == 20
        values = exp[cfg["features"]].to_numpy(float)
        assert np.isfinite(values).all() and (np.linalg.norm(values, axis=1) > 0).all()
        tn, tp = dl.collect_terminals(lineage, set(cfg["valid"]) & set(exp.index))
        cfg["terms"] = list(zip(tn, tp))
        assert len(tn) == len(set(tn))
    shared = set.intersection(*(set(c["terms"]) for c in configs.values()))
    common = [e for e in configs["ce_protein"]["terms"] if e in shared]
    assert len(common) == 188, f"Inputs changed: expected 188 shared edges, got {len(common)}"
    assert set(configs["ce_rna"]["terms"]) <= set(configs["ce_protein"]["terms"])
    types = build_type_map([c for c, _ in set.union(*(set(x["terms"]) for x in configs.values()))])
    annotation = pd.read_csv(ROOT / "data/c_briggsae/science.adu8249/shared/CellTable_20241027.txt",
                             sep="\t", keep_default_na=False).set_index("Lineage")
    manifest = []
    for config, cfg in configs.items():
        for child, parent in cfg["terms"]:
            region = next((r for r in ("ABa", "ABp", "MS", "E", "C", "D", "P4")
                           if child.startswith(r)), "other")
            manifest.append(dict(config=config, terminal=child, natural_parent=parent,
                matched=(child, parent) in shared, cell_type=types[child], region=region,
                annotation_parent=annotation.loc[child, "Parent"] if child in annotation.index else "absent",
                child_annotation_missing=annotation.loc[child, "Missing"] if child in annotation.index else "absent",
                parent_annotation_missing=annotation.loc[parent, "Missing"] if parent in annotation.index else "absent"))
    frame = pd.DataFrame(manifest)
    frame.to_csv(out / "cell_manifest.csv", index=False)
    for group in ("cell_type", "region"):
        summary = frame.groupby(["config", group], sort=False).agg(
            available=("matched", "size"), matched=("matched", "sum")).reset_index()
        summary["retained_fraction"] = summary.matched / summary.available
        summary.to_csv(out / f"coverage_by_{group}.csv", index=False)
    mismatches = frame.query("annotation_parent != 'absent' and annotation_parent != natural_parent")
    assert mismatches.empty, mismatches.to_string()
    audit = []
    for key, cfg in configs.items():
        nodes = sorted({x for edge in common for x in edge})
        positions = np.asarray([cfg["xyz"][x] for x in nodes], float)
        audit.append(dict(config=key, source=str(Path(cfg["path"]).relative_to(ROOT)),
            cutoff=cfg["cutoff"], available_edges=len(cfg["terms"]), matched_edges=len(common),
            retained_fraction=len(common)/len(cfg["terms"]),
            matched_z_extent=float(np.ptp(positions[:, 0])),
            matched_horizontal_extent_1=float(np.ptp(positions[:, 1])),
            matched_horizontal_extent_2=float(np.ptp(positions[:, 2])),
            preprocessing="stored z-scored proteins" if key == "ce_protein" else "stored RNA values; no additional transform",
            cutoff_basis="predefined cutoff; developmental alignment not independently revalidated"))
    pd.DataFrame(audit).to_csv(out / "input_audit.csv", index=False)
    return configs, common, pe.build_grandparent_map(lineage)


def costs(cfg, terms, geometry):
    tn, tp = zip(*terms)
    p = np.asarray([cfg["xyz"][x] for x in tp], float)
    c = np.asarray([cfg["xyz"][x] for x in tn], float)
    if geometry == "xy":
        p, c = p[:, 1:], c[:, 1:]
    ep = cfg["exp"].loc[list(tp), cfg["features"]].to_numpy(float)
    ec = cfg["exp"].loc[list(tn), cfg["features"]].to_numpy(float)
    xm, em = cdist(p, c), cdist(ep, ec, metric="cosine")
    assert np.isfinite(xm).all() and np.isfinite(em).all()
    # Check the vectorized matrix against the established scalar convention.
    for i, j in [(0, 0), (0, len(tn)-1), (len(tn)-1, 0)]:
        assert np.isclose(xm[i, j], np.linalg.norm(p[i]-c[j]))
        assert np.isclose(em[i, j], cosine(ep[i], ec[j]))
    return xm, em


def endpoint_assignment(primary, secondary):
    """Resolve endpoint ties by the other objective, with numerical verification.

    Test diminishing perturbations in matrices scaled to unit maximum. Accept
    only assignments retaining the primary minimum to floating-point summation
    tolerance, with the same secondary minimum at the next smaller perturbation.
    This is a checked numerical tie-break, not an exhaustive tie enumeration.
    """
    n = len(primary)
    rows, ordinary = linear_sum_assignment(primary)
    p = primary / max(float(np.max(np.abs(primary))), 1e-15)
    s = secondary / max(float(np.max(np.abs(secondary))), 1e-15)
    best = np.sum(p[rows, ordinary], dtype=np.longdouble)
    tol = 32*np.finfo(float).eps*n
    previous = None
    for power in range(4, 14):
        eps = 10.0**(-power)
        _, chosen = linear_sum_assignment(p + eps*s)
        delta = float(np.sum(p[rows, chosen], dtype=np.longdouble)-best)
        value = float(s[rows, chosen].sum())
        if abs(delta) <= tol:
            if previous is not None and abs(value-previous[0]) <= tol:
                return chosen, dict(primary_gap=float(primary[rows, chosen].sum()-primary[rows, ordinary].sum()),
                    secondary_saving=float(secondary[rows, ordinary].sum()-secondary[rows, chosen].sum()),
                    perturbation=eps, primary_tolerance_scaled=tol)
            previous = (value, chosen)
        else:
            previous = None
    raise AssertionError("Endpoint tie-break did not stabilize while preserving the primary optimum")


def analyze(cfg, terms, geometry, gp, intervals, draws):
    tn, tp = map(np.asarray, zip(*terms))
    xm, em = costs(cfg, terms, geometry)
    groups = pe.build_cousin_groups(tn, gp)
    mx, me, vx, ve, cov = exact_cousin_stats(xm, em, groups)
    assert vx > 0 and ve > 0
    sx, se = np.sqrt(vx), np.sqrt(ve)
    lx, le = np.trace(xm), np.trace(em)
    alpha = np.linspace(0, 1, intervals+1)
    perms, xs, es, ers = [], [], [], []
    rows = np.arange(len(tn))
    state_endpoint, state_tie = endpoint_assignment(em, xm)
    travel_endpoint, travel_tie = endpoint_assignment(xm, em)
    for a in alpha:
        ri, ci = linear_sum_assignment(a*xm/sx + (1-a)*em/se)
        if a == 0:
            ci = state_endpoint
        elif a == 1:
            ci = travel_endpoint
        assert np.array_equal(ri, rows) and np.array_equal(np.sort(ci), rows)
        perms.append(ci)
        xs.append(xm[ri, ci].sum())
        es.append(em[ri, ci].sum())
        ers.append(np.mean(tp[ri] == tp[ci]))
    x, e, er = map(np.asarray, (xs, es, ers))
    assert np.all(np.diff(x) <= 1e-7) and np.all(np.diff(e) >= -1e-9)
    assert np.isclose(x[-1], xm[linear_sum_assignment(xm)].sum())
    assert np.isclose(e[0], em[linear_sum_assignment(em)].sum())
    assert pe.lineage_edge_ratio(rows, rows, tp) == 1
    landmarks = compute_front_landmarks((x-lx)/sx, (e-le)/se, er,
                                        (mx-lx)/sx, (me-le)/se)
    assert landmarks["endpoint_ok"] and landmarks["u_monotone"]
    dx, de = x[0]-x[-1], e[-1]-e[0]
    X, Y = (x-x[-1])/dx, (e-e[0])/de
    L = np.array([(lx-x[-1])/dx, (le-e[0])/de])
    pidx = int(np.argmin(np.hypot(X-L[0], Y-L[1])))
    assert np.isclose(np.hypot(X[pidx]-L[0], Y[pidx]-L[1]), landmarks["d_lp"])
    # Polyline u only parameterizes actual solutions; no interpolated solution is scored.
    seg = np.hypot(np.diff(X[::-1]), np.diff(Y[::-1]))
    u = (np.r_[0, np.cumsum(seg)] / seg.sum())[::-1]
    assert np.isclose(u[pidx], landmarks["u_lineage_lp"])
    max_indices = np.flatnonzero(er == er.max())
    midx = int(max_indices[np.argmin(np.abs(u[max_indices]-u[pidx]))])
    rx, re = pe.compute_group_shuffle_costs(xm, em, groups, n_random=draws, seed=42)
    curve = pd.DataFrame(dict(sweep_index=np.arange(intervals+1), alpha_travel=alpha,
        travel=x, cell_state=e, travel_sigma=(x-lx)/sx, cell_state_sigma=(e-le)/se,
        D1=X, D2=Y, u=u, edge_retention=er,
        nearest_assignment=np.arange(intervals+1)==pidx,
        maximum_retention=np.arange(intervals+1)==midx))
    null = pd.DataFrame(dict(travel=rx, cell_state=re,
        travel_sigma=(rx-lx)/sx, cell_state_sigma=(re-le)/se))
    metrics = dict(n=len(tn), intervals=intervals, distinct_cost_points=landmarks["n_front_points"],
        distinct_biological_assignments=len({tuple(tp[np.argsort(p)]) for p in perms}),
        d_LP=landmarks["d_lp"], d_NP=landmarks["d_np"], u_L=landmarks["u_lineage_lp"],
        max_edge_retention=float(er.max()), retention_at_nearest=float(er[pidx]),
        retention_travel=float(er[-1]), retention_state=float(er[0]),
        natural_D1=L[0], natural_D2=L[1], natural_travel=lx, natural_cell_state=le,
        null_mean_travel=mx, null_mean_state=me, null_variance_travel=vx,
        null_variance_state=ve, null_covariance=cov, travel_sigma=sx, state_sigma=se,
        travel_span=dx, state_span=de,
        endpoint_ties_checked=True,
        **{f"state_endpoint_{k}": v for k, v in state_tie.items()},
        **{f"travel_endpoint_{k}": v for k, v in travel_tie.items()})
    return dict(curve=curve, null=null, metrics=metrics, permutations=np.asarray(perms),
                terminals=tn, parents=tp, xm=xm, em=em)


def draw_comparison(results, geometry):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.0, 6.6))
    fig.subplots_adjust(left=.07, right=.975, top=.815, bottom=.32, wspace=.25)
    fig.suptitle("Terminal parentage across protein and RNA representations", fontsize=14, y=.98)
    fig.text(.5, .932, f"PILOT · 188 matched terminal edges · {GEOMETRIES[geometry]}",
             ha="center", fontsize=9, color="#555555")
    handles, data = [], []
    points = []
    for key, (label, color, ls) in STYLES.items():
        r = results[("matched", geometry, key)]
        f, m = r["curve"], r["metrics"]
        near = f.loc[f.nearest_assignment].iloc[0]
        maximum = f.loc[f.maximum_retention].iloc[0]
        ax.plot(f.D1, f.D2, color=color, ls=ls, lw=1.7)
        ax.scatter(m["natural_D1"], m["natural_D2"], marker="X", s=66,
                   color=color, edgecolors="white", linewidths=.6, zorder=6)
        ax.plot([m["natural_D1"], near.D1], [m["natural_D2"], near.D2], color=color, lw=1.1)
        ax.scatter(near.D1, near.D2, s=28, facecolors="white", edgecolors=color, zorder=5)
        bx.plot(f.u, f.edge_retention, color=color, ls=ls, lw=1.6)
        bx.scatter(near.u, near.edge_retention, s=35, facecolors="white", edgecolors=color, zorder=5)
        bx.scatter(maximum.u, maximum.edge_retention, s=32, marker="D", color=color,
                   edgecolors="white", linewidths=.4, zorder=5)
        handles.append(Line2D([], [], color=color, ls=ls, label=label))
        data.append([label, f'{m["d_LP"]:.3f}', f'{m["d_NP"]:.3f}',
                     f'{m["retention_at_nearest"]:.3f}', f'{m["max_edge_retention"]:.3f}'])
        points.extend([(m["natural_D1"], m["natural_D2"]), (near.D1, near.D2)])
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.52, .915),
               ncol=2, fontsize=8, columnspacing=2)
    ax.set(xlabel="Travel distance (endpoint normalized)",
           ylabel="Cell-state distance (endpoint normalized)", title="A   Pareto geometry")
    ax.set_xlim(-.025, 1.04)
    ax.set_ylim(-.025, 1.04)
    bx.set(xlabel="Canonical front position, u", ylabel="Natural-edge retention",
           title="B   Parentage retained along the sampled front", xlim=(-.025, 1.025), ylim=(0, 1))
    for a in (ax, bx):
        a.grid(alpha=.3)
    # A zoom makes the natural-lineage distances legible without stretching the main axes.
    z = ax.inset_axes([.40, .42, .56, .50])
    pts = np.asarray(points)
    low, high = pts.min(axis=0), pts.max(axis=0)
    pad = np.maximum((high-low)*.25, .012)
    for key, (_, color, ls) in STYLES.items():
        r = results[("matched", geometry, key)]
        f, m = r["curve"], r["metrics"]
        n = f.loc[f.nearest_assignment].iloc[0]
        z.plot(f.D1, f.D2, color=color, ls=ls, lw=1.3)
        z.plot([m["natural_D1"], n.D1], [m["natural_D2"], n.D2], color=color)
        z.scatter(m["natural_D1"], m["natural_D2"], marker="X", s=52,
                  color=color, edgecolors="white", linewidths=.5, zorder=5)
        z.scatter(n.D1, n.D2, s=24, facecolors="white", edgecolors=color, zorder=4)
    z.set(xlim=(low[0]-pad[0], high[0]+pad[0]), ylim=(low[1]-pad[1], high[1]+pad[1]))
    z.set_title("Near the natural lineages", fontsize=7)
    z.tick_params(labelsize=6)
    z.grid(alpha=.25)
    marker_handles = [Line2D([], [], marker="X", color="#333333", ls="", label="Natural lineage"),
        Line2D([], [], marker="o", markerfacecolor="white", color="#333333", ls="", label="Nearest sampled assignment"),
        Line2D([], [], marker="D", color="#333333", ls="", label="Maximum observed retention")]
    fig.legend(handles=marker_handles, loc="lower center", bbox_to_anchor=(.5, .238), ncol=3, fontsize=8)
    tax = fig.add_axes([.14, .063, .74, .16])
    tax.axis("off")
    table = tax.table(cellText=data,
        colLabels=["Configuration", r"$d_{LP}$", r"$d_{NP}$", "Retention at nearest", "Maximum retention"],
        colWidths=[.36, .12, .12, .21, .19], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(.3)
        cell.set_edgecolor("#dddddd")
        if row == 0:
            cell.set_facecolor("#f3f3f3")
        if row and col == 0:
            cell.get_text().set_color(list(STYLES.values())[row-1][1])
    fig.text(.5, .025, "Endpoint spans are configuration-specific. Lines connect sampled assignments; intervening lineages are not evaluated.",
             ha="center", fontsize=7, color="#555555")
    return fig


def draw_nulls(results):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.subplots_adjust(top=.77, bottom=.15, wspace=.30)
    fig.suptitle("First-cousin-null scaling · matched terminal pilot", fontsize=13)
    for ax, geometry in zip(axes, GEOMETRIES):
        for key, (label, color, ls) in STYLES.items():
            r = results[("matched", geometry, key)]
            f, n, m = r["curve"], r["null"], r["metrics"]
            ax.scatter(n.travel_sigma, n.cell_state_sigma, s=3, color=color, alpha=.06, rasterized=True)
            ax.plot(f.travel_sigma, f.cell_state_sigma, color=color, ls=ls, label=label)
            ax.scatter((m["null_mean_travel"]-m["natural_travel"])/m["travel_sigma"],
                       (m["null_mean_state"]-m["natural_cell_state"])/m["state_sigma"],
                       marker="+", color=color, s=42, zorder=5)
        ax.scatter(0, 0, marker="X", color="black", s=65, zorder=6)
        ax.axhline(0, color="gray", ls=":", lw=.6)
        ax.axvline(0, color="gray", ls=":", lw=.6)
        ax.set(xlabel="Travel distance (null SD; natural lineage = 0)",
               ylabel="Cell-state distance (null SD; natural lineage = 0)",
               title="Existing 3D (axial scale unverified)" if geometry == "raw3d" else "XY only")
        ax.grid(alpha=.3)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(.5, .92), ncol=2)
    fig.text(.5, .02, "Colored clouds: 1,000 first-cousin shuffles per curve. Plus signs: exact null means. Black cross: natural lineage.",
             ha="center", fontsize=7)
    return fig


def draw_coverage(out, results):
    fig = plt.figure(figsize=(11, 6.6))
    fig.suptitle("Pilot checks · coverage and dependence on the terminal set", fontsize=13, y=.96)
    ax = fig.add_axes([.17, .26, .26, .57])
    coverage = pd.read_csv(out / "coverage_by_cell_type.csv")
    df = coverage.query("config == 'ce_protein'").sort_values("available")
    ax.barh(df.cell_type.str.replace("_", " "), df.available, color="#dddddd", label="Available in protein analysis")
    ax.barh(df.cell_type.str.replace("_", " "), df.matched, color="#0072B2", label="In the common 188-edge set")
    for y, (_, row) in enumerate(df.iterrows()):
        ax.text(row.available+1, y, f"{row.matched}/{row.available}", va="center", fontsize=7)
    ax.set(xlim=(0, df.available.max()*1.25), xlabel="Terminal edges",
           title="A   Which cells survive matching?")
    ax.legend(loc="upper left", bbox_to_anchor=(-.4, -.15), fontsize=7)
    tax = fig.add_axes([.49, .28, .49, .54])
    tax.axis("off")
    rows = []
    for cohort in ("available", "ce231", "matched"):
        for key, (label, _, _) in STYLES.items():
            if (cohort, "raw3d", key) not in results:
                continue
            a, b = (results[(cohort, g, key)]["metrics"] for g in GEOMETRIES)
            short = {"ce_protein": "CE protein", "ce_rna": "CE RNA",
                     "cb_af16_p1": "CB AF16 p1", "cb_af16_p5": "CB AF16 p5"}[key]
            rows.append([short, str(a["n"]), f'{a["max_edge_retention"]:.3f}',
                         f'{b["max_edge_retention"]:.3f}', f'{a["d_LP"]:.3f}', f'{b["d_LP"]:.3f}'])
    table = tax.table(cellText=rows,
        colLabels=["Configuration", "n", "Max ER\n3D*", "Max ER\nXY", r"$d_{LP}$"+"\n3D*", r"$d_{LP}$"+"\nXY"],
        colWidths=[.28, .10, .17, .15, .15, .15], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(.4)
        cell.set_edgecolor("#dddddd")
        cell.set_facecolor("#eeeeee" if row == 0 else ("#f2f7fa" if row >= 7 else "white"))
    tax.set_title("B   Full available, CE-only matched, and common sets", fontsize=9, pad=12)
    fig.text(.50, .225, "Blue rows: the same 188 edges across all configurations.\n"
             "CE-only matched set: the same 231 edges across the two CE modalities.", fontsize=8)
    fig.text(.10, .08, "*3D uses the existing coordinate convention; C. briggsae axial calibration remains unresolved.\n"
             "Coverage is selective: only 1 of 32 programmed-death cells and 7 of 20 glial cells in the protein set remain.\n"
             "The RNA matrices are fixed across tracking embryos. These are descriptive sensitivities, not population uncertainty.", fontsize=8)
    return fig


def render(out, results):
    figs = [("comparison_existing_3d", draw_comparison(results, "raw3d")),
            ("comparison_xy", draw_comparison(results, "xy")),
            ("null_sd_companion", draw_nulls(results)),
            ("coverage_and_set_sensitivity", draw_coverage(out, results))]
    with PdfPages(out / "pilot_review.pdf") as review:
        for stem, fig in figs:
            fig.savefig(out / f"{stem}.pdf")
            fig.savefig(out / f"{stem}.png", dpi=180)
            review.savefig(fig)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--intervals", type=int, default=300)
    parser.add_argument("--dense-intervals", type=int, default=1200)
    parser.add_argument("--render-only", action="store_true", help="Render existing caches without rerunning optimization")
    args = parser.parse_args()
    assert args.intervals > 0 and args.dense_intervals > args.intervals
    assert args.dense_intervals % args.intervals == 0
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    ps.configure()
    if args.render_only:
        metrics = pd.read_csv(out / "metrics.csv")
        fronts = pd.read_csv(out / "fronts.csv")
        nulls = pd.read_csv(out / "null_clouds.csv")
        results = {}
        for _, m in metrics.query("cohort != 'matched_dense'").iterrows():
            key = (m.cohort, m.geometry, m.config)
            f = fronts[(fronts.cohort == m.cohort) & (fronts.geometry == m.geometry) & (fronts.config == m.config)]
            n = nulls[(nulls.cohort == m.cohort) & (nulls.geometry == m.geometry) & (nulls.config == m.config)]
            results[key] = dict(metrics=m.to_dict(), curve=f, null=n)
        render(out, results)
        return
    configs, common, gp = load_inputs(out)
    results, metrics, fronts, clouds, arrays, convergence = {}, [], [], [], {}, []
    tasks = [("matched", key, common) for key in configs]
    tasks += [("available", key, cfg["terms"]) for key, cfg in configs.items()]
    tasks += [("ce231", key, configs["ce_rna"]["terms"]) for key in ("ce_protein", "ce_rna")]
    for cohort, key, terms in tasks:
        for geometry in GEOMETRIES:
            r = analyze(configs[key], terms, geometry, gp, args.intervals, 1000)
            results[(cohort, geometry, key)] = r
            ident = dict(cohort=cohort, geometry=geometry, config=key)
            metrics.append(dict(**ident, **r["metrics"]))
            fronts.append(r["curve"].assign(**ident))
            clouds.append(r["null"].assign(**ident))
            prefix = f"{cohort}__{geometry}__{key}"
            for name in ("permutations", "terminals", "parents", "xm", "em"):
                arrays[f"{prefix}__{name}"] = r[name]
            print(prefix, {k: round(r["metrics"][k], 4) for k in ("d_LP", "d_NP", "max_edge_retention")}, flush=True)
            if cohort == "matched":
                dense = analyze(configs[key], terms, geometry, gp, args.dense_intervals, 1000)
                # Nested grids must contain the same objective costs at corresponding weights.
                stride = args.dense_intervals // args.intervals
                assert np.allclose(dense["curve"].travel.to_numpy()[::stride], r["curve"].travel)
                assert np.allclose(dense["curve"].cell_state.to_numpy()[::stride], r["curve"].cell_state)
                dense_ident = dict(cohort="matched_dense", geometry=geometry, config=key)
                metrics.append(dict(**dense_ident, **dense["metrics"]))
                fronts.append(dense["curve"].assign(**dense_ident))
                for name in ("permutations", "terminals", "parents"):
                    arrays[f"matched_dense__{geometry}__{key}__{name}"] = dense[name]
                convergence.append(dict(config=key, geometry=geometry,
                    base_intervals=args.intervals, dense_intervals=args.dense_intervals,
                    **{f"delta_{k}": dense["metrics"][k]-r["metrics"][k]
                       for k in ("d_LP", "u_L", "max_edge_retention", "d_NP")}))
    # Holding geometry fixed across modalities must hold the entire travel matrix fixed.
    for geometry in GEOMETRIES:
        assert np.array_equal(results[("matched", geometry, "ce_protein")]["xm"],
                              results[("matched", geometry, "ce_rna")]["xm"])
        assert np.array_equal(results[("matched", geometry, "cb_af16_p1")]["em"],
                              results[("matched", geometry, "cb_af16_p5")]["em"])
    for key in configs:
        a, b = (results[("matched", g, key)] for g in GEOMETRIES)
        assert np.array_equal(a["em"], b["em"])
        # State optimum value is geometry-invariant; secondary travel tie-break
        # may legitimately select different assignments within that optimum.
        assert np.isclose(a["metrics"]["natural_cell_state"], b["metrics"]["natural_cell_state"])
        assert np.isclose(a["curve"].cell_state.iloc[0], b["curve"].cell_state.iloc[0])
    pd.DataFrame(metrics).to_csv(out / "metrics.csv", index=False)
    pd.concat(fronts, ignore_index=True).to_csv(out / "fronts.csv", index=False)
    pd.concat(clouds, ignore_index=True).to_csv(out / "null_clouds.csv", index=False)
    pd.DataFrame(convergence).to_csv(out / "convergence.csv", index=False)
    np.savez_compressed(out / "assignments_and_costs.npz", **arrays)
    sources = [Path(__file__), ROOT / "terminal_pareto/data_loader.py",
        ROOT / "terminal_pareto/pareto_engine.py", ROOT / "terminal_pareto/subtree_analysis.py",
        ROOT / "terminal_pareto/fig5_table1_ce_canonical_metrics.py",
        ROOT / "data/cell_lineage.json", ROOT / "data/protein/aggregated_all/s3_zscore.csv",
        ROOT / "data/c_briggsae/science.adu8249/c_elegans_tf.csv",
        ROOT / "data/c_briggsae/science.adu8249/c_briggsae_tf.csv",
        ROOT / "expression_embedding/results/elegans_protein_linear_baseline/top20_protein_names.csv",
        ROOT / "expression_embedding/results/cross_species_rna_linear/rna_selected_features.tsv",
        ROOT / "data/c_briggsae/yiming/metadata.csv",
        ROOT / "data/c_briggsae/science.adu8249/shared/CellTable_20241027.txt"]
    sources += [Path(c["path"]) for c in configs.values()]
    provenance = dict(intervals=args.intervals, dense_intervals=args.dense_intervals,
        seed=42, null_draws=1000, matched_edges=len(common),
        python=sys.version, numpy=np.__version__,
        features={k: c["features"] for k, c in configs.items()},
        sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        limitations=["C. briggsae z-to-xy calibration unresolved; raw3d is a legacy convention",
            "Predefined cutoffs and RNA measurement/imputation provenance require publication audit",
            "Protein and RNA feature panels differ; this does not isolate modality",
            "Shared lineage annotation is not independent species-specific parentage validation",
            "RNA expression is fixed across tracking replicates; no RNA replicate uncertainty",
            "Sampled weighted optima only; interior ties are solver-selected, not exhaustively enumerated",
            "Endpoints minimize the other objective within numerical primary-optimum tolerance"])
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    render(out, results)
    print("Convergence (dense minus base):", flush=True)
    print(pd.DataFrame(convergence).to_string(index=False), flush=True)
    print(f"Pilot complete: {out / 'pilot_review.pdf'}", flush=True)


if __name__ == "__main__":
    main()
