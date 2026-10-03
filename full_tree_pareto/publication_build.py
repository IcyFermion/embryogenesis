"""Build Figure 7 and its supplement from isolated, validated pooled caches.

python -m full_tree_pareto.publication_build [--layout-only] [--publish]
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np

from full_tree_pareto import pooled_analysis as pa
from full_tree_pareto import cousin_references as cr
from terminal_pareto import plot_style as ps
from terminal_pareto.front_coordinates import EndpointTransform

# Single source for the sweep/draw settings used by builds, checkpoints and
# captions. pooled_analysis.py is hash-pinned in cache identities, so these live
# here rather than there; changing them requires a new run directory.
INTERVALS = 300
DRAWS = 10000
SEED = 42
DEFAULT_WORKERS = 24
DISPLAY_DRAWS = 1000
PANEL_COLUMNS = 4
NUMBER_WORDS = dict(enumerate(("Zero", "One", "Two", "Three", "Four", "Five", "Six",
                               "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve")))
WRAPPERS = ("fig7_ce_full_tree_layerwise", "figs_ce_full_tree_heuristics")
METHOD_COLORS = dict(zip(pa.METHODS, ("#0072B2", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")))
NULL_COLORS = dict(zip(pa.NULLS, (ps.SEMANTIC_COLORS["first_cousin_null"], "#8195A5", "#9E6BA0", "#6F3B5C", "#6B6ECF")))
NULL_COLORS.update(zip(cr.REFERENCES, (ps.NULL_MODEL_COLORS["second_cousin"],
                                      ps.NULL_MODEL_COLORS["third_cousin"])))
MAIN_REFERENCES = (pa.NULLS[0], *cr.REFERENCES)
SUPPLEMENT_REFERENCES = (*MAIN_REFERENCES, *pa.NULLS[1:])
X_LABEL = r"Pooled travel distance, $D_1$"
Y_LABEL = r"Cell-state distance, $D_2$"
RETENTION_CMAP = LinearSegmentedColormap.from_list(
    "full_tree_retention", ["#17365D", ps.COLORS["blue"], "#72C7EC"])


def endpoint(frame, scope):
    a = frame.loc[frame.weight_index.idxmax()]
    b = frame.loc[frame.weight_index.idxmin()]
    return EndpointTransform.from_endpoints(
        reference_analysis_id=f"{pa.PROFILE}:layerwise:{scope}",
        travel_optimum_assignment_id=f"layerwise:{scope}:{int(a.weight_index)}",
        state_optimum_assignment_id=f"layerwise:{scope}:{int(b.weight_index)}",
        travel_optimum_costs=(a.travel, a.state),
        state_optimum_costs=(b.travel, b.state))


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def retention_front(ax, frame, transform, natural, *, endpoints=False):
    x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
    order = np.argsort(x)
    ax.plot(x[order], y[order], color=METHOD_COLORS[pa.METHODS[0]], lw=.7, alpha=.65, zorder=3)
    scatter = ax.scatter(x, y, c=frame.retention, cmap=RETENTION_CMAP, norm=Normalize(0, 1), s=9, linewidths=0, zorder=4)
    nx_, ny_ = transform.transform(*np.asarray(natural))
    ax.scatter(nx_, ny_, marker="X", s=42, color=ps.COLORS["black"], edgecolors="white", lw=.5, zorder=7)
    best = np.flatnonzero(frame.retention.to_numpy() == frame.retention.max())
    i = best[len(best)//2]
    ax.scatter(x[i], y[i], marker="o", s=28, c=[frame.retention.iloc[i]], cmap=RETENTION_CMAP, norm=Normalize(0, 1), edgecolors=ps.COLORS["black"], lw=.8, zorder=6)
    if endpoints:
        ax.annotate("A", (0, 1), xytext=(-8, 5), textcoords="offset points", fontsize=7)
        ax.annotate("B", (1, 0), xytext=(5, -3), textcoords="offset points", fontsize=7)
    ax.grid(True, alpha=.3)
    ax.margins(.13)
    return scatter


def plot_rounds(ctx, fronts, layers, out):
    rows = math.ceil(len(ctx.layers) / PANEL_COLUMNS)
    fig, axes = plt.subplots(rows, PANEL_COLUMNS, figsize=(7.15, 1.925 * rows),
                             layout="constrained", squeeze=False)
    for ax in axes.flat[len(ctx.layers):]:
        ax.remove()
    axes = axes.flat[:len(ctx.layers)]
    transforms = {}
    for i, ax in enumerate(axes, 1):
        scope = f"round_{i}"
        frame = fronts[(fronts.method == pa.METHODS[0]) & (fronts.scope == scope)]
        transform = endpoint(frame, scope)
        row = layers.set_index("scope").loc[scope]
        scatter = retention_front(ax, frame, transform, [row.natural_travel, row.natural_state])
        ax.set_title(f"Round {i}" + (" (bottom)" if i == 1 else " (top)" if i == len(ctx.layers) else "") + f"\n{int(row.edges)} edges", fontsize=8)
        ax.tick_params(labelsize=6.5)
        transforms[scope] = transform.metadata(axis_limits={"x": ax.get_xlim(), "y": ax.get_ylim()})
    cb = fig.colorbar(scatter, ax=list(axes), fraction=.025, pad=.015, aspect=32)
    cb.set_label("Natural edges retained", fontsize=8)
    fig.supxlabel(X_LABEL, fontsize=8.5)
    fig.supylabel(Y_LABEL, fontsize=8.5)
    fig.text(.003, 1.015, "A", fontweight="bold", fontsize=11)
    save(fig, out, "fig7A_ce_full_tree_layerwise_rounds")
    return transforms


def plot_comparison(ctx, fronts, nulls, transform, out, *, supplement=False):
    fig, ax = plt.subplots(figsize=(7.15, 4.0 if not supplement else 5.2))
    fig.subplots_adjust(left=.10, right=.97, bottom=.14, top=.97)
    models = pa.METHODS if supplement else pa.METHODS[:1]
    references = SUPPLEMENT_REFERENCES if supplement else MAIN_REFERENCES
    handles = []
    for method in models:
        frame = fronts[(fronts.method == method) & (fronts.scope == "aggregate") & fronts.nondominated]
        x, y = transform.transform(frame.travel.to_numpy(), frame.state.to_numpy())
        order = np.argsort(x)
        line, = ax.plot(x[order], y[order], color=METHOD_COLORS[method],
                        lw=1.8 if method == pa.METHODS[0] else 1.2,
                        ls="-" if method in pa.METHODS[:2] else "--", zorder=4, label=method)
        handles.append(line)
    for name in references:
        values = nulls[name]
        subset = np.linspace(0, len(values)-1, min(DISPLAY_DRAWS, len(values)), dtype=int)
        x, y = transform.transform(values[subset, 0], values[subset, 1])
        ax.scatter(x, y, s=7, alpha=.13, color=NULL_COLORS[name], edgecolors="none", rasterized=True, zorder=1)
        mean = values.mean(axis=0)
        mx, my = transform.transform(*mean)
        marker = ax.scatter(mx, my, marker="+", s=55, lw=1.3, color=NULL_COLORS[name], label=name, zorder=5)
        handles.append(marker)
    x, y = transform.transform(*ctx.natural)
    natural = ax.scatter(x, y, marker="X", s=65, c=ps.COLORS["black"], edgecolors="white", lw=.55, zorder=7, label="Natural lineage")
    handles.append(natural)
    inset_limits = None
    if not supplement:
        # Match the terminal presentation: broad randomization is contextual,
        # with the exact same endpoint transform but separately scaled axes.
        name = "Random rebuild"
        values = nulls[name]
        subset = np.linspace(0, len(values)-1, min(DISPLAY_DRAWS, len(values)), dtype=int)
        x, y = transform.transform(values[subset, 0], values[subset, 1])
        inset = ax.inset_axes([.66, .70, .31, .25])
        inset.scatter(x, y, s=6, alpha=.18, color=NULL_COLORS[name],
                      edgecolors="none", rasterized=True)
        inset.scatter(*transform.transform(*values.mean(axis=0)), marker="+", s=35,
                      lw=1.1, color=NULL_COLORS[name])
        inset.set_title(name, fontsize=7, pad=2)
        inset.set_xlabel("Pooled travel distance", fontsize=6, labelpad=1)
        inset.set_ylabel("Cell-state distance", fontsize=6, labelpad=1)
        inset.tick_params(labelsize=5.5, length=2, pad=1)
        inset.xaxis.set_major_locator(plt.MaxNLocator(4))
        inset.yaxis.set_major_locator(plt.MaxNLocator(3))
        inset.grid(True, alpha=.25)
        inset.margins(.10)
        inset_limits = {"x": inset.get_xlim(), "y": inset.get_ylim()}
        handles.append(Line2D([0], [0], marker="+", ls="", color=NULL_COLORS[name],
                              label="Random rebuild (inset)"))
    ax.grid(True, alpha=.25)
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.margins(.07)
    legend_options = {} if supplement else dict(bbox_to_anchor=(.66, .59), borderaxespad=0)
    ax.legend(handles=handles, loc="upper left", fontsize=6.6, frameon=True,
              facecolor="white", edgecolor="none", framealpha=.9, **legend_options)
    if not supplement:
        ax.text(-.10, 1.015, "B", transform=ax.transAxes, fontweight="bold", fontsize=11)
    metadata = transform.metadata(axis_limits={"x": ax.get_xlim(), "y": ax.get_ylim()})
    metadata["reference_models"] = list(references)
    if inset_limits:
        metadata["random_rebuild_inset"] = dict(axis_limits=inset_limits,
                                                shares_endpoint_transform=True)
    save(fig, out, "figs_ce_full_tree_heuristics_panel" if supplement else "fig7B_ce_full_tree_collective")
    return metadata


PREAMBLE = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=12mm]{geometry}
\usepackage{graphicx,caption,booktabs}
\pagestyle{empty}
\captionsetup{font=footnotesize,labelfont=bf}
\begin{document}
"""


def write_wrappers(ctx, out):
    n, leaves, internal = len(ctx.edges), len(ctx.leaves), len(ctx.internal)
    rounds = ", ".join(str(len(r)) for r in ctx.layers)
    main = PREAMBLE + r"""\setcounter{figure}{6}
\begin{figure}[p]
\centering
\includegraphics[width=0.95\textwidth]{fig7A_ce_full_tree_layerwise_rounds.pdf}
\vspace{1mm}
\includegraphics[width=0.95\textwidth]{fig7B_ce_full_tree_collective.pdf}
\caption{\textbf{Full-tree trade-offs under layerwise assignment.}
The shared \emph{C. elegans} protein/tracking cohort contains INTERNAL internal
cells, LEAVES measured leaves and four roots, with EDGES scored edges.
Travel is the equal mean of pairwise displacements in three embryos (cutoffs
255, 247 and 225), each divided by its natural full-cohort travel total;
coordinates are not averaged. Cell-state distance is Euclidean distance
across the selected 20 z-scored protein features.
(A) NROUNDS asynchronous bottom-up contraction rounds (ROUNDS edges)
each reassign children to their natural parent slots using exact linear
assignment across NWEIGHTS shared weights after global first-cousin-null SD scaling.
Color indicates biological-parent
retention; outlined circles mark maximum-retention sampled assignments.
Rounds are not uniform developmental depths.
(B) Costs summed across rounds at each weight form the aggregate layerwise
front, compared with first-, second- and third-cousin shuffles. These permute
measured leaf identities within groups sharing an ancestor two, three or four
canonical generations back, respectively, including cutoff/coverage leaves;
all internal states, parent-slot multiplicities and four root identities stay
fixed. Position and protein state move together. All EDGES edges are scored,
including unchanged internal edges. Random rebuild appears in an inset, using
the same endpoint transform but separately scaled axes, as in the terminal plot.
Dots show NSHOWN of NDRAWS draws; plus signs show full-draw means. Black crosses
mark natural lineage. Coordinates use endpoint normalization: the travel
optimum is \((0,1)\) and the state optimum \((1,0)\), separately per round in A
and using aggregate layerwise endpoints in B. Observations are not clipped.
These are descriptive randomization references, not calibrated significance
tests; separation does not establish biological efficiency.}
\label{fig:ce_full_tree_layerwise}
\end{figure}
\end{document}
"""
    counts = dict(NROUNDS=NUMBER_WORDS.get(len(ctx.layers), str(len(ctx.layers))),
                  NWEIGHTS=f"{INTERVALS+1:,}", NSHOWN=f"{min(DISPLAY_DRAWS, DRAWS):,}",
                  NDRAWS=f"{DRAWS:,}")
    for key, value in (("INTERNAL", internal), ("LEAVES", leaves), ("EDGES", n),
                       ("NROUNDS", counts["NROUNDS"]), ("ROUNDS", rounds),
                       ("NWEIGHTS", counts["NWEIGHTS"]), ("NSHOWN", counts["NSHOWN"]),
                       ("NDRAWS", counts["NDRAWS"])):
        main = main.replace(key, str(value))
    supplement = PREAMBLE + r"""\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{figure}{3}
\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figs_ce_full_tree_heuristics_panel.pdf}
\caption{\textbf{Full-tree reconstruction strategies and reference models.}
Five strategies are regenerated on the same pooled-travel cohort and protein
representation as Figure 7, using NWEIGHTS endpoint-inclusive weights. Each curve
shows its own sampled non-dominated cost pairs, not a cross-method global
optimum. All curves, natural lineage and reference draws use the \emph{same}
aggregate layerwise endpoints from Figure 7B. These common display anchors do
not make the methods' feasible sets equivalent. All EDGES edges and all
measured nodes are scored once in each reconstruction; stored parent arrays
are replay-validated for binary degree, four components and absence of cycles.
The three cousin shuffles permute measured leaves within canonical ancestor
groups two, three or four generations back, with internal states and roots fixed;
internal-layer shuffle permutes internal identities within canonical depths
and holds leaves fixed; full assignment shuffle permutes internal and terminal
identities separately. The latter two permit root identities to move.
Random rebuild fixes the four roots and samples a binary topology over measured
states. The Gaussian reference, retained only here, fits separate-clock
covariance to all observed edges, freely simulates internal states from four
observed roots, then clamps every measured leaf to its observed position and
protein state. Independent spatial processes are fitted per embryo and their
costs pooled; protein variance is per canonical transition. Terminal edges
reconnect free internal states to observed leaves. This is not joint endpoint
conditioning or a calibrated significance test; Gaussian internal states are
not restricted to the measured-state inventory used by reconstruction methods.
Dots show deterministic NSHOWN-draw subsets; plus signs use all NDRAWS draws.
Axes are continuous and observations are not clipped.}
\label{fig:ce_full_tree_heuristics_supplement}
\vspace{4mm}
\footnotesize
\begin{tabular}{@{}p{0.30\textwidth}p{0.66\textwidth}@{}}
\toprule
Strategy & Constraints and algorithm \\
\midrule
Layerwise assignment & Exact assignment within fixed contraction rounds;
natural parent-slot multiplicities and four root identities preserved. \\
Degree-constrained forest & Greedy Kruskal internal forest with one fixed root
per component and binary capacity; terminal Hungarian assignment. \\
Top-down rebuild & Greedy attachment of internal states from four fixed roots;
terminal Hungarian assignment to remaining binary openings. \\
Bottom-up by layer & Fixed topology; internal identities reassigned by
depth from the bottom up, including movable root identities. Leaves fixed. \\
Paired bottom-up & Minimum-weight daughter pairing followed by parent
assignment; four root identities emerge, rather than being fixed. \\
Terminal-only (excluded) & Historical constructed-state algorithm has unresolved
scope/scoring problems; no pooled curve generated or displayed. \\
\bottomrule
\end{tabular}
\end{figure}
\end{document}
"""
    for key, value in (("EDGES", n), ("NWEIGHTS", counts["NWEIGHTS"]),
                       ("NSHOWN", counts["NSHOWN"]), ("NDRAWS", counts["NDRAWS"])):
        supplement = supplement.replace(key, str(value))
    for stem, source in zip(WRAPPERS, (main, supplement)):
        (out / f"{stem}.tex").write_text(source)


def compile_wrappers(out):
    compiler = shutil.which("tectonic")
    if not compiler:
        raise RuntimeError("Activate dev: tectonic is required")
    for stem in WRAPPERS:
        result = subprocess.run([compiler, "--keep-logs", f"{stem}.tex"], cwd=out, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (out / f"{stem}.build.log").write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(result.stdout)
        info = subprocess.check_output(["pdfinfo", str(out / f"{stem}.pdf")], text=True)
        pages = int(next(line.split(":")[1] for line in info.splitlines() if line.startswith("Pages:")))
        if pages != 1:
            raise ValueError(f"Expected one-page figure: {stem} has {pages}")


def archive_working_publication(run):
    """Retain a hash-verified layout snapshot before overwriting review assets."""
    run = Path(run)
    source = run / "publication"
    if not source.exists() or not any(source.iterdir()):
        return None
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = run / "layout_history" / stamp
    hashes = {str(p.relative_to(source)): pa.digest(p)
              for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source, archive / "publication")
    for name, sha in hashes.items():
        if pa.digest(archive / "publication" / name) != sha:
            raise ValueError("Working layout archive verification failed")
    pa.write_json(archive / "manifest.json", dict(files=hashes,
                  analysis_files={str(p.relative_to(run / "analysis")): pa.digest(p)
                                  for p in (run / "analysis").rglob("*") if p.is_file()}))
    print(f"Preserved preceding layout: {archive}", flush=True)
    return archive


def render(ctx, fronts, layers, nulls, run):
    ps.configure()
    archive_working_publication(run)
    out = Path(run) / "publication"
    out.mkdir(parents=True, exist_ok=True)
    transforms = plot_rounds(ctx, fronts, layers, out)
    frame = fronts[(fronts.method == pa.METHODS[0]) & (fronts.scope == "aggregate")]
    transform = endpoint(frame, "aggregate")
    transforms["main_collective"] = plot_comparison(ctx, fronts, nulls, transform, out)
    transforms["supplement"] = plot_comparison(ctx, fronts, nulls, transform, out, supplement=True)
    pa.write_json(out / "display_transforms.json", transforms)
    write_wrappers(ctx, out)
    compile_wrappers(out)
    return out


def verify_release(publication):
    publication = Path(publication)
    record = json.loads((publication / "release_manifest.json").read_text())
    for filename, sha in record["files"].items():
        if pa.digest(publication / filename) != sha:
            raise ValueError(f"Changed publication file: {filename}")
    run = Path(record["run"])
    for filename, sha in record["analysis_files"].items():
        if pa.digest(run / "analysis" / filename) != sha:
            raise ValueError(f"Changed analysis file: {filename}")
    return record


def sweep_settings():
    return dict(intervals=INTERVALS, draws=DRAWS, seed=SEED)


def promote(run):
    """Explicit publication replacement; archive/hash all former assets first."""
    run = Path(run).resolve()
    source = run / "publication"
    root = pa.ROOT / "full_tree_pareto/output"
    target = root / "publication"
    if (target / "cross_species_release_manifest.json").exists():
        raise ValueError("Production includes Figures 10/11. Coordinate an additive Figure 7 release "
                         "that preserves them and updates their preserved-file hashes.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    # Revalidate source/config and replay parent arrays before promotion.
    ctx, _, _, _ = pa.build(run=run, layout_only=True, **sweep_settings())
    cr.build(ctx, run, draws=DRAWS, seed=SEED, layout_only=True)
    for stem in WRAPPERS:
        info = subprocess.check_output(["pdfinfo", str(source / f"{stem}.pdf")], text=True)
        if int(next(x.split(":")[1] for x in info.splitlines() if x.startswith("Pages:"))) != 1:
            raise ValueError("Invalid wrapper page count")
    # Stage and verify the complete replacement before moving the live bundle.
    # Every failure path removes the stage and any partial archive.
    stage = Path(tempfile.mkdtemp(prefix=".pooled-stage-", dir=root))
    archive = root / "legacy/releases" / stamp
    previous = None
    created_archive = swapped = False
    try:
        for p in source.iterdir():
            if p.suffix in (".pdf", ".png", ".tex", ".json"):
                shutil.copy2(p, stage / p.name)
        record = dict(profile=pa.PROFILE, run=str(run), published_at=stamp,
                      files={p.name: pa.digest(p) for p in stage.iterdir() if p.is_file()},
                      analysis_files={str(p.relative_to(run / "analysis")): pa.digest(p)
                                      for p in (run / "analysis").rglob("*") if p.is_file()})
        pa.write_json(stage / "release_manifest.json", record)
        verify_release(stage)
        if target.exists():
            archive.mkdir(parents=True)
            created_archive = True
            hashes = {str(p.relative_to(target)): pa.digest(p) for p in sorted(target.rglob("*")) if p.is_file()}
            previous = archive / "publication"
            target.rename(previous)
            pa.write_json(archive / "manifest.json", hashes)
            for name, sha in hashes.items():
                if pa.digest(previous / name) != sha:
                    raise ValueError("Legacy archive verification failed")
        stage.rename(target)
        swapped = True
        verify_release(target)
    except Exception:
        if swapped:
            shutil.rmtree(target)
        if previous is not None and previous.exists():
            previous.rename(target)
        if created_archive:
            shutil.rmtree(archive, ignore_errors=True)
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return target


def positive_workers(value):
    """Validate the shared CLI worker option without affecting cache identity."""
    try:
        workers = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("workers must be a positive integer") from exc
    if workers < 1:
        raise argparse.ArgumentTypeError("workers must be a positive integer")
    return workers


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=pa.DEFAULT_RUN)
    parser.add_argument("--workers", type=positive_workers, default=DEFAULT_WORKERS,
                        help=f"Parallel solver processes (default: {DEFAULT_WORKERS}; use 1 for serial execution).")
    parser.add_argument("--layout-only", action="store_true")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--verify", action="store_true")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.verify:
        verify_release(pa.ROOT / "full_tree_pareto/output/publication")
        print("Verified pooled full-tree publication and analysis hashes")
        return
    result = pa.build(run=args.run, workers=args.workers, layout_only=args.layout_only,
                      **sweep_settings())
    result[3].update(cr.build(result[0], args.run, draws=DRAWS, seed=SEED,
                              layout_only=args.layout_only))
    out = render(*result, args.run)
    if args.publish:
        out = promote(args.run)
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main()
