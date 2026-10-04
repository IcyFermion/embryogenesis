"""Publication tables. Table 1: per-subtree canonical and null-relative statistics.

Formatting moved from ``write_stats_table`` in
``terminal_pareto/fig5_table1_ce_canonical_metrics.py``; that module is
hash-pinned by the tracking-metric audit cache and stays byte-identical, so
its writer remains for the legacy full pipeline. This version formats the
validated canonical-metric CSV and takes every symbol from the notation
registry.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from publication import notation as nt

TABLE1_STEM = "table1_ce_subtree_statistics"
REGION_LABELS = {"ABa": "ABa", "ABp": "ABp", "P1": "P1", "root": "Root"}
REGION_ORDER = {"ABa": 0, "ABp": 1, "P1": 2, "root": 3}


def _rows(metrics: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    df = metrics[metrics["endpoint_ok"].astype(bool)].copy()
    df = df.merge(summary[["subtree", "entropy_nats"]], on="subtree", how="left")
    df["_r"] = df["region"].map(REGION_ORDER)
    return df.sort_values(["_r", "n"], ascending=[True, False])


def subtree_statistics_tex(metrics: pd.DataFrame, summary: pd.DataFrame, min_cells: int) -> str:
    """Standalone booktabs Table 1 from validated canonical metrics and the subtree summary."""
    u, lp, cp = (nt.symbol(k) for k in ("canonical_position", "natural_front_distance", "null_front_distance"))
    null_mean, closest, r = nt.symbol("null_mean"), nt.symbol("closest_point"), nt.symbol("cousin_relative")
    body, last_region = [], None
    for _, row in _rows(metrics, summary).iterrows():
        if last_region is not None and row["region"] != last_region:
            body.append(r"\midrule")
        last_region = row["region"]
        on_front = "yes" if bool(row["exactly_on_front"]) else "no"
        body.append(
            f"{row['subtree']} & {REGION_LABELS[row['region']]} & "
            f"{int(row['n'])} & {on_front} & "
            f"{row['u_lineage_lp']:.3f} & {row['d_lp']:.3f} & "
            f"{row['d_np']:.3f} & {row['relative_distance']:.3f} & "
            f"{row['u_max_er']:.3f} & {row['max_er']:.3f} & "
            f"{row['entropy_nats']:.3f} \\\\")
    rows_tex = "\n".join(body)
    return rf"""\documentclass[10pt]{{article}}
\usepackage[margin=14mm]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{caption}}
\usepackage{{graphicx}}
\pagestyle{{empty}}
\captionsetup{{font=small,labelfont=bf}}

\begin{{document}}
\begin{{table}}[p]
\centering
\caption{{\textbf{{Per-subtree statistics for the investigated lineages.}}
%
For every subtree containing at least {min_cells} usable terminal cells: the
number \(n\) of usable terminal cells within the subtree; whether the natural terminal
assignment is exactly Pareto-optimal at some objective weighting (on
front); \({u}_L\), the position along the subtree's own Pareto front
(normalized arc length, travel optimum 0 to cell-state optimum 1) of the
closest sampled Pareto assignment in endpoint-normalized geometry;
\({lp}\), the endpoint-normalized distance from the natural terminal
assignment to the closest sampled Pareto assignment along the subtree's own
Pareto front;
\({cp}\), the endpoint-normalized distance from the first-cousin-null
mean \({null_mean}\) to that same sampled Pareto assignment \({closest}\);
\({r}\), the relative distance under the original first-cousin-null-SD
geometry (null reference distance \({r}=1\); exactly-on-front subtrees have
\({r}=0\)); because \({r}\) uses a different coordinate system, it is not the
ratio \({lp}/{cp}\); \({u}_{{ER}}\), the normalized front position of the
maximum-edge-retention assignment (choosing the tied vertex nearest \({u}_L\)
when multiple assignments attain the maximum); maximum edge retention, the
largest fraction of natural terminal parent--child edges preserved anywhere
along the front; and fate diversity, the Shannon entropy
\(H=-\sum_k p_k\ln p_k\) for the terminal cell-type proportions in nats
(0 for a single fate; larger values indicate a more even mixture of fates).
Subtrees are grouped by first-order lineage region. The endpoint-normalized
coordinates are those used in the canonical subtree summary (Figure 5);
\({r}\) is retained as a null-model sensitivity measure.}}
\label{{tab:ce_subtree_statistics}}
\resizebox{{\textwidth}}{{!}}{{
\begin{{tabular}}{{llrrrrrrrrr}}
\toprule
Subtree & Region & \(n\) & On front & \({u}_L\) & \({lp}\) & \({cp}\) &
\({r}\) & \({u}_{{ER}}\) & Max ER & Fate diversity (nats) \\
{rows_tex}
\bottomrule
\end{{tabular}}
}}
\end{{table}}
\end{{document}}
"""


def write_subtree_statistics(out_dir: Path, metrics, summary, min_cells: int) -> Path:
    path = Path(out_dir) / f"{TABLE1_STEM}.tex"
    path.write_text(subtree_statistics_tex(metrics, summary, min_cells))
    return path
