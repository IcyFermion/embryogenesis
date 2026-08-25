"""Canonical coordinates and metrics for the subtree Pareto analysis.

This is the authoritative metric builder for Figure 5.  For each
investigated subtree it records (1) the compromise position ``u`` along its
own Pareto front and (2) the null-independent endpoint-normalized distance
``d_LP`` from the natural lineage to the nearest *attainable* Pareto
assignment.  The older first-cousin-null-relative distance ``r`` is retained
in the metrics table as a sensitivity measure and is shown only in the
supplementary summary rendered by ``fig5_figs1_ce_canonical_summary.py``.

Design rules:

- never combine travel and cell-state distance into a weighted biological
  cost; the axes are normalized separately and any Euclidean distance is
  used only in the explicitly defined display geometry;
- primary proximity is measured to the **closest sampled Pareto assignment**
  after separately normalizing travel and cell-state distances by the two
  single-objective endpoint spans. Interpolated segments are not attainable
  lineage assignments and never determine ``d_LP``;
- the established cousin-null ``r`` is computed independently in its
  original null-SD geometry and is never used to normalize ``d_LP``;
- nested subtrees are descriptive observations, not independent replicates;
- preserve maximum-edge-retention ties in diagnostics, but use one
  deterministic table representative: the actual tied max-ER vertex nearest
  the primary ``u_L``.

Usage::

    python terminal_pareto/fig5_table1_ce_canonical_metrics.py [--min-cells 12]
        [--iteration 300] [--selftest]

The metric table is written before plotting. All numerical helpers are pure
and independently testable (``--selftest`` runs the validation checklist).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from terminal_pareto import data_loader as dl
from terminal_pareto import lineage_metrics as lm
from terminal_pareto import pareto_engine as pe
from terminal_pareto import plot_style as ps
from terminal_pareto.subtree_analysis import first_cousin_null_summary

OUT = Path(__file__).resolve().parent / "output" / "publication"
ANALYSIS_OUT = Path(__file__).resolve().parent / "output"

# ── Analysis parameters (match the completed terminal pipeline) ──
ITERATION = 300          # Pareto sweep resolution (publication default)
MIN_CELLS = 12           # publication subtree threshold
U_GRID = np.linspace(0.0, 1.0, 201)   # common arc-length grid for the ensemble
D_LINEAGE_TOL = 1e-6     # d_lineage below this counts as "on the front"
D_NULL_MIN = 1e-10       # below this, r is flagged rather than divided


# ═══════════════════════════════════════════════════════════════
# Pure numerical helpers
# ═══════════════════════════════════════════════════════════════

def dedup_front(x, y, edge, rtol=1e-10, atol=1e-12):
    """Remove consecutive duplicate front coordinates.

    Consecutive alpha solutions frequently return the same assignment (and
    hence the same cost point). Within a run of identical coordinates the
    maximum edge retention is kept, so max-ER information survives dedup.
    Returns (x, y, edge) as float arrays.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    edge = np.asarray(edge, dtype=float)
    if len(x) < 2:
        return x, y, edge
    xs, ys, es = [x[0]], [y[0]], [edge[0]]
    for i in range(1, len(x)):
        if (np.isclose(x[i], xs[-1], rtol=rtol, atol=atol)
                and np.isclose(y[i], ys[-1], rtol=rtol, atol=atol)):
            es[-1] = max(es[-1], edge[i])
        else:
            xs.append(x[i])
            ys.append(y[i])
            es.append(edge[i])
    return np.array(xs), np.array(ys), np.array(es)


def endpoint_normalize(x, y, a_idx, b_idx, atol=1e-9):
    """Map the travel optimum A to (0, 1) and cell-state optimum B to (1, 0).

    ``x`` is the travel axis, ``y`` the cell-state axis. The front is assumed
    to run from A (index ``a_idx``) to B (index ``b_idx``):

        X = (x - x(A)) / (x(B) - x(A))
        Y = (y - y(B)) / (y(A) - y(B))

    Returns ``(X, Y, denom_x, denom_y)`` or ``None`` when either denominator
    is non-positive or numerically degenerate (caller reports and skips).
    """
    xa, ya = float(x[a_idx]), float(y[a_idx])
    xb, yb = float(x[b_idx]), float(y[b_idx])
    denom_x = xb - xa
    denom_y = ya - yb
    if not (denom_x > atol and denom_y > atol):
        return None
    X = (x - xa) / denom_x
    Y = (y - yb) / denom_y
    return X, Y, denom_x, denom_y


def arc_length_parameterize(x, y):
    """Cumulative arc length of the polyline, normalized to [0, 1].

    Returns ``(u, cum)`` where ``u[i]`` is the arc-length coordinate of
    vertex ``i`` (``u[0] = 0``, ``u[-1] = 1``). ``cum`` holds raw cumulative
    lengths. If the polyline has zero total length, ``u`` is all zeros.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    seg = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = cum[-1]
    if total <= 0:
        return np.zeros_like(cum), cum
    return cum / total, cum


def project_point_to_polyline(q, x, y, u):
    """Project ``q`` onto the polyline, including projection onto segments.

    Returns ``(u_star, p, distance)``: the arc-length coordinate of the
    closest point, the projected point, and the Euclidean distance. Falls
    back to the first vertex with infinite distance if the polyline has no
    segments.
    """
    q = np.asarray(q, dtype=float)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    best_dist, best_p, best_u = np.inf, None, np.nan
    for i in range(len(x) - 1):
        p1 = np.array([x[i], y[i]])
        p2 = np.array([x[i + 1], y[i + 1]])
        d = p2 - p1
        l2 = float(d @ d)
        if l2 <= 0:
            continue
        t = np.clip((q - p1) @ d / l2, 0.0, 1.0)
        p = p1 + t * d
        dist = float(np.linalg.norm(q - p))
        if dist < best_dist:
            best_dist, best_p, best_u = dist, p, u[i] + t * (u[i + 1] - u[i])
    if best_p is None:
        return 0.0, np.array([x[0], y[0]]), np.inf
    return best_u, best_p, best_dist


def max_er_spans(edge, u, rtol=1e-6, atol=1e-8):
    """All front positions attaining the maximum edge retention.

    Returns ``(m, mask, spans)``: the max value, a boolean mask over the
    (deduped) front, and the contiguous arc-length spans as a list of
    ``(u_start, u_end)`` tuples. The spans may be disjoint; every span is
    preserved so ties remain visible (no single representative midpoint).
    """
    edge = np.asarray(edge, dtype=float)
    u = np.asarray(u, dtype=float)
    m = float(np.max(edge))
    mask = np.isclose(edge, m, rtol=rtol, atol=atol)
    idx = np.flatnonzero(mask)
    spans = []
    if len(idx):
        start = prev = int(idx[0])
        for i in idx[1:]:
            i = int(i)
            if i > prev + 1:
                spans.append((float(u[start]), float(u[prev])))
                start = i
            prev = i
        spans.append((float(u[start]), float(u[prev])))
    return m, mask, spans


def resample_front(x, y, u, u_grid=U_GRID):
    """Piecewise-linear resampling of a normalized front onto ``u_grid``."""
    return np.interp(u_grid, u, x), np.interp(u_grid, u, y)


def ensemble_curves(resampled, u_grid=U_GRID, weights=None):
    """Ensemble median curve and a normal-deviation band.

    ``resampled`` is a list of ``(X, Y)`` arrays on ``u_grid``. The median
    is the coordinate-wise median. The band is the 25th--75th percentile of
    the per-front *signed deviation measured along the median curve's
    normal*, so the band curves delimit the middle 50% of subtree FRONTS --
    a valid two-dimensional envelope. (Marginal X/Y quartiles would be
    invalid: their paired coordinates need not belong to the same
    subtrees.) ``weights`` (optional, e.g. subtree sizes) selects a
    size-weighted median for comparison only -- used as a diagnostic for
    large-subtree dominance.
    """
    M = np.stack([f[0] for f in resampled])   # (n_fronts, n_grid)
    N = np.stack([f[1] for f in resampled])
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        w = w / w.sum()
        order = np.argsort(M, axis=0)
        cumw = np.cumsum(w[order], axis=0)
        med_x = np.array([M[order[:, j], j][np.searchsorted(cumw[:, j], 0.5)]
                          for j in range(M.shape[1])])
        order = np.argsort(N, axis=0)
        cumw = np.cumsum(w[order], axis=0)
        med_y = np.array([N[order[:, j], j][np.searchsorted(cumw[:, j], 0.5)]
                          for j in range(N.shape[1])])
    else:
        med_x = np.median(M, axis=0)
        med_y = np.median(N, axis=0)
    # Normal direction of the median curve at each u.
    dx = np.gradient(med_x, u_grid)
    dy = np.gradient(med_y, u_grid)
    nrm = np.hypot(dx, dy)
    nrm[nrm <= 0] = 1.0
    nx = -dy / nrm
    ny = dx / nrm
    # Signed deviation of each front along the normal.
    dev = ((M - med_x[None, :]) * nx[None, :]
           + (N - med_y[None, :]) * ny[None, :])
    lo = np.percentile(dev, 25, axis=0)
    hi = np.percentile(dev, 75, axis=0)
    return dict(
        med_X=med_x, med_Y=med_y,
        band_lo_X=med_x + lo * nx, band_lo_Y=med_y + lo * ny,
        band_hi_X=med_x + hi * nx, band_hi_Y=med_y + hi * ny,
    )


def medoid_front(resampled, u_grid=U_GRID):
    """Index of the observed resampled front minimizing summed curve distance."""
    M = np.stack([f[0] for f in resampled])
    N = np.stack([f[1] for f in resampled])
    n = len(resampled)
    dist = np.zeros((n, n))
    for i in range(n):
        dx = M - M[i]
        dy = N - N[i]
        dist[i] = np.sqrt(dx ** 2 + dy ** 2).sum(axis=1)
    return int(np.argmin(dist.sum(axis=1))), dist


def compute_front_landmarks(xa, ea, edge, nx, ny, atol=1e-9):
    """All per-subtree quantities for the canonical map.

    ``xa``, ``ea``: natural-lineage-centred, null-SD-scaled front arrays
    from ``pe.compute_std_scaled_pareto`` (cell-state optimum B at index 0,
    travel optimum A at the last index). ``nx``, ``ny``: first-cousin null
    mean in the same display coordinates.

    Proximity uses the **closest sampled Pareto assignment** in the raw
    null-SD coordinates (the established convention that defines
    ``rel_pareto_dist``); ``relative_distance`` below is that established
    ratio. Interpolated positions between two optimal assignments are NOT
    attainable lineages, so segment projection is kept only as a named
    diagnostic (``*_segment`` fields) and never drives the display.

    Returns a flat dict of scalars (plus the normalized front arrays
    ``X``, ``Y``, ``u``, ``edge``, ``max_er_mask``, ``max_er_spans``) with
    validation flags. Degenerate subtrees (zero or negative endpoint
    scales) are reported, not silently dropped.
    """
    x, y, e = dedup_front(xa, ea, edge)
    # Reverse so traversal runs A (travel optimum) -> B (cell-state optimum).
    x = x[::-1].copy()
    y = y[::-1].copy()
    e = e[::-1].copy()
    denx = float(x[-1] - x[0])     # x(B) - x(A)
    deny = float(y[0] - y[-1])     # y(A) - y(B)
    out = dict(
        n_front_raw=len(np.asarray(xa, dtype=float)),
        n_front_points=len(x),
        x_scale=denx, y_scale=deny,
        endpoint_ok=bool(denx > atol and deny > atol),
        max_er=float(e.max()) if len(e) else np.nan,
    )
    if not out["endpoint_ok"]:
        out.update(dict(
            X=np.array([np.nan]), Y=np.array([np.nan]), u=np.array([np.nan]),
            u_lineage=np.nan, u_max_er=np.nan, span_lo=np.nan, span_hi=np.nan,
            max_er_mask=np.zeros(len(e), dtype=bool), max_er_spans=[],
            d_lineage=np.nan, d_null=np.nan, relative_distance=np.nan,
            D1_lineage=np.nan, D2_lineage=np.nan,
            D1_null=np.nan, D2_null=np.nan,
            d_lp=np.nan, d_np=np.nan, u_lineage_lp=np.nan,
            P_lp_norm=(np.nan, np.nan), N_norm=(np.nan, np.nan),
            d_lineage_segment=np.nan, u_lineage_segment=np.nan,
            relative_distance_segment=np.nan,
            natural_on_grid=False, exactly_on_front=False,
            u_monotone=False, d_null_ok=False, M_raw=(np.nan, np.nan),
        ))
        return out

    X = (x - x[0]) / denx
    Y = (y - y[-1]) / deny
    u, _cum = arc_length_parameterize(X, Y)
    u_monotone = bool(len(u) > 1 and np.all(np.diff(u) > 0))

    # ── Proximity (established convention): closest SAMPLED front point in
    #    the raw null-SD coordinates. Interpolated points between two
    #    optimal assignments are not attainable lineages and are excluded.
    dists = np.sqrt(x ** 2 + y ** 2)
    p_idx = int(np.argmin(dists))
    d_lineage = float(dists[p_idx])
    d_null = float(np.hypot(nx - x[p_idx], ny - y[p_idx]))
    d_null_ok = bool(d_null > D_NULL_MIN)
    r = d_lineage / d_null if d_null_ok else np.nan
    u_L = float(u[p_idx])
    # The lineage's own normalized position (used only by the overlay).
    lx = (0.0 - x[0]) / denx
    ly = (0.0 - y[-1]) / deny

    # ── Null-independent endpoint-normalized proximity. In this geometry T
    #    is (0, 1), S is (1, 0), and the natural lineage is (lx, ly). The
    #    closest point must be an ACTUAL sampled Pareto assignment: line
    #    segments between assignments are drawing aids, not feasible lineage
    #    configurations. This is the primary Figure 5 distance metric.
    d_lp_vertices = np.hypot(X - lx, Y - ly)
    p_lp_idx = int(np.argmin(d_lp_vertices))
    d_lp = float(d_lp_vertices[p_lp_idx])
    u_lineage_lp = float(u[p_lp_idx])
    # Place the first-cousin-null mean in the same endpoint-normalized
    # geometry. d_NP uses the SAME sampled P* selected above for d_LP, making
    # the two absolute distances directly comparable. It is deliberately not
    # combined with the older cousin-SD ratio r, which uses another geometry.
    nlx = (nx - x[0]) / denx
    nly = (ny - y[-1]) / deny
    d_np = float(np.hypot(nlx - X[p_lp_idx], nly - Y[p_lp_idx]))

    # ── Segment-projection diagnostic (NOT used for the display). The
    #    polyline connects discrete optimal assignments; a small distance to
    #    an interpolated point does NOT mean the lineage is near an
    #    attainable assignment. Kept, quantified, and named separately.
    u_L_seg, P_seg, d_lineage_seg = project_point_to_polyline(
        np.array([lx, ly]), X, Y, u)
    d_null_seg = float(np.hypot(nlx - P_seg[0], nly - P_seg[1]))
    r_seg = (d_lineage_seg / d_null_seg
             if d_null_seg > D_NULL_MIN else np.nan)

    # ── Max-ER ties: all spans remain in diagnostics, but the publication
    #    representative is the ACTUAL tied vertex nearest the primary
    #    endpoint-normalized u_L. This retains
    #    maximum ER while minimizing the displayed movement from the natural
    #    lineage's closest-front compromise.
    m, mask, spans = max_er_spans(e, u)
    span_lo = min((s[0] for s in spans), default=np.nan)
    span_hi = max((s[1] for s in spans), default=np.nan)
    if spans:
        idx = np.flatnonzero(mask)
        max_idx = int(idx[np.argmin(np.abs(u[idx] - u_lineage_lp))])
        u_max_er = float(u[max_idx])
        m_raw = (float(x[max_idx]), float(y[max_idx]))
    else:
        u_max_er = np.nan
        m_raw = (np.nan, np.nan)
    # Is the natural lineage itself a sampled front point (grid convention)?
    natural_on_grid = bool(np.any(
        np.isclose(x, 0.0, rtol=1e-7, atol=1e-9)
        & np.isclose(y, 0.0, rtol=1e-7, atol=1e-9)))
    out.update(dict(
        X=X, Y=Y, u=u, edge=e,
        u_lineage=u_L, u_max_er=float(u_max_er),
        span_lo=float(span_lo), span_hi=float(span_hi),
        max_er=m, max_er_mask=mask, max_er_spans=spans,
        d_lineage=d_lineage, d_null=float(d_null),
        relative_distance=(float(r) if d_null_ok else np.nan),
        D1_lineage=float(lx), D2_lineage=float(ly),
        D1_null=float(nlx), D2_null=float(nly),
        d_lp=d_lp, d_np=d_np, u_lineage_lp=u_lineage_lp,
        P_lp_norm=(float(X[p_lp_idx]), float(Y[p_lp_idx])),
        N_norm=(float(nlx), float(nly)),
        d_lineage_segment=float(d_lineage_seg),
        u_lineage_segment=float(u_L_seg),
        relative_distance_segment=(float(r_seg)
                                   if np.isfinite(r_seg) else np.nan),
        natural_on_grid=natural_on_grid,
        exactly_on_front=bool(d_lineage <= D_LINEAGE_TOL),
        u_monotone=u_monotone, d_null_ok=d_null_ok,
        L_norm=(float(lx), float(ly)),
        P_norm=(float(X[p_idx]), float(Y[p_idx])),
        M_raw=m_raw,
    ))
    return out


# ═══════════════════════════════════════════════════════════════
# Analysis pipeline
# ═══════════════════════════════════════════════════════════════

def load_analysis():
    """Primary publication configuration (same inputs as subtree_analysis)."""
    lineage = dl.load_json(dl._REPO_ROOT + "/data/cell_lineage.json")
    tree_index = lm.build_lineage_tree_index(lineage)
    xyz_ce, valid_ce = dl.load_elegans_tracking(dl.T_CE)
    protein_exp = dl.load_protein_expression()
    prot_sel = dl.load_prot_sel()
    v_prot = [n for n in valid_ce if n in protein_exp.index]
    gp_map = pe.build_grandparent_map(lineage)
    return dict(lineage=lineage, tree_index=tree_index, xyz_ce=xyz_ce,
                protein_exp=protein_exp, prot_sel=prot_sel, v_prot=v_prot,
                gp_map=gp_map)


def analyze_subtree(name, terms, ctx, iteration=ITERATION, seed=42):
    """Front + null + landmarks for one subtree. Returns a record dict."""
    tn = [c for c, _p in terms]
    tp = [_p for c, _p in terms]
    xm, em, _ = pe.build_cost_matrices(tn, tp, ctx["xyz_ce"],
                                       ctx["protein_exp"], ctx["prot_sel"])
    null = first_cousin_null_summary(xm, em, tn, ctx["gp_map"], seed=seed)
    xa, ea, edge, _kp = pe.compute_std_scaled_pareto(
        xm, em, tp, null["_raw"], iteration=iteration)
    rec = compute_front_landmarks(
        xa, ea, edge,
        nx=null["null_xyz_mean_sigma"], ny=null["null_exp_mean_sigma"])
    rec.update(subtree=name, n=len(tn),
               region=assign_regions(name, ctx["tree_index"]))
    # Private keys for the selftest: raw front and null mean in display
    # coordinates (excluded from the published CSV by explicit column list).
    rec["_front"] = (np.asarray(xa, float), np.asarray(ea, float),
                      np.asarray(edge, float))
    rec["_null"] = (null["null_xyz_mean_sigma"],
                     null["null_exp_mean_sigma"])
    raw_null = null["_raw"]
    rec["_null_cloud"] = (
        (np.asarray(raw_null["random_xyz"], float)
         - raw_null["lineage_xyz"]) / raw_null["xyz_std"],
        (np.asarray(raw_null["random_exp"], float)
         - raw_null["lineage_exp"]) / raw_null["exp_std"],
    )
    return rec


def assign_regions(name, tree_index):
    """First-order lineage region (ABa, ABp, P1, or root) for a subtree."""
    curr = name
    while curr not in {"ABa", "ABp", "P1"}:
        p = tree_index[curr]["parent"]
        if p is None:
            break
        curr = p
    return curr if curr in {"ABa", "ABp", "P1"} else "root"


def analyze_all(min_cells=MIN_CELLS, iteration=ITERATION, ctx=None):
    """Landmarks for every subtree with >= min_cells usable terminals."""
    if ctx is None:
        ctx = load_analysis()
    subtrees = dl.collect_all_subtrees(ctx["lineage"], ctx["v_prot"],
                                       min_cells=min_cells)
    records = []
    skipped = []
    for name, terms in subtrees:
        rec = analyze_subtree(name, terms, ctx, iteration=iteration)
        if not rec["endpoint_ok"]:
            skipped.append((name, rec["x_scale"], rec["y_scale"],
                            rec["n_front_points"]))
        records.append(rec)
        print(f"  {name:10s} n={rec['n']:3d}  region={rec['region']:5s} "
              f"u_L={rec['u_lineage']:.3f}  r={rec['relative_distance']:.3f} "
              f"max_er={rec['max_er']:.3f}  on_front={int(rec['exactly_on_front'])}")
    print(f"\n{len(records)} subtrees; {len(skipped)} skipped (degenerate "
          f"endpoint scales): {skipped}")
    return records


def metrics_table(records):
    """Flat per-subtree table for the CSV (written before plotting).

    ``relative_distance`` is the established ``rel_pareto_dist`` convention
    (closest sampled Pareto assignment in null-SD coordinates, first-cousin
    null as the reference). The ``*_segment`` columns are diagnostics only:
    distances to interpolated polyline positions, which are not attainable
    lineage assignments. Every max-ER span is preserved in
    ``max_er_spans`` (semicolon-separated ``lo-hi`` pairs).
    """
    rows = []
    for r in records:
        rows.append(dict(
            subtree=r["subtree"], n=r["n"], region=r["region"],
            n_front_raw=r["n_front_raw"], n_front_points=r["n_front_points"],
            x_scale=r["x_scale"], y_scale=r["y_scale"],
            u_lineage=r["u_lineage"], u_max_er=r["u_max_er"],
            n_spans=len(r["max_er_spans"]),
            max_er_spans=";".join(f"{a:.6g}-{b:.6g}"
                                  for a, b in r["max_er_spans"]),
            span_lo=r["span_lo"], span_hi=r["span_hi"],
            max_er=r["max_er"], d_lineage=r["d_lineage"], d_null=r["d_null"],
            relative_distance=r["relative_distance"],
            D1_lineage=r["D1_lineage"], D2_lineage=r["D2_lineage"],
            D1_null=r["D1_null"], D2_null=r["D2_null"],
            d_lp=r["d_lp"], d_np=r["d_np"],
            u_lineage_lp=r["u_lineage_lp"],
            d_lineage_segment=r["d_lineage_segment"],
            u_lineage_segment=r["u_lineage_segment"],
            relative_distance_segment=r["relative_distance_segment"],
            exactly_on_front=r["exactly_on_front"],
            natural_on_grid=r["natural_on_grid"],
            endpoint_ok=r["endpoint_ok"], u_monotone=r["u_monotone"],
            d_null_ok=r["d_null_ok"],
        ))
    return pd.DataFrame(rows)


# Region labels used by the Table 1 supplement.
REGION_LABELS = {"ABa": "ABa", "ABp": "ABp",
                 "P1": "P1", "root": "Root"}


def write_stats_table(records, min_cells, out_dir=OUT, ctx=None):
    """LaTeX (booktabs) Table 1 for the canonical subtree analysis.

    Built from the canonical-map records (``records``) so the table carries
    the coordinates a reader needs to identify unlabeled map points: ``u_L``
    (position of the closest sampled Pareto assignment in endpoint-normalized
    geometry), null-independent distance ``d_LP``, endpoint-normalized
    first-cousin-null distance ``d_NP``, cousin-null-relative distance ``r``,
    ``u_ER`` (the actual tied maximum-retention vertex nearest ``u_L``), and
    maximum edge retention. Fate diversity is merged from the authoritative
    summary CSV. Rows are grouped by first-order lineage region; the tabular
    is resized to ``\\textwidth`` so it fits at final manuscript width.
    """
    rows = []
    for r in records:
        if not r["endpoint_ok"]:
            continue
        rows.append(dict(
            subtree=r["subtree"], n=r["n"], region=r["region"],
            on_front="yes" if r["exactly_on_front"] else "no",
            u_L=r["u_lineage_lp"],
            d_lp=r["d_lp"],
            d_np=r["d_np"],
            rel_dist=r["relative_distance"],
            u_ER=r["u_max_er"],
            max_er=r["max_er"],
        ))
    df = pd.DataFrame(rows)
    summary = (Path(__file__).resolve().parent / "output"
               / f"subtree_summary_min{min_cells}.csv")
    summ = pd.read_csv(summary)[["subtree", "entropy_nats"]]
    df = df.merge(summ, on="subtree", how="left")
    order = {"ABa": 0, "ABp": 1, "P1": 2, "root": 3}
    df["_r"] = df["region"].map(order)
    df = df.sort_values(["_r", "n"], ascending=[True, False])

    body = []
    last_reg = None
    for _, row in df.iterrows():
        if last_reg is not None and row["region"] != last_reg:
            body.append(r"\midrule")
        last_reg = row["region"]
        body.append(
            f"{row['subtree']} & {REGION_LABELS[row['region']]} & "
            f"{int(row['n'])} & {row['on_front']} & "
            f"{row['u_L']:.3f} & {row['d_lp']:.3f} & "
            f"{row['d_np']:.3f} & {row['rel_dist']:.3f} & "
            f"{row['u_ER']:.3f} & {row['max_er']:.3f} & "
            f"{row['entropy_nats']:.3f} \\\\")
    rows_tex = "\n".join(body)

    doc = rf"""\documentclass[10pt]{{article}}
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
front); \(u_L\), the position along the subtree's own Pareto front
(normalized arc length, travel optimum 0 to cell-state optimum 1) of the
closest sampled Pareto assignment in endpoint-normalized geometry;
\(d_{{LP}}\), the endpoint-normalized distance from the natural terminal
assignment to the closest sampled Pareto assignment along the subtree's own
Pareto front;
\(d_{{NP}}\), the endpoint-normalized distance from the first-cousin-null
mean \(N\) to that same sampled Pareto assignment \(P^*\);
\(r\), the relative distance under the original first-cousin-null-SD
geometry (null reference distance \(r=1\); exactly-on-front subtrees have
\(r=0\)); because \(r\) uses a different coordinate system, it is not the
ratio \(d_{{LP}}/d_{{NP}}\); \(u_{{ER}}\), the normalized front position of the
maximum-edge-retention assignment (choosing the tied vertex nearest \(u_L\)
when multiple assignments attain the maximum); maximum edge retention, the
largest fraction of natural terminal parent--child edges preserved anywhere
along the front; and fate diversity, the Shannon entropy
\(H=-\sum_k p_k\ln p_k\) for the terminal cell-type proportions in nats
(0 for a single fate; larger values indicate a more even mixture of fates).
Subtrees are grouped by first-order lineage region. The endpoint-normalized
coordinates are those used in the canonical subtree summary (Figure 5);
\(r\) is retained as a null-model sensitivity measure.}}
\label{{tab:ce_subtree_statistics}}
\resizebox{{\textwidth}}{{!}}{{
\begin{{tabular}}{{llrrrrrrrrr}}
\toprule
Subtree & Region & \(n\) & On front & \(u_L\) & \(d_{{LP}}\) & \(d_{{NP}}\) &
\(r\) & \(u_{{ER}}\) & Max ER & Fate diversity (nats) \\
{rows_tex}
\bottomrule
\end{{tabular}}
}}
\end{{table}}
\end{{document}}
"""
    path = out_dir / "table1_ce_subtree_statistics.tex"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(doc)
    print(f"Wrote {path.relative_to(Path.cwd())} ({len(df)} rows)")
    return path


def write_curves_npz(records, resampled, ens, medoid_idx,
                     out_dir=ANALYSIS_OUT):
    """Common u grid, normalized individual fronts, median, and band."""
    names = [r["subtree"] for r in records if r["endpoint_ok"]]
    fX = np.stack([f[0] for f in resampled])
    fY = np.stack([f[1] for f in resampled])
    path = out_dir / "ce_subtree_canonical_curves.npz"
    np.savez(path, u_grid=U_GRID, subtree_names=np.array(names),
             fronts_X=fX, fronts_Y=fY,
             med_X=ens["med_X"], med_Y=ens["med_Y"],
             band_lo_X=ens["band_lo_X"], band_lo_Y=ens["band_lo_Y"],
             band_hi_X=ens["band_hi_X"], band_hi_Y=ens["band_hi_Y"],
             medoid_idx=medoid_idx)
    return path


# ═══════════════════════════════════════════════════════════════
# Validation (--selftest)
# ═══════════════════════════════════════════════════════════════

def _check(checks, name, ok, detail=""):
    checks.append((name, bool(ok), str(detail)))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")


def _projection_unit_tests():
    """Item 3: segment projection on endpoints, a line segment, on-front."""
    x = np.array([0.0, 1.0, 1.0])       # A -> diagonal segment -> B
    y = np.array([1.0, 0.5, 0.0])
    u = np.array([0.0, 0.5, 1.0])
    results = []
    uq, pq, dq = project_point_to_polyline(np.array([0.0, 1.0]), x, y, u)
    results.append(("vertex A", uq == 0.0 and dq < 1e-12, (uq, pq, dq)))
    uq, pq, dq = project_point_to_polyline(np.array([1.0, 0.0]), x, y, u)
    results.append(("vertex B", uq == 1.0 and dq < 1e-12, (uq, pq, dq)))
    # Interior point of the first segment: (0.5, 0.75).
    uq, pq, dq = project_point_to_polyline(np.array([0.5, 0.75]), x, y, u)
    results.append(("segment interior", dq < 1e-12 and abs(uq - 0.25) < 1e-12,
                    (uq, pq, dq)))
    # Point beyond endpoint A: projection must clamp to A (u = 0).
    uq, pq, dq = project_point_to_polyline(np.array([-1.0, 1.0]), x, y, u)
    results.append(("beyond endpoint (clamps)", abs(uq - 0.0) < 1e-12
                    and np.allclose(pq, [0.0, 1.0]) and abs(dq - 1.0) < 1e-12,
                    (uq, pq, dq)))
    return results


def _selftest(args):
    print("=" * 78)
    print("fig5_table1_ce_canonical_metrics --selftest (validation checklist)")
    print("=" * 78)
    ctx = load_analysis()
    checks = []

    # ── Item 3: projection unit tests ──
    print("\n[3] Projection unit tests")
    for label, ok, detail in _projection_unit_tests():
        _check(checks, f"projection/{label}", ok, detail)

    # ── Full run at the publication threshold ──
    print(f"\nFull analysis (min_cells={args.min_cells}, "
          f"iteration={args.iteration})")
    records = analyze_all(min_cells=args.min_cells,
                          iteration=args.iteration, ctx=ctx)
    valid = [r for r in records if r["endpoint_ok"]]

    # ── Item 1: endpoint mapping exact ──
    print("\n[1] Endpoint mapping A=(0,1), B=(1,0)")
    bad = [r["subtree"] for r in valid
           if not (np.isclose(r["X"][0], 0.0, atol=1e-9)
                   and np.isclose(r["Y"][0], 1.0, atol=1e-9)
                   and np.isclose(r["X"][-1], 1.0, atol=1e-9)
                   and np.isclose(r["Y"][-1], 0.0, atol=1e-9))]
    _check(checks, "endpoint mapping", not bad, f"{len(valid)} valid, bad={bad[:5]}")

    # ── Item 2: u monotone, 0 -> 1 ──
    print("\n[2] Arc-length parameterization")
    bad = [r["subtree"] for r in valid if not r["u_monotone"]]
    ends_ok = all(abs(r["u"][0]) < 1e-12 and abs(r["u"][-1] - 1.0) < 1e-12
                  for r in valid)
    _check(checks, "u monotone/0..1", not bad and ends_ok,
           f"non-monotone={bad[:5]}")

    # ── Item 4: r vs the established rel_pareto_dist (same convention) ──
    print("\n[4] relative_distance vs summary-table rel_pareto_dist")
    summary_path = OUT.parent / f"subtree_summary_min{args.min_cells}.csv"
    if summary_path.exists():
        summ = pd.read_csv(summary_path)
        mismatches = []
        dlineage_mismatch = []
        for r in valid:
            row = summ[summ["subtree"] == r["subtree"]]
            if not len(row):
                continue
            ref = float(row.iloc[0]["rel_pareto_dist"])
            if not np.isclose(r["relative_distance"], ref,
                              rtol=1e-6, atol=1e-6):
                mismatches.append((r["subtree"], r["relative_distance"],
                                   ref))
            # d_lineage must equal the closest-sampled-point distance in the
            # RAW front (the established convention), not a segment distance.
            xa, ea = r["_front"][0], r["_front"][1]
            dmin = float(np.sqrt(xa ** 2 + ea ** 2).min())
            if not np.isclose(r["d_lineage"], dmin, rtol=1e-9, atol=1e-12):
                dlineage_mismatch.append((r["subtree"], r["d_lineage"], dmin))
        _check(checks, "relative_distance matches rel_pareto_dist",
               not mismatches, f"mismatches={mismatches[:5]}")
        _check(checks, "d_lineage = closest sampled point (not segment)",
               not dlineage_mismatch,
               f"mismatches={dlineage_mismatch[:5]}")
    else:
        _check(checks, "relative_distance matches rel_pareto_dist", False,
               f"summary CSV not found: {summary_path}")

    # Primary d_LP must likewise use an attainable vertex, but in the
    # endpoint-normalized geometry rather than null-SD coordinates.
    dlp_mismatch = []
    for r in valid:
        lineage = np.asarray(r["L_norm"])
        vertices = np.column_stack([r["X"], r["Y"]])
        distance = np.linalg.norm(vertices - lineage, axis=1)
        closest = int(np.argmin(distance))
        if (not np.isclose(r["d_lp"], distance[closest], atol=1e-12)
                or not np.isclose(r["u_lineage_lp"], r["u"][closest],
                                  atol=1e-12)):
            dlp_mismatch.append(r["subtree"])
    _check(checks, "d_LP = closest endpoint-normalized sampled point",
           not dlp_mismatch, f"mismatches={dlp_mismatch[:5]}")

    # d_NP uses the same endpoint-normalized P* selected for d_LP.
    dnp_mismatch = []
    for r in valid:
        p_star = np.asarray(r["P_lp_norm"])
        null = np.asarray(r["N_norm"])
        expected = float(np.linalg.norm(null - p_star))
        if not np.isclose(r["d_np"], expected, atol=1e-12):
            dnp_mismatch.append(r["subtree"])
    _check(checks, "d_NP uses the same endpoint-normalized P* as d_LP",
           not dnp_mismatch, f"mismatches={dnp_mismatch[:5]}")

    # ── Item 5: exactly-on-front subtrees ──
    print("\n[5] Exactly-on-front subtrees")
    on = [r for r in valid if r["exactly_on_front"]]
    # (a) segment convention consistent: grid-on-front implies d_lineage ~ 0.
    grid_on = [r for r in valid if r["natural_on_grid"]]
    inconsistent = [r["subtree"] for r in grid_on
                    if not r["exactly_on_front"]]
    # (b) cross-check the grid convention against the summary table.
    grid_mismatch = []
    if summary_path.exists():
        summ = pd.read_csv(summary_path)
        for r in valid:
            row = summ[summ["subtree"] == r["subtree"]]
            if not len(row):
                continue
            ref = bool(row.iloc[0]["natural_on_front"])
            if ref != r["natural_on_grid"]:
                grid_mismatch.append((r["subtree"], ref,
                                      r["natural_on_grid"]))
    _check(checks, "on-front consistency",
           not inconsistent and not grid_mismatch,
           f"{len(on)}/{len(valid)} on front; grid-vs-summary mismatches="
           f"{grid_mismatch[:5]}")

    # ── Item 6: max-ER markers correspond to max(edge_retention) ──
    print("\n[6] Max-ER marker / span validity")
    bad6 = []
    for r in valid:
        if not len(r["max_er_spans"]):
            bad6.append((r["subtree"], "no spans"))
            continue
        e = r["edge"]
        u = r["u"]
        for ua, ub in r["max_er_spans"]:
            sel = (u >= ua - 1e-12) & (u <= ub + 1e-12)
            if not np.all(np.isclose(e[sel], r["max_er"], rtol=1e-6,
                                     atol=1e-8)):
                bad6.append((r["subtree"], "span contains sub-max ER"))
                break
        # The representative circle must sit on an ACTUAL max-ER vertex.
        i_rep = int(np.argmin(np.abs(u - r["u_max_er"])))
        if not np.isclose(e[i_rep], r["max_er"], rtol=1e-6, atol=1e-8):
            bad6.append((r["subtree"], "u_max_er not on a max-ER point"))
    _check(checks, "max-ER spans + representative valid", not bad6,
           f"bad={bad6[:5]}")

    # ── Item 7 (part 1): invariance to duplicate alpha solutions ──
    print("\n[7a] Invariance to duplicated alpha solutions (dedup)")
    bad7 = []
    for r in valid[:8]:
        # Interleave duplicates of every second point (consecutive runs).
        xi, yi, ei = [], [], []
        for i in range(len(r["_front"][0])):
            xi += [r["_front"][0][i]] * (2 if i % 2 else 1)
            yi += [r["_front"][1][i]] * (2 if i % 2 else 1)
            ei += [r["_front"][2][i]] * (2 if i % 2 else 1)
        rec2 = compute_front_landmarks(np.array(xi), np.array(yi),
                                       np.array(ei),
                                       nx=r["_null"][0], ny=r["_null"][1])
        if (abs(rec2["u_lineage"] - r["u_lineage"]) > 1e-9
                or abs(rec2["u_lineage_lp"] - r["u_lineage_lp"]) > 1e-9
                or abs(rec2["d_lp"] - r["d_lp"]) > 1e-9
                or abs(rec2["u_max_er"] - r["u_max_er"]) > 1e-9
                or abs(rec2["relative_distance"]
                       - r["relative_distance"]) > 1e-9):
            bad7.append(r["subtree"])
    _check(checks, "dedup invariance", not bad7, f"bad={bad7[:5]}")

    # ── Item 7 (part 2): sweep 300 -> 1000 stability on a size spread ──
    print("\n[7b] Sweep 300 -> 1000 stability (size-spread subset)")
    subtrees = dl.collect_all_subtrees(ctx["lineage"], ctx["v_prot"],
                                       min_cells=args.min_cells)
    sizes = sorted(range(len(subtrees)),
                   key=lambda i: -len(subtrees[i][1]))[::6]
    subset = [subtrees[i] for i in sizes[:8]]
    deltas = []
    for name, terms in subset:
        r300 = analyze_subtree(name, terms, ctx, iteration=300)
        r1000 = analyze_subtree(name, terms, ctx, iteration=1000)
        if not (r300["endpoint_ok"] and r1000["endpoint_ok"]):
            continue
        deltas.append((name,
                       abs(r1000["u_lineage_lp"] - r300["u_lineage_lp"]),
                       abs(r1000["u_max_er"] - r300["u_max_er"]),
                       abs(r1000["relative_distance"]
                           - r300["relative_distance"]),
                       abs(r1000["d_lp"] - r300["d_lp"])))
    dmax_u = max((d[1] for d in deltas), default=0.0)
    dmax_er = max((d[2] for d in deltas), default=0.0)
    dmax_r = max((d[3] for d in deltas), default=0.0)
    dmax_lp = max((d[4] for d in deltas), default=0.0)
    _check(checks, "sweep stability",
           dmax_u < 0.01 and dmax_r < 0.02 and dmax_lp < 0.02,
           f"{len(deltas)} subtrees: max|du_L|={dmax_u:.4f} "
           f"max|du_ER|={dmax_er:.4f} max|dr|={dmax_r:.4f} "
           f"max|d(d_LP)|={dmax_lp:.4f}")

    # ── Item 8: median vs medoid ──
    print("\n[8] Median curve vs unweighted medoid front")
    resampled = [resample_front(r["X"], r["Y"], r["u"]) for r in valid]
    med_idx, _dist = medoid_front(resampled)
    ens = ensemble_curves(resampled)
    medoid_X, medoid_Y = resampled[med_idx]
    diff = np.sqrt((ens["med_X"] - medoid_X) ** 2
                   + (ens["med_Y"] - medoid_Y) ** 2).max()
    _check(checks, "median close to medoid", diff < 0.05,
           f"medoid={valid[med_idx]['subtree']} max|median-medoid|={diff:.4f}")

    # ── Item 9: min_cells sensitivity (shared values must agree) ──
    print("\n[9] min_cells sensitivity (shared subtree values)")
    if args.min_cells == 12:
        rec10 = analyze_all(min_cells=10, iteration=args.iteration, ctx=ctx)
        rec15 = analyze_all(min_cells=15, iteration=args.iteration, ctx=ctx)
        tbl = metrics_table(records)
        tbl10 = metrics_table(rec10)
        tbl15 = metrics_table(rec15)
        shared_cols = ["n", "x_scale", "y_scale", "u_lineage",
                       "u_lineage_lp", "d_lp", "d_np", "u_max_er", "max_er",
                       "d_lineage", "d_null", "relative_distance"]
        m10 = tbl10.set_index("subtree")[shared_cols]
        m15 = tbl15.set_index("subtree")[shared_cols]
        m12 = tbl.set_index("subtree")[shared_cols]
        shared = m12.index.intersection(m15.index).intersection(m10.index)
        ok10 = np.allclose(m10.reindex(shared).values,
                           m12.reindex(shared).values,
                           rtol=1e-9, atol=1e-12)
        ok15 = np.allclose(m15.reindex(shared).values,
                           m12.reindex(shared).values,
                           rtol=1e-9, atol=1e-12)
        _check(checks, "shared values agree across min_cells",
               ok10 and ok15,
               f"{len(shared)} shared subtrees; min10={len(m10)} rows, "
               f"min15={len(m15)} rows")
    else:
        _check(checks, "shared values agree across min_cells", True,
               "skipped (run with --min-cells 12)")

    # ── Item 10: compile + output regeneration ──
    print("\n[10] Outputs")
    tbl = metrics_table(records)
    csv_path = ANALYSIS_OUT / "ce_subtree_canonical_metrics.csv"
    tbl.to_csv(csv_path, index=False)
    ens = ensemble_curves(resampled)
    med_idx, _ = medoid_front(resampled)
    npz_path = write_curves_npz(valid, resampled, ens, med_idx)
    table_path = write_stats_table(records, args.min_cells, out_dir=OUT,
                                  ctx=ctx)
    _check(checks, "outputs written",
           csv_path.exists() and npz_path.exists() and table_path.exists(),
           f"{csv_path.name}, {npz_path.name}, {table_path.name}")

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    print(f"\n{'=' * 78}\n{len(checks) - n_fail}/{len(checks)} checks passed")
    if n_fail:
        raise SystemExit(1)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-cells", type=int, default=MIN_CELLS)
    parser.add_argument("--iteration", type=int, default=ITERATION)
    parser.add_argument("--selftest", action="store_true",
                        help="Run the validation checklist and exit.")
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args(argv)

    ps.configure()

    if args.selftest:
        return _selftest(args)

    ctx = load_analysis()
    print(f"Canonical subtree Pareto map (min_cells={args.min_cells}, "
          f"iteration={args.iteration})")
    records = analyze_all(min_cells=args.min_cells,
                          iteration=args.iteration, ctx=ctx)
    valid = [r for r in records if r["endpoint_ok"]]
    print(f"{len(valid)}/{len(records)} subtrees valid for the canonical map.")

    # Metric table BEFORE plotting.
    tbl = metrics_table(records)
    csv_path = ANALYSIS_OUT / "ce_subtree_canonical_metrics.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tbl.to_csv(csv_path, index=False)
    print(f"Wrote {len(tbl)} rows to {csv_path}")

    # Table 1 (replaces the retired statistics
    # panel): canonical-map coordinates + summary metrics, resized to
    # textwidth for final manuscript width.
    write_stats_table(records, args.min_cells, out_dir=args.out, ctx=ctx)

    # Ensemble + curves npz.
    resampled = [resample_front(r["X"], r["Y"], r["u"]) for r in valid]
    ens = ensemble_curves(resampled)
    med_idx, dist = medoid_front(resampled)
    npz_path = write_curves_npz(valid, resampled, ens, med_idx)
    print(f"Wrote {npz_path.name} (medoid = {valid[med_idx]['subtree']})")

    # Size-weighted median diagnostic: does any subtree dominate the curve?
    wmed = ensemble_curves(resampled, weights=[r["n"] for r in valid])
    wdiff = np.sqrt((ens["med_X"] - wmed["med_X"]) ** 2
                    + (ens["med_Y"] - wmed["med_Y"]) ** 2).max()
    print(f"Max |unweighted - size-weighted median| = {wdiff:.4f} "
          f"(small => large subtrees do not dominate the curve)")

    # Figure 5 and Figure S1 are rendered from this metric cache by
    # fig5_figs1_ce_canonical_summary.py.
    return records


if __name__ == "__main__":
    main()
