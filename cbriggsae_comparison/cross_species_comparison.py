"""
Cross-species comparison: C. elegans vs C. briggsae embryo tracking.
Compares 5 datasets using pairwise distance correlations and travel metrics.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

REPO = Path("/home/bingran/code/embryogenesis")
OUTPUT_DIR = REPO / "cbriggsae_comparison" / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

sns.set_style("whitegrid")


# ================================================================
# Data loading
# ================================================================

def parse_sbd(filepath):
    with open(filepath) as f:
        content = f.read()
    blocks = content.split("---")
    cells = []
    for block in blocks:
        lines = [l.strip() for l in block.strip().split("\n") if l.strip()]
        if len(lines) < 5:
            continue
        p0 = lines[0].split()
        if len(p0) < 5:
            continue
        active = p0[0] == "1"
        name = p0[-1]
        n_obs_str = lines[3].split()[0]
        try:
            n_obs = int(n_obs_str)
        except ValueError:
            continue
        if n_obs == 0:
            continue
        if 4 + n_obs > len(lines):
            n_obs = len(lines) - 4
        frames, xs, ys = [], [], []
        for j in range(4, 4 + n_obs):
            op = lines[j].split()
            if len(op) < 4:
                continue
            frames.append(int(op[0]))
            xs.append(float(op[1]))
            ys.append(float(op[2]))
        cells.append({
            "cell": name, "active": active,
            "frame": np.array(frames), "x": np.array(xs), "y": np.array(ys),
        })
    return cells


# C. briggsae
cb_csv_raw = pd.read_csv(REPO / "data/c_briggsae/CD140715HLH1cbp1.csv")
cb_csv = cb_csv_raw[["cell", "time", "x", "y"]].copy()  # 2D only for fair comparison
cb_sbd2_raw = parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd")
cb_sbd3_raw = parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_3a.sbd")
cb_sbd2_active = [c for c in cb_sbd2_raw if c["active"]]
cb_sbd3_active = [c for c in cb_sbd3_raw if c["active"]]

# C. elegans
ce_tracks_raw = pd.read_csv(REPO / "data/embryo1/tracks.txt", sep="\t")
# Columns: t, z, y, x — x,y are already named correctly
ce_tracks = ce_tracks_raw[["name", "t", "x", "y"]].copy()
ce_tracks = ce_tracks.rename(columns={"name": "cell", "t": "time"})
ce_sbd_raw = parse_sbd(REPO / "data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd")
ce_sbd_active = [c for c in ce_sbd_raw if c["active"]]


# ================================================================
# Utility functions
# ================================================================

def get_csv_positions(df, timepoint):
    """Get (x, y) positions from a CSV-like DataFrame at a given timepoint."""
    subset = df[df.time == timepoint]
    return {row["cell"]: (row["x"], row["y"]) for _, row in subset.iterrows()}


def interpolate_sbd_positions(cells, target_frame):
    """Linearly interpolate (x,y) for all cells at target_frame."""
    positions = {}
    for c in cells:
        frames = c["frame"]
        if len(frames) < 2:
            continue
        if not (frames.min() <= target_frame <= frames.max()):
            continue
        x = np.interp(target_frame, frames, c["x"])
        y = np.interp(target_frame, frames, c["y"])
        positions[c["cell"]] = (x, y)
    return positions


def compute_distance_matrix(positions, cell_list):
    coords = np.array([positions[c] for c in cell_list])
    n = len(cell_list)
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(np.sum((coords[i] - coords[j]) ** 2))
            dm[i, j] = d
            dm[j, i] = d
    return dm


def find_best_frame(cells, csv_cells, max_f, step=5):
    best_f, best_c = 0, 0
    for frame in range(10, max_f, step):
        at_frame = {c["cell"] for c in cells
                     if c["frame"].min() <= frame <= c["frame"].max()}
        common = len(at_frame & csv_cells)
        if common > best_c:
            best_c = common
            best_f = frame
    return best_f, best_c


def correlate_pairwise(pos_a, pos_b):
    """Compute pairwise distance Spearman correlation between two position dicts."""
    common = sorted(set(pos_a.keys()) & set(pos_b.keys()))
    if len(common) < 10:
        return {"n_cells": len(common), "n_pairs": 0, "spearman_r": np.nan}
    dm_a = compute_distance_matrix(pos_a, common)
    dm_b = compute_distance_matrix(pos_b, common)
    triu_idx = np.triu_indices_from(dm_a, k=1)
    r, _ = spearmanr(dm_a[triu_idx], dm_b[triu_idx])
    return {"n_cells": len(common), "n_pairs": len(common) * (len(common) - 1) // 2,
            "spearman_r": r, "cells": common}


def compute_travel_2d(cells_2d_or_df, is_df=False, cell_col="cell", x_col="x", y_col="y",
                       time_col="time"):
    """Compute travel metrics for cells."""
    travel = {}
    if is_df:
        for cell, grp in cells_2d_or_df.groupby(cell_col):
            grp = grp.sort_values(time_col)
            if len(grp) < 2:
                continue
            xs, ys = grp[x_col].values, grp[y_col].values
            dx, dy = np.diff(xs), np.diff(ys)
            steps = np.sqrt(dx**2 + dy**2)
            travel[cell] = {
                "total": steps.sum(), "mean_step": steps.mean(),
                "n_steps": len(steps),
                "net_disp": np.sqrt((xs[-1] - xs[0])**2 + (ys[-1] - ys[0])**2),
            }
    else:
        for c in cells_2d_or_df:
            if len(c["x"]) < 2:
                continue
            dx, dy = np.diff(c["x"]), np.diff(c["y"])
            steps = np.sqrt(dx**2 + dy**2)
            travel[c["cell"]] = {
                "total": steps.sum(), "mean_step": steps.mean(),
                "n_steps": len(steps),
                "net_disp": np.sqrt(
                    (c["x"][-1] - c["x"][0])**2 + (c["y"][-1] - c["y"][0])**2),
            }
    return travel


def correlate_travel_metric(ta, tb, common, metric):
    va = np.array([ta[c][metric] for c in common])
    vb = np.array([tb[c][metric] for c in common])
    r, _ = spearmanr(va, vb)
    return r


# ================================================================
# Timepoint matching
# ================================================================

def find_comparable_timepoints(df_a, cells_b, max_f_b, is_a_csv=True, is_b_sbd=True,
                                late_fraction=0.7):
    """Find a matching timepoint pair that maximizes cell overlap at late stage."""
    if is_a_csv:
        t_min, t_max = int(df_a.time.min()), int(df_a.time.max())
        late_start = int(t_min + late_fraction * (t_max - t_min))
        best_t, best_f, best_c = 0, 0, 0
        for t in range(late_start, t_max + 1):
            pos_a = get_csv_positions(df_a, t)
            cells_a = set(pos_a.keys())
            bf, bc = find_best_frame(cells_b, cells_a, max_f_b)
            if bc > best_c:
                best_c = bc
                best_t = t
                best_f = bf
        return best_t, best_f, best_c
    else:
        # Both are SBD-like — pick comparable frames from each
        # Just use 80% coverage frame for both
        def endpoint_80(cells):
            fc = {}
            for c in cells:
                for f in c["frame"]:
                    fc[f] = fc.get(f, 0) + 1
            peak = max(fc.values())
            return max(f for f, cnt in fc.items() if cnt >= peak * 0.8)

        f_a = endpoint_80(df_a)
        f_b = endpoint_80(cells_b)
        pos_a = interpolate_sbd_positions(df_a, f_a)
        common = len(set(pos_a.keys()) & {c["cell"] for c in cells_b
                                          if c["frame"].min() <= f_b <= c["frame"].max()})
        return f_a, f_b, common


# ================================================================
# Run comparisons
# ================================================================

# Pre-compute travel metrics for all datasets
print("Computing travel metrics...")
cb_csv_travel = compute_travel_2d(cb_csv, is_df=True, cell_col="cell")
cb_sbd2_travel = compute_travel_2d(cb_sbd2_active, is_df=False)
cb_sbd3_travel = compute_travel_2d(cb_sbd3_active, is_df=False)
ce_tracks_travel = compute_travel_2d(ce_tracks, is_df=True, cell_col="cell")
ce_sbd_travel = compute_travel_2d(ce_sbd_active, is_df=False)

# Comparison pairs: (label, type_a, data_a, type_b, data_b)
# type: "csv" = evenly-spaced timepoints, "sbd" = irregular frames
comparisons = [
    # (label, type_a, data_a, type_b, data_b, group)
    # Within C. briggsae (reference, already known)
    ("Cb-CSV vs Cb-SBD2", "csv", cb_csv, "sbd", cb_sbd2_active, "C. briggsae"),
    ("Cb-CSV vs Cb-SBD3", "csv", cb_csv, "sbd", cb_sbd3_active, "C. briggsae"),
    # Within C. elegans (new)
    ("Ce-tracks vs Ce-SBD", "csv", ce_tracks, "sbd", ce_sbd_active, "C. elegans"),
    # Cross-species: gold-standard to gold-standard
    ("Cb-CSV vs Ce-tracks", "csv", cb_csv, "csv", ce_tracks, "Cross-species"),
    # Cross-species: SBD to SBD
    ("Cb-SBD2 vs Ce-SBD", "sbd", cb_sbd2_active, "sbd", ce_sbd_active, "Cross-species"),
    # Cross-species: mixed
    ("Cb-CSV vs Ce-SBD", "csv", cb_csv, "sbd", ce_sbd_active, "Cross-species"),
    ("Ce-tracks vs Cb-SBD2", "csv", ce_tracks, "sbd", cb_sbd2_active, "Cross-species"),
]

# For CSV data, get max time; for SBD data, get max frame
def get_max(data, dtype):
    if dtype == "csv":
        return int(data.time.max())
    else:
        return int(max(c["frame"].max() for c in data))

results = []

for label, type_a, data_a, type_b, data_b, group in comparisons:
    print(f"\n{'='*60}")
    print(f"{label} [{group}]")
    print(f"{'='*60}")

    max_a = get_max(data_a, type_a)
    max_b = get_max(data_b, type_b)

    # ---- Pairwise distance correlation ----
    if type_a == "csv" and type_b == "csv":
        # Both CSV-like: try matching at late timepoints
        t_a_min, t_a_max = int(data_a.time.min()), int(data_a.time.max())
        t_b_min, t_b_max = int(data_b.time.min()), int(data_b.time.max())
        # Use same fraction of development
        for frac in [0.65, 0.70, 0.75, 0.80]:
            t_a = int(t_a_min + frac * (t_a_max - t_a_min))
            t_b = int(t_b_min + frac * (t_b_max - t_b_min))
            pos_a = get_csv_positions(data_a, t_a)
            pos_b = get_csv_positions(data_b, t_b)
            common = len(set(pos_a.keys()) & set(pos_b.keys()))
            if common > 50:
                break
        pw_result = correlate_pairwise(pos_a, pos_b)
        print(f"  Pairwise: t_a={t_a}, t_b={t_b}, "
              f"n_cells={pw_result['n_cells']}, n_pairs={pw_result['n_pairs']:,}, "
              f"Spearman r={pw_result['spearman_r']:.4f}")

    elif type_a == "sbd" and type_b == "sbd":
        # Both SBD: use 80% coverage frames
        def endpoint_80(cells):
            fc = {}
            for c in cells:
                for f in c["frame"]:
                    fc[f] = fc.get(f, 0) + 1
            peak = max(fc.values())
            return max(f for f, cnt in fc.items() if cnt >= peak * 0.8)

        f_a = endpoint_80(data_a)
        f_b = endpoint_80(data_b)
        pos_a = interpolate_sbd_positions(data_a, f_a)
        pos_b = interpolate_sbd_positions(data_b, f_b)
        pw_result = correlate_pairwise(pos_a, pos_b)
        print(f"  Pairwise: f_a={f_a}, f_b={f_b}, "
              f"n_cells={pw_result['n_cells']}, n_pairs={pw_result['n_pairs']:,}, "
              f"Spearman r={pw_result['spearman_r']:.4f}")

    else:
        # Mixed CSV + SBD
        if type_a == "csv":
            df_csv, cells_sbd = data_a, data_b
            max_f = max_b
            csv_label, sbd_label = "a", "b"
        else:
            df_csv, cells_sbd = data_b, data_a
            max_f = max_a
            csv_label, sbd_label = "b", "a"

        # Find best matching late timepoint
        t_min, t_max = int(df_csv.time.min()), int(df_csv.time.max())
        best_t, best_f, best_c = 0, 0, 0
        for t in range(int(t_min + 0.6 * (t_max - t_min)), t_max + 1, 5):
            pos_csv = get_csv_positions(df_csv, t)
            cells_set = set(pos_csv.keys())
            bf, bc = find_best_frame(cells_sbd, cells_set, max_f)
            if bc > best_c:
                best_c = bc
                best_t = t
                best_f = bf

        pos_csv = get_csv_positions(df_csv, best_t)
        pos_sbd = interpolate_sbd_positions(cells_sbd, best_f)
        pw_result = correlate_pairwise(pos_csv, pos_sbd)
        print(f"  Pairwise: t={best_t}, f={best_f}, "
              f"n_cells={pw_result['n_cells']}, n_pairs={pw_result['n_pairs']:,}, "
              f"Spearman r={pw_result['spearman_r']:.4f}")

    # ---- Travel distance correlations ----
    travel_maps = {
        "Cb-CSV": cb_csv_travel,
        "Cb-SBD2": cb_sbd2_travel,
        "Cb-SBD3": cb_sbd3_travel,
        "Ce-tracks": ce_tracks_travel,
        "Ce-SBD": ce_sbd_travel,
    }
    ta_name = label.split(" vs ")[0]
    tb_name = label.split(" vs ")[1]
    ta = travel_maps[ta_name]
    tb = travel_maps[tb_name]
    common_t = sorted(set(ta.keys()) & set(tb.keys()))

    travel_rs = {}
    for metric in ["total", "mean_step", "net_disp", "n_steps"]:
        travel_rs[metric] = correlate_travel_metric(ta, tb, common_t, metric)

    print(f"  Travel: n={len(common_t)} cells")
    print(f"    Total={travel_rs['total']:.4f}, MeanStep={travel_rs['mean_step']:.4f}, "
          f"NetDisp={travel_rs['net_disp']:.4f}, NSteps={travel_rs['n_steps']:.4f}")

    results.append({
        "label": label, "group": group,
        "pw_spearman_r": pw_result["spearman_r"],
        "pw_n_cells": pw_result["n_cells"],
        "pw_n_pairs": pw_result["n_pairs"],
        "travel_n": len(common_t),
        "travel_total": travel_rs["total"],
        "travel_mean_step": travel_rs["mean_step"],
        "travel_net_disp": travel_rs["net_disp"],
        "travel_n_steps": travel_rs["n_steps"],
    })


# ================================================================
# Summary table and visualization
# ================================================================
print("\n\n" + "=" * 80)
print("SUMMARY: Cross-species comparison")
print("=" * 80)

print(f"\n{'Comparison':<28} {'Group':<16} {'PW_r':>8} {'PW_cells':>9} {'PW_pairs':>10} "
      f"{'Tr_MeanStep':>12} {'Tr_Total':>10} {'Tr_NetDisp':>10} {'Tr_cells':>9}")
print("-" * 112)
for r in results:
    print(f"{r['label']:<28} {r['group']:<16} {r['pw_spearman_r']:>8.4f} "
          f"{r['pw_n_cells']:>9} {r['pw_n_pairs']:>10,} "
          f"{r['travel_mean_step']:>12.4f} {r['travel_total']:>10.4f} "
          f"{r['travel_net_disp']:>10.4f} {r['travel_n']:>9}")

# Save results
df_results = pd.DataFrame(results)
df_results.to_csv(OUTPUT_DIR / "cross_species_results.csv", index=False)
print(f"\nResults saved to {OUTPUT_DIR / 'cross_species_results.csv'}")


# ================================================================
# Visualization
# ================================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel 1: Pairwise distance correlations
ax = axes[0]
groups_order = ["C. briggsae", "C. elegans", "Cross-species"]
colors_map = {"C. briggsae": "#2196F3", "C. elegans": "#4CAF50", "Cross-species": "#FF5722"}
x_pos = []
labels_pw = []
values_pw = []
colors_pw = []

for i, r in enumerate(results):
    x_pos.append(i)
    short_label = r["label"].replace("Cb-", "C.b ").replace("Ce-", "C.e ").replace(" vs ", "\nvs\n")
    labels_pw.append(short_label)
    values_pw.append(r["pw_spearman_r"])
    colors_pw.append(colors_map[r["group"]])

bars = ax.bar(range(len(results)), values_pw, color=colors_pw, alpha=0.85, edgecolor='white')
for bar, val in zip(bars, values_pw):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{val:.3f}', ha='center', fontsize=9, fontweight='bold')
ax.set_xticks(range(len(results)))
ax.set_xticklabels(labels_pw, fontsize=7)
ax.set_ylabel("Spearman Rank Correlation")
ax.set_title("Pairwise Distance Correlation (static spatial arrangement)")
ax.set_ylim(0, 1.1)
ax.grid(True, alpha=0.3, axis='y')
# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors_map[g], label=g) for g in groups_order]
ax.legend(handles=legend_elements, fontsize=9)

# Panel 2: Travel mean_step correlations (best metric) vs pairwise
ax = axes[1]
x = np.arange(len(results))
width = 0.35
bars1 = ax.bar(x - width/2, [r["pw_spearman_r"] for r in results], width,
               label="Pairwise distance (static)", color="#607D8B", alpha=0.85)
bars2 = ax.bar(x + width/2, [r["travel_mean_step"] for r in results], width,
               label="Mean step (dynamic)", color="#FF9800", alpha=0.85)
ax.set_xticks(x)
ax.set_xticklabels(labels_pw, fontsize=7)
ax.set_ylabel("Spearman Rank Correlation")
ax.set_title("Static (pairwise distance) vs Dynamic (mean step) correlation")
ax.set_ylim(0, 1.1)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "cross_species_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved cross_species_comparison.png")
