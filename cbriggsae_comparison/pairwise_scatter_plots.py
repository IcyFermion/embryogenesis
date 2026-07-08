"""
Generate pairwise distance scatter plots for all 7 comparisons.
Each plot shows: for every cell pair, distance in dataset A vs distance in dataset B.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

REPO = Path("/home/bingran/code/embryogenesis")
OUTPUT_DIR = REPO / "cbriggsae_comparison" / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


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


print("Loading data...")
cb_csv_raw = pd.read_csv(REPO / "data/c_briggsae/CD140715HLH1cbp1.csv")
cb_csv = cb_csv_raw[["cell", "time", "x", "y"]].copy()
cb_sbd2 = [c for c in parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd") if c["active"]]
cb_sbd3 = [c for c in parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_3a.sbd") if c["active"]]
ce_tracks_raw = pd.read_csv(REPO / "data/embryo1/tracks.txt", sep="\t")
ce_tracks = ce_tracks_raw[["name", "t", "x", "y"]].copy()
ce_tracks = ce_tracks.rename(columns={"name": "cell", "t": "time"})
ce_sbd = [c for c in parse_sbd(REPO / "data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd") if c["active"]]


# ================================================================
# Helper functions
# ================================================================
def get_csv_positions(df, timepoint):
    subset = df[df.time == timepoint]
    return {row["cell"]: (row["x"], row["y"]) for _, row in subset.iterrows()}


def interpolate_sbd_positions(cells, target_frame):
    positions = {}
    for c in cells:
        if len(c["frame"]) < 2:
            continue
        if not (c["frame"].min() <= target_frame <= c["frame"].max()):
            continue
        x = np.interp(target_frame, c["frame"], c["x"])
        y = np.interp(target_frame, c["frame"], c["y"])
        positions[c["cell"]] = (x, y)
    return positions


def compute_pairwise_flat(positions, cell_list):
    """Return flattened upper-triangle pairwise distances."""
    coords = np.array([positions[c] for c in cell_list])
    n = len(cell_list)
    flats = []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(np.sum((coords[i] - coords[j]) ** 2))
            flats.append(d)
    return np.array(flats)


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


def endpoint_80(cells):
    fc = {}
    for c in cells:
        for f in c["frame"]:
            fc[f] = fc.get(f, 0) + 1
    peak = max(fc.values())
    return max(f for f, cnt in fc.items() if cnt >= peak * 0.8)


# ================================================================
# Compute all pairwise distance vectors for all 7 comparisons
# ================================================================

comparisons_config = [
    # (label, group, type_a, data_a, type_b, data_b)
    ("Cb-CSV vs Cb-SBD2",    "C. briggsae",   "csv", cb_csv, "sbd", cb_sbd2),
    ("Cb-CSV vs Cb-SBD3",    "C. briggsae",   "csv", cb_csv, "sbd", cb_sbd3),
    ("Ce-tracks vs Ce-SBD",  "C. elegans",    "csv", ce_tracks, "sbd", ce_sbd),
    ("Cb-CSV vs Ce-tracks",  "Cross-species", "csv", cb_csv, "csv", ce_tracks),
    ("Cb-SBD2 vs Ce-SBD",    "Cross-species", "sbd", cb_sbd2, "sbd", ce_sbd),
    ("Cb-CSV vs Ce-SBD",     "Cross-species", "csv", cb_csv, "sbd", ce_sbd),
    ("Ce-tracks vs Cb-SBD2", "Cross-species", "csv", ce_tracks, "sbd", cb_sbd2),
]

plot_data = []

for label, group, type_a, data_a, type_b, data_b in comparisons_config:
    print(f"\n{label} [{group}]")

    max_a = int(data_a.time.max()) if type_a == "csv" else int(max(c["frame"].max() for c in data_a))
    max_b = int(data_b.time.max()) if type_b == "csv" else int(max(c["frame"].max() for c in data_b))

    # ---- Determine timepoint matching ----
    if type_a == "csv" and type_b == "csv":
        t_a_min, t_a_max = int(data_a.time.min()), int(data_a.time.max())
        t_b_min, t_b_max = int(data_b.time.min()), int(data_b.time.max())
        for frac in [0.65, 0.70, 0.75]:
            t_a = int(t_a_min + frac * (t_a_max - t_a_min))
            t_b = int(t_b_min + frac * (t_b_max - t_b_min))
            pos_a = get_csv_positions(data_a, t_a)
            pos_b = get_csv_positions(data_b, t_b)
            if len(set(pos_a.keys()) & set(pos_b.keys())) > 50:
                break
        tp_label = f"t_a={t_a}, t_b={t_b}"

    elif type_a == "sbd" and type_b == "sbd":
        f_a = endpoint_80(data_a)
        f_b = endpoint_80(data_b)
        pos_a = interpolate_sbd_positions(data_a, f_a)
        pos_b = interpolate_sbd_positions(data_b, f_b)
        t_a, t_b = f_a, f_b
        tp_label = f"f_a={f_a}, f_b={f_b}"

    else:
        if type_a == "csv":
            df_csv, cells_sbd = data_a, data_b
            max_f = max_b
        else:
            df_csv, cells_sbd = data_b, data_a
            max_f = max_a

        t_min, t_max = int(df_csv.time.min()), int(df_csv.time.max())
        best_t, best_f, best_c = 0, 0, 0
        for t in range(int(t_min + 0.6 * (t_max - t_min)), t_max + 1, 5):
            pos_csv = get_csv_positions(df_csv, t)
            bf, bc = find_best_frame(cells_sbd, set(pos_csv.keys()), max_f)
            if bc > best_c:
                best_c = bc
                best_t = t
                best_f = bf

        pos_csv = get_csv_positions(df_csv, best_t)
        pos_sbd = interpolate_sbd_positions(cells_sbd, best_f)

        if type_a == "csv":
            pos_a, pos_b = pos_csv, pos_sbd
            t_a, t_b = best_t, best_f
            tp_label = f"t={best_t}, f={best_f}"
        else:
            pos_a, pos_b = pos_sbd, pos_csv
            t_a, t_b = best_f, best_t
            tp_label = f"f={best_f}, t={best_t}"

    # ---- Compute pairwise distances ----
    common = sorted(set(pos_a.keys()) & set(pos_b.keys()))
    flat_a = compute_pairwise_flat(pos_a, common)
    flat_b = compute_pairwise_flat(pos_b, common)
    sr, _ = spearmanr(flat_a, flat_b)

    print(f"  {tp_label}: n_cells={len(common)}, n_pairs={len(flat_a):,}, Spearman r={sr:.4f}")

    plot_data.append({
        "label": label, "group": group, "tp_label": tp_label,
        "n_cells": len(common), "n_pairs": len(flat_a),
        "spearman_r": sr, "flat_a": flat_a, "flat_b": flat_b,
        "name_a": label.split(" vs ")[0], "name_b": label.split(" vs ")[1],
    })


# ================================================================
# Generate scatter plots
# ================================================================
print("\nGenerating plots...")

# ---- Panel 1: Within-species (3 plots) ----
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
within = [d for d in plot_data if d["group"] != "Cross-species"]
for idx, d in enumerate(within):
    ax = axes[idx]
    # Subsample for plotting if too many pairs (use 50k max)
    n = len(d["flat_a"])
    if n > 50000:
        idx_sub = np.random.RandomState(42).choice(n, 50000, replace=False)
        fa, fb = d["flat_a"][idx_sub], d["flat_b"][idx_sub]
        sub_label = f" (showing 50k of {n:,} pairs)"
    else:
        fa, fb = d["flat_a"], d["flat_b"]
        sub_label = ""

    hb = ax.hexbin(fa, fb, gridsize=60, cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax, label="log10(count)")
    max_val = max(fa.max(), fb.max())
    ax.plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="y=x")
    ax.set_xlabel(f"{d['name_a']} pairwise distance (px)")
    ax.set_ylabel(f"{d['name_b']} pairwise distance (px)")
    ax.set_title(f"{d['label']}\n{d['tp_label']}, n={d['n_cells']} cells, "
                 f"Spearman r={d['spearman_r']:.4f}{sub_label}")
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "pairwise_scatter_within_species.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved pairwise_scatter_within_species.png")

# ---- Panel 2: Cross-species (4 plots) ----
fig, axes = plt.subplots(2, 2, figsize=(14, 13))
cross = [d for d in plot_data if d["group"] == "Cross-species"]
for idx, d in enumerate(cross):
    ax = axes[idx // 2, idx % 2]
    n = len(d["flat_a"])
    if n > 50000:
        idx_sub = np.random.RandomState(42).choice(n, 50000, replace=False)
        fa, fb = d["flat_a"][idx_sub], d["flat_b"][idx_sub]
        sub_label = f" (showing 50k of {n:,} pairs)"
    else:
        fa, fb = d["flat_a"], d["flat_b"]
        sub_label = ""

    hb = ax.hexbin(fa, fb, gridsize=60, cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax, label="log10(count)")
    max_val = max(fa.max(), fb.max())
    ax.plot([0, max_val], [0, max_val], "r--", alpha=0.5, label="y=x")
    ax.set_xlabel(f"{d['name_a']} pairwise distance (px)")
    ax.set_ylabel(f"{d['name_b']} pairwise distance (px)")
    ax.set_title(f"{d['label']}\n{d['tp_label']}, n={d['n_cells']} cells, "
                 f"Spearman r={d['spearman_r']:.4f}{sub_label}")
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "pairwise_scatter_cross_species.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved pairwise_scatter_cross_species.png")

# ---- Panel 3: All 7 in one grid ----
fig, axes = plt.subplots(2, 4, figsize=(22, 11))
axes_flat = axes.flatten()
for idx, d in enumerate(plot_data):
    ax = axes_flat[idx]
    n = len(d["flat_a"])
    if n > 30000:
        idx_sub = np.random.RandomState(42).choice(n, 30000, replace=False)
        fa, fb = d["flat_a"][idx_sub], d["flat_b"][idx_sub]
    else:
        fa, fb = d["flat_a"], d["flat_b"]

    hb = ax.hexbin(fa, fb, gridsize=50, cmap="viridis", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax, label="log10")
    max_val = max(fa.max(), fb.max())
    ax.plot([0, max_val], [0, max_val], "r--", alpha=0.4, linewidth=1, label="y=x")

    # Color border by group
    border_colors = {"C. briggsae": "#2196F3", "C. elegans": "#4CAF50", "Cross-species": "#FF5722"}
    for spine in ax.spines.values():
        spine.set_edgecolor(border_colors.get(d["group"], "gray"))
        spine.set_linewidth(2.5)

    ax.set_xlabel(d["name_a"], fontsize=8)
    ax.set_ylabel(d["name_b"], fontsize=8)
    ax.set_title(f"{d['label']}\n{d['tp_label']}, n={d['n_cells']} cells, r={d['spearman_r']:.3f}",
                 fontsize=9)
    ax.legend(fontsize=6)

# Hide the 8th subplot
axes_flat[7].set_visible(False)

# Add group legend at the bottom
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#2196F3", label="C. briggsae"),
    Patch(facecolor="#4CAF50", label="C. elegans"),
    Patch(facecolor="#FF5722", label="Cross-species"),
]
fig.legend(handles=legend_elements, loc="lower right", fontsize=10,
           title="Border color = group")

plt.suptitle("Pairwise Distance Correlation: All Comparisons", fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "pairwise_scatter_all_7.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved pairwise_scatter_all_7.png")

print("\nDone. All scatter plots generated.")
