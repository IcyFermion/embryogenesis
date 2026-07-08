"""Quick exploration of all 5 datasets for the cross-species comparison."""
import numpy as np
import pandas as pd
from pathlib import Path

REPO = Path("/home/bingran/code/embryogenesis")


def parse_sbd(filepath):
    """Parse SIMI*BIOCELL .sbd file."""
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


# ---- Load all 5 datasets ----
# 1. C. briggsae CSV
cb_csv = pd.read_csv(REPO / "data/c_briggsae/CD140715HLH1cbp1.csv")
cb_csv_cells = set(cb_csv["cell"].unique())

# 2 & 3. C. briggsae SBD
cb_sbd2_raw = parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd")
cb_sbd3_raw = parse_sbd(REPO / "data/c_briggsae/nadin/NM_C_briggsae_3a.sbd")
cb_sbd2 = [c for c in cb_sbd2_raw if c["active"]]
cb_sbd3 = [c for c in cb_sbd3_raw if c["active"]]
cb_sbd2_names = {c["cell"] for c in cb_sbd2}
cb_sbd3_names = {c["cell"] for c in cb_sbd3}

# 4. C. elegans tracks.txt
ce_tracks = pd.read_csv(REPO / "data/embryo1/tracks.txt", sep="\t")
# NOTE: columns are t, z, y, x — remap to standard x, y, z
ce_tracks = ce_tracks.rename(columns={"z": "z_coord", "y": "y_coord", "x": "x_coord"})
# Actually: the columns ARE t, z, y, x so z→z, y→y, x→x is correct
# Just rename for clarity
ce_tracks = ce_tracks.rename(columns={"z": "z_val", "x": "x_val"})
# Wait, let me just keep original: t, z, y, x — where z,y,x are the 3 coords
ce_tracks_cells = set(ce_tracks["name"].unique())

# 5. C. elegans SBD
ce_sbd_raw = parse_sbd(REPO / "data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd")
ce_sbd = [c for c in ce_sbd_raw if c["active"]]
ce_sbd_names = {c["cell"] for c in ce_sbd}

# ---- Dataset summaries ----
print("=" * 70)
print("DATASET SUMMARIES")
print("=" * 70)
print(f"C. briggsae CSV:     {len(cb_csv):>6} rows, {len(cb_csv_cells):>4} cells, t=[{cb_csv.time.min()},{cb_csv.time.max()}]")
print(f"C. briggsae SBD2:    {len(cb_sbd2):>4} active cells, frames=[{min(c['frame'].min() for c in cb_sbd2)},{max(c['frame'].max() for c in cb_sbd2)}]")
print(f"C. briggsae SBD3:    {len(cb_sbd3):>4} active cells, frames=[{min(c['frame'].min() for c in cb_sbd3)},{max(c['frame'].max() for c in cb_sbd3)}]")
print(f"C. elegans tracks:   {len(ce_tracks):>6} rows, {len(ce_tracks_cells):>4} cells, t=[{ce_tracks.t.min()},{ce_tracks.t.max()}]")
print(f"C. elegans SBD:      {len(ce_sbd):>4} active cells, frames=[{min(c['frame'].min() for c in ce_sbd)},{max(c['frame'].max() for c in ce_sbd)}]")

# ---- Coverage profiles ----
print("\n" + "=" * 70)
print("TEMPORAL COVERAGE")
print("=" * 70)

# C. elegans tracks: cells per timepoint
ce_t_cov = ce_tracks.groupby("t")["name"].nunique()
print(f"C. elegans tracks: max {ce_t_cov.max()} cells at t={ce_t_cov.idxmax()}, "
      f"cells at t={ce_tracks.t.max()}: {ce_t_cov.get(ce_tracks.t.max(), 0)}")

# C. elegans SBD: cells per frame
ce_sbd_fcov = {}
for c in ce_sbd:
    for f in c["frame"]:
        ce_sbd_fcov[f] = ce_sbd_fcov.get(f, 0) + 1
ce_peak = max(ce_sbd_fcov.values())
ce_peak_f = [f for f, c in ce_sbd_fcov.items() if c == ce_peak][0]
ce_80pct_f = max(f for f, c in ce_sbd_fcov.items() if c >= ce_peak * 0.8)
print(f"C. elegans SBD: peak {ce_peak} cells at frame {ce_peak_f}, "
      f"80% cutoff at frame {ce_80pct_f} ({ce_sbd_fcov[ce_80pct_f]} cells)")

# ---- Cell name overlaps ----
print("\n" + "=" * 70)
print("CELL NAME OVERLAP MATRIX")
print("=" * 70)

datasets = {
    "Cb-CSV": cb_csv_cells,
    "Cb-SBD2": cb_sbd2_names,
    "Cb-SBD3": cb_sbd3_names,
    "Ce-tracks": ce_tracks_cells,
    "Ce-SBD": ce_sbd_names,
}

names = list(datasets.keys())
print(f"{'':>12}", end="")
for n in names:
    print(f"{n:>10}", end="")
print()
for n1 in names:
    print(f"{n1:>12}", end="")
    for n2 in names:
        overlap = len(datasets[n1] & datasets[n2])
        print(f"{overlap:>10}", end="")
    print()

# ---- Key cross-species numbers ----
print("\n" + "=" * 70)
print("KEY CROSS-SPECIES OVERLAPS")
print("=" * 70)
print(f"Cb-CSV  ∩ Ce-tracks: {len(cb_csv_cells & ce_tracks_cells):>4} cells")
print(f"Cb-CSV  ∩ Ce-SBD:    {len(cb_csv_cells & ce_sbd_names):>4} cells")
print(f"Cb-SBD2 ∩ Ce-tracks: {len(cb_sbd2_names & ce_tracks_cells):>4} cells")
print(f"Cb-SBD2 ∩ Ce-SBD:   {len(cb_sbd2_names & ce_sbd_names):>4} cells")
print(f"Ce-tracks ∩ Ce-SBD:  {len(ce_tracks_cells & ce_sbd_names):>4} cells")

# Sample cell names from C. elegans datasets
print(f"\nC. elegans tracks sample cells: {sorted(ce_tracks_cells)[:30]}")
print(f"C. elegans SBD sample cells:    {sorted(ce_sbd_names)[:30]}")
