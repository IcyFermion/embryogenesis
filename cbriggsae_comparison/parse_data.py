"""
Parse C. briggsae embryo tracking data from two sources:
1. CD140715HLH1cbp1.csv - well-structured CSV with evenly-spaced timepoints
2. NM_C_briggsae_2.sbd / NM_C_briggsae_3a.sbd - SIMI*BIOCELL format with
   irregular frame intervals
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

DATA_DIR = Path(__file__).parent.parent / "data" / "c_briggsae"


def parse_csv_data(csv_path=None):
    """Parse the CSV tracking data (Source 1).

    Returns DataFrame with columns: cell, time, x, y, z
    """
    if csv_path is None:
        csv_path = DATA_DIR / "CD140715HLH1cbp1.csv"

    df = pd.read_csv(csv_path)
    # Keep only the core columns
    df = df[["cell", "time", "x", "y", "z"]].copy()
    return df


def parse_sbd_file(filepath):
    """Parse a single SIMI*BIOCELL .sbd file (Source 2).

    Each cell block is separated by '---'.
    Block structure:
      Line 0: "1 1 0 0 <cell_name>"  (active) or "0 0 -1 -1 <cell_name>" (terminal)
      Line 1: "<start_time> 3 -1 -1 [flag]"
      Line 2: "<cell_id> <something> -1 -1 -1 <flag> <cell_name>"
      Line 3: "<n_observations> [extra text]"
      Lines 4..4+n: "<frame> <x> <y> <diameter> -1 -1 -1"

    Returns list of dicts, one per active cell with tracking data.
    """
    with open(filepath) as f:
        content = f.read()

    blocks = content.split("---")
    cells = []

    for block in blocks:
        lines = [l.strip() for l in block.strip().split("\n") if l.strip()]
        if len(lines) < 5:
            continue

        # Parse header line
        parts0 = lines[0].split()
        if len(parts0) < 5:
            continue
        active = parts0[0] == "1"
        cell_name = parts0[-1]

        # Parse time line
        start_time = int(lines[1].split()[0])

        # Parse cell ID line
        cell_id = int(lines[2].split()[0])

        # Parse n_observations (may have trailing text like "excr cell")
        n_obs_str = lines[3].split()[0]
        try:
            n_obs = int(n_obs_str)
        except ValueError:
            continue

        if n_obs == 0:
            continue
        if 4 + n_obs > len(lines):
            n_obs = len(lines) - 4

        # Parse observation lines
        frames, xs, ys, diams = [], [], [], []
        for j in range(4, 4 + n_obs):
            obs_parts = lines[j].split()
            if len(obs_parts) < 4:
                continue
            frames.append(int(obs_parts[0]))
            xs.append(float(obs_parts[1]))
            ys.append(float(obs_parts[2]))
            diams.append(float(obs_parts[3]))

        cells.append({
            "cell": cell_name,
            "active": active,
            "start_time": start_time,
            "cell_id": cell_id,
            "n_obs": len(frames),
            "frame": np.array(frames),
            "x": np.array(xs),
            "y": np.array(ys),
            "diameter": np.array(diams),
        })

    return cells


def sbd_cells_to_dataframe(cells):
    """Convert parsed SBD cell list to a long-form DataFrame.

    Each row is one observation: cell, frame, x, y
    """
    rows = []
    for c in cells:
        for i in range(c["n_obs"]):
            rows.append({
                "cell": c["cell"],
                "frame": c["frame"][i],
                "x": c["x"][i],
                "y": c["y"][i],
            })
    return pd.DataFrame(rows)


def load_all_data():
    """Load all three datasets.

    Returns:
        csv_df: DataFrame for CSV data
        sbd2_cells: list of cell dicts for NM_C_briggsae_2.sbd
        sbd3_cells: list of cell dicts for NM_C_briggsae_3a.sbd
    """
    csv_df = parse_csv_data()

    sbd2_path = DATA_DIR / "nadin" / "NM_C_briggsae_2.sbd"
    sbd3_path = DATA_DIR / "nadin" / "NM_C_briggsae_3a.sbd"

    sbd2_cells = parse_sbd_file(sbd2_path)
    sbd3_cells = parse_sbd_file(sbd3_path)

    return csv_df, sbd2_cells, sbd3_cells


if __name__ == "__main__":
    csv_df, sbd2_cells, sbd3_cells = load_all_data()

    print(f"CSV: {len(csv_df)} rows, {csv_df.cell.nunique()} cells, "
          f"time range [{csv_df.time.min()}, {csv_df.time.max()}]")

    active2 = [c for c in sbd2_cells if c["active"]]
    active3 = [c for c in sbd3_cells if c["active"]]
    print(f"SBD2: {len(active2)} active cells, "
          f"frame range [{min(c['frame'].min() for c in active2)}, "
          f"{max(c['frame'].max() for c in active2)}]")
    print(f"SBD3: {len(active3)} active cells, "
          f"frame range [{min(c['frame'].min() for c in active3)}, "
          f"{max(c['frame'].max() for c in active3)}]")

    # Common cell names
    csv_cells = set(csv_df.cell.unique())
    sbd2_names = {c["cell"] for c in active2}
    sbd3_names = {c["cell"] for c in active3}
    print(f"\nCell name overlap CSV∩SBD2: {len(csv_cells & sbd2_names)}")
    print(f"Cell name overlap CSV∩SBD3: {len(csv_cells & sbd3_names)}")
    print(f"SBD2∩SBD3: {len(sbd2_names & sbd3_names)}")
