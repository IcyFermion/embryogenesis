"""Verify SBD 4th field interpretation by cross-referencing with CSV z."""
import pandas as pd
import numpy as np
from pathlib import Path

REPO = Path("/home/bingran/code/embryogenesis")

def parse_one_sbd_cell(filepath, target_cell):
    with open(filepath) as f:
        content = f.read()
    for block in content.split("---"):
        lines = [l.strip() for l in block.strip().split("\n") if l.strip()]
        if len(lines) < 5: continue
        p0 = lines[0].split()
        if len(p0) < 5: continue
        name = p0[-1]
        if name != target_cell: continue
        if p0[0] != "1": continue
        n_obs = int(lines[3].split()[0])
        rows = []
        for j in range(4, 4 + n_obs):
            op = lines[j].split()
            if len(op) >= 4:
                rows.append([int(op[0]), float(op[1]), float(op[2]), float(op[3])])
        arr = np.array(rows)
        birth_frame = int(lines[2].split()[0])
        start_time = int(lines[1].split()[0])
        return arr, birth_frame, start_time
    return None, None, None

# 1. SBD cell ABalpaaaa
arr, birth, start = parse_one_sbd_cell(
    REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd", "ABalpaaaa")
print("=== SBD (NM_C_briggsae_2) ABalpaaaa ===")
print(f"  Birth frame (line 3, field 1): {birth}")
print(f"  Start time (line 2, field 1): {start}")
print(f"  Frame range: {arr[:,0].min():.0f} - {arr[:,0].max():.0f}")
print(f"  X range: {arr[:,1].min():.0f} - {arr[:,1].max():.0f}")
print(f"  Y range: {arr[:,2].min():.0f} - {arr[:,2].max():.0f}")
print(f"  4th field range: {arr[:,3].min():.0f} - {arr[:,3].max():.0f}")
print(f"  4th field values: {arr[:,3].tolist()}")

# 2. Same cell in CSV
csv = pd.read_csv(REPO / "data/c_briggsae/CD140715HLH1cbp1.csv")
sub = csv[csv["cell"] == "ABalpaaaa"]
if len(sub) > 0:
    print(f"\n=== CSV ABalpaaaa ===")
    print(f"  Time range: {sub.time.min()} - {sub.time.max()}")
    print(f"  X range: {sub.x.min():.0f} - {sub.x.max():.0f}")
    print(f"  Y range: {sub.y.min():.0f} - {sub.y.max():.0f}")
    print(f"  Z range: {sub.z.min():.1f} - {sub.z.max():.1f}")
    print(f"  Z values: {sub.z.values.tolist()}")
else:
    print("\nABalpaaaa NOT found in CSV")

# 3. Same cell in C. elegans SBD
arr2, birth2, start2 = parse_one_sbd_cell(
    REPO / "data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd", "ABalpaaaa")
if arr2 is not None:
    print(f"\n=== C. elegans SBD ABalpaaaa ===")
    print(f"  Birth frame: {birth2}")
    print(f"  4th field values: {arr2[:,3].tolist()}")

# 4. 4th field stats across all SBD files
def all_4th_vals(filepath):
    with open(filepath) as f: content = f.read()
    vals = []
    for block in content.split("---"):
        lines = [l.strip() for l in block.strip().split("\n") if l.strip()]
        if len(lines) < 5: continue
        try: n_obs = int(lines[3].split()[0])
        except: continue
        for j in range(4, 4 + min(n_obs, len(lines)-4)):
            op = lines[j].split()
            if len(op) >= 4: vals.append(float(op[3]))
    return np.array(vals)

print("\n=== 4th field stats across all SBD files ===")
for fname in ["NM_C_briggsae_2.sbd", "NM_C_briggsae_3a.sbd",
              "IB+RS_N2_1_intestine.sbd"]:
    vals = all_4th_vals(REPO / "data/c_briggsae/nadin" / fname)
    print(f"  {fname}: min={vals.min():.1f}, max={vals.max():.1f}, "
          f"median={np.median(vals):.1f}, mean={vals.mean():.1f}")

# 5. Compare: CSV z vs SBD 4th field for a few shared cells
print("\n=== CSV z vs SBD 4th field for shared cells ===")
def get_one_sbd_cell(filepath, cell_name):
    arr, _, _ = parse_one_sbd_cell(filepath, cell_name)
    return arr

csv_cells = set(csv["cell"].unique())
for cell in ["AB", "ABa", "ABal", "ABala", "P1", "EMS"]:
    arr_sbd = get_one_sbd_cell(REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd", cell)
    sub_csv = csv[csv["cell"] == cell]
    if arr_sbd is not None and len(sub_csv) > 0:
        sbd_4th_mean = arr_sbd[:,3].mean()
        csv_z_mean = sub_csv["z"].mean()
        print(f"  {cell}: SBD 4th mean={sbd_4th_mean:.1f}, CSV z mean={csv_z_mean:.1f}")
