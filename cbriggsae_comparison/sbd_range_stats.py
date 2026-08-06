"""Quick stats: x, y, z range for all three SBD files."""
import numpy as np
from pathlib import Path

REPO = Path("/home/bingran/code/embryogenesis")

files = [
    ("C.b SBD2", REPO / "data/c_briggsae/nadin/NM_C_briggsae_2.sbd"),
    ("C.b SBD3", REPO / "data/c_briggsae/nadin/NM_C_briggsae_3a.sbd"),
    ("C.e SBD",  REPO / "data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd"),
]

def extract_coords(filepath):
    with open(filepath) as f:
        content = f.read()
    xs, ys, zs = [], [], []
    for block in content.split("---"):
        lines = [l.strip() for l in block.strip().split("\n") if l.strip()]
        if len(lines) < 5:
            continue
        p0 = lines[0].split()
        if len(p0) < 5 or p0[0] != "1":
            continue
        try:
            n_obs = int(lines[3].split()[0])
        except ValueError:
            continue
        for j in range(4, 4 + min(n_obs, len(lines) - 4)):
            op = lines[j].split()
            if len(op) >= 4:
                xs.append(float(op[1]))
                ys.append(float(op[2]))
                zs.append(float(op[3]))
    return np.array(xs), np.array(ys), np.array(zs)

print(f"{'File':<16} {'Dim':>3} {'Min':>8} {'Max':>8} {'Range':>8} {'Std':>8}")
print("-" * 58)
for label, fp in files:
    xs, ys, zs = extract_coords(fp)
    for dim, arr in [("x", xs), ("y", ys), ("z", zs)]:
        print(f"{label:<16} {dim:>3} {arr.min():>8.1f} {arr.max():>8.1f} "
              f"{arr.max()-arr.min():>8.1f} {arr.std():>8.1f}")
    print()
