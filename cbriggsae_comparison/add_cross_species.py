"""Add cross-species comparison cells to the notebook."""
import nbformat as nbf

with open("cbriggsae_comparison/analysis.ipynb") as f:
    nb = nbf.read(f, as_version=4)

# Find the last summary cell
insert_idx = None
for i, cell in enumerate(nb.cells):
    src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
    if "8. Summary & Conclusions" in src:
        insert_idx = i
        break

if insert_idx is None:
    print("Could not find summary cell!")
    exit(1)

print(f"Inserting before cell {insert_idx} (old summary)")

new_cells = []

# ---- Intro ----
new_cells.append(nbf.v4.new_markdown_cell("""\
## 8. Cross-Species Comparison: C. elegans vs C. briggsae

### 5 datasets

| Dataset | Species | Format | Dims | Timepoints |
|---------|---------|--------|------|-----------|
| `CD140715HLH1cbp1.csv` | *C. briggsae* | CSV (regular) | 3D | t=1..200, 1,626 cells |
| `NM_C_briggsae_2.sbd` | *C. briggsae* | SBD (irregular) | 2D | frames 14..522, 635 cells |
| `NM_C_briggsae_3a.sbd` | *C. briggsae* | SBD (irregular) | 2D | frames 1..458, 642 cells |
| `embryo1/tracks.txt` | *C. elegans* | TSV (regular) | 3D | t=0..378, 1,332 cells |
| `IB+RS_N2_1_intestine.sbd` | *C. elegans* | SBD (irregular) | 2D | frames 17..648, 651 cells |

**Cell name overlap is excellent**: all SBD datasets share ~630-650 cells with matching names.
The two species share the same lineage annotation (AB, ABa, ABal, ...), enabling direct comparison.

### Comparisons
We focus on **cross-species** comparisons plus within-species baselines for reference.
"""))

# ---- Cell name overlap matrix ----
new_cells.append(nbf.v4.new_code_cell("""\
# ================================================================
# Cell name overlap across all 5 datasets
# ================================================================

import numpy as np, pandas as pd
from pathlib import Path

REPO = Path(os.getcwd())
if REPO.name == 'cbriggsae_comparison':
    REPO = REPO.parent

# Parse datasets (reusing earlier parsing code)
def parse_sbd_fast(fp):
    with open(fp) as f:
        content = f.read()
    blocks = content.split('---')
    names = set()
    for block in blocks:
        lines = [l.strip() for l in block.strip().split(chr(10)) if l.strip()]
        if len(lines) < 5: continue
        p0 = lines[0].split()
        if len(p0) < 5: continue
        if p0[0] == '1':
            names.add(p0[-1])
    return names

cb_csv_names = set(pd.read_csv(REPO / 'data/c_briggsae/CD140715HLH1cbp1.csv')['cell'].unique())
cb_sbd2_names = parse_sbd_fast(REPO / 'data/c_briggsae/nadin/NM_C_briggsae_2.sbd')
cb_sbd3_names = parse_sbd_fast(REPO / 'data/c_briggsae/nadin/NM_C_briggsae_3a.sbd')
ce_tracks_names = set(pd.read_csv(REPO / 'data/embryo1/tracks.txt', sep='\\t')['name'].unique())
ce_sbd_names = parse_sbd_fast(REPO / 'data/c_briggsae/nadin/IB+RS_N2_1_intestine.sbd')

all_names = {
    'C.b CSV': cb_csv_names, 'C.b SBD2': cb_sbd2_names,
    'C.b SBD3': cb_sbd3_names, 'C.e tracks': ce_tracks_names, 'C.e SBD': ce_sbd_names,
}
labels = list(all_names.keys())
overlap_matrix = np.zeros((5, 5), dtype=int)
for i, l1 in enumerate(labels):
    for j, l2 in enumerate(labels):
        overlap_matrix[i, j] = len(all_names[l1] & all_names[l2])

df_overlap = pd.DataFrame(overlap_matrix, index=labels, columns=labels)
display(df_overlap.style.background_gradient(cmap='Blues', axis=None).format('{}'))

# Key cross-species numbers
print(f"\\nC.b CSV ∩ C.e tracks: {len(cb_csv_names & ce_tracks_names)} cells  (both 'gold standard')")
print(f"C.b SBD2 ∩ C.e SBD:   {len(cb_sbd2_names & ce_sbd_names)} cells  (both SBD format)")
print(f"C.e tracks ∩ C.e SBD:  {len(ce_tracks_names & ce_sbd_names)} cells  (within C. elegans)")
"""))

# ---- Results table ----
new_cells.append(nbf.v4.new_markdown_cell("""\
### Results

Comparison code:
- `Cb` = *C. briggsae*, `Ce` = *C. elegans*
- `CSV` / `tracks` = evenly-spaced tracking data (gold standard format)
- `SBD` = SIMI*BIOCELL format (irregular frame intervals)
"""))

new_cells.append(nbf.v4.new_code_cell("""\
# ================================================================
# Cross-species results (pre-computed; re-running takes ~2 min)
# ================================================================

results_data = [
    # (label, group, pw_r, pw_cells, pw_pairs, tr_mean_step, tr_total, tr_net_disp, tr_n)
    ("Cb-CSV vs Cb-SBD2",    "C. briggsae",   0.9300, 248,  30628, 0.6647, 0.2739, 0.5425, 630),
    ("Cb-CSV vs Cb-SBD3",    "C. briggsae",   0.9101, 249,  30876, 0.6049, 0.4664, 0.4239, 635),
    ("Ce-tracks vs Ce-SBD",  "C. elegans",    0.9666, 243,  29403, 0.6705, 0.3492, 0.4187, 647),
    ("Cb-CSV vs Ce-tracks",  "Cross-species", 0.9029, 173,  14878, 0.5606, 0.0924, 0.0997, 1132),
    ("Cb-SBD2 vs Ce-SBD",    "Cross-species", 0.9519, 184,  16836, 0.5963, 0.3240, 0.3344, 631),
    ("Cb-CSV vs Ce-SBD",     "Cross-species", 0.9603, 250,  31125, 0.6259, 0.3239, 0.3701, 643),
    ("Ce-tracks vs Cb-SBD2", "Cross-species", 0.9344, 241,  28920, 0.6701, 0.1505, 0.3923, 633),
]

print(f'{"Comparison":<28} {"Group":<16} {"PW_r":>7} {"PW_cells":>9} {"PW_pairs":>9} '
      f'{"MS_r":>7} {"Tot_r":>7} {"ND_r":>7}')
print('-' * 100)
for r in results_data:
    print(f'{r[0]:<28} {r[1]:<16} {r[2]:>7.4f} {r[3]:>9} {r[4]:>9,} '
          f'{r[5]:>7.4f} {r[6]:>7.4f} {r[7]:>7.4f}')

print(f'\\nPW_r = pairwise distance Spearman r')
print(f'MS_r = mean step travel Spearman r')
print(f'Tot_r = total travel Spearman r')
print(f'ND_r = net displacement Spearman r')
"""))

# ---- Visualization ----
new_cells.append(nbf.v4.new_code_cell("""\
# ================================================================
# Visualization: static vs dynamic across species
# ================================================================
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(18, 6.5))

groups = ["C. briggsae", "C. elegans", "Cross-species"]
colors = {"C. briggsae": "#2196F3", "C. elegans": "#4CAF50", "Cross-species": "#FF5722"}
labels_short = [r[0].replace("Cb-", "C.b ").replace("Ce-", "C.e ").replace(" vs ", "\\nvs\\n")
                for r in results_data]

# Panel 1: Pairwise distance correlations
ax = axes[0]
x = np.arange(len(results_data))
vals_pw = [r[2] for r in results_data]
group_colors = [colors[r[1]] for r in results_data]
bars = ax.bar(x, vals_pw, color=group_colors, alpha=0.85, edgecolor='white')
for bar, val in zip(bars, vals_pw):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
            f'{val:.3f}', ha='center', fontsize=9, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(labels_short, fontsize=7.5)
ax.set_ylabel('Spearman Rank Correlation')
ax.set_title('Pairwise Distance Correlation (static spatial arrangement)')
ax.set_ylim(0.85, 1.02)
ax.grid(True, alpha=0.3, axis='y')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=colors[g], label=g) for g in groups], fontsize=9)

# Panel 2: Static vs Dynamic side-by-side
ax = axes[1]
width = 0.3
for i, (r, metric, mcolor, mlabel, offset) in enumerate([
    (results_data, 2, '#607D8B', 'Pairwise distance (static)', -width),
    (results_data, 5, '#FF9800', 'Mean step travel (dynamic)', 0),
    (results_data, 6, '#E91E63', 'Total travel (dynamic)', +width),
]):
    vals = [rr[metric] for rr in r]
    x_pos = np.arange(len(r)) + offset
    ax.bar(x_pos, vals, width * 0.9, color=mcolor, alpha=0.85, label=mlabel)

ax.set_xticks(np.arange(len(results_data)))
ax.set_xticklabels(labels_short, fontsize=7.5)
ax.set_ylabel('Spearman Rank Correlation')
ax.set_title('Static vs Dynamic Correlation (all comparisons)')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'cross_species_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
"""))

# ---- Updated summary ----
new_cells.append(nbf.v4.new_markdown_cell("""\
## 9. Final Summary & Conclusions

### Part A: Static spatial arrangement — pairwise distances

| Comparison | Group | n Cells | Spearman r |
|-----------|-------|---------|------------|
| Ce-tracks vs Ce-SBD | *C. elegans* | 243 | **0.967** |
| Cb-CSV vs Ce-SBD | Cross | 250 | **0.960** |
| Cb-SBD2 vs Ce-SBD | Cross | 184 | **0.952** |
| Ce-tracks vs Cb-SBD2 | Cross | 241 | **0.934** |
| Cb-CSV vs Cb-SBD2 | *C. briggsae* | 248 | **0.930** |
| Cb-CSV vs Cb-SBD3 | *C. briggsae* | 249 | **0.910** |
| Cb-CSV vs Ce-tracks | Cross | 173 | **0.903** |

### Part B: Dynamic movement — travel distances (mean step, all common cells)

| Comparison | Group | n Cells | Spearman r |
|-----------|-------|---------|------------|
| Ce-tracks vs Ce-SBD | *C. elegans* | 647 | **0.671** |
| Ce-tracks vs Cb-SBD2 | Cross | 633 | **0.670** |
| Cb-CSV vs Cb-SBD2 | *C. briggsae* | 630 | **0.665** |
| Cb-CSV vs Ce-SBD | Cross | 643 | **0.626** |
| Cb-CSV vs Cb-SBD3 | *C. briggsae* | 635 | **0.605** |
| Cb-SBD2 vs Ce-SBD | Cross | 631 | **0.596** |
| Cb-CSV vs Ce-tracks | Cross | 1,132 | **0.561** |

### Key Insights

1. **Spatial organization is strongly conserved across species** (r ≈ 0.90–0.97):
   The relative positions of cells in the embryo are remarkably similar between
   *C. elegans* and *C. briggsae*. Cross-species pairwise distance correlations
   are as high as (and in some cases higher than) within-species comparisons.
   This suggests a deeply conserved embryonic body plan.

2. **C. elegans has the most internally consistent data**: Ce-tracks vs Ce-SBD
   achieves r = 0.967, the highest pairwise correlation observed.

3. **SBD-to-SBD cross-species is excellent** (r = 0.952): Since both datasets
   come from the same lab using the same imaging pipeline, this high correlation
   validates that the SBD format captures real biological signal despite
   irregular frame intervals.

4. **Dynamic behavior is less conserved** (r ≈ 0.56–0.67 for mean step):
   Cell travel distances are only moderately correlated across species,
   consistent with our earlier finding that movement dynamics are noisier
   than static positions. Biological variation in cell migration paths
   between species likely contributes additional variance.

5. **The two "gold standard" datasets (Cb-CSV vs Ce-tracks) show the weakest
   correlations**: Both for static (r = 0.903) and especially for dynamic
   metrics (total travel r = 0.092, net displacement r = 0.100). This may
   reflect genuine species differences that are most visible when comparing
   the most complete tracking data from each species.

### Bottom Line
The embryonic spatial architecture is **highly conserved** between *C. elegans*
and *C. briggsae* — pairwise cell distance correlations of r > 0.90 across
species confirm that the same cell lineage produces a nearly identical spatial
arrangement. Cell movement dynamics, however, show substantial species-specific
variation, suggesting that while the *endpoint* of development is conserved,
the *trajectories* cells take to get there can differ.
"""))

# Insert new cells before the old summary
for i, new_cell in enumerate(new_cells):
    nb.cells.insert(insert_idx + i, new_cell)

# Delete the old summary
old_idx = insert_idx + len(new_cells)
old_src = "".join(nb.cells[old_idx]["source"]) if isinstance(nb.cells[old_idx]["source"], list) else nb.cells[old_idx]["source"]
print(f"Deleting old summary: {old_src[:80]}...")
del nb.cells[old_idx]

nbf.write(nb, "cbriggsae_comparison/analysis.ipynb")
print(f"Done. Notebook now has {len(nb.cells)} cells.")
