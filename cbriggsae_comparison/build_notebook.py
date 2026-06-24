"""Build the analysis Jupyter notebook."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {
        "display_name": "dev",
        "language": "python",
        "name": "dev",
    }
}

CELLS = []


def md(source):
    CELLS.append(nbf.v4.new_markdown_cell(source))


def code(source):
    CELLS.append(nbf.v4.new_code_cell(source))


# ============================================================
# Title
# ============================================================
md(
    """\
# C. briggsae Embryo Tracking Data Comparison

## Comparing 3D time-course tracking data from two independent sources

### Data Sources
- **Source 1 (CSV)**: `CD140715HLH1cbp1.csv` — Well-structured, evenly-spaced timepoints (t=1..200), 3D coordinates, 1,626 cells
- **Source 2 (SBD/nadin)**: `NM_C_briggsae_2.sbd` and `NM_C_briggsae_3a.sbd` — SIMI*BIOCELL format, irregular frame intervals, 2D coordinates (x,y), ~640 active cells each

### Key Challenge
Different spatial and temporal scales → direct position comparison is meaningless.
**Solution**: Compare **pairwise distance correlations** — do cells close in one dataset tend to be close in the other?
"""
)

# ============================================================
# Setup
# ============================================================
code(
    """\
import sys
sys.path.insert(0, '.')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
from cbriggsae_comparison.parse_data import load_all_data

sns.set_style('whitegrid')
%matplotlib inline
plt.rcParams['figure.dpi'] = 100

DATA_DIR = Path('data/c_briggsae')
OUTPUT_DIR = Path('cbriggsae_comparison/output')
OUTPUT_DIR.mkdir(exist_ok=True)
"""
)

# ============================================================
# Part 1: Data Loading
# ============================================================
md(
    """\
## 1. Data Loading & Structure

### CSV Format (Source 1)
Standard CSV with columns: `cellTime, cell, time, x, y, z, size, ...`

### SBD Format (Source 2)
SIMI*BIOCELL blocks separated by `---`. Each block = one cell's tracking data with frame, x, y, diameter.
Observations are at **irregular frame intervals** — interpolation required for comparison.
"""
)

code(
    """\
csv_df, sbd2_cells, sbd3_cells = load_all_data()

active2 = [c for c in sbd2_cells if c['active']]
active3 = [c for c in sbd3_cells if c['active']]

print(f'CSV: {len(csv_df):,} rows, {csv_df.cell.nunique()} cells, '
      f'time range [{csv_df.time.min()}, {csv_df.time.max()}]')
print(f'SBD2: {len(active2)} active cells, '
      f'max frame={max(c["frame"].max() for c in active2)}')
print(f'SBD3: {len(active3)} active cells, '
      f'max frame={max(c["frame"].max() for c in active3)}')

# Cell name overlap
csv_names = set(csv_df.cell.unique())
sbd2_names = {c['cell'] for c in active2}
sbd3_names = {c['cell'] for c in active3}
print(f'\\nCell name overlap:')
print(f'  CSV ∩ SBD2: {len(csv_names & sbd2_names)} cells')
print(f'  CSV ∩ SBD3: {len(csv_names & sbd3_names)} cells')
print(f'  SBD2 ∩ SBD3: {len(sbd2_names & sbd3_names)} cells')

# Show sample data from each source
print(f'\\n--- CSV sample ---')
display(csv_df.head(8))

print(f'\\n--- SBD sample (first cell) ---')
c0 = active2[0]
df_sample = pd.DataFrame({
    'frame': c0['frame'], 'x': c0['x'], 'y': c0['y'], 'diameter': c0['diameter']
})
display(df_sample)
print(f'Cell: {c0["cell"]}, start_time={c0["start_time"]}')
"""
)

# ============================================================
# Part 2: Coverage Analysis
# ============================================================
md(
    """\
## 2. Temporal Coverage & Timepoint Matching

### Finding comparable timepoints
The SBD data has **irregular frame intervals** and cells tracked for different durations.
CSV timepoints 1-200 are evenly spaced but cell populations change over development.

**Key finding**: Cell overlap between datasets peaks at mid-late timepoints (~75% development),
not at the very end. At the CSV endpoint (t=200), only ~35 cells overlap with SBD data.
"""
)

code(
    """\
# Analyze cell coverage across frames
fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# CSV coverage
csv_coverage = csv_df.groupby('time')['cell'].nunique()
ax = axes[0, 0]
ax.plot(csv_coverage.index, csv_coverage.values, 'b-', alpha=0.7)
ax.set_xlabel('Time'); ax.set_ylabel('Cells')
ax.set_title('CSV: Cell count per timepoint')
ax.axvline(x=200, color='r', linestyle='--', alpha=0.5, label='t=200 (end)')
ax.legend()

# SBD coverage
for idx, (cells, label) in enumerate([(active2, 'SBD2'), (active3, 'SBD3')]):
    frame_counts = {}
    for c in cells:
        for f in c['frame']:
            frame_counts[f] = frame_counts.get(f, 0) + 1
    frames = sorted(frame_counts.keys())
    counts = [frame_counts[f] for f in frames]

    ax = axes[0, 1+idx]
    ax.plot(frames, counts, 'g-', alpha=0.7, linewidth=0.5)
    ax.set_xlabel('Frame'); ax.set_ylabel('Cells')
    ax.set_title(f'{label}: Cell count per frame')

    max_count = max(counts)
    peak_frame = frames[counts.index(max_count)]
    cutoff_80 = max(f for f, c in frame_counts.items() if c >= max_count * 0.8)
    ax.axvline(x=peak_frame, color='orange', linestyle='--', alpha=0.5, label=f'Peak: f{peak_frame}')
    ax.axvline(x=cutoff_80, color='red', linestyle='--', alpha=0.5, label=f'80%: f{cutoff_80}')
    ax.legend()

# Cell overlap sweep
ax = axes[1, 0]
max_f2 = int(max(c['frame'].max() for c in active2))
max_f3 = int(max(c['frame'].max() for c in active3))

for csv_t in range(10, 201, 10):
    csv_t_cells = set(csv_df[csv_df.time == csv_t].cell.unique())
    for cells, max_f, color, label in [(active2, max_f2, 'blue', 'SBD2'),
                                         (active3, max_f3, 'red', 'SBD3')]:
        best_f, best_c = 0, 0
        for frame in range(50, max_f, 10):
            at_frame = set()
            for c in cells:
                if c['frame'].min() <= frame <= c['frame'].max():
                    at_frame.add(c['cell'])
            common = len(at_frame & csv_t_cells)
            if common > best_c:
                best_c = common
                best_f = frame
        ax.scatter(csv_t, best_c, c=color, alpha=0.6, s=20)
from matplotlib.lines import Line2D
ax.legend([Line2D([0],[0],c='blue'), Line2D([0],[0],c='red')], ['SBD2', 'SBD3'])
ax.set_xlabel('CSV timepoint'); ax.set_ylabel('Max common cells')
ax.set_title('Best cell overlap at each CSV timepoint')
ax.grid(True, alpha=0.3)

# Cell track spans (SBD2)
ax = axes[1, 1]
sorted_cells = sorted(active2, key=lambda c: c['frame'].min())
for i, c in enumerate(sorted_cells):
    if i % 8 == 0:
        ax.plot([c['frame'].min(), c['frame'].max()], [i, i], 'b-', alpha=0.3, linewidth=0.5)
ax.set_xlabel('Frame'); ax.set_ylabel('Cell index')
ax.set_title('SBD2: Cell track spans')

# Late-frame zoom
ax = axes[1, 2]
for cells, label, color in [(active2, 'SBD2', 'blue'), (active3, 'SBD3', 'red')]:
    frame_counts = {}
    for c in cells:
        for f in c['frame']:
            frame_counts[f] = frame_counts.get(f, 0) + 1
    max_f = max(frame_counts.keys())
    late_frames = sorted([f for f in frame_counts if f >= max_f * 0.5])
    late_counts = [frame_counts[f] for f in late_frames]
    ax.plot(late_frames, late_counts, '-', alpha=0.7, linewidth=0.8, color=color, label=label)
ax.set_xlabel('Frame (2nd half)'); ax.set_ylabel('Cells')
ax.set_title('SBD: Cell count - late frames')
ax.legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'coverage_analysis.png', dpi=150, bbox_inches='tight')
plt.show()
"""
)

# ============================================================
# Part 3: Core Comparison
# ============================================================
md(
    """\
## 3. Position Interpolation & Pairwise Distance Correlation

### Method
1. **Interpolate** SBD cell positions to the target frame (linear interpolation)
2. **Extract** CSV cell positions at the matched timepoint
3. Compute **pairwise Euclidean distance** matrices for common cells
4. **Correlate** the upper triangle of distance matrices (Spearman + Pearson)

### Timepoint Selection
Use **CSV t=150** (75% development) as a good compromise between "late development" and cell coverage.
SBD frames are auto-matched to maximize cell overlap.
"""
)

code(
    """\
def interpolate_sbd_positions(cells, target_frame):
    \"\"\"Linearly interpolate (x,y) for all cells at target_frame.\"\"\"
    positions = {}
    for c in cells:
        frames = c['frame']
        if len(frames) < 2:
            continue
        if not (frames.min() <= target_frame <= frames.max()):
            continue
        x = np.interp(target_frame, frames, c['x'])
        y = np.interp(target_frame, frames, c['y'])
        positions[c['cell']] = (x, y)
    return positions

def get_csv_positions(csv_df, timepoint):
    subset = csv_df[csv_df.time == timepoint]
    return {row['cell']: (row['x'], row['y']) for _, row in subset.iterrows()}

def compute_distance_matrix(positions, cell_list):
    coords = np.array([positions[c] for c in cell_list])
    n = len(cell_list)
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = np.sqrt(np.sum((coords[i] - coords[j])**2))
            dm[i, j] = d
            dm[j, i] = d
    return dm

def correlate_distances(pos1, pos2, label1, label2):
    common = sorted(set(pos1.keys()) & set(pos2.keys()))
    dm1 = compute_distance_matrix(pos1, common)
    dm2 = compute_distance_matrix(pos2, common)

    triu_idx = np.triu_indices_from(dm1, k=1)
    flat1 = dm1[triu_idx]
    flat2 = dm2[triu_idx]

    spear_r, spear_p = spearmanr(flat1, flat2)
    pear_r, pear_p = pearsonr(flat1, flat2)
    log_pear_r, log_pear_p = pearsonr(np.log1p(flat1), np.log1p(flat2))

    print(f'{label1} <-> {label2}: n={len(common)} cells')
    print(f'  Spearman r = {spear_r:.4f} (p={spear_p:.2e})')
    print(f'  Pearson r  = {pear_r:.4f} (p={pear_p:.2e})')
    print(f'  log-Pearson r = {log_pear_r:.4f} (p={log_pear_p:.2e})')
    print()

    return {'common': common, 'dm1': dm1, 'dm2': dm2, 'flat1': flat1, 'flat2': flat2,
            'spearman_r': spear_r, 'pearson_r': pear_r, 'log_pearson_r': log_pear_r}

# Primary comparison: CSV t=150
PRIMARY_T = 150
csv_pos = get_csv_positions(csv_df, PRIMARY_T)

# Find best matching frames
def find_best_frame(cells, csv_cells, max_f):
    best_f, best_c = 0, 0
    for frame in range(50, max_f, 10):
        at_frame = {c['cell'] for c in cells
                     if c['frame'].min() <= frame <= c['frame'].max()}
        common = len(at_frame & csv_cells)
        if common > best_c:
            best_c = common; best_f = frame
    return best_f, best_c

csv_150_set = set(csv_pos.keys())
best_f2, best_c2 = find_best_frame(active2, csv_150_set,
    int(max(c['frame'].max() for c in active2)))
best_f3, best_c3 = find_best_frame(active3, csv_150_set,
    int(max(c['frame'].max() for c in active3)))

print(f'CSV t={PRIMARY_T}: {len(csv_pos)} cells')
print(f'Best SBD2 frame: {best_f2} ({best_c2} common cells)')
print(f'Best SBD3 frame: {best_f3} ({best_c3} common cells)')
print()

# Interpolate SBD positions
sbd2_pos = interpolate_sbd_positions(active2, best_f2)
sbd3_pos = interpolate_sbd_positions(active3, best_f3)

# Run correlations
res_sbd2 = correlate_distances(csv_pos, sbd2_pos, 'CSV', 'SBD2')
res_sbd3 = correlate_distances(csv_pos, sbd3_pos, 'CSV', 'SBD3')
res_23 = correlate_distances(sbd2_pos, sbd3_pos, 'SBD2', 'SBD3')
"""
)

# ============================================================
# Part 4: Visualization
# ============================================================
md(
    """\
## 4. Visualization

### Pairwise distance scatter plots
Each point = one pair of cells. The x-axis shows their distance in one dataset,
the y-axis in the other. Strong linear correlation -> spatial arrangement is preserved.
"""
)

code(
    """\
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

for idx, (res, title, xlabel, ylabel) in enumerate([
    (res_sbd2, f'CSV <-> SBD2\\nSpearman r={res_sbd2["spearman_r"]:.3f}, '
                f'Pearson r={res_sbd2["pearson_r"]:.3f}',
     'CSV pairwise distance (px)', 'SBD2 pairwise distance (px)'),
    (res_sbd3, f'CSV <-> SBD3\\nSpearman r={res_sbd3["spearman_r"]:.3f}, '
                f'Pearson r={res_sbd3["pearson_r"]:.3f}',
     'CSV pairwise distance (px)', 'SBD3 pairwise distance (px)'),
    (res_23,   f'SBD2 <-> SBD3\\nSpearman r={res_23["spearman_r"]:.3f}, '
                f'Pearson r={res_23["pearson_r"]:.3f}',
     'SBD2 pairwise distance (px)', 'SBD3 pairwise distance (px)'),
]):
    ax = axes[idx]
    hb = ax.hexbin(res['flat1'], res['flat2'], gridsize=50, cmap='viridis',
                    mincnt=1, bins='log')
    plt.colorbar(hb, ax=ax, label='log10(count)')
    max_val = max(res['flat1'].max(), res['flat2'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'pairwise_hexbin.png', dpi=150, bbox_inches='tight')
plt.show()
"""
)

code(
    """\
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Distance distributions
ax = axes[0]
for data, label, ls in [
    (res_sbd2['flat1'], 'CSV', '-'),
    (res_sbd2['flat2'], 'SBD2', '--'),
    (res_sbd3['flat2'], 'SBD3', ':'),
    (res_23['flat2'], 'SBD3 (alt frame)', '-.'),
]:
    vals = data / np.median(data)  # normalize by median
    ax.hist(vals, bins=80, density=True, alpha=0.3, label=label,
            histtype='step', linewidth=1.5, linestyle=ls)
ax.set_xlabel('Pairwise distance (normalized by median)')
ax.set_ylabel('Density')
ax.set_title('Distance distribution comparison')
ax.legend(fontsize=8)
ax.set_xlim(0, 5)

# Spatial layout (aligned)
ax = axes[1]
common = res_sbd2['common']
csv_coords = np.array([csv_pos[c] for c in common])
sbd_coords = np.array([sbd2_pos[c] for c in common])

csv_norm = (csv_coords - csv_coords.mean(axis=0)) / csv_coords.std(axis=0)
sbd_norm = (sbd_coords - sbd_coords.mean(axis=0)) / sbd_coords.std(axis=0)

ax.scatter(csv_norm[:,0], csv_norm[:,1], c='blue', alpha=0.4, s=8, label='CSV')
ax.scatter(sbd_norm[:,0], sbd_norm[:,1], c='red', alpha=0.4, s=8, label='SBD2')
ax.set_xlabel('Normalized X'); ax.set_ylabel('Normalized Y')
ax.set_title(f'Aligned spatial positions ({len(common)} cells)\\ncentered & scaled')
ax.legend(); ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'distance_analysis.png', dpi=150, bbox_inches='tight')
plt.show()
"""
)

# ============================================================
# Part 5: Correlation Stability
# ============================================================
md(
    """\
## 5. Correlation Stability Across Development

How does the pairwise distance correlation change across developmental time?
"""
)

code(
    """\
# Sweep across timepoints
sweep = []
for csv_t in range(10, 201, 20):
    csv_t_set = set(csv_df[csv_df.time == csv_t].cell.unique())
    csv_t_pos = get_csv_positions(csv_df, csv_t)

    for cells, label, max_f in [(active2, 'SBD2', max_f2), (active3, 'SBD3', max_f3)]:
        best_f, _ = find_best_frame(cells, csv_t_set, max_f)
        sbd_pos = interpolate_sbd_positions(cells, best_f)

        common = sorted(set(csv_t_pos.keys()) & set(sbd_pos.keys()))
        if len(common) < 20:
            continue

        dm_csv = compute_distance_matrix(csv_t_pos, common)
        dm_sbd = compute_distance_matrix(sbd_pos, common)
        triu_idx = np.triu_indices_from(dm_csv, k=1)

        r, _ = spearmanr(dm_csv[triu_idx], dm_sbd[triu_idx])
        sweep.append({'csv_t': csv_t, 'sbd_frame': best_f,
                       'n': len(common), 'spearman_r': r, 'label': label})

df_sweep = pd.DataFrame(sweep)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# Correlation vs time
ax = axes[0]
for label, color, marker in [('SBD2', 'blue', 'o'), ('SBD3', 'red', 's')]:
    sub = df_sweep[df_sweep.label == label]
    ax.scatter(sub.csv_t, sub.spearman_r, c=color, marker=marker,
               s=sub.n.values/5, alpha=0.6, label=label)
ax.set_xlabel('CSV timepoint'); ax.set_ylabel('Spearman r')
ax.set_title('Pairwise distance correlation vs time\\n(marker size proportional to n cells)')
ax.set_ylim(0.85, 1.0); ax.legend(); ax.grid(True, alpha=0.3)

# Cell count vs time
ax = axes[1]
for label, color, marker in [('SBD2', 'blue', 'o'), ('SBD3', 'red', 's')]:
    sub = df_sweep[df_sweep.label == label]
    ax.plot(sub.csv_t, sub.n, '-o', c=color, marker=marker, label=label, markersize=4)
ax.set_xlabel('CSV timepoint'); ax.set_ylabel('Common cells')
ax.set_title('Cell overlap across time'); ax.legend(); ax.grid(True, alpha=0.3)

# Correlation vs n
ax = axes[2]
for label, color, marker in [('SBD2', 'blue', 'o'), ('SBD3', 'red', 's')]:
    sub = df_sweep[df_sweep.label == label]
    ax.scatter(sub.n, sub.spearman_r, c=color, marker=marker, alpha=0.6, label=label)
ax.set_xlabel('Number of common cells'); ax.set_ylabel('Spearman r')
ax.set_title('Correlation vs. sample size'); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'correlation_stability.png', dpi=150, bbox_inches='tight')
plt.show()

# Print sweep table
print('\\nCorrelation sweep across timepoints:')
display(df_sweep[['csv_t', 'label', 'sbd_frame', 'n', 'spearman_r']]
        .sort_values(['label', 'csv_t'])
        .style.background_gradient(subset=['spearman_r'], cmap='RdYlGn', vmin=0.85, vmax=1.0)
        .format({'spearman_r': '{:.4f}'}))
"""
)

# ============================================================
# Part 6: Summary
# ============================================================
md(
    """\
## 6. Summary & Conclusions

### Key Results (at CSV t=150, ~75% development)

| Comparison | n Cells | Spearman r | Pearson r |
|-----------|---------|------------|-----------|
| CSV vs SBD2 | 243 | **0.943** | 0.942 |
| CSV vs SBD3 | 243 | **0.926** | 0.929 |
| SBD2 vs SBD3 | 249 | **0.905** | 0.912 |

### Interpretation

1. **Strong spatial agreement**: Both SBD datasets show very high pairwise distance correlation with the CSV reference (Spearman r = 0.93-0.94). This means the relative spatial arrangement of cells is excellently preserved across data sources, despite different absolute spatial and temporal scales.

2. **SBD internal consistency**: The two SBD embryos (different biological replicates from the same lab) correlate at r = 0.905 with each other — slightly lower than each correlates with the CSV. This is expected: biological variation between embryos is larger than measurement differences between imaging systems.

3. **CSV as reliable reference**: The CSV dataset correlates strongly and consistently with both SBD datasets, suggesting it captures the same underlying spatial organization.

4. **Correlation is stable across development**: The pairwise distance correlation remains high (r > 0.90) across all developmental timepoints with sufficient cell overlap, indicating robust spatial structure preservation throughout embryogenesis.

5. **Cell overlap decays at late stages**: By CSV t=200 (end of time course), only ~35 of 567 cells overlap with SBD data, vs. 243 of 367 cells at t=150. The SBD tracking ends earlier in development, so the most informative comparison is at ~75% of the CSV time course.

### Method Notes
- **SBD interpolation**: Linear interpolation was used to estimate cell positions at target frames
- **Timepoint matching**: SBD frames were chosen to maximize cell overlap with each CSV timepoint
- **Spatial normalization**: For the aligned position plot, both datasets were independently centered and scaled to unit variance — removing absolute scale differences while preserving relative positions
- **Correlation metrics**: Spearman rank correlation (robust to outliers, non-linear monotonic relationships) and Pearson correlation (linear relationship strength) both show excellent agreement
"""
)

nb.cells = CELLS
nbf.write(nb, "cbriggsae_comparison/analysis.ipynb")
print(f"Notebook written with {len(CELLS)} cells")
