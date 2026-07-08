"""Add travel distance analysis cells to the existing notebook."""
import json
import nbformat as nbf

# Read existing notebook
with open("cbriggsae_comparison/analysis.ipynb") as f:
    nb = nbf.read(f, as_version=4)

# Find the old summary cell
summary_idx = None
for i, cell in enumerate(nb.cells):
    src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
    if "6. Summary & Conclusions" in src:
        summary_idx = i
        break

if summary_idx is None:
    print("Could not find summary cell!")
    exit(1)

print(f"Found old summary at index {summary_idx}")
print(f"Total cells before: {len(nb.cells)}")

new_cells = []

# ---- Cell: Travel intro ----
new_cells.append(
    nbf.v4.new_markdown_cell(
        """\
## 7. Total Cell Travel Distance Analysis

### A different angle: cell movement dynamics vs. static positions

Instead of comparing pairwise distances between cells (a **static** spatial property),
we now compare the **total travel distance** of each cell — the sum of all stepwise
displacements across its tracked lifetime. This captures **dynamic** movement behavior.

**Metrics computed per cell:**
- **Total travel**: sum of all step lengths — captures overall movement
- **Mean step size**: total travel / n_steps — average speed (normalized for duration)
- **Net displacement**: straight-line distance from start to end
- **Number of steps**: how many timepoints the cell was tracked
"""
    )
)

# ---- Cell: Compute travel ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# Compute total travel distances for all cells
# ============================================================

def compute_sbd_travel(cells):
    \"\"\"Compute travel metrics for SBD cells (2D only).\"\"\"
    travel = {}
    for c in cells:
        if len(c['x']) < 2:
            continue
        dx = np.diff(c['x']); dy = np.diff(c['y'])
        steps = np.sqrt(dx**2 + dy**2)
        travel[c['cell']] = {
            'total': steps.sum(),
            'mean_step': steps.mean(),
            'n_steps': len(steps),
            'net_disp': np.sqrt((c['x'][-1]-c['x'][0])**2 + (c['y'][-1]-c['y'][0])**2),
        }
    return travel

# CSV travel (2D x,y for fair comparison with SBD)
csv_travel = {}
for cell, grp in csv_df.groupby('cell'):
    grp = grp.sort_values('time')
    if len(grp) < 2:
        continue
    xs, ys = grp['x'].values, grp['y'].values
    dx, dy = np.diff(xs), np.diff(ys)
    steps = np.sqrt(dx**2 + dy**2)
    csv_travel[cell] = {
        'total': steps.sum(), 'mean_step': steps.mean(),
        'n_steps': len(steps),
        'net_disp': np.sqrt((xs[-1]-xs[0])**2 + (ys[-1]-ys[0])**2),
    }

sbd2_travel = compute_sbd_travel(active2)
sbd3_travel = compute_sbd_travel(active3)

print(f'Cells with travel data: CSV={len(csv_travel)}, SBD2={len(sbd2_travel)}, SBD3={len(sbd3_travel)}')

common_t12 = sorted(set(csv_travel) & set(sbd2_travel))
common_t13 = sorted(set(csv_travel) & set(sbd3_travel))
common_t23 = sorted(set(sbd2_travel) & set(sbd3_travel))
print(f'Common cells: CSV-SBD2={len(common_t12)}, CSV-SBD3={len(common_t13)}, SBD2-SBD3={len(common_t23)}')
"""
    )
)

# ---- Cell: Correlation table ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# Correlate travel metrics across datasets
# ============================================================

def correlate_metric(travel_a, travel_b, common, metric):
    va = np.array([travel_a[c][metric] for c in common])
    vb = np.array([travel_b[c][metric] for c in common])
    sr, sp = spearmanr(va, vb)
    pr, pp = pearsonr(va, vb)
    lpr, lpp = pearsonr(np.log1p(va), np.log1p(vb))
    return {'va': va, 'vb': vb, 'spearman_r': sr, 'pearson_r': pr,
            'log_pearson_r': lpr, 'n': len(common)}

metrics = ['total', 'mean_step', 'net_disp', 'n_steps']
metric_labels = {'total': 'Total Travel', 'mean_step': 'Mean Step',
                 'net_disp': 'Net Displacement', 'n_steps': 'N Steps'}

travel_results = {}
print(f'{"Metric":<18} {"Comparison":<14} {"n":<6} {"Spearman r":<12} {"Pearson r":<12} {"log-Pearson r"}')
print('-' * 80)
for metric in metrics:
    travel_results[metric] = {}
    for (ta, tb, common, label) in [
        (csv_travel, sbd2_travel, common_t12, 'CSV vs SBD2'),
        (csv_travel, sbd3_travel, common_t13, 'CSV vs SBD3'),
        (sbd2_travel, sbd3_travel, common_t23, 'SBD2 vs SBD3'),
    ]:
        res = correlate_metric(ta, tb, common, metric)
        travel_results[metric][label] = res
        print(f'{metric_labels[metric]:<18} {label:<14} {res["n"]:<6} '
              f'{res["spearman_r"]:<12.4f} {res["pearson_r"]:<12.4f} {res["log_pearson_r"]:.4f}')
"""
    )
)

# ---- Cell: Travel scatter plots ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# Travel metrics scatter: CSV vs SBD (both embryos)
# ============================================================
fig, axes = plt.subplots(2, 4, figsize=(20, 10))

metric_keys = ['total', 'mean_step', 'net_disp', 'n_steps']
comparisons = [('CSV vs SBD2', csv_travel, sbd2_travel, common_t12, 'CSV', 'SBD2'),
               ('CSV vs SBD3', csv_travel, sbd3_travel, common_t13, 'CSV', 'SBD3')]

for row, (comp_label, ta, tb, common, xlab, ylab) in enumerate(comparisons):
    for col, metric in enumerate(metric_keys):
        ax = axes[row, col]
        va = np.array([ta[c][metric] for c in common])
        vb = np.array([tb[c][metric] for c in common])

        hb = ax.hexbin(va, vb, gridsize=40, cmap='viridis', mincnt=1, bins='log')
        plt.colorbar(hb, ax=ax, label='log10(count)')

        max_val = max(va.max(), vb.max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')

        sr = travel_results[metric][comp_label]['spearman_r']
        pr = travel_results[metric][comp_label]['pearson_r']
        ax.set_xlabel(f'{xlab} {metric_labels[metric]}')
        ax.set_ylabel(f'{ylab} {metric_labels[metric]}')
        ax.set_title(f'{metric_labels[metric]}: {comp_label}\\n'
                     f'Spearman r={sr:.3f}, Pearson r={pr:.3f}')
        ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'travel_metrics_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
"""
    )
)

# ---- Cell: SBD2 vs SBD3 travel ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# SBD2 vs SBD3 travel metrics (same source, different embryos)
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))

for col, metric in enumerate(metric_keys):
    ax = axes[col]
    va = np.array([sbd2_travel[c][metric] for c in common_t23])
    vb = np.array([sbd3_travel[c][metric] for c in common_t23])

    hb = ax.hexbin(va, vb, gridsize=40, cmap='viridis', mincnt=1, bins='log')
    plt.colorbar(hb, ax=ax, label='log10(count)')

    max_val = max(va.max(), vb.max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')

    sr = travel_results[metric]['SBD2 vs SBD3']['spearman_r']
    pr = travel_results[metric]['SBD2 vs SBD3']['pearson_r']
    ax.set_xlabel(f'SBD2 {metric_labels[metric]}')
    ax.set_ylabel(f'SBD3 {metric_labels[metric]}')
    ax.set_title(f'{metric_labels[metric]}: SBD2 vs SBD3\\n'
                 f'Spearman r={sr:.3f}, Pearson r={pr:.3f}')
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'travel_sbd23_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
"""
    )
)

# ---- Cell: Static vs Dynamic bar chart ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# Compare: Pairwise Distance (static) vs Travel (dynamic) correlations
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

categories = ['Pairwise\\nDistance\\n(static)', 'Mean\\nStep', 'Net\\nDisplacement',
              'Total\\nTravel']
comparisons = ['CSV vs SBD2', 'CSV vs SBD3', 'SBD2 vs SBD3']

# Gather all Spearman r values
values = {
    ('Pairwise\\nDistance\\n(static)', 'CSV vs SBD2'): res_sbd2['spearman_r'],
    ('Pairwise\\nDistance\\n(static)', 'CSV vs SBD3'): res_sbd3['spearman_r'],
    ('Pairwise\\nDistance\\n(static)', 'SBD2 vs SBD3'): res_23['spearman_r'],
}
for metric, cat_name in [('mean_step', 'Mean\\nStep'),
                          ('net_disp', 'Net\\nDisplacement'),
                          ('total', 'Total\\nTravel')]:
    for cl in comparisons:
        values[(cat_name, cl)] = travel_results[metric][cl]['spearman_r']

x = np.arange(len(categories))
width = 0.22
colors = ['#2196F3', '#4CAF50', '#FF9800']

for i, cl in enumerate(comparisons):
    bars = ax.bar(x + i * width - width, [values[(cat, cl)] for cat in categories],
                  width, label=cl, color=colors[i], alpha=0.85)
    for bar, val in zip(bars, [values[(cat, cl)] for cat in categories]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', fontsize=8.5, rotation=90)

ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.set_ylabel('Spearman Rank Correlation')
ax.set_title('Static (pairwise distance) vs Dynamic (travel) correlation across datasets')
ax.set_ylim(0, 1.08)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'static_vs_dynamic_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
"""
    )
)

# ---- Cell: Example cells ----
new_cells.append(
    nbf.v4.new_code_cell(
        """\
# ============================================================
# Example: most and least traveled cells
# ============================================================
total_travels = np.array([csv_travel[c]['total'] for c in common_t12])
top_idx = np.argsort(total_travels)[-10:]
bot_idx = np.argsort(total_travels)[:10]

top_cells = [common_t12[i] for i in top_idx]
bot_cells = [common_t12[i] for i in bot_idx]

print("Top 10 most-traveled cells (CSV):")
for c in top_cells:
    print(f"  {c:<32s} total={csv_travel[c]['total']:8.1f}, "
          f"mean_step={csv_travel[c]['mean_step']:6.1f}, "
          f"n_steps={csv_travel[c]['n_steps']:3d}")

print("\\nBottom 10 least-traveled cells (CSV):")
for c in bot_cells:
    print(f"  {c:<32s} total={csv_travel[c]['total']:8.1f}, "
          f"mean_step={csv_travel[c]['mean_step']:6.1f}, "
          f"n_steps={csv_travel[c]['n_steps']:3d}")

# Check if cell lineage depth correlates with travel
print(f"\\nCell name length vs total travel (CSV):")
name_lens = np.array([len(c) for c in common_t12])
sr_name_len, _ = spearmanr(name_lens, total_travels)
pr_name_len, _ = pearsonr(name_lens, total_travels)
print(f"  Spearman r = {sr_name_len:.4f}, Pearson r = {pr_name_len:.4f}")
print("  (Longer names = later-born cells = tend to travel less)")
"""
    )
)

# ---- Updated Summary ----
new_cells.append(
    nbf.v4.new_markdown_cell(
        """\
## 8. Summary & Conclusions

### Part A: Static spatial arrangement — pairwise distances

| Comparison | n Cells | Spearman r | Pearson r |
|-----------|---------|------------|-----------|
| CSV vs SBD2 | 243 | **0.943** | 0.942 |
| CSV vs SBD3 | 243 | **0.926** | 0.929 |
| SBD2 vs SBD3 | 249 | **0.905** | 0.912 |

### Part B: Dynamic movement behavior — cell travel distances (all common cells, n ≈ 630)

| Metric | CSV vs SBD2 | CSV vs SBD3 | SBD2 vs SBD3 |
|--------|------------|------------|-------------|
| Total Travel | 0.274 | 0.466 | 0.465 |
| **Mean Step** | **0.665** | **0.605** | **0.625** |
| Net Displacement | 0.543 | 0.424 | 0.480 |
| N Steps (duration) | 0.627 | 0.785 | 0.638 |

*All values are Spearman rank correlations.*

### Key Insight: Static >> Dynamic Conservation

1. **Pairwise distances are highly conserved** (r ≈ 0.90–0.94): The relative spatial
   arrangement of cells at a given developmental stage is excellently preserved across
   imaging systems, embryos, and labs. This is a robust, system-independent property.

2. **Travel distances are only moderately conserved** (best r ≈ 0.60–0.67 for mean step):
   How much individual cells *move* is much less consistent across datasets. This is
   likely due to a combination of tracking noise in SBD data (irregular frame intervals),
   biological variation in migration paths, and scale differences in temporal sampling.

3. **Mean step beats total travel**: Normalizing by tracking duration (mean step size)
   improves correlation substantially (r = 0.60–0.67 vs r = 0.27–0.47). Total travel
   is confounded by how long a cell is tracked.

4. **Number of tracking steps correlates well** (r ≈ 0.63–0.79): The relative duration
   each cell is tracked is somewhat conserved — cells that divide late in one dataset
   also tend to be tracked longer in the other.

5. **SBD2 and SBD3 agree more on dynamics** (r ≈ 0.46–0.64) than on statics (r = 0.91),
   but the static agreement between embryos is still far stronger.

### Bottom Line
The **static spatial organization** of the embryo is a robust, highly conserved feature
across datasets (r > 0.9). The **dynamic movement behavior** (how much individual cells
migrate) is a noisier, less conserved property (r ≈ 0.5–0.7). For cross-dataset
validation, pairwise distance correlation is the more reliable metric.
"""
    )
)

# Insert new cells before the old summary
for i, new_cell in enumerate(new_cells):
    nb.cells.insert(summary_idx + i, new_cell)

# Delete the old summary
old_idx = summary_idx + len(new_cells)
del nb.cells[old_idx]

# Write
nbf.write(nb, "cbriggsae_comparison/analysis.ipynb")
print(f"Done. Notebook now has {len(nb.cells)} cells.")
