# A-P Axis Determination for C. elegans Embryo Tracking Data

This subdirectory contains the analysis and utility code for determining the **anterior-posterior (A-P) axis** from 3D cell tracking microscopy data of *C. elegans* embryos.

## Key Finding

**The A-P axis is aligned with the microscopy X axis. Anterior is in the -X direction (decreasing X).**

This means:
- AB lineage cells (anterior) have **lower X** coordinates
- P lineage cells (posterior) have **higher X** coordinates
- The A-P position can be computed as `ap = -x` (flipping X so anterior is positive)

This holds across all three available embryos (embryo1, embryo2, embryo3) and is temporally stable across the entire pre-twitching period (t=0–255).

## Biological Background

In *C. elegans* embryogenesis, the A-P axis is established at the first cell division (t≈7 in this dataset):
- **AB** (larger, radius ~30μm): anterior blastomere, divides first
- **P1** (smaller, radius ~18μm): posterior blastomere, divides later

All AB descendants form the anterior of the embryo; P1 descendants (P lineage) form the posterior.

## Methods

Three methods are implemented, all converging on the same result:

### 1. PCA of Nuclear Positions
Principal component analysis of all nuclear (x, y, z) positions. PC1 (direction of maximum variance) captures the embryo's long axis, which is the A-P axis.

**Result:** PC1 aligns with X (component magnitude = 0.999, R² ≈ 0.64–0.95)

### 2. Convex Hull + Inertia Tensor (Insley & Shaham 2018)
Literature method: compute convex hull of all nuclei → fill interior uniformly → compute inertia tensor → eigenvector of smallest principal moment = long axis.

**Reference:** Insley & Shaham (2018) *PLoS ONE*, PMC5874040.

**Result:** Long axis aligns with X (component magnitude = 0.9998)

### 3. Lineage Centroid Separation (Biological Ground Truth)
Direct measurement: AB lineage centroid minus P lineage centroid defines the A-P vector.

**Result:** Δx ≈ -100 to -140 μm across timepoints, Δy < 30 μm, Δz < 20 μm

### Temporal Stability

| Method | Mean Angular Dev | Within 5° | Within 10° |
|--------|-----------------|-----------|------------|
| PCA    | 2.02° ± 0.93°   | 99.5%     | 100.0%     |
| Hull   | 1.17° ± 0.78°   | 100.0%    | 100.0%     |

## Usage

### Quick Check

```python
from ap_axis import quick_check

result = quick_check("data/embryo1/tracks.txt", time_cutoff=255)
# {'ap_axis': 'X', 'anterior_direction': '-X', 'consistent': True, ...}
```

### Full Analysis

```python
from ap_axis import compute_ap_axis, get_ap_position, transform_to_ap
import pandas as pd

df = pd.read_csv("data/embryo1/tracks.txt", sep="\t")

# Compute A-P axis
result = compute_ap_axis(df, method="pca", time_cutoff=255)
print(result)  # APAaxisResult(axis=X, anterior=-X, method=pca, ...)

# Get A-P position for all cells
df["ap"] = get_ap_position(df)

# Or transform the full DataFrame
df_ap = transform_to_ap(df)
```

### Coordinate Transform

To align the A-P axis with a standard coordinate system (anterior = +X):

| Original | Transformed | Meaning |
|----------|-------------|---------|
| x        | ap = -x     | A-P position (anterior positive) |
| y        | dv = y      | Dorsal-ventral (placeholder) |
| z        | lr = z      | Left-right (placeholder) |

Note: D-V and L-R axes have NOT been definitively determined from this analysis. Only the A-P axis is established with high confidence.

## Files

| File | Purpose |
|------|---------|
| `ap_axis.py` | Core utility module (importable) |
| `analyze_ap_axis.py` | Full analysis script with all methods |
| `visualize.py` | Generate diagnostic figures |
| `results/` | Output directory (JSON results + PNG figures) |

## Running

```bash
# Full analysis
python3 analyze_ap_axis.py

# Visualizations
python3 visualize.py

# Quick test
python3 -c "from ap_axis import quick_check; print(quick_check('../data/embryo1/tracks.txt'))"
```

## References

- Insley, P. & Shaham, S. (2018). "Automated C. elegans embryo alignments reveal brain neuropil position invariance despite lax cell body placement." *PLoS ONE*. [PMC5874040](https://pmc.ncbi.nlm.nih.gov/articles/PMC5874040/)
- Sulston, J.E., Schierenberg, E., White, J.G. & Thomson, J.N. (1983). "The embryonic cell lineage of the nematode Caenorhabditis elegans." *Developmental Biology*, 100(1), 64-119.
- Christensen, R.P. et al. (2015). "Untwisting the Caenorhabditis elegans embryo." *eLife*, 4, e10070.
