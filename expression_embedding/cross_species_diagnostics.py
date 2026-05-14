"""
Cross-species RNA-seq diagnostic: C. elegans vs C. briggsae.

Assesses whether the raw mRNA data is already cross-species comparable, or whether
a strong species (domain) gap exists that an encoder would need to correct.

Outputs five diagnostic sections:
  1. PCA (species separation vs biological variation)
  2. UMAP (species coloring + terminal type coloring)
  3. Silhouette score by species
  4. Per-TF distribution comparison (mean / std scatter)
  5. Cross-species centroid distance table
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
ELEGANS_RNA_PATH = DATA_DIR / "c_briggsae" / "science.adu8249" / "c_elegans_tf.csv"
BRIGGSAE_RNA_PATH = DATA_DIR / "c_briggsae" / "science.adu8249" / "c_briggsae_tf.csv"
LINEAGE_PATH = DATA_DIR / "cell_lineage.json"
OUTPUT_DIR = Path(__file__).resolve().parent / "results" / "cross_species_diagnostics"

# Avoid loading heavy GUI backend on headless systems
import matplotlib as mpl

mpl.use("Agg")


# ==============================================================================
# Data loading & preprocessing
# ==============================================================================


def load_raw_data():
    """Load elegans and briggsae RNA CSVs, filter to shared TF genes.

    Returns
    -------
    ele : pd.DataFrame  (cells x genes)
    bri : pd.DataFrame  (cells x genes)
    shared_genes : list[str]
    """
    ele_raw = pd.read_csv(ELEGANS_RNA_PATH, index_col=0).T  # cells × genes
    bri_raw = pd.read_csv(BRIGGSAE_RNA_PATH, index_col=0).T  # cells × genes

    ele_genes = set(ele_raw.columns)
    bri_genes = set(bri_raw.columns)
    shared_genes = sorted(ele_genes & bri_genes)

    print(f"Elegans raw: {ele_raw.shape[0]} cells × {ele_raw.shape[1]} genes")
    print(f"Briggsae raw: {bri_raw.shape[0]} cells × {bri_raw.shape[1]} genes")
    print(f"Shared TF orthologs: {len(shared_genes)}")

    ele = ele_raw[shared_genes].copy()
    bri = bri_raw[shared_genes].copy()

    return ele, bri, shared_genes


def preprocess(ele, bri, cross_species_zscore=True):
    """Apply standard preprocessing and return combined matrix with species labels.

    Steps: log1p → per-cell L2-norm → z-score per gene (shared or per-species)

    Parameters
    ----------
    cross_species_zscore : bool
        If True, z-score across the combined matrix (what the encoder will see).
        If False, z-score each species separately (shows raw modality gap).

    Returns
    -------
    combined : np.ndarray  (n_total_cells, n_genes)
    species  : np.ndarray  (n_total_cells,)  — 0=elegans, 1=briggsae
    cell_names : np.ndarray (n_total_cells,) — original cell names
    gene_names : np.ndarray (n_genes,)
    """
    ele_vals = ele.values.astype(np.float64)
    bri_vals = bri.values.astype(np.float64)

    # Step 1: log1p
    ele_vals = np.log1p(ele_vals)
    bri_vals = np.log1p(bri_vals)

    # Step 2: per-cell L2 normalization
    ele_norms = np.linalg.norm(ele_vals, axis=1, keepdims=True)
    ele_norms[ele_norms == 0] = 1.0
    ele_vals = ele_vals / ele_norms

    bri_norms = np.linalg.norm(bri_vals, axis=1, keepdims=True)
    bri_norms[bri_norms == 0] = 1.0
    bri_vals = bri_vals / bri_norms

    # Step 3: z-score
    if cross_species_zscore:
        combined = np.vstack([ele_vals, bri_vals])
        scaler = StandardScaler()
        combined = scaler.fit_transform(combined)
    else:
        scaler_ele = StandardScaler()
        scaler_bri = StandardScaler()
        ele_scaled = scaler_ele.fit_transform(ele_vals)
        bri_scaled = scaler_bri.fit_transform(bri_vals)
        combined = np.vstack([ele_scaled, bri_scaled])

    n_ele = ele_vals.shape[0]
    n_bri = bri_vals.shape[0]
    species = np.array([0] * n_ele + [1] * n_bri, dtype=int)
    cell_names = np.array(list(ele.index) + list(bri.index))

    return combined, species, cell_names, ele.columns.values


# ==============================================================================
# Cell type assignment
# ==============================================================================


def _flatten_lineage(node, out=None):
    """Recursively flatten lineage tree keyed by did."""
    if out is None:
        out = {}
    did = node.get("did", "")
    if did:
        out[did] = node
    for child in node.get("children", []):
        _flatten_lineage(child, out)
    return out


def assign_cell_types(cell_names):
    """Assign a dominant terminal cell type to every RNA cell.

    Matches cells to lineage nodes via the ``did`` field.  Terminal nodes get
    their annotated type; internal nodes get a soft label from their descendant
    terminals (one-hot sum → normalised → argmax).

    Returns
    -------
    dominant_type : dict[str, str]   cell_name → type string (or 'unknown')
    type_list      : list[str]       sorted canonical type names
    """
    with open(LINEAGE_PATH) as f:
        root = json.load(f)

    all_nodes = _flatten_lineage(root)

    # Collect terminal types
    terminal_nodes = []
    for did, node in all_nodes.items():
        if not node.get("children"):
            terminal_nodes.append(did)

    # Descendant lists for all nodes
    descendant_list = defaultdict(list)

    def _collect_descendants(node, ancestors):
        did = node.get("did", "")
        if not did:
            return
        children = node.get("children", [])
        if not children:
            for anc in ancestors:
                descendant_list[anc].append(did)
        else:
            for child in children:
                _collect_descendants(child, ancestors + [did])

    _collect_descendants(root, [])

    # Terminal type map
    terminal_type = {}
    for did in terminal_nodes:
        node = all_nodes[did]
        t = (node.get("data") or {}).get("type", "")
        terminal_type[did] = t if t else "unknown"

    # Build one-hot encoding for terminal types
    unique_types = sorted(set(terminal_type.values()) - {"unknown"})
    unique_types.append("unknown")
    type_to_idx = {t: i for i, t in enumerate(unique_types)}

    terminal_one_hot = {}
    for did, t in terminal_type.items():
        vec = np.zeros(len(unique_types))
        vec[type_to_idx[t]] = 1.0
        terminal_one_hot[did] = vec

    # Soft labels for all nodes
    soft_labels = {}
    for did in all_nodes:
        desc = descendant_list.get(did, [])
        if not desc:
            soft_labels[did] = terminal_one_hot.get(did, np.zeros(len(unique_types)))
        else:
            summed = np.sum([terminal_one_hot.get(d, np.zeros(len(unique_types))) for d in desc], axis=0)
            total = summed.sum()
            soft_labels[did] = summed / total if total > 0 else summed

    # Map cell names to dominant type
    dominant_type = {}
    for name in cell_names:
        node = all_nodes.get(name)
        if node is None:
            dominant_type[name] = "unknown"
        else:
            sl = soft_labels.get(name)
            if sl is None or sl.sum() == 0:
                dominant_type[name] = "unknown"
            else:
                dominant_type[name] = unique_types[int(np.argmax(sl))]

    return dominant_type, unique_types


# ==============================================================================
# 1. PCA diagnostics
# ==============================================================================


def run_pca_diagnostics(combined, species, cell_names, suffix="", output_dir=OUTPUT_DIR):
    """PCA scatter colored by species; report per-species stats."""
    pca = PCA(n_components=min(50, combined.shape[1]))
    pca.fit(combined)
    coords = pca.transform(combined)

    pc1, pc2 = coords[:, 0], coords[:, 1]
    var1, var2 = pca.explained_variance_ratio_[:2]

    ele_mask = species == 0
    bri_mask = species == 1

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pc1[ele_mask], pc2[ele_mask], s=2, alpha=0.5, label="C. elegans", c="#1f77b4")
    ax.scatter(pc1[bri_mask], pc2[bri_mask], s=2, alpha=0.5, label="C. briggsae", c="#d62728")
    ax.set_xlabel(f"PC1 ({var1:.2%} var)")
    ax.set_ylabel(f"PC2 ({var2:.2%} var)")
    ax.set_title(f"PCA — Combined RNA-seq{suffix}")
    ax.legend(markerscale=5)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"pca_species{suffix.replace(' ', '_')}.png", dpi=150)
    plt.close(fig)

    # Per-species stats
    stats = {}
    for label, mask in [("elegans", ele_mask), ("briggsae", bri_mask)]:
        stats[label] = {
            "pc1_mean": float(np.mean(pc1[mask])),
            "pc1_std": float(np.std(pc1[mask])),
            "pc2_mean": float(np.mean(pc2[mask])),
            "pc2_std": float(np.std(pc2[mask])),
        }

    print(f"\n─── PCA Diagnostics{suffix} ───")
    print(f"PC1 variance explained: {var1:.4f}  ({var1:.2%})")
    print(f"PC2 variance explained: {var2:.4f}  ({var2:.2%})")
    for sp, s in stats.items():
        print(f"  {sp:10s}  PC1: {s['pc1_mean']:+.4f} ± {s['pc1_std']:.4f}   PC2: {s['pc2_mean']:+.4f} ± {s['pc2_std']:.4f}")

    return coords, stats


# ==============================================================================
# 2. UMAP diagnostics
# ==============================================================================


def run_umap_diagnostics(combined, species, cell_names, cell_types, suffix="", output_dir=OUTPUT_DIR):
    """Side-by-side UMAP colored by species and by terminal type."""
    import umap

    reducer = umap.UMAP(n_components=2, random_state=42, n_jobs=1, verbose=False)
    emb = reducer.fit_transform(combined)

    ele_mask = species == 0
    bri_mask = species == 1

    # Map dominant types to colors
    type_arr = np.array([cell_types.get(n, "unknown") for n in cell_names])
    all_types = sorted(set(type_arr))
    # Use a qualitative colormap with enough colors
    type_to_color = {}
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(all_types), 20)))
    for i, t in enumerate(all_types):
        type_to_color[t] = palette[i % len(palette)]
    type_colors = np.array([type_to_color[t] for t in type_arr])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(15, 6.5))

    # Species plot
    ax0.scatter(emb[ele_mask, 0], emb[ele_mask, 1], s=2, alpha=0.5, label="C. elegans", c="#1f77b4")
    ax0.scatter(emb[bri_mask, 0], emb[bri_mask, 1], s=2, alpha=0.5, label="C. briggsae", c="#d62728")
    ax0.set_title(f"UMAP by Species{suffix}")
    ax0.legend(markerscale=5)

    # Type plot
    for t in all_types:
        mask = type_arr == t
        ax1.scatter(emb[mask, 0], emb[mask, 1], s=2, alpha=0.5, label=t, c=[type_to_color[t]])
    ax1.set_title(f"UMAP by Terminal Type{suffix}")
    ax1.legend(markerscale=4, fontsize=5, ncol=2, loc="lower right")

    fig.tight_layout()
    fig.savefig(output_dir / f"umap{suffix.replace(' ', '_')}.png", dpi=150)
    plt.close(fig)

    print(f"  UMAP saved to umap{suffix.replace(' ', '_')}.png")

    return emb


# ==============================================================================
# 3. Silhouette score
# ==============================================================================


def compute_silhouette(combined, species, suffix="", output_dir=OUTPUT_DIR):
    """Silhouette score by species label on the combined PCA space."""
    pca = PCA(n_components=min(10, combined.shape[1]))
    coords = pca.fit_transform(combined)

    score = silhouette_score(coords, species)
    print(f"\n─── Silhouette Score (species){suffix} ───")
    print(f"  Silhouette: {score:.4f}")

    # Interpretation
    if score < 0.1:
        level = "adapter is optional"
    elif score < 0.3:
        level = "adapter recommended"
    elif score < 0.5:
        level = "adapter is essential"
    else:
        level = "stronger correction needed"
    print(f"  Interpretation: {level}")

    return score


# ==============================================================================
# 4. Per-TF distribution comparison
# ==============================================================================


def per_tf_distribution_comparison(ele_l2, bri_l2, suffix="", output_dir=OUTPUT_DIR):
    """Per-TF mean and std scatter across species.

    Parameters
    ----------
    ele_l2, bri_l2 : pd.DataFrame
        Already-preprocessed data (log1p + per-cell L2, without z-scoring).
    """
    ele_mean = ele_l2.values.mean(axis=0)
    bri_mean = bri_l2.values.mean(axis=0)
    ele_std = ele_l2.values.std(axis=0)
    bri_std = bri_l2.values.std(axis=0)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13, 6))

    # Mean scatter
    lim_mean = max(ele_mean.max(), bri_mean.max()) * 1.05
    ax0.scatter(ele_mean, bri_mean, s=8, alpha=0.7)
    ax0.plot([0, lim_mean], [0, lim_mean], "k--", linewidth=0.8)
    ax0.set_xlabel("C. elegans mean")
    ax0.set_ylabel("C. briggsae mean")
    ax0.set_title(f"Per-TF Mean Expression{suffix}")
    # R² annotation
    r2_mean = np.corrcoef(ele_mean, bri_mean)[0, 1] ** 2
    ax0.text(0.95, 0.05, f"R² = {r2_mean:.3f}", transform=ax0.transAxes, ha="right", va="bottom")

    # Std scatter
    lim_std = max(ele_std.max(), bri_std.max()) * 1.05
    ax1.scatter(ele_std, bri_std, s=8, alpha=0.7)
    ax1.plot([0, lim_std], [0, lim_std], "k--", linewidth=0.8)
    ax1.set_xlabel("C. elegans std")
    ax1.set_ylabel("C. briggsae std")
    ax1.set_title(f"Per-TF Std Expression{suffix}")
    r2_std = np.corrcoef(ele_std, bri_std)[0, 1] ** 2
    ax1.text(0.95, 0.05, f"R² = {r2_std:.3f}", transform=ax1.transAxes, ha="right", va="bottom")

    fig.tight_layout()
    fig.savefig(output_dir / f"per_tf_distribution{suffix.replace(' ', '_')}.png", dpi=150)
    plt.close(fig)

    # Quick summary
    pearson_mean = np.corrcoef(ele_mean, bri_mean)[0, 1]
    pearson_std = np.corrcoef(ele_std, bri_std)[0, 1]
    print(f"\n─── Per-TF Distribution{suffix} ───")
    print(f"  Mean Pearson r: {pearson_mean:.4f}  (R² = {r2_mean:.4f})")
    print(f"  Std  Pearson r: {pearson_std:.4f}  (R² = {r2_std:.4f})")

    return ele_mean, bri_mean, ele_std, bri_std


# ==============================================================================
# 5. Cross-species centroid distances
# ==============================================================================


def centroid_distance_analysis(combined, species, cell_names, cell_types, suffix="", output_dir=OUTPUT_DIR):
    """For each shared terminal type, compare same-type cross-species distance
    to different-type within-species distance."""
    ele_mask = species == 0
    bri_mask = species == 1

    type_arr = np.array([cell_types.get(n, "unknown") for n in cell_names])

    # Find types with enough cells in both species
    shared_types = []
    for t in sorted(set(type_arr) - {"unknown", "programmed_death"}):
        n_ele = ((type_arr == t) & ele_mask).sum()
        n_bri = ((type_arr == t) & bri_mask).sum()
        if n_ele >= 3 and n_bri >= 3:
            shared_types.append(t)

    if len(shared_types) < 2:
        print("  Not enough shared types for centroid analysis")
        return

    # Compute centroids
    ele_centroids = {}
    bri_centroids = {}
    for t in shared_types:
        ele_centroids[t] = combined[(type_arr == t) & ele_mask].mean(axis=0)
        bri_centroids[t] = combined[(type_arr == t) & bri_mask].mean(axis=0)

    rows = []
    for t in shared_types:
        # d_same: same type, different species
        d_same = float(np.linalg.norm(ele_centroids[t] - bri_centroids[t]))

        # d_diff_within_species: different type, same species
        d_ele_diff = []
        for t2 in shared_types:
            if t2 != t:
                d_ele_diff.append(float(np.linalg.norm(ele_centroids[t] - ele_centroids[t2])))
        d_ele_diff_mean = np.mean(d_ele_diff) if d_ele_diff else float("nan")

        d_bri_diff = []
        for t2 in shared_types:
            if t2 != t:
                d_bri_diff.append(float(np.linalg.norm(bri_centroids[t] - bri_centroids[t2])))
        d_bri_diff_mean = np.mean(d_bri_diff) if d_bri_diff else float("nan")

        rows.append(
            {
                "type": t,
                "d_same": d_same,
                "d_diff_ele": d_ele_diff_mean,
                "d_diff_bri": d_bri_diff_mean,
                "ratio_ele": d_same / d_ele_diff_mean if d_ele_diff_mean > 0 else float("nan"),
                "ratio_bri": d_same / d_bri_diff_mean if d_bri_diff_mean > 0 else float("nan"),
            }
        )

    df = pd.DataFrame(rows).set_index("type")
    print(f"\n─── Cross-Species Centroid Distances{suffix} ───")
    print(df.to_string(float_format=lambda x: f"{x:.4f}"))

    # Summary stat: fraction of types where d_same < d_diff_within_species
    better_ele = sum(1 for r in rows if r["ratio_ele"] < 1.0)
    better_bri = sum(1 for r in rows if r["ratio_bri"] < 1.0)
    n = len(rows)
    print(f"\n  Types where d_same < d_diff_within_species:")
    print(f"    vs elegans diffs:  {better_ele}/{n} ({better_ele/n:.0%})")
    print(f"    vs briggsae diffs: {better_bri}/{n} ({better_bri/n:.0%})")

    df.to_csv(output_dir / f"centroid_distances{suffix.replace(' ', '_')}.csv")
    return df


# ==============================================================================
# Main
# ==============================================================================


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- Load ----
    ele, bri, shared_genes = load_raw_data()
    all_cell_names = np.array(list(ele.index) + list(bri.index))

    # ---- Cell types (computed once from lineage, same for both species) ----
    cell_types, type_list = assign_cell_types(all_cell_names)
    print(f"\nCell types found: {len(type_list)}")
    print(f"Type distribution: { {t: sum(1 for v in cell_types.values() if v == t) for t in type_list[:8]} } ...")

    # ---- Run diagnostics for BOTH preprocessing variants ----
    for cross_zscore, label in [(False, " (without cross-species z-score)"), (True, " (with cross-species z-score)")]:
        print(f"\n{'='*60}")
        print(f"Variant:{label}")
        print(f"{'='*60}")

        combined, species, cell_names_arr, gene_names = preprocess(ele, bri, cross_species_zscore=cross_zscore)
        suffix = label

        # 1. PCA
        run_pca_diagnostics(combined, species, cell_names_arr, suffix=suffix, output_dir=OUTPUT_DIR)

        # 2. UMAP
        run_umap_diagnostics(combined, species, cell_names_arr, cell_types, suffix=suffix, output_dir=OUTPUT_DIR)

        # 3. Silhouette
        compute_silhouette(combined, species, suffix=suffix, output_dir=OUTPUT_DIR)

        # 4. Per-TF distribution (on shared genes, using pre-z-score data)
        # z-scoring removes mean/std differences, so use the pre-z-scored data for this
        # Re-prepare without z-score for distribution comparison
        ele_logl2, bri_logl2 = _prep_for_dist(ele, bri)
        per_tf_distribution_comparison(ele_logl2, bri_logl2, suffix=suffix, output_dir=OUTPUT_DIR)

        # 5. Centroid distances (on z-scored combined matrix as the encoder would see it)
        centroid_distance_analysis(combined, species, cell_names_arr, cell_types, suffix=suffix, output_dir=OUTPUT_DIR)

    print(f"\nAll diagnostics saved to {OUTPUT_DIR}")


def _prep_for_dist(ele, bri):
    """Helper: log1p + L2 normalize each species separately, return as DataFrames for
    per_tf_distribution_comparison. Z-scoring is skipped so we see raw distribution overlap."""
    ele_vals = ele.values.astype(np.float64)
    bri_vals = bri.values.astype(np.float64)

    ele_log = np.log1p(ele_vals)
    bri_log = np.log1p(bri_vals)

    ele_norms = np.linalg.norm(ele_log, axis=1, keepdims=True)
    ele_norms[ele_norms == 0] = 1.0
    bri_norms = np.linalg.norm(bri_log, axis=1, keepdims=True)
    bri_norms[bri_norms == 0] = 1.0

    ele_l2 = ele_log / ele_norms
    bri_l2 = bri_log / bri_norms

    return (
        pd.DataFrame(ele_l2, index=ele.index, columns=ele.columns),
        pd.DataFrame(bri_l2, index=bri.index, columns=bri.columns),
    )


if __name__ == "__main__":
    main()
