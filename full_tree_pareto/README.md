# Full-tree Pareto analysis

This module contains the experiment runners and exploratory notebook for Pareto analysis of the complete embryonic lineage tree. The analysis balances travel distance against cell-state distance while reconstructing assignments throughout the lineage, rather than only among terminal cells.

The optimization implementation remains in `../pareto_core.py`. Module entry points resolve source data from the repository root and store expensive numerical results under `output/internal_opt/` inside this directory.

## Entry points

- `run_pipeline_multi.py` is the primary parameterized runner. Pass one of `embryo1_prot`, `embryo2_prot`, `embryo3_prot`, `elegans_rna`, or `briggsae_rna`.
- `run_pipeline.py` reproduces the original single-configuration protein experiment and compares it with legacy cached results.
- `full_tree_pareto.ipynb` loads cached results for exploration, comparison, and figure development.

From the repository root, for example:

```bash
python full_tree_pareto/run_pipeline_multi.py embryo1_prot
```

## Cached results

Null-model and heuristic results are cached under:

```text
full_tree_pareto/output/internal_opt/<data-config>/<embryo>/
├── null/
└── heuristics/
```

The full-tree optimization is computationally intensive. Runners load an existing `.npz` cache when available and recompute a result only when its expected cache file is absent. Keep the cache directory with the module when relocating or archiving experiments. These generated files are intentionally excluded from Git by the repository-level `output/` ignore rule.
