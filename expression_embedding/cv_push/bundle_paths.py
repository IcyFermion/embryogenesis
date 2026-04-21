from pathlib import Path
import sys

BUNDLE_DIR = Path(__file__).resolve().parent
EXPERIMENTS_DIR = BUNDLE_DIR.parent
ROOT_DIR = EXPERIMENTS_DIR.parent
RESULTS_DIR = EXPERIMENTS_DIR / "results"
DATA_DIR = ROOT_DIR / "data"
PROTEIN_DATA_DIR = DATA_DIR / "protein" / "aggregated_all"

CELL_LINEAGE_PATH = DATA_DIR / "cell_lineage.json"
CELL_TYPE_CSV = DATA_DIR / "2023-06-29_entropy_cell_key_V2.csv"
EMBRYO1_TRACKS = DATA_DIR / "embryo1" / "tracks.txt"
VISCELLO_EXPR_RDS = DATA_DIR / "viscello" / "lin_sc_expr_190602.rds"

S3_CSV = PROTEIN_DATA_DIR / "s3.csv"
S3_ZSCORE_CSV = RESULTS_DIR / "s3_zscore.csv"
S3_PCA_10D_CSV = PROTEIN_DATA_DIR / "s3_pca_10d.csv"
LINEAGE_BINARY_EXPRESSION_CSV = PROTEIN_DATA_DIR / "lineage_binary_expression.csv"

AR_EMBEDDINGS_CSV = RESULTS_DIR / "ar_embeddings_32d.csv"
AR_FEATURE_IMPORTANCE_TSV = RESULTS_DIR / "ar_feature_importance.tsv"
AR_MODEL_CHECKPOINT_PT = RESULTS_DIR / "ar_model_checkpoint.pt"

PHASE2_ALTERNATIVE_RESULTS_CSV = RESULTS_DIR / "phase2_alternative_results.csv"
PHASE2_PUSH_RESULTS_CSV = RESULTS_DIR / "phase2_push_results.csv"
PHASE2_CV_SUMMARY_CSV = RESULTS_DIR / "phase2_cv_summary.csv"
PHASE2_OPTIMIZATION_SUMMARY_CSV = RESULTS_DIR / "phase2_optimization_summary.csv"
PHASE2_STAGE2_PLOT = RESULTS_DIR / "phase2_stage2_results.png"

STAGE1_RESULTS_CSV = RESULTS_DIR / "stage1_results.csv"
STAGE1_RESULTS_PKL = RESULTS_DIR / "stage1_results.pkl"
STAGE2_RESULTS_CSV = RESULTS_DIR / "stage2_results.csv"
STAGE2_RESULTS_PKL = RESULTS_DIR / "stage2_results.pkl"

FEATURE_SELECT_EMBEDDINGS_CSV = RESULTS_DIR / "embeddings_32d.csv"
FEATURE_SELECT_SELECTED_PROTEINS_TSV = RESULTS_DIR / "nn_selected_proteins_rev.tsv"

print(BUNDLE_DIR, EXPERIMENTS_DIR, ROOT_DIR)
def ensure_repo_on_path() -> None:
    root_str = str(ROOT_DIR)
    if root_str not in sys.path:
        sys.path.insert(0, root_str)


def ensure_results_dir() -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    return RESULTS_DIR
