"""Ridge regression baseline for cross-species RNA cell type classification.

Trains an L2-regularised multinomial logistic regression on the same data, split,
and preprocessing used by the neural network encoder.  Reports held-out terminal
hard-label accuracy per species for direct comparison.

Usage:
    cd /home/bingran/code/embryogenesis
    conda run -n dev python3 expression_embedding/cross_species_rna/linear_baseline.py
"""

import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold

# Load the same data through the NN pipeline
_here = Path(__file__).resolve().parents[0]
_bundle = _here.parent
if str(_bundle) not in sys.path:
    sys.path.insert(0, str(_bundle))

from cross_species_rna.config import Config
from cross_species_rna.data_loader import load_cross_species_data


def main():
    config = Config()

    print("Loading data (same preprocessing + split as NN)...")
    data = load_cross_species_data(config)

    X_train = data["X_train"].numpy()
    y_train = data["y_train"].numpy()
    hm_train = data["hard_mask_train"]
    species_train = data["species_train"]

    X_val = data["X_val"].numpy()
    y_val = data["y_val"].numpy()
    hm_val = data["hard_mask_val"]
    species_val = data["species_val"]

    y_hard_train = y_train.argmax(axis=1)
    y_hard_val = y_val.argmax(axis=1)

    class_names = data["class_names"]

    # ---- 1. Cross-validation on training hard-labeled cells ----
    X_term_train = X_train[hm_train]
    y_term_train = y_hard_train[hm_train]

    print(f"\nCV data: {X_term_train.shape[0]} hard-labeled training cells, "
          f"{X_term_train.shape[1]} features, {len(class_names)} classes")

    # Per-species CV
    for sp, sp_label in [(0, "elegans"), (1, "briggsae")]:
        sp_mask = (species_train == sp) & hm_train
        X_sp = X_train[sp_mask]
        y_sp = y_hard_train[sp_mask]
        n_sp = X_sp.shape[0]
        n_splits = min(5, max(2, n_sp // 10))
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        fold_accs = []
        for tr_idx, val_idx in cv.split(X_sp, y_sp):
            m = LogisticRegression(C=1.0, max_iter=5000, tol=1e-4, random_state=42)
            m.fit(X_sp[tr_idx], y_sp[tr_idx])
            fold_accs.append(accuracy_score(y_sp[val_idx], m.predict(X_sp[val_idx])))
        print(f"  [{sp_label}] L2 Logistic CV: {np.mean(fold_accs):.4f} +/- {np.std(fold_accs):.4f}  ({n_splits}-fold)")

    # ---- 2. Joint CV (both species) ----
    n_splits_joint = min(5, max(2, X_term_train.shape[0] // 10))
    cv_joint = StratifiedKFold(n_splits=n_splits_joint, shuffle=True, random_state=42)
    fold_accs = []
    for tr_idx, val_idx in cv_joint.split(X_term_train, y_term_train):
        m = LogisticRegression(C=1.0, max_iter=5000, tol=1e-4, random_state=42)
        m.fit(X_term_train[tr_idx], y_term_train[tr_idx])
        fold_accs.append(accuracy_score(y_term_train[val_idx], m.predict(X_term_train[val_idx])))
    cv_mean = np.mean(fold_accs)
    cv_std = np.std(fold_accs)
    print(f"  [joint]  L2 Logistic CV: {cv_mean:.4f} +/- {cv_std:.4f}  ({n_splits_joint}-fold)")

    # ---- 3. Train on full training set, evaluate on held-out val ----
    model = LogisticRegression(C=1.0, max_iter=5000, tol=1e-4, random_state=42)
    model.fit(X_train[hm_train], y_hard_train[hm_train])

    print(f"\n{'='*60}")
    print(f"Linear (Ridge) vs NN -- Held-out Hard-Label Accuracy")
    print(f"{'='*60}")
    print(f"{'Species':12s} {'Linear':>10s} {'NN':>10s}  {'n_val':>6s}")
    print(f"{'-'*45}")

    nn_ele_acc = 0.9647
    nn_bri_acc = 0.9765

    for sp, sp_label, nn_acc in [(0, "elegans", nn_ele_acc), (1, "briggsae", nn_bri_acc)]:
        sp_mask = (species_val == sp) & hm_val
        n_sp = sp_mask.sum()
        if n_sp > 0:
            lin_acc = accuracy_score(y_hard_val[sp_mask], model.predict(X_val[sp_mask]))
            print(f"{sp_label:12s} {lin_acc:10.4f} {nn_acc:10.4f}  {n_sp:6d}")

    # Joint
    lin_joint = accuracy_score(y_hard_val[hm_val], model.predict(X_val[hm_val]))
    nn_joint = 0.9706
    print(f"{'joint':12s} {lin_joint:10.4f} {nn_joint:10.4f}  {hm_val.sum():6d}")

    # ---- 4. C-value sweep ----
    print(f"\n{'='*60}")
    print("Regularisation strength sweep (C values)")
    print(f"{'='*60}")
    for c_val in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        fold_accs_c = []
        for tr_idx, val_idx in cv_joint.split(X_term_train, y_term_train):
            m = LogisticRegression(C=c_val, max_iter=5000, tol=1e-4, random_state=42)
            m.fit(X_term_train[tr_idx], y_term_train[tr_idx])
            fold_accs_c.append(accuracy_score(y_term_train[val_idx], m.predict(X_term_train[val_idx])))
        print(f"  C={c_val:5.2f}  CV: {np.mean(fold_accs_c):.4f} +/- {np.std(fold_accs_c):.4f}")

    print(f"\nSummary: Linear (ridge) CV = {cv_mean:.4f} +/- {cv_std:.4f}")
    print(f"         NN best joint val acc   = {nn_joint:.4f}")
    print(f"         Linear held-out joint   = {lin_joint:.4f}")
    delta = lin_joint - nn_joint
    print(f"         Delta (linear - NN)     = {delta:+.4f}")


if __name__ == "__main__":
    main()
