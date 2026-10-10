#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier


FEATURES = [
    "mhcflurry_affinity",
    "mhcflurry_affinity_percentile",
    "WT_Binding_Affinity",
    "WT_Binding_Percentile",
    "affinity_fold_improvement",
    "percentile_fold_improvement",
    "Peptide_Length",
    "Mutation_Position_in_Peptide",
]


PRIORITY_WEIGHTS = {
    "Priority 1": 3,
    "Priority 2": 2,
    "Priority 3": 1,
    "Below stringent threshold": 0,
}


def main() -> int:

    root = Path(__file__).resolve().parents[1]

    input_file = (
        root
        / "example_results"
        / "prioritized_candidates.csv"
    )

    output_file = (
        root
        / "example_results"
        / "ai_ranked_candidates.csv"
    )

    feature_file = (
        root
        / "example_results"
        / "feature_importance.csv"
    )

    patient_file = (
        root
        / "example_results"
        / "patient_ai_summary.csv"
    )

    hla_file = (
        root
        / "example_results"
        / "hla_ai_summary.csv"
    )

    gene_file = (
        root
        / "example_results"
        / "gene_ai_summary.csv"
    )

    qc_file = (
        root
        / "example_results"
        / "ai_ml_qc.txt"
    )

    if not input_file.exists():
        print(f"ERROR: input file not found: {input_file}")
        return 1

    df = pd.read_csv(input_file)

    required = FEATURES + [
        "Priority",
        "Priority1",
        "Candidate_ID",
        "Sample_ID",
        "Hugo_Symbol",
        "allele",
    ]

    missing = [
        col
        for col in required
        if col not in df.columns
    ]

    if missing:
        print("ERROR: missing required columns:")
        for col in missing:
            print(f"  - {col}")
        return 1

    # --------------------------------------------------------------
    # Prepare feature matrix.
    # --------------------------------------------------------------

    X = df[FEATURES].copy()

    for col in FEATURES:
        X[col] = pd.to_numeric(
            X[col],
            errors="coerce",
        )

    X = X.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    if X.isna().any().any():
        print(
            "ERROR: missing/non-finite ML feature values detected."
        )
        return 1

    y = df["Priority1"].astype(int)

    if y.nunique() < 2:
        print(
            "ERROR: supervised ML requires at least two classes."
        )
        return 1

    # --------------------------------------------------------------
    # Random Forest
    #
    # Exact PI configuration.
    # --------------------------------------------------------------

    rf = RandomForestClassifier(
        n_estimators=500,
        random_state=42,
        class_weight="balanced",
        min_samples_leaf=2,
        n_jobs=-1,
    )

    rf.fit(X, y)

    df["RF_score"] = rf.predict_proba(X)[:, 1]

    feature_importance = pd.DataFrame(
        {
            "Feature": FEATURES,
            "RF_importance": rf.feature_importances_,
        }
    ).sort_values(
        "RF_importance",
        ascending=False,
    ).reset_index(drop=True)

    # --------------------------------------------------------------
    # Logistic Regression
    #
    # Exact PI configuration.
    # --------------------------------------------------------------

    lr = Pipeline(
        [
            (
                "scaler",
                StandardScaler(),
            ),
            (
                "classifier",
                LogisticRegression(
                    max_iter=3000,
                    class_weight="balanced",
                ),
            ),
        ]
    )

    lr.fit(X, y)

    df["LR_score"] = lr.predict_proba(X)[:, 1]

    # --------------------------------------------------------------
    # XGBoost
    #
    # Exact PI configuration.
    # --------------------------------------------------------------

    xgb = XGBClassifier(
        n_estimators=300,
        max_depth=3,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        eval_metric="logloss",
        random_state=42,
    )

    xgb.fit(X, y)

    df["XGB_score"] = xgb.predict_proba(X)[:, 1].astype(np.float64)

    # --------------------------------------------------------------
    # Ensemble AI score.
    # --------------------------------------------------------------

    df["AI_score"] = (
        0.40 * df["RF_score"]
        + 0.30 * df["XGB_score"]
        + 0.30 * df["LR_score"]
    )

    # --------------------------------------------------------------
    # Biological priority weighting.
    # --------------------------------------------------------------

    df["Priority_weight"] = (
        df["Priority"]
        .map(PRIORITY_WEIGHTS)
    )

    if df["Priority_weight"].isna().any():
        print(
            "ERROR: unexpected priority label detected."
        )
        return 1

    df["Final_rank_score"] = (
        df["Priority_weight"] * 100
        + df["AI_score"] * 100
    )

    # Higher final score = higher priority.
    # AI score breaks ties within biological priority classes.
    df = df.sort_values(
        "Final_rank_score",
        ascending=False,
    ).reset_index(drop=True)

    df["Final_rank"] = (
        np.arange(len(df)) + 1
    )

    # --------------------------------------------------------------
    # Summary tables.
    # --------------------------------------------------------------

    patient_summary = (
        df.groupby("Sample_ID")
        .agg(
            total_candidates=(
                "Candidate_ID",
                "size",
            ),
            priority1_candidates=(
                "Priority1",
                "sum",
            ),
            mean_mutant_affinity=(
                "mhcflurry_affinity",
                "mean",
            ),
            mean_mutant_percentile=(
                "mhcflurry_affinity_percentile",
                "mean",
            ),
            maximum_AI_score=(
                "AI_score",
                "max",
            ),
        )
        .reset_index()
    )

    hla_summary = (
        df.groupby("allele")
        .agg(
            candidates=(
                "Candidate_ID",
                "size",
            ),
            priority1=(
                "Priority1",
                "sum",
            ),
            mean_affinity=(
                "mhcflurry_affinity",
                "mean",
            ),
            mean_percentile=(
                "mhcflurry_affinity_percentile",
                "mean",
            ),
        )
        .reset_index()
    )

    gene_summary = (
        df.groupby("Hugo_Symbol")
        .agg(
            candidates=(
                "Candidate_ID",
                "size",
            ),
            priority1=(
                "Priority1",
                "sum",
            ),
            best_AI_score=(
                "AI_score",
                "max",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "priority1",
                "best_AI_score",
            ],
            ascending=False,
        )
        .reset_index(drop=True)
    )

    # --------------------------------------------------------------
    # Save outputs.
    # --------------------------------------------------------------

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df.to_csv(
        output_file,
        index=False,
    )

    feature_importance.to_csv(
        feature_file,
        index=False,
    )

    patient_summary.to_csv(
        patient_file,
        index=False,
    )

    hla_summary.to_csv(
        hla_file,
        index=False,
    )

    gene_summary.to_csv(
        gene_file,
        index=False,
    )

    # --------------------------------------------------------------
    # QC.
    # --------------------------------------------------------------

    if len(df) != len(X):
        print(
            "ERROR: ML ranking changed row count."
        )
        return 1

    if df["AI_score"].isna().any():
        print(
            "ERROR: missing AI scores detected."
        )
        return 1

    if not df["AI_score"].between(
        0,
        1,
    ).all():
        print(
            "ERROR: AI score outside [0,1]."
        )
        return 1

    if not np.isfinite(
        df["Final_rank_score"]
    ).all():
        print(
            "ERROR: non-finite final rank score."
        )
        return 1

    if df["Final_rank"].duplicated().any():
        print(
            "ERROR: duplicate final ranks."
        )
        return 1

    with open(qc_file, "w") as handle:

        handle.write(
            "AI/ML RANKING QC\n"
        )

        handle.write(
            "=" * 72 + "\n"
        )

        handle.write(
            f"Input rows: {len(X)}\n"
        )

        handle.write(
            f"Features: {len(FEATURES)}\n"
        )

        handle.write(
            f"Priority 1: {int(y.sum())}\n"
        )

        handle.write(
            f"Non-Priority 1: "
            f"{int(len(y) - y.sum())}\n\n"
        )

        handle.write(
            "Features used:\n"
        )

        for feature in FEATURES:
            handle.write(
                f"  {feature}\n"
            )

        handle.write(
            "\nModels:\n"
        )

        handle.write(
            "  Random Forest: "
            "500 trees, balanced, min_leaf=2, "
            "random_state=42\n"
        )

        handle.write(
            "  Logistic Regression: "
            "StandardScaler, max_iter=3000, "
            "balanced\n"
        )

        handle.write(
            "  XGBoost: 300 trees, max_depth=3, "
            "learning_rate=0.05, subsample=0.8, "
            "colsample_bytree=0.8, random_state=42\n"
        )

        handle.write(
            "\nEnsemble:\n"
        )

        handle.write(
            "  AI_score = "
            "0.40 RF + 0.30 XGB + 0.30 LR\n"
        )

        handle.write(
            "\nFinal ranking:\n"
        )

        handle.write(
            "  Final_rank_score = "
            "Priority_weight*100 + AI_score*100\n"
        )

        handle.write(
            "\nQC checks: PASSED\n"
        )

    print("=" * 72)
    print("AI/ML-ASSISTED NEOANTIGEN RANKING")
    print("=" * 72)

    print(
        f"Input candidates:            {len(X)}"
    )

    print(
        f"Features used:               {len(FEATURES)}"
    )

    print(
        f"Priority 1 candidates:       {int(y.sum())}"
    )

    print(
        f"Non-Priority 1 candidates:   "
        f"{int(len(y) - y.sum())}"
    )

    print()
    print("Models:")
    print("  Random Forest")
    print("  Logistic Regression")
    print("  XGBoost")

    print()
    print(
        "AI score range: "
        f"{df['AI_score'].min():.6f} - "
        f"{df['AI_score'].max():.6f}"
    )

    print(
        "Missing AI scores: "
        f"{int(df['AI_score'].isna().sum())}"
    )

    print()
    print("Top RF features:")

    for _, row in feature_importance.head(8).iterrows():
        print(
            f"  {row['Feature']:35s} "
            f"{row['RF_importance']:.6f}"
        )

    print()
    print(f"Ranked output: {output_file}")
    print(f"Feature importance: {feature_file}")
    print(f"Patient summary: {patient_file}")
    print(f"HLA summary: {hla_file}")
    print(f"Gene summary: {gene_file}")
    print(f"QC log: {qc_file}")

    print()
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
