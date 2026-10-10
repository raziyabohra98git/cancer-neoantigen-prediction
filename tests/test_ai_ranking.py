from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "example_results"

EXPECTED_FEATURES = {
    "mhcflurry_affinity",
    "mhcflurry_affinity_percentile",
    "WT_Binding_Affinity",
    "WT_Binding_Percentile",
    "affinity_fold_improvement",
    "percentile_fold_improvement",
    "Peptide_Length",
    "Mutation_Position_in_Peptide",
}

PRIORITY_WEIGHTS = {
    "Priority 1": 3,
    "Priority 2": 2,
    "Priority 3": 1,
    "Below stringent threshold": 0,
}


def load_ranked_results():
    path = RESULTS / "ai_ranked_candidates.csv"
    assert path.exists(), (
        "Run ./run_demo.sh before running the tests."
    )
    return pd.read_csv(path)


def test_ai_scores_are_complete_and_valid():
    df = load_ranked_results()

    assert len(df) == 600
    assert df["AI_score"].notna().all()
    assert np.isfinite(df["AI_score"]).all()
    assert df["AI_score"].between(0, 1).all()

    assert df["RF_score"].between(0, 1).all()
    assert df["XGB_score"].between(0, 1).all()
    assert df["LR_score"].between(0, 1).all()


def test_ensemble_score_uses_documented_weights():
    df = load_ranked_results()

    expected = (
        0.40 * df["RF_score"]
        + 0.30 * df["XGB_score"]
        + 0.30 * df["LR_score"]
    )

    np.testing.assert_allclose(
        df["AI_score"],
        expected,
        rtol=1e-10,
        atol=1e-10,
    )


def test_final_rank_score_uses_priority_weights():
    df = load_ranked_results()

    weights = df["Priority"].map(PRIORITY_WEIGHTS)

    assert weights.notna().all()

    expected = weights * 100 + df["AI_score"] * 100

    np.testing.assert_allclose(
        df["Final_rank_score"],
        expected,
        rtol=1e-10,
        atol=1e-10,
    )


def test_final_ranks_are_unique_and_ordered():
    df = load_ranked_results()

    expected_ranks = list(range(1, len(df) + 1))

    assert sorted(df["Final_rank"].tolist()) == expected_ranks

    assert df["Final_rank_score"].is_monotonic_decreasing

    priority_order = {
        "Priority 1": 1,
        "Priority 2": 2,
        "Priority 3": 3,
        "Below stringent threshold": 4,
    }

    observed_order = df["Priority"].map(priority_order)

    assert observed_order.notna().all()
    assert observed_order.is_monotonic_increasing


def test_feature_importance_is_valid():
    path = RESULTS / "feature_importance.csv"
    assert path.exists(), (
        "Run ./run_demo.sh before running the tests."
    )

    df = pd.read_csv(path)

    assert set(df["Feature"]) == EXPECTED_FEATURES
    assert len(df) == len(EXPECTED_FEATURES)
    assert df["RF_importance"].notna().all()
    assert np.isfinite(df["RF_importance"]).all()
    assert (df["RF_importance"] >= 0).all()

    assert np.isclose(
        df["RF_importance"].sum(),
        1.0,
        atol=1e-6,
    )
