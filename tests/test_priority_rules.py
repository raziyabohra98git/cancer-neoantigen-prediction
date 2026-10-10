from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from priority_classification import classify_priority


def test_priority_distribution():

    df = pd.read_csv(
        ROOT
        / "example_results"
        / "mutant_wt_comparison.csv"
    )

    result = classify_priority(df)

    counts = (
        result["Priority"]
        .value_counts()
        .to_dict()
    )

    assert counts["Priority 1"] == 40
    assert counts["Priority 2"] == 80
    assert counts["Priority 3"] == 120
    assert counts["Below stringent threshold"] == 360


def test_p1_rules_are_all_satisfied():

    df = pd.read_csv(
        ROOT
        / "example_results"
        / "mutant_wt_comparison.csv"
    )

    result = classify_priority(df)

    p1 = result[
        result["Priority"] == "Priority 1"
    ]

    assert (
        p1["mhcflurry_affinity"] <= 50
    ).all()

    assert (
        p1["mhcflurry_affinity_percentile"] <= 0.5
    ).all()

    assert (
        p1["affinity_fold_improvement"] >= 2
    ).all()

    assert (
        p1["percentile_fold_improvement"] >= 2
    ).all()

    assert p1[
        "mutant_stronger_affinity"
    ].all()

    assert p1[
        "mutant_better_percentile"
    ].all()


def test_no_priority_mismatch_with_public_fixture():

    df = pd.read_csv(
        ROOT
        / "example_results"
        / "mutant_wt_comparison.csv"
    )

    result = classify_priority(df)

    assert (
        result["Priority"]
        == result["Synthetic_Target_Class"]
    ).all()


def test_peptide_lengths_are_project_compatible():

    df = pd.read_csv(
        ROOT
        / "example_data"
        / "synthetic_candidates.csv"
    )

    assert set(
        df["Peptide_Length"].unique()
    ) == {8, 9, 10, 11}

    assert (
        df["Mutant_Peptide"].str.len()
        == df["Peptide_Length"]
    ).all()

    assert (
        df["WT_Peptide"].str.len()
        == df["Peptide_Length"]
    ).all()
