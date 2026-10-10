#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


def main() -> int:
    root = Path(__file__).resolve().parents[1]

    candidate_file = (
        root
        / "example_results"
        / "primary_binding_candidates.csv"
    )

    wt_file = (
        root
        / "example_data"
        / "synthetic_wt_predictions.csv"
    )

    output_file = (
        root
        / "example_results"
        / "wt_aware_candidates.csv"
    )

    if not candidate_file.exists():
        print(f"ERROR: missing candidate file: {candidate_file}")
        return 1

    if not wt_file.exists():
        print(f"ERROR: missing WT prediction file: {wt_file}")
        return 1

    candidates = pd.read_csv(candidate_file)
    wt = pd.read_csv(wt_file)

    candidate_required = [
        "Candidate_ID",
        "WT_Peptide",
        "allele",
    ]

    wt_required = [
        "WT_Peptide",
        "allele",
        "WT_Binding_Affinity",
        "WT_Binding_Percentile",
    ]

    missing_candidates = [
        col for col in candidate_required
        if col not in candidates.columns
    ]

    missing_wt = [
        col for col in wt_required
        if col not in wt.columns
    ]

    if missing_candidates:
        print("ERROR: candidate table is missing:")
        for col in missing_candidates:
            print(f"  - {col}")
        return 1

    if missing_wt:
        print("ERROR: WT table is missing:")
        for col in missing_wt:
            print(f"  - {col}")
        return 1

    # WT prediction table must contain one record per
    # WT peptide + HLA prediction key.
    duplicate_wt_keys = wt.duplicated(
        subset=["WT_Peptide", "allele"]
    ).sum()

    if duplicate_wt_keys:
        print(
            "ERROR: duplicate WT peptide/HLA keys detected: "
            f"{duplicate_wt_keys}"
        )
        return 1

    for col in [
        "WT_Binding_Affinity",
        "WT_Binding_Percentile",
    ]:
        wt[col] = pd.to_numeric(
            wt[col],
            errors="coerce",
        )

    if wt[
        [
            "WT_Binding_Affinity",
            "WT_Binding_Percentile",
        ]
    ].isna().any().any():
        print("ERROR: missing/non-numeric WT values detected.")
        return 1

    merged = candidates.merge(
        wt,
        on=["WT_Peptide", "allele"],
        how="left",
        validate="many_to_one",
    )

    missing_wt_values = merged[
        [
            "WT_Binding_Affinity",
            "WT_Binding_Percentile",
        ]
    ].isna().any(axis=1).sum()

    if missing_wt_values:
        print(
            "ERROR: candidates without WT predictions: "
            f"{missing_wt_values}"
        )
        return 1

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    merged.to_csv(
        output_file,
        index=False,
    )

    unique_candidate_keys = (
        candidates[
            ["WT_Peptide", "allele"]
        ]
        .drop_duplicates()
        .shape[0]
    )

    print("=" * 72)
    print("WT PREDICTION MATCHING")
    print("=" * 72)
    print(f"Candidate rows:             {len(candidates)}")
    print(f"Unique WT peptide/HLA keys: {unique_candidate_keys}")
    print(f"WT prediction rows:         {len(wt)}")
    print(f"Merged rows:                {len(merged)}")
    print(f"Missing WT values:          {missing_wt_values}")
    print(
        f"WT rows reused across candidates: "
        f"{len(candidates) - unique_candidate_keys}"
    )
    print()
    print(f"Output: {output_file}")
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
