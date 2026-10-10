#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


AFFINITY_THRESHOLD = 500.0
PERCENTILE_THRESHOLD = 2.0


def main() -> int:
    root = Path(__file__).resolve().parents[1]

    candidate_file = (
        root
        / "example_data"
        / "synthetic_candidates.csv"
    )

    mutant_file = (
        root
        / "example_data"
        / "synthetic_mutant_predictions.csv"
    )

    output_file = (
        root
        / "example_results"
        / "primary_binding_candidates.csv"
    )

    if not candidate_file.exists():
        print(f"ERROR: candidate file not found: {candidate_file}")
        return 1

    if not mutant_file.exists():
        print(f"ERROR: mutant prediction file not found: {mutant_file}")
        return 1

    candidates = pd.read_csv(candidate_file)
    mutant = pd.read_csv(mutant_file)

    candidate_required = [
        "Candidate_ID",
        "Sample_ID",
        "Hugo_Symbol",
        "Protein_ID",
        "HGVSp",
        "Mutant_Peptide",
        "WT_Peptide",
        "Peptide_Length",
        "Mutation_Position_in_Peptide",
        "allele",
    ]

    mutant_required = [
        "Candidate_ID",
        "Mutant_Peptide",
        "allele",
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
    ]

    missing_candidates = [
        col
        for col in candidate_required
        if col not in candidates.columns
    ]

    missing_mutant = [
        col
        for col in mutant_required
        if col not in mutant.columns
    ]

    if missing_candidates:
        print("ERROR: candidate metadata missing:")
        for col in missing_candidates:
            print(f"  - {col}")
        return 1

    if missing_mutant:
        print("ERROR: mutant prediction columns missing:")
        for col in missing_mutant:
            print(f"  - {col}")
        return 1

    # Candidate IDs must uniquely identify the candidate metadata.
    if candidates["Candidate_ID"].duplicated().any():
        print("ERROR: duplicate Candidate_ID values in candidates.")
        return 1

    # One mutant prediction per candidate.
    if mutant["Candidate_ID"].duplicated().any():
        print(
            "ERROR: duplicate Candidate_ID values in "
            "mutant predictions."
        )
        return 1

    # Merge candidate metadata with mutant binding predictions.
    merged = candidates.merge(
        mutant[
            [
                "Candidate_ID",
                "Mutant_Peptide",
                "allele",
                "mhcflurry_affinity",
                "mhcflurry_affinity_percentile",
            ]
        ],
        on=["Candidate_ID", "Mutant_Peptide", "allele"],
        how="inner",
        validate="one_to_one",
    )

    if len(merged) != len(candidates):
        print(
            "ERROR: not every candidate received exactly one "
            "mutant prediction."
        )
        print(
            f"Candidate rows: {len(candidates)}"
        )
        print(
            f"Merged rows:    {len(merged)}"
        )
        return 1

    # Explicit numeric conversion.
    for col in [
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
    ]:
        merged[col] = pd.to_numeric(
            merged[col],
            errors="coerce",
        )

    if merged[
        [
            "mhcflurry_affinity",
            "mhcflurry_affinity_percentile",
        ]
    ].isna().any().any():
        print(
            "ERROR: missing/non-numeric mutant prediction values."
        )
        return 1

    # Primary MHC-I binding filter.
    primary = merged[
        (merged["mhcflurry_affinity"] <= AFFINITY_THRESHOLD)
        & (
            merged["mhcflurry_affinity_percentile"]
            <= PERCENTILE_THRESHOLD
        )
    ].copy()

    primary = primary.sort_values(
        [
            "mhcflurry_affinity",
            "mhcflurry_affinity_percentile",
        ]
    ).reset_index(drop=True)

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    primary.to_csv(
        output_file,
        index=False,
    )

    print("=" * 72)
    print("PRIMARY MHC-I BINDING FILTER")
    print("=" * 72)
    print(f"Candidate metadata rows:     {len(candidates)}")
    print(f"Mutant prediction rows:      {len(mutant)}")
    print(f"Merged candidate rows:       {len(merged)}")
    print()
    print(
        f"Affinity threshold:          "
        f"<= {AFFINITY_THRESHOLD} nM"
    )
    print(
        f"Percentile threshold:        "
        f"<= {PERCENTILE_THRESHOLD}"
    )
    print(f"Primary candidates:          {len(primary)}")
    print(
        f"Retention:                   "
        f"{len(primary) / len(merged) * 100:.2f}%"
    )
    print()
    print("Required downstream columns present:")

    downstream_columns = [
        "WT_Peptide",
        "allele",
        "Sample_ID",
        "Hugo_Symbol",
        "HGVSp",
        "Peptide_Length",
        "Mutation_Position_in_Peptide",
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
    ]

    for col in downstream_columns:
        status = "YES" if col in primary.columns else "NO"
        print(f"  {col:35s} {status}")

    print()
    print(f"Output: {output_file}")
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
