#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


REQUIRED_COLUMNS = [
    "Candidate_ID",
    "Mutant_Peptide",
    "allele",
    "mhcflurry_affinity",
    "mhcflurry_affinity_percentile",
]


def main() -> int:
    root = Path(__file__).resolve().parents[1]

    input_file = (
        root
        / "example_data"
        / "synthetic_mutant_predictions.csv"
    )

    if not input_file.exists():
        print(f"ERROR: input file not found: {input_file}")
        return 1

    df = pd.read_csv(input_file)

    missing = [
        col for col in REQUIRED_COLUMNS
        if col not in df.columns
    ]

    if missing:
        print("ERROR: missing required columns:")
        for col in missing:
            print(f"  - {col}")
        return 1

    # Convert prediction columns explicitly to numeric.
    for col in [
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
    ]:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    checks = {}

    checks["rows"] = len(df)

    checks["missing_affinity"] = int(
        df["mhcflurry_affinity"].isna().sum()
    )

    checks["missing_percentile"] = int(
        df["mhcflurry_affinity_percentile"].isna().sum()
    )

    checks["nonpositive_affinity"] = int(
        (df["mhcflurry_affinity"] <= 0).sum()
    )

    checks["nonpositive_percentile"] = int(
        (df["mhcflurry_affinity_percentile"] <= 0).sum()
    )

    checks["duplicate_candidate_ids"] = int(
        df["Candidate_ID"].duplicated().sum()
    )

    checks["duplicate_candidate_peptide_hla"] = int(
        df[
            ["Candidate_ID", "Mutant_Peptide", "allele"]
        ].duplicated().sum()
    )

    # Synthetic fixture should contain only the intended
    # MHC-I prediction values and valid positive numbers.
    errors = []

    if checks["rows"] == 0:
        errors.append("Prediction table is empty.")

    if checks["missing_affinity"] > 0:
        errors.append("Missing mutant affinity values detected.")

    if checks["missing_percentile"] > 0:
        errors.append(
            "Missing mutant percentile values detected."
        )

    if checks["nonpositive_affinity"] > 0:
        errors.append(
            "Non-positive mutant affinity detected."
        )

    if checks["nonpositive_percentile"] > 0:
        errors.append(
            "Non-positive mutant percentile detected."
        )

    if checks["duplicate_candidate_ids"] > 0:
        errors.append(
            "Duplicate Candidate_ID values detected."
        )

    if checks["duplicate_candidate_peptide_hla"] > 0:
        errors.append(
            "Duplicate candidate/peptide/HLA prediction keys detected."
        )

    print("=" * 72)
    print("MUTANT BINDING PREDICTION QC")
    print("=" * 72)
    print(f"Input: {input_file}")
    print(f"Rows: {checks['rows']}")
    print()

    print("QC checks:")
    print(
        f"  Missing affinity:              "
        f"{checks['missing_affinity']}"
    )
    print(
        f"  Missing percentile:            "
        f"{checks['missing_percentile']}"
    )
    print(
        f"  Non-positive affinity:         "
        f"{checks['nonpositive_affinity']}"
    )
    print(
        f"  Non-positive percentile:       "
        f"{checks['nonpositive_percentile']}"
    )
    print(
        f"  Duplicate Candidate_ID:        "
        f"{checks['duplicate_candidate_ids']}"
    )
    print(
        f"  Duplicate prediction keys:    "
        f"{checks['duplicate_candidate_peptide_hla']}"
    )

    print()

    if errors:
        print("QC STATUS: FAILED")
        for error in errors:
            print(f"  ERROR: {error}")
        return 1

    print("QC STATUS: PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
