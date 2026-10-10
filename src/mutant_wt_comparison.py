#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


def main() -> int:
    root = Path(__file__).resolve().parents[1]

    input_file = (
        root
        / "example_results"
        / "wt_aware_candidates.csv"
    )

    output_file = (
        root
        / "example_results"
        / "mutant_wt_comparison.csv"
    )

    if not input_file.exists():
        print(f"ERROR: input file not found: {input_file}")
        return 1

    df = pd.read_csv(input_file)

    required = [
        "Candidate_ID",
        "Sample_ID",
        "Hugo_Symbol",
        "HGVSp",
        "Mutant_Peptide",
        "WT_Peptide",
        "Peptide_Length",
        "Mutation_Position_in_Peptide",
        "allele",
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
        "WT_Binding_Affinity",
        "WT_Binding_Percentile",
    ]

    missing = [
        col for col in required
        if col not in df.columns
    ]

    if missing:
        print("ERROR: missing required columns:")
        for col in missing:
            print(f"  - {col}")
        return 1

    # ----------------------------------------------------------
    # Explicit numeric conversion.
    # ----------------------------------------------------------

    numeric_columns = [
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
        "WT_Binding_Affinity",
        "WT_Binding_Percentile",
    ]

    for col in numeric_columns:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    if df[numeric_columns].isna().any().any():
        print(
            "ERROR: missing/non-numeric mutant or WT values."
        )
        return 1

    if (df[numeric_columns] <= 0).any().any():
        print(
            "ERROR: non-positive affinity/percentile value detected."
        )
        return 1

    # ----------------------------------------------------------
    # Mutant vs WT comparison.
    #
    # Lower affinity = stronger predicted binding.
    # Lower percentile = stronger predicted binding.
    # ----------------------------------------------------------

    df["affinity_fold_improvement"] = (
        df["WT_Binding_Affinity"]
        / df["mhcflurry_affinity"]
    )

    df["percentile_fold_improvement"] = (
        df["WT_Binding_Percentile"]
        / df["mhcflurry_affinity_percentile"]
    )

    df["mutant_stronger_affinity"] = (
        df["mhcflurry_affinity"]
        < df["WT_Binding_Affinity"]
    )

    df["mutant_better_percentile"] = (
        df["mhcflurry_affinity_percentile"]
        < df["WT_Binding_Percentile"]
    )

    df["mutant_stronger_both"] = (
        df["mutant_stronger_affinity"]
        & df["mutant_better_percentile"]
    )

    # ----------------------------------------------------------
    # Stringent mutant binding threshold used by the final
    # prioritization stage.
    # ----------------------------------------------------------

    df["stringent_binding"] = (
        (df["mhcflurry_affinity"] <= 50.0)
        & (
            df["mhcflurry_affinity_percentile"]
            <= 0.5
        )
    )

    # ----------------------------------------------------------
    # Basic QC.
    # ----------------------------------------------------------

    if not df["affinity_fold_improvement"].gt(0).all():
        print(
            "ERROR: invalid affinity fold-improvement values."
        )
        return 1

    if not df["percentile_fold_improvement"].gt(0).all():
        print(
            "ERROR: invalid percentile fold-improvement values."
        )
        return 1

    output_columns = [
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
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
        "WT_Binding_Affinity",
        "WT_Binding_Percentile",
        "affinity_fold_improvement",
        "percentile_fold_improvement",
        "mutant_stronger_affinity",
        "mutant_better_percentile",
        "mutant_stronger_both",
        "stringent_binding",
        "Synthetic_Target_Class",
    ]

    output_columns = [
        col for col in output_columns
        if col in df.columns
    ]

    result = df[output_columns].copy()

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result.to_csv(
        output_file,
        index=False,
    )

    print("=" * 72)
    print("MUTANT vs WT COMPARISON")
    print("=" * 72)

    print(f"Input rows:                  {len(df)}")
    print(f"Output rows:                 {len(result)}")
    print()

    print("Stringent mutant binding:")
    print(
        f"  Passed:                    "
        f"{int(df['stringent_binding'].sum())}"
    )
    print(
        f"  Failed:                    "
        f"{int((~df['stringent_binding']).sum())}"
    )
    print()

    print("Mutant stronger than WT:")
    print(
        f"  Affinity:                  "
        f"{int(df['mutant_stronger_affinity'].sum())}"
    )
    print(
        f"  Percentile:                "
        f"{int(df['mutant_better_percentile'].sum())}"
    )
    print(
        f"  Both:                      "
        f"{int(df['mutant_stronger_both'].sum())}"
    )
    print()

    print("Fold-improvement ranges:")
    print(
        f"  Affinity:                  "
        f"{df['affinity_fold_improvement'].min():.4f}"
        f" - "
        f"{df['affinity_fold_improvement'].max():.4f}"
    )
    print(
        f"  Percentile:                "
        f"{df['percentile_fold_improvement'].min():.4f}"
        f" - "
        f"{df['percentile_fold_improvement'].max():.4f}"
    )
    print()

    print(f"Output: {output_file}")
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
