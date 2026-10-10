#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


# ------------------------------------------------------------------
# Reference set used in the completed project.
#
# These are public gene symbols; their presence here does NOT mean
# they were observed in the synthetic demonstration dataset.
# ------------------------------------------------------------------

LUAD_CANCER_GENES = {
    "ALK",
    "ARID1A",
    "BRAF",
    "CDKN2A",
    "EGFR",
    "ERBB2",
    "KEAP1",
    "KRAS",
    "MET",
    "MGA",
    "NF1",
    "NTRK1",
    "NTRK2",
    "NTRK3",
    "PIK3CA",
    "RB1",
    "RBM10",
    "RET",
    "ROS1",
    "SETD2",
    "SMARCA4",
    "STK11",
    "TP53",
}


def main() -> int:

    root = Path(__file__).resolve().parents[1]

    input_file = (
        root
        / "example_results"
        / "ai_ranked_candidates.csv"
    )

    output_dir = (
        root
        / "example_results"
    )

    required = [
        "Candidate_ID",
        "Sample_ID",
        "Hugo_Symbol",
        "HGVSp",
        "Mutant_Peptide",
        "allele",
        "Priority",
        "Priority1",
        "AI_score",
    ]

    if not input_file.exists():
        print(
            f"ERROR: input file not found: {input_file}"
        )
        return 1

    df = pd.read_csv(input_file)

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
    # Biological analysis must use the final P1 subset.
    # --------------------------------------------------------------

    p1 = df[
        df["Priority"] == "Priority 1"
    ].copy()

    if p1.empty:
        print(
            "ERROR: no Priority 1 candidates found."
        )
        return 1

    # --------------------------------------------------------------
    # Unique genes.
    # --------------------------------------------------------------

    unique_genes = (
        p1["Hugo_Symbol"]
        .dropna()
        .drop_duplicates()
        .sort_values()
        .reset_index(drop=True)
    )

    unique_genes_df = pd.DataFrame(
        {
            "Hugo_Symbol": unique_genes
        }
    )

    # --------------------------------------------------------------
    # Exact mutation definition:
    #
    # gene + HGVSp
    #
    # Recurrence is evaluated across distinct patients.
    # --------------------------------------------------------------

    p1["Mutation_Key"] = (
        p1["Hugo_Symbol"].astype(str)
        + ":"
        + p1["HGVSp"].astype(str)
    )

    p1["Patient_Mutation_Key"] = (
        p1["Sample_ID"].astype(str)
        + "|"
        + p1["Mutation_Key"]
    )

    unique_mutations = (
        p1["Mutation_Key"]
        .drop_duplicates()
        .sort_values()
        .reset_index(drop=True)
    )

    unique_mutations_df = pd.DataFrame(
        {
            "Mutation_Key": unique_mutations
        }
    )

    # --------------------------------------------------------------
    # Mutation recurrence across patients.
    # --------------------------------------------------------------

    mutation_recurrence = (
        p1.groupby(
            "Mutation_Key"
        )
        .agg(
            patients=(
                "Sample_ID",
                "nunique",
            ),
            candidate_rows=(
                "Candidate_ID",
                "size",
            ),
            gene=(
                "Hugo_Symbol",
                "first",
            ),
            HGVSp=(
                "HGVSp",
                "first",
            ),
        )
        .reset_index()
    )

    mutation_recurrence[
        "recurrent_across_patients"
    ] = mutation_recurrence[
        "patients"
    ] > 1

    mutation_recurrence = mutation_recurrence.sort_values(
        [
            "recurrent_across_patients",
            "patients",
            "candidate_rows",
        ],
        ascending=[
            False,
            False,
            False,
        ],
    ).reset_index(drop=True)

    # --------------------------------------------------------------
    # Patient-level biological summary.
    # --------------------------------------------------------------

    patient_summary = (
        p1.groupby("Sample_ID")
        .agg(
            P1_candidates=(
                "Candidate_ID",
                "size",
            ),
            unique_genes=(
                "Hugo_Symbol",
                "nunique",
            ),
            unique_mutations=(
                "Mutation_Key",
                "nunique",
            ),
            unique_mutant_peptides=(
                "Mutant_Peptide",
                "nunique",
            ),
            unique_HLA_alleles=(
                "allele",
                "nunique",
            ),
            maximum_AI_score=(
                "AI_score",
                "max",
            ),
        )
        .reset_index()
    )

    # --------------------------------------------------------------
    # Gene contribution.
    # --------------------------------------------------------------

    gene_contribution = (
        p1.groupby("Hugo_Symbol")
        .agg(
            P1_candidates=(
                "Candidate_ID",
                "size",
            ),
            patients=(
                "Sample_ID",
                "nunique",
            ),
            unique_mutations=(
                "Mutation_Key",
                "nunique",
            ),
            unique_peptides=(
                "Mutant_Peptide",
                "nunique",
            ),
            unique_HLA_alleles=(
                "allele",
                "nunique",
            ),
            best_AI_score=(
                "AI_score",
                "max",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "P1_candidates",
                "best_AI_score",
            ],
            ascending=[
                False,
                False,
            ],
        )
        .reset_index(drop=True)
    )

    # --------------------------------------------------------------
    # HLA contribution.
    # --------------------------------------------------------------

    hla_contribution = (
        p1.groupby("allele")
        .agg(
            P1_candidates=(
                "Candidate_ID",
                "size",
            ),
            patients=(
                "Sample_ID",
                "nunique",
            ),
            unique_genes=(
                "Hugo_Symbol",
                "nunique",
            ),
            mean_AI_score=(
                "AI_score",
                "mean",
            ),
            maximum_AI_score=(
                "AI_score",
                "max",
            ),
        )
        .reset_index()
        .sort_values(
            "P1_candidates",
            ascending=False,
        )
        .reset_index(drop=True)
    )

    # --------------------------------------------------------------
    # Cancer-gene reference overlap.
    # --------------------------------------------------------------

    observed_genes = set(
        unique_genes.astype(str)
    )

    overlap = sorted(
        observed_genes & LUAD_CANCER_GENES
    )

    cancer_gene_representation = pd.DataFrame(
        {
            "Reference_Gene": sorted(
                LUAD_CANCER_GENES
            ),
        }
    )

    cancer_gene_representation[
        "Observed_in_P1"
    ] = cancer_gene_representation[
        "Reference_Gene"
    ].isin(observed_genes)

    cancer_gene_representation[
        "P1_candidate_rows"
    ] = cancer_gene_representation[
        "Reference_Gene"
    ].map(
        p1["Hugo_Symbol"].value_counts()
    ).fillna(0).astype(int)

    # --------------------------------------------------------------
    # Save outputs.
    # --------------------------------------------------------------

    unique_genes_df.to_csv(
        output_dir
        / "unique_p1_genes.csv",
        index=False,
    )

    unique_mutations_df.to_csv(
        output_dir
        / "unique_p1_mutations.csv",
        index=False,
    )

    patient_summary.to_csv(
        output_dir
        / "patient_biological_summary.csv",
        index=False,
    )

    gene_contribution.to_csv(
        output_dir
        / "gene_contribution.csv",
        index=False,
    )

    mutation_recurrence.to_csv(
        output_dir
        / "mutation_recurrence.csv",
        index=False,
    )

    hla_contribution.to_csv(
        output_dir
        / "hla_contribution.csv",
        index=False,
    )

    cancer_gene_representation.to_csv(
        output_dir
        / "cancer_gene_representation.csv",
        index=False,
    )

    # --------------------------------------------------------------
    # QC.
    # --------------------------------------------------------------

    recurrent_count = int(
        mutation_recurrence[
            "recurrent_across_patients"
        ].sum()
    )

    qc_file = (
        output_dir
        / "biological_analysis_qc.txt"
    )

    with open(qc_file, "w") as handle:

        handle.write(
            "BIOLOGICAL / GENE ANALYSIS QC\n"
        )
        handle.write(
            "=" * 72 + "\n"
        )

        handle.write(
            f"P1 candidate rows: {len(p1)}\n"
        )

        handle.write(
            f"Unique P1 genes: "
            f"{len(unique_genes_df)}\n"
        )

        handle.write(
            f"Unique P1 mutations: "
            f"{len(unique_mutations_df)}\n"
        )

        handle.write(
            f"Unique mutant peptides: "
            f"{p1['Mutant_Peptide'].nunique()}\n"
        )

        handle.write(
            f"Unique HLA alleles: "
            f"{p1['allele'].nunique()}\n"
        )

        handle.write(
            f"Patients: "
            f"{p1['Sample_ID'].nunique()}\n"
        )

        handle.write(
            f"Exact mutations recurrent across patients: "
            f"{recurrent_count}\n"
        )

        handle.write(
            f"Reference genes: "
            f"{len(LUAD_CANCER_GENES)}\n"
        )

        handle.write(
            f"Reference overlap: "
            f"{len(overlap)}\n"
        )

        handle.write(
            "\nQC checks: PASSED\n"
        )

    # --------------------------------------------------------------
    # Console summary.
    # --------------------------------------------------------------

    print("=" * 72)
    print("BIOLOGICAL / CANCER-GENE ANALYSIS")
    print("=" * 72)

    print(
        f"P1 candidates:              {len(p1)}"
    )

    print(
        f"Unique genes:               "
        f"{p1['Hugo_Symbol'].nunique()}"
    )

    print(
        f"Unique mutations:           "
        f"{p1['Mutation_Key'].nunique()}"
    )

    print(
        f"Unique mutant peptides:     "
        f"{p1['Mutant_Peptide'].nunique()}"
    )

    print(
        f"Unique HLA alleles:         "
        f"{p1['allele'].nunique()}"
    )

    print(
        f"Patients:                   "
        f"{p1['Sample_ID'].nunique()}"
    )

    print()

    print(
        "Exact mutations recurrent across patients: "
        f"{recurrent_count}"
    )

    print(
        "Cancer-reference genes:     "
        f"{len(LUAD_CANCER_GENES)}"
    )

    print(
        "Reference overlap:          "
        f"{len(overlap)}"
    )

    if overlap:
        print()
        print("Overlapping reference genes:")
        for gene in overlap:
            print(f"  {gene}")

    print()

    print("Outputs:")

    for filename in [
        "unique_p1_genes.csv",
        "unique_p1_mutations.csv",
        "patient_biological_summary.csv",
        "gene_contribution.csv",
        "mutation_recurrence.csv",
        "hla_contribution.csv",
        "cancer_gene_representation.csv",
        "biological_analysis_qc.txt",
    ]:
        print(
            f"  {output_dir / filename}"
        )

    print()
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
