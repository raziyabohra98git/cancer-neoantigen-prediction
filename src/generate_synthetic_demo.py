#!/usr/bin/env python3

from pathlib import Path
import random

import numpy as np
import pandas as pd


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "example_data"
OUT.mkdir(parents=True, exist_ok=True)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

PATIENTS = [
    "SYN-001",
    "SYN-002",
    "SYN-003",
    "SYN-004",
    "SYN-005",
]

HLA_ALLELES = [
    "A*01:01",
    "A*02:01",
    "A*03:01",
    "A*24:02",
    "B*07:02",
    "B*08:01",
    "B*35:01",
    "B*44:03",
]


def random_peptide(length: int) -> str:
    return "".join(
        random.choice(AMINO_ACIDS)
        for _ in range(length)
    )


def mutate_peptide(wt: str, position: int, replacement_index: int) -> str:
    """
    Create a mutant peptide with one amino-acid substitution.

    position is 1-based.
    replacement_index is used to make the synthetic candidates
    deterministic and distinct within a WT peptide group.
    """
    idx = position - 1
    original = wt[idx]

    alternatives = [
        aa for aa in AMINO_ACIDS
        if aa != original
    ]

    replacement = alternatives[
        replacement_index % len(alternatives)
    ]

    chars = list(wt)
    chars[idx] = replacement

    return "".join(chars)


# ------------------------------------------------------------------
# Define 60 reusable WT peptide/HLA keys.
#
# Crucially, keys are partitioned by priority class:
#
#   P1      :  8 keys x 5 candidates  = 40
#   P2      : 16 keys x 5 candidates  = 80
#   P3      : 24 keys x 5 candidates  = 120
#   Below   : 12 keys x 30 candidates = 360
#
# Thus WT prediction reuse occurs, but no WT key crosses priority
# classes.
# ------------------------------------------------------------------

GROUPS = []

for i in range(8):
    GROUPS.append(("Priority 1", i, 5))

for i in range(16):
    GROUPS.append(("Priority 2", i + 8, 5))

for i in range(24):
    GROUPS.append(("Priority 3", i + 24, 5))

for i in range(12):
    GROUPS.append(("Below stringent threshold", i + 48, 30))


rows = []
mutant_predictions = []
wt_predictions = []

candidate_number = 1


for priority_class, group_id, group_size in GROUPS:

    # One WT peptide/HLA key is shared by all candidates in this group.
    peptide_length = [8, 9, 10, 11][group_id % 4]

    wt_peptide = random_peptide(peptide_length)
    allele = HLA_ALLELES[
        group_id % len(HLA_ALLELES)
    ]

    # --------------------------------------------------------------
    # Construct ONE WT prediction for this group.
    # Every candidate using this key will match this prediction.
    # --------------------------------------------------------------

    if priority_class == "Priority 1":

        mutant_affinity = random.uniform(10.0, 40.0)
        mutant_percentile = random.uniform(0.05, 0.30)

        affinity_fold = random.uniform(2.5, 6.0)
        percentile_fold = random.uniform(2.5, 6.0)

        wt_affinity = mutant_affinity * affinity_fold
        wt_percentile = mutant_percentile * percentile_fold

        # All candidates in the group receive these identical
        # mutant values, so every candidate satisfies P1.

        mutant_values = [
            (mutant_affinity, mutant_percentile)
            for _ in range(group_size)
        ]

    elif priority_class == "Priority 2":

        mutant_affinity = random.uniform(10.0, 40.0)
        mutant_percentile = random.uniform(0.05, 0.30)

        # Stronger mutant on both measurements, but affinity
        # fold improvement deliberately remains below 2.
        affinity_fold = random.uniform(1.20, 1.70)
        percentile_fold = random.uniform(2.50, 5.00)

        wt_affinity = mutant_affinity * affinity_fold
        wt_percentile = mutant_percentile * percentile_fold

        mutant_values = [
            (mutant_affinity, mutant_percentile)
            for _ in range(group_size)
        ]

    elif priority_class == "Priority 3":

        mutant_affinity = random.uniform(10.0, 40.0)
        mutant_percentile = random.uniform(0.05, 0.30)

        # Deliberately make mutant NOT stronger on both measures.
        wt_affinity = mutant_affinity * random.uniform(0.50, 0.90)
        wt_percentile = mutant_percentile * random.uniform(0.50, 0.90)

        mutant_values = [
            (mutant_affinity, mutant_percentile)
            for _ in range(group_size)
        ]

    else:

        mode = group_id % 3

        if mode == 0:
            # Fails affinity threshold.
            mutant_affinity = random.uniform(50.1, 500.0)
            mutant_percentile = random.uniform(0.05, 0.49)

        elif mode == 1:
            # Fails percentile threshold.
            mutant_affinity = random.uniform(10.0, 40.0)
            mutant_percentile = random.uniform(0.51, 2.0)

        else:
            # Fails both stringent thresholds.
            mutant_affinity = random.uniform(50.1, 500.0)
            mutant_percentile = random.uniform(0.51, 2.0)

        wt_affinity = random.uniform(25.0, 700.0)
        wt_percentile = random.uniform(0.1, 4.0)

        mutant_values = [
            (mutant_affinity, mutant_percentile)
            for _ in range(group_size)
        ]

    # Store exactly one WT prediction per WT peptide + HLA key.
    wt_predictions.append(
        {
            "WT_Peptide": wt_peptide,
            "allele": allele,
            "WT_Binding_Affinity": wt_affinity,
            "WT_Binding_Percentile": wt_percentile,
        }
    )

    # --------------------------------------------------------------
    # Generate candidate records sharing that WT key.
    # --------------------------------------------------------------

    for local_index in range(group_size):

        patient = PATIENTS[
            (candidate_number - 1) % len(PATIENTS)
        ]

        mutation_position = (
            local_index % peptide_length
        ) + 1

        mutant_peptide = mutate_peptide(
            wt_peptide,
            mutation_position,
            local_index,
        )

        gene = (
            f"GENE{((candidate_number - 1) % 120) + 1:03d}"
        )

        protein_id = (
            f"PROT{((candidate_number - 1) % 160) + 1:04d}"
        )

        hgvsp = (
            f"p.Synthetic{((candidate_number - 1) % 900) + 1}"
        )

        candidate_id = f"CAND-{candidate_number:04d}"

        mut_affinity, mut_percentile = mutant_values[
            local_index
        ]

        rows.append(
            {
                "Candidate_ID": candidate_id,
                "Sample_ID": patient,
                "Hugo_Symbol": gene,
                "Protein_ID": protein_id,
                "HGVSp": hgvsp,
                "Mutant_Peptide": mutant_peptide,
                "WT_Peptide": wt_peptide,
                "Peptide_Length": peptide_length,
                "Mutation_Position_in_Peptide": mutation_position,
                "allele": allele,
                "Synthetic_Target_Class": priority_class,
            }
        )

        mutant_predictions.append(
            {
                "Candidate_ID": candidate_id,
                "Mutant_Peptide": mutant_peptide,
                "allele": allele,
                "mhcflurry_affinity": mut_affinity,
                "mhcflurry_affinity_percentile": mut_percentile,
            }
        )

        candidate_number += 1


candidates = pd.DataFrame(rows)
mutant = pd.DataFrame(mutant_predictions)
wt = pd.DataFrame(wt_predictions)


manifest = pd.DataFrame(
    {
        "dataset": [
            "synthetic_candidates.csv",
            "synthetic_mutant_predictions.csv",
            "synthetic_wt_predictions.csv",
        ],
        "description": [
            "Synthetic candidate metadata and 8-11 aa peptide windows.",
            "Synthetic MHC-I mutant binding predictions; not real MHCflurry output.",
            "Synthetic WT binding predictions keyed by WT peptide and HLA allele.",
        ],
    }
)


candidates.to_csv(
    OUT / "synthetic_candidates.csv",
    index=False,
)

mutant.to_csv(
    OUT / "synthetic_mutant_predictions.csv",
    index=False,
)

wt.to_csv(
    OUT / "synthetic_wt_predictions.csv",
    index=False,
)

manifest.to_csv(
    OUT / "synthetic_dataset_manifest.csv",
    index=False
)


# ------------------------------------------------------------------
# Reconstruct the biological classification from the actual PI rules.
# ------------------------------------------------------------------

check = (
    candidates[
        [
            "Candidate_ID",
            "WT_Peptide",
            "allele",
            "Synthetic_Target_Class",
        ]
    ]
    .merge(
        mutant[
            [
                "Candidate_ID",
                "mhcflurry_affinity",
                "mhcflurry_affinity_percentile",
            ]
        ],
        on="Candidate_ID",
        how="left",
        validate="one_to_one",
    )
    .merge(
        wt,
        on=["WT_Peptide", "allele"],
        how="left",
        validate="many_to_one",
    )
)


check["affinity_fold_improvement"] = (
    check["WT_Binding_Affinity"]
    / check["mhcflurry_affinity"]
)

check["percentile_fold_improvement"] = (
    check["WT_Binding_Percentile"]
    / check["mhcflurry_affinity_percentile"]
)

check["mutant_stronger_affinity"] = (
    check["mhcflurry_affinity"]
    < check["WT_Binding_Affinity"]
)

check["mutant_better_percentile"] = (
    check["mhcflurry_affinity_percentile"]
    < check["WT_Binding_Percentile"]
)

check["mutant_stronger_both"] = (
    check["mutant_stronger_affinity"]
    & check["mutant_better_percentile"]
)

check["stringent_binding"] = (
    (check["mhcflurry_affinity"] <= 50)
    & (check["mhcflurry_affinity_percentile"] <= 0.5)
)

check["Priority1"] = (
    (check["mhcflurry_affinity"] <= 50)
    & (check["mhcflurry_affinity_percentile"] <= 0.5)
    & (check["affinity_fold_improvement"] >= 2)
    & (check["percentile_fold_improvement"] >= 2)
    & check["mutant_stronger_affinity"]
    & check["mutant_better_percentile"]
)


def classify(row):

    if row["Priority1"]:
        return "Priority 1"

    if (
        row["stringent_binding"]
        and row["mutant_stronger_both"]
    ):
        return "Priority 2"

    if row["stringent_binding"]:
        return "Priority 3"

    return "Below stringent threshold"


check["Calculated_Priority"] = check.apply(
    classify,
    axis=1,
)


expected = pd.Series(
    {
        "Below stringent threshold": 360,
        "Priority 1": 40,
        "Priority 2": 80,
        "Priority 3": 120,
    }
)

observed = (
    check["Calculated_Priority"]
    .value_counts()
    .reindex(expected.index, fill_value=0)
)

# ------------------------------------------------------------------
# Hard validation checks.
# ------------------------------------------------------------------

if not observed.equals(expected.astype(int)):
    raise RuntimeError(
        "Synthetic fixture validation failed.\n\n"
        f"Expected:\n{expected}\n\n"
        f"Observed:\n{observed}"
    )

if len(candidates) != 600:
    raise RuntimeError("Expected exactly 600 candidates.")

if len(mutant) != 600:
    raise RuntimeError("Expected exactly 600 mutant predictions.")

if len(wt) != 60:
    raise RuntimeError(
        "Expected exactly 60 unique WT peptide/HLA predictions."
    )

if not candidates["Peptide_Length"].isin(
    [8, 9, 10, 11]
).all():
    raise RuntimeError(
        "Unexpected peptide length detected."
    )

if not (
    candidates["Mutant_Peptide"].str.len()
    == candidates["Peptide_Length"]
).all():
    raise RuntimeError(
        "Mutant peptide length mismatch."
    )

if not (
    candidates["WT_Peptide"].str.len()
    == candidates["Peptide_Length"]
).all():
    raise RuntimeError(
        "WT peptide length mismatch."
    )

if (
    candidates["Mutant_Peptide"]
    == candidates["WT_Peptide"]
).any():
    raise RuntimeError(
        "Mutant peptide equals WT peptide."
    )


print("=" * 72)
print("SYNTHETIC DEMONSTRATION DATA GENERATED + VALIDATED")
print("=" * 72)

print(f"Candidate rows:            {len(candidates)}")
print(f"Mutant prediction rows:    {len(mutant)}")
print(f"Unique WT predictions:     {len(wt)}")

print()
print("Peptide lengths:")
print(
    candidates["Peptide_Length"]
    .value_counts()
    .sort_index()
)

print()
print("Validated priority distribution:")
print(observed.to_string())

print()
print("WT-key reuse:")
print(
    f"{len(candidates)} candidate rows "
    f"share {len(wt)} unique WT peptide/HLA keys"
)

print()
print("Files:")

for path in sorted(
    OUT.glob("synthetic_*.csv")
):
    print(f"  {path}")
