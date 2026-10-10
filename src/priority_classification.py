#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd


# ------------------------------------------------------------------
# PI biological thresholds
# ------------------------------------------------------------------

STRINGENT_AFFINITY = 50.0
STRINGENT_PERCENTILE = 0.5
FOLD_IMPROVEMENT_THRESHOLD = 2.0


def classify_priority(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply the mutually exclusive priority rules used in the
    verified PI final classification.

    Because lower affinity and lower percentile indicate stronger
    predicted MHC-I binding:

    P1 requires:
        mutant affinity <= 50 nM
        mutant percentile <= 0.5
        affinity fold improvement >= 2
        percentile fold improvement >= 2
        mutant affinity < WT affinity
        mutant percentile < WT percentile

    P2:
        stringent binding + mutant stronger on both measures,
        excluding P1.

    P3:
        stringent binding, excluding P1 and P2.

    Below stringent:
        all remaining primary candidates.
    """

    result = df.copy()

    result["Priority1"] = (
        (result["mhcflurry_affinity"] <= STRINGENT_AFFINITY)
        & (
            result["mhcflurry_affinity_percentile"]
            <= STRINGENT_PERCENTILE
        )
        & (
            result["affinity_fold_improvement"]
            >= FOLD_IMPROVEMENT_THRESHOLD
        )
        & (
            result["percentile_fold_improvement"]
            >= FOLD_IMPROVEMENT_THRESHOLD
        )
        & result["mutant_stronger_affinity"]
        & result["mutant_better_percentile"]
    )

    result["Priority"] = "Below stringent threshold"

    priority2 = (
        ~result["Priority1"]
        & result["stringent_binding"]
        & result["mutant_stronger_both"]
    )

    priority3 = (
        ~result["Priority1"]
        & ~priority2
        & result["stringent_binding"]
    )

    result.loc[
        priority3,
        "Priority"
    ] = "Priority 3"

    result.loc[
        priority2,
        "Priority"
    ] = "Priority 2"

    result.loc[
        result["Priority1"],
        "Priority"
    ] = "Priority 1"

    priority_order = {
        "Priority 1": 1,
        "Priority 2": 2,
        "Priority 3": 3,
        "Below stringent threshold": 4,
    }

    result["Priority_Order"] = (
        result["Priority"].map(priority_order)
    )

    return result


def main() -> int:

    root = Path(__file__).resolve().parents[1]

    input_file = (
        root
        / "example_results"
        / "mutant_wt_comparison.csv"
    )

    output_file = (
        root
        / "example_results"
        / "prioritized_candidates.csv"
    )

    summary_file = (
        root
        / "example_results"
        / "priority_summary.csv"
    )

    qc_file = (
        root
        / "example_results"
        / "priority_qc.txt"
    )

    if not input_file.exists():
        print(f"ERROR: input file not found: {input_file}")
        return 1

    df = pd.read_csv(input_file)

    required = [
        "Candidate_ID",
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

    result = classify_priority(df)

    # --------------------------------------------------------------
    # Validate that exactly one priority class was assigned to
    # every candidate.
    # --------------------------------------------------------------

    valid_priorities = {
        "Priority 1",
        "Priority 2",
        "Priority 3",
        "Below stringent threshold",
    }

    unexpected = set(result["Priority"]) - valid_priorities

    if unexpected:
        print(
            "ERROR: unexpected priority classes detected:"
        )
        for value in sorted(unexpected):
            print(f"  - {value}")
        return 1

    if result["Priority"].isna().any():
        print("ERROR: missing priority assignments.")
        return 1

    # The classifier must preserve every input row.
    if len(result) != len(df):
        print(
            "ERROR: priority classification changed "
            "the number of rows."
        )
        return 1

    # --------------------------------------------------------------
    # Independent logical consistency checks.
    # --------------------------------------------------------------

    p1 = result["Priority"] == "Priority 1"
    p2 = result["Priority"] == "Priority 2"
    p3 = result["Priority"] == "Priority 3"
    below = (
        result["Priority"]
        == "Below stringent threshold"
    )

    # Every P1 must satisfy every P1 requirement.
    p1_valid = (
        result.loc[p1, "mhcflurry_affinity"]
        .le(STRINGENT_AFFINITY)
        & result.loc[p1, "mhcflurry_affinity_percentile"]
        .le(STRINGENT_PERCENTILE)
        & result.loc[p1, "affinity_fold_improvement"]
        .ge(FOLD_IMPROVEMENT_THRESHOLD)
        & result.loc[p1, "percentile_fold_improvement"]
        .ge(FOLD_IMPROVEMENT_THRESHOLD)
        & result.loc[p1, "mutant_stronger_affinity"]
        & result.loc[p1, "mutant_better_percentile"]
    )

    if not p1_valid.all():
        print(
            "ERROR: one or more P1 records violate the "
            "P1 criteria."
        )
        return 1

    # P2 must be stringent + stronger on both, but not P1.
    p2_valid = (
        result.loc[p2, "stringent_binding"]
        & result.loc[p2, "mutant_stronger_both"]
        & ~result.loc[p2, "Priority1"]
    )

    if not p2_valid.all():
        print(
            "ERROR: one or more P2 records violate "
            "the P2 criteria."
        )
        return 1

    # P3 must be stringent but neither P1 nor P2.
    p3_valid = (
        result.loc[p3, "stringent_binding"]
        & ~result.loc[p3, "Priority1"]
        & ~result.loc[p3, "mutant_stronger_both"]
    )

    if not p3_valid.all():
        print(
            "ERROR: one or more P3 records violate "
            "the P3 criteria."
        )
        return 1

    # Below-stringent records must not satisfy stringent binding.
    below_valid = ~result.loc[below, "stringent_binding"]

    if not below_valid.all():
        print(
            "ERROR: below-stringent record passed "
            "the stringent threshold."
        )
        return 1

    # --------------------------------------------------------------
    # Sort according to biological priority.
    # Within each class, stronger predicted mutant binding is
    # placed first.
    # --------------------------------------------------------------

    result = result.sort_values(
        [
            "Priority_Order",
            "mhcflurry_affinity",
            "mhcflurry_affinity_percentile",
        ],
        ascending=[
            True,
            True,
            True,
        ],
    ).reset_index(drop=True)

    result["Priority_Rank"] = (
        range(1, len(result) + 1)
    )

    # --------------------------------------------------------------
    # Summary table.
    # --------------------------------------------------------------

    summary = (
        result.groupby(
            "Priority",
            sort=False,
        )
        .agg(
            candidates=("Candidate_ID", "size"),
            mean_mutant_affinity=(
                "mhcflurry_affinity",
                "mean",
            ),
            mean_mutant_percentile=(
                "mhcflurry_affinity_percentile",
                "mean",
            ),
            mean_affinity_fold=(
                "affinity_fold_improvement",
                "mean",
            ),
            mean_percentile_fold=(
                "percentile_fold_improvement",
                "mean",
            ),
        )
        .reset_index()
    )

    summary["Priority_Order"] = (
        summary["Priority"].map(
            {
                "Priority 1": 1,
                "Priority 2": 2,
                "Priority 3": 3,
                "Below stringent threshold": 4,
            }
        )
    )

    summary = summary.sort_values(
        "Priority_Order"
    ).reset_index(drop=True)

    # --------------------------------------------------------------
    # Save outputs.
    # --------------------------------------------------------------

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result.to_csv(
        output_file,
        index=False,
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    with open(qc_file, "w") as handle:
        handle.write(
            "PRIORITY CLASSIFICATION QC\n"
        )
        handle.write("=" * 72 + "\n")
        handle.write(
            f"Input rows: {len(df)}\n"
        )
        handle.write(
            f"Output rows: {len(result)}\n\n"
        )
        handle.write(
            "PI classification thresholds:\n"
        )
        handle.write(
            f"  Mutant affinity <= "
            f"{STRINGENT_AFFINITY} nM\n"
        )
        handle.write(
            f"  Mutant percentile <= "
            f"{STRINGENT_PERCENTILE}\n"
        )
        handle.write(
            f"  Fold improvement >= "
            f"{FOLD_IMPROVEMENT_THRESHOLD}\n\n"
        )
        handle.write(
            "Priority counts:\n"
        )

        for _, row in summary.iterrows():
            handle.write(
                f"  {row['Priority']}: "
                f"{int(row['candidates'])}\n"
            )

        handle.write(
            "\nLogical validation: PASSED\n"
        )

    print("=" * 72)
    print("PRIORITY CLASSIFICATION")
    print("=" * 72)

    print(
        f"Input rows:                  {len(df)}"
    )
    print(
        f"Output rows:                 {len(result)}"
    )

    print()
    print("Priority distribution:")

    counts = (
        result["Priority"]
        .value_counts()
        .reindex(
            [
                "Priority 1",
                "Priority 2",
                "Priority 3",
                "Below stringent threshold",
            ],
            fill_value=0,
        )
    )

    for priority, count in counts.items():
        print(
            f"  {priority:28s} {int(count)}"
        )

    print()
    print("Logical QC:")
    print("  P1 criteria: PASSED")
    print("  P2 criteria: PASSED")
    print("  P3 criteria: PASSED")
    print("  Below-stringent criteria: PASSED")

    # --------------------------------------------------------------
    # Synthetic fixture validation.
    #
    # This is ONLY a test of our public synthetic fixture.
    # The real internship data are not being loaded here.
    # --------------------------------------------------------------

    if "Synthetic_Target_Class" in result.columns:

        mismatch = (
            result["Priority"]
            != result["Synthetic_Target_Class"]
        ).sum()

        print()
        print(
            "Synthetic fixture classification check:"
        )
        print(
            f"  Mismatched rows:            {mismatch}"
        )

        if mismatch != 0:
            print(
                "ERROR: calculated priorities do not "
                "match the synthetic validation labels."
            )
            return 1

    print()
    print(f"Output: {output_file}")
    print(f"Summary: {summary_file}")
    print(f"QC log: {qc_file}")
    print()
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
