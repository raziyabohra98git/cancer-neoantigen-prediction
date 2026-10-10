#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys

import pandas as pd


def read_gene_list(path: Path, column: str) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    df = pd.read_csv(path)

    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found in {path}"
        )

    genes = (
        df[column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    return sorted(
        {
            gene
            for gene in genes
            if gene
        }
    )


def is_synthetic_demo(genes: list[str]) -> bool:
    """
    Public demo protection.

    The demo uses placeholder symbols such as GENE001.
    They must not be sent to g:Profiler and interpreted
    as biological findings.
    """
    return bool(genes) and all(
        gene.startswith("GENE")
        for gene in genes
    )


def run_gprofiler(
    query_genes: list[str],
    background_genes: list[str],
) -> pd.DataFrame:

    from gprofiler import GProfiler

    gp = GProfiler(
        return_dataframe=True
    )

    result = gp.profile(
        organism="hsapiens",
        query=query_genes,
        sources=[
            "GO:BP",
            "REAC",
        ],
        user_threshold=0.05,
        significance_threshold_method="fdr",
        domain_scope="custom",
        background=background_genes,
        all_results=True,
        no_iea=True,
        no_evidences=False,
    )

    if result is None:
        return pd.DataFrame()

    return result


def main() -> int:

    parser = argparse.ArgumentParser(
        description=(
            "GO Biological Process and Reactome "
            "enrichment using a custom background."
        )
    )

    parser.add_argument(
        "--query",
        required=True,
        help="CSV containing query genes.",
    )

    parser.add_argument(
        "--query-column",
        default="Hugo_Symbol",
        help="Query gene column.",
    )

    parser.add_argument(
        "--background",
        required=True,
        help="CSV containing background genes.",
    )

    parser.add_argument(
        "--background-column",
        default="Hugo_Symbol",
        help="Background gene column.",
    )

    parser.add_argument(
        "--output",
        default="example_results/pathway_enrichment.csv",
    )

    parser.add_argument(
        "--qc-output",
        default="example_results/pathway_enrichment_qc.txt",
    )

    args = parser.parse_args()

    query_path = Path(args.query)
    background_path = Path(args.background)
    output_path = Path(args.output)
    qc_path = Path(args.qc_output)

    try:
        query_genes = read_gene_list(
            query_path,
            args.query_column,
        )

        background_genes = read_gene_list(
            background_path,
            args.background_column,
        )

    except (
        FileNotFoundError,
        ValueError,
    ) as exc:

        print(f"ERROR: {exc}")
        return 1

    print("=" * 72)
    print("PATHWAY ENRICHMENT")
    print("=" * 72)

    print(
        f"Query genes:                {len(query_genes)}"
    )

    print(
        f"Background genes:           {len(background_genes)}"
    )

    if not query_genes:
        print("ERROR: query gene set is empty.")
        return 1

    if not background_genes:
        print("ERROR: background gene set is empty.")
        return 1

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    qc_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------------
    # Synthetic public demo: do not perform biological enrichment.
    # --------------------------------------------------------------

    if is_synthetic_demo(query_genes):

        empty = pd.DataFrame(
            columns=[
                "source",
                "native",
                "name",
                "p_value",
                "significant",
            ]
        )

        empty.to_csv(
            output_path,
            index=False,
        )

        with open(qc_path, "w") as handle:

            handle.write(
                "PATHWAY ENRICHMENT QC\n"
            )
            handle.write("=" * 72 + "\n")
            handle.write(
                f"Query genes: {len(query_genes)}\n"
            )
            handle.write(
                f"Background genes: "
                f"{len(background_genes)}\n"
            )
            handle.write(
                "Status: SKIPPED\n"
            )
            handle.write(
                "Reason: synthetic placeholder gene symbols\n"
            )
            handle.write(
                "No biological enrichment was interpreted.\n"
            )

        print()
        print("DEMO STATUS: ENRICHMENT SKIPPED")
        print(
            "Reason: public demo uses synthetic "
            "placeholder gene symbols."
        )
        print(
            "No biological pathway interpretation "
            "has been claimed."
        )
        print()
        print(f"Output: {output_path}")
        print(f"QC log: {qc_path}")

        return 0

    # --------------------------------------------------------------
    # Real-data pathway analysis.
    # --------------------------------------------------------------

    try:
        result = run_gprofiler(
            query_genes,
            background_genes,
        )

    except Exception as exc:

        print(
            "ERROR: g:Profiler analysis failed."
        )
        print(f"Details: {exc}")
        return 1

    result.to_csv(
        output_path,
        index=False,
    )

    significant_count = 0

    if (
        "significant" in result.columns
        and len(result) > 0
    ):
        significant_count = int(
            result["significant"]
            .fillna(False)
            .astype(bool)
            .sum()
        )

    source_counts = {}

    if (
        "source" in result.columns
        and len(result) > 0
    ):
        source_counts = (
            result["source"]
            .value_counts()
            .to_dict()
        )

    with open(qc_path, "w") as handle:

        handle.write(
            "PATHWAY ENRICHMENT QC\n"
        )
        handle.write("=" * 72 + "\n")
        handle.write(
            f"Query genes: {len(query_genes)}\n"
        )
        handle.write(
            f"Background genes: "
            f"{len(background_genes)}\n"
        )
        handle.write(
            "Organism: hsapiens\n"
        )
        handle.write(
            "Sources: GO:BP + Reactome\n"
        )
        handle.write(
            "Custom background: YES\n"
        )
        handle.write(
            "Correction: FDR\n"
        )
        handle.write(
            "Threshold: 0.05\n"
        )
        handle.write(
            "all_results=True\n"
        )
        handle.write(
            "no_iea=True\n"
        )
        handle.write(
            "no_evidences=False\n"
        )
        handle.write(
            f"Returned annotations: {len(result)}\n"
        )
        handle.write(
            f"Significant annotations: "
            f"{significant_count}\n"
        )

        handle.write(
            "\nSource counts:\n"
        )

        for source, count in sorted(
            source_counts.items()
        ):
            handle.write(
                f"  {source}: {count}\n"
            )

        handle.write(
            "\nQC STATUS: PASSED\n"
        )

    print()
    print(
        f"Returned annotations:       {len(result)}"
    )
    print(
        f"FDR-significant annotations: "
        f"{significant_count}"
    )
    print()
    print(f"Output: {output_path}")
    print(f"QC log: {qc_path}")
    print("QC STATUS: PASSED")

    return 0


if __name__ == "__main__":
    sys.exit(main())
