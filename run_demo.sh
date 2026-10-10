#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$ROOT"

echo
echo "======================================================================"
echo " CANCER NEOANTIGEN PREDICTION — PUBLIC DEMONSTRATION"
echo "======================================================================"
echo

echo "[1/9] Generating synthetic demonstration data..."
python src/generate_synthetic_demo.py

echo
echo "[2/9] Running mutant binding QC..."
python src/binding_qc.py

echo
echo "[3/9] Applying primary MHC-I binding filter..."
python src/binding_filter.py

echo
echo "[4/9] Matching WT predictions..."
python src/wt_matching.py

echo
echo "[5/9] Calculating mutant-vs-WT comparison..."
python src/mutant_wt_comparison.py

echo
echo "[6/9] Assigning biological priorities..."
python src/priority_classification.py

echo
echo "[7/9] Running AI-assisted ranking..."
python src/ml_ranking.py

echo
echo "[8/9] Running biological/cancer-gene analysis..."
python src/cancer_gene_analysis.py

echo
echo "[9/9] Running pathway-analysis interface..."
python src/pathway_enrichment.py \
    --query example_results/unique_p1_genes.csv \
    --query-column Hugo_Symbol \
    --background example_results/unique_p1_genes.csv \
    --background-column Hugo_Symbol

echo
echo "======================================================================"
echo " DEMONSTRATION COMPLETED"
echo "======================================================================"
echo
echo "Generated outputs:"
echo "  example_results/primary_binding_candidates.csv"
echo "  example_results/wt_aware_candidates.csv"
echo "  example_results/mutant_wt_comparison.csv"
echo "  example_results/prioritized_candidates.csv"
echo "  example_results/ai_ranked_candidates.csv"
echo "  example_results/feature_importance.csv"
echo "  example_results/patient_ai_summary.csv"
echo "  example_results/hla_ai_summary.csv"
echo "  example_results/gene_ai_summary.csv"
echo "  example_results/unique_p1_genes.csv"
echo "  example_results/unique_p1_mutations.csv"
echo "  example_results/patient_biological_summary.csv"
echo "  example_results/gene_contribution.csv"
echo "  example_results/mutation_recurrence.csv"
echo "  example_results/hla_contribution.csv"
echo "  example_results/cancer_gene_representation.csv"
echo
echo "Synthetic pathway enrichment was intentionally skipped."
echo "See src/pathway_enrichment.py for real-data usage."
echo
