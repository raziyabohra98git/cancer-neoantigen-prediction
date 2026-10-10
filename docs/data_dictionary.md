# Data Dictionary

## 1. Important distinction

All public demonstration records are synthetic. Sample identifiers, gene
symbols, peptides, and prediction values in the demo are not real patient
results.

## 2. Input data

| Field | Meaning |
|---|---|
| Candidate_ID | Unique identifier for a candidate record |
| Sample_ID | Synthetic sample identifier in the public demo |
| Hugo_Symbol | Gene symbol associated with the candidate |
| HGVSp | Protein-level variant description |
| WT_Peptide | Wild-type peptide sequence |
| allele | HLA allele used for the binding prediction |
| Peptide_Length | Number of amino acids in the peptide |
| Mutation_Position_in_Peptide | Position of the mutation within the peptide |
| mhcflurry_affinity | Predicted mutant peptide-HLA binding affinity, in nM |
| mhcflurry_affinity_percentile | Predicted mutant binding percentile |
| WT_Binding_Affinity | Predicted wild-type peptide-HLA affinity, in nM |
| WT_Binding_Percentile | Predicted wild-type binding percentile |

## 3. Derived measurements

| Field | Meaning |
|---|---|
| affinity_fold_improvement | WT affinity divided by mutant affinity |
| percentile_fold_improvement | WT percentile divided by mutant percentile |
| Priority | Category assigned using explicit decision rules |
| RF_score | Random Forest model score |
| XGB_score | XGBoost model score |
| LR_score | Logistic Regression model score |
| AI_score | Weighted ensemble of the three model scores |
| Final_Rank_Score | Combined priority-weighted and AI-assisted ranking score |

Column names may differ between intermediate files. The executable scripts
in `src/` define the actual input and output schemas.

## 4. Prediction matching

WT predictions are matched using the WT peptide and HLA allele. Identical
peptide-HLA keys can be reused across multiple candidate records. Therefore,
the number of candidate rows can exceed the number of unique WT predictions.

## 5. Output files

| File | Purpose |
|---|---|
| primary_binding_candidates.csv | Candidates passing the primary binding filter |
| wt_aware_candidates.csv | Candidate rows with matched WT predictions |
| mutant_wt_comparison.csv | Mutant-versus-WT measurements and ratios |
| prioritized_candidates.csv | Candidates with assigned priority categories |
| ai_ranked_candidates.csv | AI-assisted candidate ranking and model scores |
| feature_importance.csv | Random Forest feature-importance summary |
| unique_p1_genes.csv | Unique genes among Priority 1 demo candidates |
| unique_p1_mutations.csv | Unique mutations among Priority 1 demo candidates |
| pathway_enrichment.csv | Pathway-analysis output or an explicit skip result |

## 6. Interpretation

Predicted binding, fold-improvement ratios, priority categories, model scores,
and rankings are computational outputs. They do not independently establish
T-cell recognition, immunogenicity, vaccine efficacy, or clinical utility.
