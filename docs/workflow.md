# Neoantigen Prioritization Workflow

All inputs and outputs in the public demonstration are synthetic.

```mermaid
flowchart TD
    A["Synthetic demo data"] --> B["Mutant binding QC"]
    B --> C["Primary MHC-I binding filter"]
    C --> D["Match WT peptide-HLA predictions"]
    D --> E["Compare mutant and WT binding"]
    E --> F["Rule-based priority classification"]
    F --> G["AI-assisted ranking"]
    G --> H["Candidate, gene, patient and HLA summaries"]
    H --> I{"Are gene symbols synthetic placeholders?"}
    I -->|Yes| J["Skip pathway enrichment"]
    I -->|No| K["Optional enrichment for authorized real data"]
```

This workflow demonstrates software processing and prioritization. It does not establish immunogenicity, vaccine efficacy, or clinical utility.
