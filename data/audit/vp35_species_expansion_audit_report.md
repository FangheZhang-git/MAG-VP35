# VP35 Species Expansion Audit Report

## Purpose

This dataset expands the initial Zaire ebolavirus VP35 set into a broader ebolavirus VP35 homolog set while preserving species labels. Records are processed species-by-species before combined FASTA outputs are created.

## Filtering Rules

- Protein FASTA records are retained when their headers identify VP35, polymerase cofactor VP35, or structural protein VP35.
- Broad exclusion terms such as generic protein labels or L protein labels are not used.
- Complete homolog candidates are retained within an initial 320-360 aa length window.
- Records with partial, fragment, low-quality, incomplete, stop-codon, or unknown-amino-acid signals are excluded from the cleaned protein FASTA outputs.
- Raw files are not modified.

## Species Comparison

| Species | Raw protein records | Complete VP35 records | Unique complete VP35 sequences |
|---|---:|---:|---:|
| zaire | 1008 | 110 | 4 |
| sudan | 1316 | 13 | 6 |
| bundibugyo | 280 | 7 | 2 |
| tai_forest | 42 | 2 | 1 |
| reston | 239 | 22 | 6 |
| bombali | 94 | 3 | 1 |

## Combined Dataset

- Filtered complete VP35 homolog records across all species: 157
- Unique amino-acid sequences across all species: 20

## Biological Interpretation

Within Zaire ebolavirus, amino-acid differences can be discussed as natural variation within the EBOV/Zaire dataset. Across ebolavirus species, amino-acid differences should be described as homologous amino-acid differences, not mutation frequency. The broader species dataset is intended for cross-species conservation analysis and representative VP35 ensemble selection.

## Generated Outputs

- `data/processed/by_species/{species}_vp35_filtered_sequences.fasta`
- `data/processed/by_species/{species}_vp35_filtered_metadata.csv`
- `data/processed/by_species/{species}_vp35_unique_protein_sequences.fasta`
- `data/processed/combined/vp35_all_species_filtered_sequences.fasta`
- `data/processed/combined/vp35_all_species_unique_sequences.fasta`
- `results/tables/vp35_species_expansion_audit.csv`
- `results/tables/vp35_species_haplotype_summary.csv`
- `results/tables/vp35_cross_species_length_summary.csv`
