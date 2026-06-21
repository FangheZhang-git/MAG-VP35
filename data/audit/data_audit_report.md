# VP35 Step 1 Data Audit Report

## Inputs

- Raw protein FASTA: `data/raw/raw_vp35_sequences.fasta`
- Raw metadata CSV: `data/raw/raw_vp35_metadata.csv`
- Raw coding-region FASTA: `data/raw/raw_ebola_coding_region_sequences.fasta`
- UniProt reference FASTA: `data/raw/reference_vp35_uniprot.fasta`

## Protein FASTA filtering

- Total raw protein FASTA records: 1008
- Complete VP35 protein records retained: 110
- Removed raw protein records: 898
- Unique complete VP35 amino-acid sequences: 4
- Duplicate complete VP35 protein records: 106
- Expected complete VP35 length used for quality filtering: 340 aa
- Accepted length window: 330-350 aa

### Removal reasons

| Reason | Records |
|---|---:|
| non_vp35_header | 898 |

## Metadata matching

- Metadata rows: 1010
- Metadata rows labeled as complete VP35: 111
- VP35 metadata rows matching retained FASTA records: 110
- VP35 metadata rows without matching retained FASTA sequence: 1
- Retained FASTA records without matching metadata row: 0

Metadata-only VP35 accessions: S32584

Protein-only retained FASTA accessions: none

### Metadata-only VP35 records

`S32584` is present in the metadata as a 340 aa `structural protein VP35` record, but the accession is not present in the raw protein FASTA or the raw coding-region FASTA provided with this project. It was therefore not removed by the VP35 filter; it is excluded from sequence-based analyses because no FASTA sequence is available in the raw inputs. The record is preserved in `data/audit/metadata_only_vp35_records.csv` so the metadata/sequence mismatch remains documented.

## Retained VP35 metadata summary

- Countries represented: 5
- Collection years represented, including missing: 5
- Host labels represented: 4
- Species labels represented: 1
- Records with missing collection date: 3
- Records with missing country/location: 3

## Coding-region VP35 records

- VP35 coding-region records retained: 108
- Unique VP35 coding-region sequences: 5

## Mutation summary vs UniProt reference

- Reference length: 340 aa
- Positions with at least one observed amino-acid difference: 4
- Observed mutation labels: A12V, S41N, T68M, N204D

## Generated outputs

- `data/processed/vp35_filtered_sequences.fasta`
- `data/processed/vp35_filtered_metadata.csv`
- `data/processed/vp35_unique_protein_sequences.fasta`
- `data/processed/vp35_filtered_coding_regions.fasta`
- `data/processed/vp35_unique_coding_regions.fasta`
- `data/audit/removed_records.csv`
- `data/audit/metadata_only_vp35_records.csv`
- `data/audit/sequence_length_summary.csv`
- `results/tables/vp35_haplotype_summary.csv`
- `results/tables/vp35_coding_haplotype_summary.csv`
- `results/tables/vp35_mutation_table.csv`
- `results/tables/vp35_position_conservation.csv`
- `results/tables/vp35_amino_acid_frequencies.csv`
- `results/tables/vp35_entropy_by_position.csv`
- `results/tables/vp35_record_counts_by_country.csv`
- `results/tables/vp35_record_counts_by_year.csv`
- `results/tables/vp35_record_counts_by_host.csv`

## Interpretation

The raw protein FASTA contains many non-VP35 Ebola protein records. After filtering, the current dataset contains a compact set of complete VP35 protein sequences suitable for an initial mutation and conservation analysis. Because only a small number of unique amino-acid haplotypes are present, future representative-variant selection should either use these observed haplotypes directly or expand the dataset with additional VP35 records before building a larger structural ensemble.
