# AGENT.md

## Project Name

**Mutation-Aware Multi-Objective Fragment Growing for the Design of Oral-Like Ebola VP35 Inhibitors**

## Core Mission

This project develops a computational drug-design framework to generate and prioritize oral-like small molecules predicted to bind robustly across naturally occurring Ebola VP35 variants.

The project is not trying to prove that a molecule cures Ebola. The goal is to build and evaluate a rigorous **mutation-aware molecular optimization framework** that can prioritize VP35 inhibitor candidates for future experimental testing.

## Main Research Question

Can a mutation-aware, multi-objective fragment-growing framework generate oral-like Ebola VP35 inhibitor candidates that remain more stable across diverse viral variants than molecules optimized using conventional single-structure docking?

## Scientific Logic

Ordinary docking projects often optimize molecules against only one protein structure. This can fail when the target protein mutates. Ebola VP35 naturally varies across viral samples, outbreaks, geographic regions, and species. Therefore, this project first builds a clean map of real VP35 sequence variation, then uses that variation to construct representative VP35 structural variants for docking and molecular optimization.

The project combines:

1. Ebola VP35 sequence collection and cleaning.
2. VP35 mutation and conservation analysis.
3. VP35 binding-pocket identification using PDB structures with bound ligands.
4. Representative VP35 variant ensemble construction.
5. Multi-variant docking.
6. SMILES-based fragment growing.
7. Multi-objective scoring for binding, mutation robustness, ligand efficiency, oral-like properties, and synthetic feasibility.
8. Comparison against control optimization methods.
9. Held-out VP35 variant testing.

---

# Current Stage

The project is currently in **Step 1: data collection and cleaning**.

The immediate goal is to collect raw Ebola VP35-related sequence data and metadata, then filter the raw files to keep only **polymerase cofactor VP35** records.

At this stage, do **not** start docking, SMILES generation, fragment growing, or ChEMBL screening yet.

---

# Key Data Sources and Their Roles

## 1. NCBI Virus / GenBank

**Main use:** collect many real Ebola virus protein sequences and metadata.

NCBI Virus is the main source for the raw sequence dataset. It provides protein FASTA files and metadata tables containing accession IDs, protein names, organism names, species, sequence lengths, host, geographic location, collection date, release date, isolate information, and related identifiers.

Use NCBI Virus to collect:

- Raw Ebola protein FASTA.
- Raw Ebola metadata CSV.
- Optional coding-region FASTA.

NCBI is the main data source for discovering naturally occurring VP35 variation.

## 2. UniProt

**Main use:** collect one clean reference VP35 sequence.

UniProt should not be the main source for thousands of VP35 variants. Instead, it should be used to confirm the standard reference VP35 protein identity, official name, normal length, domains, function, and cross-links.

Current target reference:

- UniProt entry: `Q05127`
- Entry name: `VP35_EBOZM`
- Protein name: `Polymerase cofactor VP35`
- Expected length: approximately 340 amino acids

Save this as:

```text
reference_vp35_uniprot.fasta
```

This reference sequence will later be used to describe mutations relative to a standard VP35 sequence, such as `S210N`.

## 3. PDB / RCSB Protein Data Bank

**Main use:** collect experimentally determined VP35 structures and bound small molecules.

PDB is not the main database for Step 1 sequence collection. It becomes important in Step 2, when the project needs VP35 3D structures to identify the small-molecule binding pocket and map sequence mutations onto structure.

Use PDB later to collect:

- VP35 experimental structures.
- VP35 structures containing bound small molecules.
- Experimental ligand poses for redocking validation.
- Pocket residues around known ligands.

## 4. ChEMBL / BindingDB / PubChem

**Main use:** collect known ligands, candidate compounds, and possible starting molecules.

These databases are not needed before VP35 sequence filtering. They become useful later when selecting starting scaffolds, known ligands, decoys, or candidate molecules for virtual screening.

## 5. Fragment Libraries

**Main use:** provide small medicinal-chemistry fragments for controlled SMILES-based fragment growing.

Do not collect fragment libraries until the VP35 structure, pocket, known ligand, and starting scaffold strategy are clearer.

---

# Current Project Folder Structure

The research project folder should look like this during the current Step 1 stage:

```text
VP35_EBOLA/
├── AGENT.md
├── raw_vp35_sequences.fasta
├── raw_vp35_metadata.csv
├── reference_vp35_uniprot.fasta
└── raw_ebola_coding_region_sequences.fasta
```

Later, after filtering and analysis, the folder may expand into:

```text
VP35_EBOLA/
├── AGENT.md
├── data/
│   ├── raw/
│   │   ├── raw_vp35_sequences.fasta
│   │   ├── raw_vp35_metadata.csv
│   │   ├── reference_vp35_uniprot.fasta
│   │   └── raw_ebola_coding_region_sequences.fasta
│   ├── processed/
│   │   ├── vp35_filtered_sequences.fasta
│   │   ├── vp35_filtered_metadata.csv
│   │   ├── vp35_unique_protein_sequences.fasta
│   │   ├── vp35_unique_metadata.csv
│   │   └── vp35_filtered_coding_regions.fasta
│   └── audit/
│       ├── data_audit_report.md
│       ├── removed_records.csv
│       └── sequence_length_summary.csv
├── scripts/
│   ├── 01_filter_vp35_records.py
│   ├── 02_sequence_quality_audit.py
│   ├── 03_align_vp35_sequences.sh
│   ├── 04_calculate_conservation.py
│   └── 05_select_representative_variants.py
├── notebooks/
│   ├── 01_data_audit.ipynb
│   ├── 02_conservation_analysis.ipynb
│   └── 03_variant_selection.ipynb
├── structures/
│   ├── pdb_raw/
│   ├── pdb_prepared/
│   └── variant_models/
├── ligands/
│   ├── known_ligands/
│   ├── seed_scaffolds/
│   └── generated_molecules/
├── docking/
│   ├── redocking_validation/
│   ├── reference_structure_docking/
│   ├── multi_variant_docking/
│   └── heldout_variant_testing/
└── results/
    ├── figures/
    ├── tables/
    └── final_candidates/
```

For now, the simple flat folder is acceptable. The expanded structure can be created later when the project becomes larger.

---

# Files Already Planned or Collected

## `raw_vp35_sequences.fasta`

**Source:** NCBI Virus protein results.

**How collected:**

1. Go to NCBI Virus.
2. Search for `Ebola virus`.
3. Select the Ebola virus result.
4. Open the **Protein** tab.
5. Click **Download All Results**.
6. Choose `Sequence Data (FASTA format)`.
7. Choose `Protein`.
8. Download all records.
9. Save as `raw_vp35_sequences.fasta`.

**Important note:** this raw file may contain all Ebola proteins, not only VP35. It may include nucleoprotein, VP40, glycoprotein, VP30, VP24, L protein, and VP35. Filtering is required.

## `raw_vp35_metadata.csv`

**Source:** NCBI Virus results table.

**How collected:**

1. Go back to the same NCBI Virus Ebola protein results page.
2. Click **Download All Results**.
3. Choose `Results Table`.
4. Choose `CSV format`.
5. Select useful columns.
6. Download all records.
7. Save as `raw_vp35_metadata.csv`.

Recommended metadata columns:

- Accession
- Protein
- Nucleotide, if available
- Organism Name
- Species
- Length
- Geo Location
- Country
- Host
- Isolate
- Collection Date
- Release Date
- GenBank Title
- Genotype, if available

This metadata file can be used for both the protein FASTA and the optional coding-region FASTA at the current stage, as long as the identifiers can be matched later.

## `reference_vp35_uniprot.fasta`

**Source:** UniProt.

**How collected:**

1. Go to UniProt.
2. Search `Q05127` or `Ebola VP35 Q05127`.
3. Open the entry for `VP35_EBOZM`, `Polymerase cofactor VP35`.
4. Choose `Download`.
5. Select `FASTA (canonical)`.
6. Save as `reference_vp35_uniprot.fasta`.

**Purpose:** one clean reference VP35 sequence for mutation naming and alignment.

## `raw_ebola_coding_region_sequences.fasta`

**Source:** NCBI Virus coding-region FASTA.

**How collected:**

1. Go to the same NCBI Virus Ebola results page.
2. Click **Download All Results**.
3. Choose `Sequence Data (FASTA format)`.
4. Choose `Coding Region`.
5. Download all records.
6. Save as `raw_ebola_coding_region_sequences.fasta`.

**Purpose:** optional nucleotide/CDS-level information. This is not the main file for docking, but it can help later if the project analyzes nucleotide mutations, codon changes, or nucleotide-level diversity.

**Metadata note:** a separate coding-region metadata file is not required right now. The same `raw_vp35_metadata.csv` can be used at this stage. If later the coding-region FASTA cannot be matched cleanly to the metadata, then download a separate coding-region metadata CSV.

---

# Important Definitions

## FASTA

FASTA is a plain-text format for biological sequences.

A FASTA entry has two parts:

```text
>ARG43219.1 |polymerase cofactor VP35 [Ebola virus]
MDSIHDWTKNITDKIDQIIHDFVDKTLPDQGDNDNWWTGWRQWIPAGIGVTGVIIA...
```

The line beginning with `>` is the header. The following lines contain the sequence.

For protein FASTA, each letter is an amino acid. For coding-region FASTA, each letter is a nucleotide.

## Protein FASTA

Contains amino-acid sequences. This is the main file for VP35 mutation and conservation analysis.

## Coding-Region FASTA

Contains nucleotide coding sequences. This is optional and useful for nucleotide-level or codon-level analysis.

## Metadata

Metadata is contextual information about each sequence, such as accession ID, species, host, collection date, country, and isolate. It is essential for quality control and outbreak/geographic analysis.

## Accession ID

An accession ID is a stable database identifier for a sequence record. It is used to connect FASTA sequences with metadata rows.

---

# Immediate Next Task

The next task is:

```text
Filter the raw NCBI protein FASTA and metadata CSV to keep only polymerase cofactor VP35 records.
```

Do not manually edit the raw FASTA file. Instead, write a reproducible script.

Suggested script name:

```text
scripts/01_filter_vp35_records.py
```

Expected inputs:

```text
raw_vp35_sequences.fasta
raw_vp35_metadata.csv
```

Expected outputs:

```text
vp35_filtered_sequences.fasta
vp35_filtered_metadata.csv
removed_non_vp35_records.csv
```

---

# VP35 Filtering Rules

A record should be kept if it is clearly VP35, usually with protein name or FASTA header containing one of:

```text
polymerase cofactor VP35
VP35
```

A record should be removed if it is clearly another Ebola protein, such as:

```text
nucleoprotein
matrix protein VP40
glycoprotein
VP30
VP24
RNA-directed RNA polymerase L protein
```

After filtering by name, perform sequence-level checks.

Recommended quality filters:

1. Keep sequences with length near the expected VP35 length, approximately 340 amino acids.
2. Remove abnormally short sequences.
3. Remove abnormally long sequences.
4. Flag sequences containing many unknown amino acids such as `X`.
5. Remove records labeled partial, fragment, low-quality, or incomplete if such information is present.
6. Deduplicate identical amino-acid sequences, but preserve how many times each unique sequence appeared.

Do not destroy raw data. Place cleaned files in `processed/` or use clearly marked output filenames.

---

# Data Audit Requirements

The project must include a clear data audit. The audit should record:

- Total raw protein records downloaded.
- Number of VP35-like records after name filtering.
- Number of records removed as non-VP35.
- Sequence length distribution.
- Number of complete VP35 sequences.
- Number of incomplete or abnormal sequences removed.
- Number of unique amino-acid sequences.
- Number of duplicate sequences.
- Number of records with missing collection date.
- Number of records with missing country/location.
- Number of records by Ebola species.
- Number of records by host.
- Number of records by country/location.
- Number of records by collection year, if dates are available.

Suggested output:

```text
data_audit_report.md
sequence_length_summary.csv
removed_records.csv
vp35_record_counts_by_country.csv
vp35_record_counts_by_year.csv
```

The data audit is important because this project does not require thousands of drug-response labels. It requires a clean and well-documented VP35 variant dataset for mutation/conservation analysis and structural variant selection.

---

# Sequence Analysis Plan

After filtering and cleaning VP35 protein sequences:

1. Align the VP35 amino-acid sequences.
2. Compare all sequences to the UniProt reference VP35.
3. Calculate mutation frequency at each amino-acid position.
4. Calculate conservation at each position.
5. Calculate Shannon entropy at each position.
6. Identify naturally observed mutations.
7. Identify common VP35 haplotypes.
8. Identify rare but potentially important mutations.
9. Identify mutations that occur in or near the future binding pocket.

Possible alignment tools:

- MAFFT
- Clustal Omega
- MUSCLE

Expected outputs:

```text
vp35_alignment.fasta
vp35_mutation_table.csv
vp35_position_conservation.csv
vp35_entropy_by_position.csv
vp35_haplotype_summary.csv
```

---

# Representative Variant Selection Plan

It is unnecessary and computationally inefficient to build structures for every VP35 sequence. Instead, select approximately 10-20 representative VP35 variants.

Representative variants may include:

- The UniProt/reference VP35 sequence.
- The most common VP35 haplotypes.
- Variants from different outbreaks.
- Variants from different geographic regions.
- Variants from different Ebola species, if included in scope.
- Variants containing mutations in or near the small-molecule binding pocket.
- Several held-out variants for final testing.

Separate the representative variants into:

```text
optimization_variants/
heldout_test_variants/
```

The held-out variants must not be used during molecular optimization.

---

# Future Structure and Docking Plan

After the sequence dataset is cleaned:

1. Download VP35 structures from PDB.
2. Prioritize structures with bound small molecules.
3. Use the observed ligand to define the binding pocket.
4. Identify pocket residues within a selected distance of the ligand.
5. Map conservation and mutations onto the 3D structure.
6. Build variant structures by introducing naturally observed mutations into the reference structure.
7. Optimize mutated side chains and remove severe steric clashes.
8. Validate docking by redocking known ligands.

Docking validation metrics:

- Top-1 pose recovery.
- Top-3 pose recovery.
- Ligand RMSD compared with experimental pose.
- Consistency across random seeds.
- Consistency across docking programs, if possible.

Do not rely on docking scores without redocking validation.

---

# Future Molecular Optimization Plan

Starting molecules may include:

- Known VP35 ligands.
- Core scaffolds extracted from known ligands.
- Promising molecules from an initial virtual screen.

The fragment-growing algorithm will use SMILES-based controlled molecular modification.

Possible modifications:

- Add a small molecular fragment.
- Replace a functional group.
- Add or remove a hydrogen-bond donor.
- Add or remove a hydrogen-bond acceptor.
- Modify a linker.
- Change a ring structure.
- Perform a bioisosteric replacement.

All generated molecules must pass chemical validity filters.

Remove molecules with:

- Invalid SMILES.
- Impossible valence.
- Highly reactive groups.
- Severe chemical liabilities.
- Unrealistic structures.
- Extremely poor oral-like properties.

---

# Multi-Objective Scoring Plan

Do not optimize docking score alone.

Each molecule should be evaluated by multiple objectives:

- Average predicted binding across VP35 variants.
- Worst-case predicted binding against any VP35 variant.
- Score variation across variants.
- Ligand efficiency.
- Interactions with conserved pocket residues.
- Molecular weight.
- Lipophilicity.
- Polar surface area.
- Hydrogen-bond donors.
- Hydrogen-bond acceptors.
- Rotatable bonds.
- Predicted solubility.
- Synthetic accessibility.
- Chemical liability alerts.
- Docking-pose consistency.

Reward molecules that perform consistently across many VP35 variants. Penalize molecules that perform well only on the reference structure.

Ligand efficiency is important because docking scores often improve artificially as molecules become larger.

---

# Control Strategies for Fair Comparison

The full mutation-aware framework must be compared with controls using the same starting molecules, fragment library, generation count, and total molecule-evaluation budget.

Control methods:

1. Random fragment growing.
2. Docking-only optimization against the reference VP35 structure.
3. Wild-type multi-objective optimization against only the reference VP35 structure.
4. Full mutation-aware multi-objective framework.

The central result should test whether the full framework produces molecules that generalize better to held-out VP35 variants.

---

# Expected Figures and Tables

## Figure 1: Complete Research Workflow

Shows:

```text
Ebola sequence collection → sequence cleaning → alignment → conservation/mutation analysis → VP35 structural ensemble → seed molecules → fragment growing → multi-objective scoring → held-out variant evaluation → final candidates
```

## Figure 2: VP35 Conservation and Mutation Map

Shows VP35 structure, small-molecule binding pocket, bound ligand, conserved residues, variable residues, mutation frequencies, and Shannon entropy.

## Figure 3: Optimization Method Comparison

Compares random fragment growing, docking-only optimization, wild-type multi-objective optimization, and full mutation-aware optimization.

Possible y-axis metrics:

- Held-out robust score.
- Worst-variant docking score.
- Ligand efficiency.
- Number of acceptable oral-like molecules.

## Figure 4: Final Candidate Robustness Heat Map and Property Table

Rows are top candidate molecules. Columns are VP35 variants. Cells show docking or robust binding scores.

The property table should include:

- SMILES.
- Average binding score.
- Worst-variant score.
- Score variation.
- Ligand efficiency.
- Molecular weight.
- cLogP.
- Polar surface area.
- Synthetic accessibility.
- Oral-like property status.

---

# Reproducibility Rules

1. Never manually edit raw data files.
2. Keep raw files unchanged.
3. All filtering and cleaning must be done by scripts.
4. Record the source database, query, download date, and file name for every dataset.
5. Every generated output should be reproducible from input files and scripts.
6. Use clear filenames with step numbers.
7. Save removed records with reasons for removal.
8. Use random seeds for any stochastic selection or optimization.
9. Keep held-out variants separated from optimization variants.
10. Do not claim biological efficacy from docking alone.

---

# Recommended Naming Conventions

Use lowercase filenames with underscores.

Examples:

```text
raw_vp35_sequences.fasta
raw_vp35_metadata.csv
reference_vp35_uniprot.fasta
raw_ebola_coding_region_sequences.fasta
vp35_filtered_sequences.fasta
vp35_filtered_metadata.csv
vp35_unique_protein_sequences.fasta
vp35_mutation_table.csv
vp35_entropy_by_position.csv
vp35_representative_variants.csv
```

Use numbered scripts:

```text
01_filter_vp35_records.py
02_sequence_quality_audit.py
03_align_vp35_sequences.sh
04_calculate_conservation.py
05_select_representative_variants.py
```

---

# Assistant / Agent Behavior Rules

When helping with this project, the agent should:

1. Be direct and step-by-step.
2. Avoid jumping ahead to docking before sequence filtering is complete.
3. Explain biological and computational terms simply when needed.
4. Keep the project focused on mutation-aware VP35 small-molecule design.
5. Distinguish clearly between required and optional data.
6. Preserve raw data and recommend reproducible scripts.
7. Avoid overstating results.
8. Avoid claiming that computational candidates are real Ebola treatments.
9. Keep the central novelty clear: mutation-aware, multi-objective optimization across naturally occurring VP35 variants.
10. Prioritize clean data audit, validation, and fair controls over flashy but unsupported claims.

---

# Current Next Step Checklist

The next concrete work should be:

```text
[ ] Confirm the project folder contains the raw FASTA, metadata CSV, UniProt reference, and optional coding-region FASTA.
[ ] Create a scripts/ folder.
[ ] Write 01_filter_vp35_records.py.
[ ] Filter protein FASTA to keep only polymerase cofactor VP35.
[ ] Filter metadata CSV to match the VP35 protein records.
[ ] Count how many VP35 records remain.
[ ] Check sequence length distribution.
[ ] Remove incomplete or abnormal VP35 records.
[ ] Deduplicate identical amino-acid sequences while preserving counts.
[ ] Write the first data audit report.
```

