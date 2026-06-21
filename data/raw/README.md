# Raw VP35 Sequence Data

`data/raw/zaire/` contains the primary Zaire ebolavirus / EBOV VP35 dataset currently used for Step 1 filtering and audit work.

`data/raw/sudan/`, `data/raw/bundibugyo/`, `data/raw/tai_forest/`, `data/raw/reston/`, and `data/raw/bombali/` are placeholders for future species-specific NCBI Virus downloads.

Each species folder should eventually contain:

- protein FASTA
- metadata CSV
- optional coding-region FASTA

Raw files should never be manually edited. Any filtering, cleaning, deduplication, or sequence analysis should be performed by reproducible scripts that write outputs under `data/processed/`, `data/audit/`, or `results/`.
