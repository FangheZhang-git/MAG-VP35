#!/usr/bin/env python3
"""Build a cross-species VP35 homolog dataset while preserving species labels."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
BY_SPECIES_DIR = PROCESSED_DIR / "by_species"
COMBINED_DIR = PROCESSED_DIR / "combined"
AUDIT_DIR = ROOT / "data" / "audit"
TABLES_DIR = ROOT / "results" / "tables"

MIN_COMPLETE_LENGTH = 320
MAX_COMPLETE_LENGTH = 360
BAD_QUALITY_TERMS = ("partial", "fragment", "low-quality", "incomplete")

SPECIES = {
    "zaire": "Zaire ebolavirus",
    "sudan": "Sudan ebolavirus",
    "bundibugyo": "Bundibugyo ebolavirus",
    "tai_forest": "Tai Forest ebolavirus",
    "reston": "Reston ebolavirus",
    "bombali": "Bombali ebolavirus",
}

VP35_TERMS = (
    "vp35",
    "polymerase cofactor vp35",
    "structural protein vp35",
)


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def accession(self) -> str:
        token = self.header.split()[0]
        return token.split(".")[0].split(":")[0]


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def species_paths(species: str) -> tuple[Path | None, Path | None, Path | None]:
    species_dir = RAW_DIR / species
    protein_path = first_existing(
        [
            species_dir / "raw_protein_sequences.fasta",
            species_dir / "raw_vp35_sequences.fasta",
        ]
    )
    metadata_path = first_existing(
        [
            species_dir / "raw_metadata.csv",
            species_dir / "raw_vp35_metadata.csv",
        ]
    )
    coding_path = first_existing(
        [
            species_dir / "raw_coding_regions.fasta",
            species_dir / "raw_ebola_coding_region_sequences.fasta",
        ]
    )
    return protein_path, metadata_path, coding_path


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    seq_lines: list[str] = []

    with path.open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(seq_lines).upper()))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.replace(" ", ""))

    if header is not None:
        records.append(FastaRecord(header, "".join(seq_lines).upper()))
    return records


def write_fasta(records: list[FastaRecord], path: Path) -> None:
    with path.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for start in range(0, len(record.sequence), 70):
                handle.write(record.sequence[start : start + 70] + "\n")


def read_metadata(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def is_vp35_text(text: str) -> bool:
    lower = text.lower()
    return any(term in lower for term in VP35_TERMS)


def sequence_filter_reason(record: FastaRecord) -> str:
    if not is_vp35_text(record.header):
        return "non_vp35_header"
    lower = record.header.lower()
    if any(term in lower for term in BAD_QUALITY_TERMS):
        return "quality_label"
    if len(record.sequence) < MIN_COMPLETE_LENGTH:
        return "too_short"
    if len(record.sequence) > MAX_COMPLETE_LENGTH:
        return "too_long"
    if "*" in record.sequence:
        return "stop_codon"
    if "X" in record.sequence:
        return "contains_unknown_amino_acid"
    return "kept"


def metadata_is_vp35(row: dict[str, str]) -> bool:
    text = " ".join([row.get("Protein", ""), row.get("GenBank_Title", "")])
    if not is_vp35_text(text):
        return False
    length = row.get("Length", "").strip()
    if not length.isdigit():
        return False
    return MIN_COMPLETE_LENGTH <= int(length) <= MAX_COMPLETE_LENGTH


def with_species_header(species: str, record: FastaRecord) -> FastaRecord:
    return FastaRecord(
        header=f"species={species} accession={record.accession} {record.header}",
        sequence=record.sequence,
    )


def summarize_lengths(species: str, records: list[FastaRecord]) -> list[dict[str, object]]:
    counts = Counter(len(record.sequence) for record in records if is_vp35_text(record.header))
    return [
        {
            "species": species,
            "length": length,
            "record_count": count,
            "within_initial_window": MIN_COMPLETE_LENGTH <= length <= MAX_COMPLETE_LENGTH,
        }
        for length, count in sorted(counts.items())
    ]


def main() -> None:
    BY_SPECIES_DIR.mkdir(parents=True, exist_ok=True)
    COMBINED_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    audit_rows: list[dict[str, object]] = []
    haplotype_rows: list[dict[str, object]] = []
    length_rows: list[dict[str, object]] = []
    all_filtered_records: list[FastaRecord] = []
    all_sequence_to_records: dict[str, list[tuple[str, FastaRecord]]] = defaultdict(list)

    for species, species_label in SPECIES.items():
        protein_path, metadata_path, coding_path = species_paths(species)
        if protein_path is None:
            raise FileNotFoundError(f"No raw protein FASTA found for {species}")
        if metadata_path is None:
            raise FileNotFoundError(f"No raw metadata CSV found for {species}")

        protein_records = read_fasta(protein_path)
        length_rows.extend(summarize_lengths(species, protein_records))

        reason_counts: Counter[str] = Counter()
        filtered_records: list[FastaRecord] = []
        for record in protein_records:
            reason = sequence_filter_reason(record)
            reason_counts[reason] += 1
            if reason == "kept":
                filtered_records.append(record)

        kept_accessions = {record.accession for record in filtered_records}
        metadata_fields, metadata_rows = read_metadata(metadata_path)
        vp35_metadata_rows = [row for row in metadata_rows if metadata_is_vp35(row)]
        filtered_metadata_rows = [
            {"species": species, "species_label": species_label, **row}
            for row in vp35_metadata_rows
            if row.get("Accession", "").split(".")[0] in kept_accessions
        ]
        metadata_only_count = len(vp35_metadata_rows) - len(filtered_metadata_rows)
        metadata_accessions = {
            row.get("Accession", "").split(".")[0] for row in filtered_metadata_rows
        }
        protein_only_count = len(kept_accessions - metadata_accessions)

        species_sequence_to_records: dict[str, list[FastaRecord]] = defaultdict(list)
        for record in filtered_records:
            species_sequence_to_records[record.sequence].append(record)
            all_sequence_to_records[record.sequence].append((species, record))

        species_filtered_records = [with_species_header(species, record) for record in filtered_records]
        all_filtered_records.extend(species_filtered_records)

        unique_records: list[FastaRecord] = []
        for index, (sequence, records) in enumerate(
            sorted(
                species_sequence_to_records.items(),
                key=lambda item: (-len(item[1]), item[1][0].accession),
            ),
            start=1,
        ):
            haplotype_id = f"{species.upper()}_VP35_H{index:03d}"
            accessions = sorted(record.accession for record in records)
            unique_records.append(
                FastaRecord(
                    header=(
                        f"{haplotype_id} species={species} species_label=\"{species_label}\" "
                        f"count={len(records)} accessions={','.join(accessions)}"
                    ),
                    sequence=sequence,
                )
            )
            haplotype_rows.append(
                {
                    "species": species,
                    "species_label": species_label,
                    "haplotype_id": haplotype_id,
                    "record_count": len(records),
                    "sequence_length": len(sequence),
                    "accessions": ";".join(accessions),
                }
            )

        write_fasta(
            species_filtered_records,
            BY_SPECIES_DIR / f"{species}_vp35_filtered_sequences.fasta",
        )
        write_fasta(
            unique_records,
            BY_SPECIES_DIR / f"{species}_vp35_unique_protein_sequences.fasta",
        )
        write_csv(
            BY_SPECIES_DIR / f"{species}_vp35_filtered_metadata.csv",
            ["species", "species_label", *metadata_fields],
            filtered_metadata_rows,
        )

        audit_rows.append(
            {
                "species": species,
                "species_label": species_label,
                "protein_fasta": protein_path.relative_to(ROOT),
                "metadata_csv": metadata_path.relative_to(ROOT),
                "coding_fasta": coding_path.relative_to(ROOT) if coding_path else "not_found",
                "raw_protein_records": len(protein_records),
                "vp35_header_records": sum(
                    count for reason, count in reason_counts.items() if reason != "non_vp35_header"
                ),
                "filtered_complete_vp35_records": len(filtered_records),
                "unique_complete_vp35_sequences": len(species_sequence_to_records),
                "metadata_rows": len(metadata_rows),
                "vp35_metadata_rows_in_length_window": len(vp35_metadata_rows),
                "filtered_metadata_rows": len(filtered_metadata_rows),
                "metadata_only_vp35_rows": metadata_only_count,
                "protein_only_vp35_records": protein_only_count,
                "too_short_records": reason_counts.get("too_short", 0),
                "too_long_records": reason_counts.get("too_long", 0),
                "quality_label_records": reason_counts.get("quality_label", 0),
                "stop_codon_records": reason_counts.get("stop_codon", 0),
                "unknown_amino_acid_records": reason_counts.get("contains_unknown_amino_acid", 0),
            }
        )

    combined_unique_records: list[FastaRecord] = []
    for index, (sequence, species_records) in enumerate(
        sorted(
            all_sequence_to_records.items(),
            key=lambda item: (-len(item[1]), sorted({species for species, _ in item[1]}), item[1][0][1].accession),
        ),
        start=1,
    ):
        haplotype_id = f"VP35_ALL_H{index:03d}"
        species_counts = Counter(species for species, _ in species_records)
        species_summary = ";".join(f"{species}:{species_counts[species]}" for species in sorted(species_counts))
        accessions = sorted(f"{species}:{record.accession}" for species, record in species_records)
        combined_unique_records.append(
            FastaRecord(
                header=(
                    f"{haplotype_id} count={len(species_records)} species_counts={species_summary} "
                    f"accessions={','.join(accessions)}"
                ),
                sequence=sequence,
            )
        )

    write_fasta(all_filtered_records, COMBINED_DIR / "vp35_all_species_filtered_sequences.fasta")
    write_fasta(combined_unique_records, COMBINED_DIR / "vp35_all_species_unique_sequences.fasta")
    write_csv(TABLES_DIR / "vp35_species_expansion_audit.csv", list(audit_rows[0].keys()), audit_rows)
    write_csv(
        TABLES_DIR / "vp35_species_haplotype_summary.csv",
        ["species", "species_label", "haplotype_id", "record_count", "sequence_length", "accessions"],
        haplotype_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_cross_species_length_summary.csv",
        ["species", "length", "record_count", "within_initial_window"],
        length_rows,
    )

    total_filtered = sum(int(row["filtered_complete_vp35_records"]) for row in audit_rows)
    total_unique = len(combined_unique_records)
    report = f"""# VP35 Species Expansion Audit Report

## Purpose

This dataset expands the initial Zaire ebolavirus VP35 set into a broader ebolavirus VP35 homolog set while preserving species labels. Records are processed species-by-species before combined FASTA outputs are created.

## Filtering Rules

- Protein FASTA records are retained when their headers identify VP35, polymerase cofactor VP35, or structural protein VP35.
- Broad exclusion terms such as generic protein labels or L protein labels are not used.
- Complete homolog candidates are retained within an initial {MIN_COMPLETE_LENGTH}-{MAX_COMPLETE_LENGTH} aa length window.
- Records with partial, fragment, low-quality, incomplete, stop-codon, or unknown-amino-acid signals are excluded from the cleaned protein FASTA outputs.
- Raw files are not modified.

## Species Comparison

| Species | Raw protein records | Complete VP35 records | Unique complete VP35 sequences |
|---|---:|---:|---:|
"""
    for row in audit_rows:
        report += (
            f"| {row['species']} | {row['raw_protein_records']} | "
            f"{row['filtered_complete_vp35_records']} | {row['unique_complete_vp35_sequences']} |\n"
        )

    report += f"""
## Combined Dataset

- Filtered complete VP35 homolog records across all species: {total_filtered}
- Unique amino-acid sequences across all species: {total_unique}

## Biological Interpretation

Within Zaire ebolavirus, amino-acid differences can be discussed as natural variation within the EBOV/Zaire dataset. Across ebolavirus species, amino-acid differences should be described as homologous amino-acid differences, not mutation frequency. The broader species dataset is intended for cross-species conservation analysis and representative VP35 ensemble selection.

## Generated Outputs

- `data/processed/by_species/{{species}}_vp35_filtered_sequences.fasta`
- `data/processed/by_species/{{species}}_vp35_filtered_metadata.csv`
- `data/processed/by_species/{{species}}_vp35_unique_protein_sequences.fasta`
- `data/processed/combined/vp35_all_species_filtered_sequences.fasta`
- `data/processed/combined/vp35_all_species_unique_sequences.fasta`
- `results/tables/vp35_species_expansion_audit.csv`
- `results/tables/vp35_species_haplotype_summary.csv`
- `results/tables/vp35_cross_species_length_summary.csv`
"""
    (AUDIT_DIR / "vp35_species_expansion_audit_report.md").write_text(report)

    print(f"Filtered complete VP35 homolog records: {total_filtered}")
    print(f"Unique cross-species VP35 homolog sequences: {total_unique}")
    print(f"Audit report: {AUDIT_DIR / 'vp35_species_expansion_audit_report.md'}")


if __name__ == "__main__":
    main()
