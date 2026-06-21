#!/usr/bin/env python3
"""Filter raw Ebola records to complete VP35 protein sequences and audit them."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
AUDIT_DIR = ROOT / "data" / "audit"
TABLES_DIR = ROOT / "results" / "tables"

RAW_PROTEIN_FASTA = RAW_DIR / "raw_vp35_sequences.fasta"
RAW_METADATA_CSV = RAW_DIR / "raw_vp35_metadata.csv"
RAW_CODING_FASTA = RAW_DIR / "raw_ebola_coding_region_sequences.fasta"
REFERENCE_FASTA = RAW_DIR / "reference_vp35_uniprot.fasta"

if not RAW_PROTEIN_FASTA.exists():
    RAW_PROTEIN_FASTA = ROOT / "raw_vp35_sequences.fasta"
if not RAW_METADATA_CSV.exists():
    RAW_METADATA_CSV = ROOT / "raw_vp35_metadata.csv"
if not RAW_CODING_FASTA.exists():
    RAW_CODING_FASTA = ROOT / "raw_ebola_coding_region_sequences.fasta"
if not REFERENCE_FASTA.exists():
    REFERENCE_FASTA = ROOT / "reference_vp35_uniprot.fasta"


KEEP_TERMS = ("polymerase cofactor vp35", " vp35", "|vp35", "structural protein vp35")
REMOVE_TERMS = (
    "nucleoprotein",
    "matrix protein vp40",
    "glycoprotein",
    "vp30",
    "vp24",
    "rna-directed rna polymerase",
    "polymerase l",
)
BAD_QUALITY_TERMS = ("partial", "fragment", "low-quality", "incomplete")
EXPECTED_LENGTH = 340
MIN_COMPLETE_LENGTH = 330
MAX_COMPLETE_LENGTH = 350


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def accession(self) -> str:
        token = self.header.split()[0]
        return token.split(".")[0].split(":")[0]


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
            seq = record.sequence
            for start in range(0, len(seq), 70):
                handle.write(seq[start : start + 70] + "\n")


def is_vp35_text(text: str) -> bool:
    lower = f" {text.lower()}"
    has_keep_term = any(term in lower for term in KEEP_TERMS)
    has_remove_term = any(term in lower for term in REMOVE_TERMS)
    return has_keep_term and not has_remove_term


def removal_reason(record: FastaRecord) -> str:
    text = record.header.lower()
    if not is_vp35_text(record.header):
        return "non_vp35_header"
    if any(term in text for term in BAD_QUALITY_TERMS):
        return "quality_label"
    if len(record.sequence) < MIN_COMPLETE_LENGTH:
        return "too_short"
    if len(record.sequence) > MAX_COMPLETE_LENGTH:
        return "too_long"
    if "*" in record.sequence:
        return "stop_codon"
    if record.sequence.count("X") > 0:
        return "contains_unknown_amino_acid"
    return "kept"


def read_metadata(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def metadata_is_vp35(row: dict[str, str]) -> bool:
    text = " ".join([row.get("Protein", ""), row.get("GenBank_Title", "")])
    return is_vp35_text(text) and row.get("Length", "").strip() == str(EXPECTED_LENGTH)


def parse_year(value: str) -> str:
    value = (value or "").strip()
    return value[:4] if len(value) >= 4 and value[:4].isdigit() else ""


def shannon_entropy(counts: Counter[str], total: int) -> float:
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy


def mutation_label(reference: str, observed: str, position: int) -> str:
    return f"{reference}{position}{observed}"


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    protein_records = read_fasta(RAW_PROTEIN_FASTA)
    coding_records = read_fasta(RAW_CODING_FASTA) if RAW_CODING_FASTA.exists() else []
    reference_records = read_fasta(REFERENCE_FASTA)
    if len(reference_records) != 1:
        raise ValueError(f"Expected exactly one reference sequence in {REFERENCE_FASTA}")
    reference = reference_records[0].sequence

    kept_records: list[FastaRecord] = []
    removed_rows: list[dict[str, object]] = []
    reason_counts: Counter[str] = Counter()
    for record in protein_records:
        reason = removal_reason(record)
        reason_counts[reason] += 1
        if reason == "kept":
            kept_records.append(record)
        else:
            removed_rows.append(
                {
                    "accession": record.accession,
                    "length": len(record.sequence),
                    "reason": reason,
                    "header": record.header,
                }
            )

    kept_accessions = {record.accession for record in kept_records}
    metadata_fields, metadata_rows = read_metadata(RAW_METADATA_CSV)
    vp35_metadata_rows = [row for row in metadata_rows if metadata_is_vp35(row)]
    filtered_metadata_rows = [
        row for row in vp35_metadata_rows if row.get("Accession", "").split(".")[0] in kept_accessions
    ]
    metadata_only_rows = [
        row
        for row in vp35_metadata_rows
        if row.get("Accession", "").split(".")[0] not in kept_accessions
    ]

    protein_only_accessions = kept_accessions - {
        row.get("Accession", "").split(".")[0] for row in filtered_metadata_rows
    }

    coding_vp35_records = [
        record for record in coding_records if is_vp35_text(record.header) and len(record.sequence) == EXPECTED_LENGTH * 3 + 3
    ]
    coding_seq_to_records: dict[str, list[FastaRecord]] = defaultdict(list)
    for record in coding_vp35_records:
        coding_seq_to_records[record.sequence].append(record)

    seq_to_records: dict[str, list[FastaRecord]] = defaultdict(list)
    for record in kept_records:
        seq_to_records[record.sequence].append(record)

    unique_records: list[FastaRecord] = []
    haplotype_rows: list[dict[str, object]] = []
    for index, (sequence, records) in enumerate(
        sorted(seq_to_records.items(), key=lambda item: (-len(item[1]), item[1][0].accession)),
        start=1,
    ):
        haplotype_id = f"VP35_H{index:03d}"
        mutations = [
            mutation_label(ref_aa, obs_aa, pos)
            for pos, (ref_aa, obs_aa) in enumerate(zip(reference, sequence), start=1)
            if ref_aa != obs_aa
        ]
        accessions = sorted(record.accession for record in records)
        unique_records.append(
            FastaRecord(
                f"{haplotype_id} count={len(records)} accessions={','.join(accessions)} mutations={','.join(mutations) or 'reference'}",
                sequence,
            )
        )
        haplotype_rows.append(
            {
                "haplotype_id": haplotype_id,
                "record_count": len(records),
                "sequence_length": len(sequence),
                "mutations_vs_reference": ";".join(mutations) if mutations else "reference",
                "accessions": ";".join(accessions),
            }
        )

    mutation_counter: Counter[tuple[int, str, str]] = Counter()
    for record in kept_records:
        for pos, (ref_aa, obs_aa) in enumerate(zip(reference, record.sequence), start=1):
            if ref_aa != obs_aa:
                mutation_counter[(pos, ref_aa, obs_aa)] += 1

    mutation_rows: list[dict[str, object]] = []
    for (pos, ref_aa, obs_aa), count in sorted(mutation_counter.items()):
        mutation_rows.append(
            {
                "position": pos,
                "reference_aa": ref_aa,
                "observed_aa": obs_aa,
                "mutation": mutation_label(ref_aa, obs_aa, pos),
                "record_count": count,
                "frequency": f"{count / len(kept_records):.6f}",
            }
        )

    conservation_rows: list[dict[str, object]] = []
    amino_acid_frequency_rows: list[dict[str, object]] = []
    entropy_rows: list[dict[str, object]] = []
    for pos in range(1, len(reference) + 1):
        counts = Counter(record.sequence[pos - 1] for record in kept_records if len(record.sequence) >= pos)
        total = sum(counts.values())
        ref_count = counts.get(reference[pos - 1], 0)
        most_common_aa, most_common_count = counts.most_common(1)[0]
        entropy = shannon_entropy(counts, total)
        conservation_rows.append(
            {
                "position": pos,
                "reference_aa": reference[pos - 1],
                "most_common_aa": most_common_aa,
                "most_common_frequency": f"{most_common_count / total:.6f}",
                "reference_frequency": f"{ref_count / total:.6f}",
                "unique_amino_acids": len(counts),
                "shannon_entropy": f"{entropy:.6f}",
                "amino_acid_counts": ";".join(f"{aa}:{counts[aa]}" for aa in sorted(counts)),
            }
        )
        entropy_rows.append(
            {
                "position": pos,
                "reference_aa": reference[pos - 1],
                "shannon_entropy": f"{entropy:.6f}",
            }
        )
        for aa in sorted(counts):
            amino_acid_frequency_rows.append(
                {
                    "position": pos,
                    "reference_aa": reference[pos - 1],
                    "amino_acid": aa,
                    "record_count": counts[aa],
                    "frequency": f"{counts[aa] / total:.6f}",
                    "is_reference_aa": aa == reference[pos - 1],
                }
            )

    length_counts = Counter(len(record.sequence) for record in protein_records)
    length_rows = [
        {"length": length, "record_count": count}
        for length, count in sorted(length_counts.items())
    ]

    country_counts = Counter((row.get("Country") or "missing").strip() or "missing" for row in filtered_metadata_rows)
    year_counts = Counter(parse_year(row.get("Collection_Date", "")) or "missing" for row in filtered_metadata_rows)
    host_counts = Counter((row.get("Host") or "missing").strip() or "missing" for row in filtered_metadata_rows)
    species_counts = Counter((row.get("Species") or "missing").strip() or "missing" for row in filtered_metadata_rows)

    write_fasta(kept_records, PROCESSED_DIR / "vp35_filtered_sequences.fasta")
    write_fasta(unique_records, PROCESSED_DIR / "vp35_unique_protein_sequences.fasta")
    write_fasta(coding_vp35_records, PROCESSED_DIR / "vp35_filtered_coding_regions.fasta")
    unique_coding_records: list[FastaRecord] = []
    coding_summary_rows: list[dict[str, object]] = []
    for index, (sequence, records) in enumerate(
        sorted(coding_seq_to_records.items(), key=lambda item: (-len(item[1]), item[1][0].accession)),
        start=1,
    ):
        coding_id = f"VP35_CDS_H{index:03d}"
        accessions = sorted(record.accession for record in records)
        unique_coding_records.append(
            FastaRecord(
                f"{coding_id} count={len(records)} accessions={','.join(accessions)}",
                sequence,
            )
        )
        coding_summary_rows.append(
            {
                "coding_haplotype_id": coding_id,
                "record_count": len(records),
                "sequence_length_nt": len(sequence),
                "accessions": ";".join(accessions),
            }
        )
    write_fasta(unique_coding_records, PROCESSED_DIR / "vp35_unique_coding_regions.fasta")
    write_csv(PROCESSED_DIR / "vp35_filtered_metadata.csv", metadata_fields, filtered_metadata_rows)
    write_csv(AUDIT_DIR / "removed_records.csv", ["accession", "length", "reason", "header"], removed_rows)
    write_csv(
        AUDIT_DIR / "metadata_only_vp35_records.csv",
        metadata_fields,
        metadata_only_rows,
    )
    write_csv(AUDIT_DIR / "sequence_length_summary.csv", ["length", "record_count"], length_rows)
    write_csv(TABLES_DIR / "vp35_haplotype_summary.csv", list(haplotype_rows[0].keys()), haplotype_rows)
    write_csv(
        TABLES_DIR / "vp35_mutation_table.csv",
        ["position", "reference_aa", "observed_aa", "mutation", "record_count", "frequency"],
        mutation_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_position_conservation.csv",
        [
            "position",
            "reference_aa",
            "most_common_aa",
            "most_common_frequency",
            "reference_frequency",
            "unique_amino_acids",
            "shannon_entropy",
            "amino_acid_counts",
        ],
        conservation_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_amino_acid_frequencies.csv",
        ["position", "reference_aa", "amino_acid", "record_count", "frequency", "is_reference_aa"],
        amino_acid_frequency_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_entropy_by_position.csv",
        ["position", "reference_aa", "shannon_entropy"],
        entropy_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_coding_haplotype_summary.csv",
        ["coding_haplotype_id", "record_count", "sequence_length_nt", "accessions"],
        coding_summary_rows,
    )
    write_csv(
        TABLES_DIR / "vp35_record_counts_by_country.csv",
        ["country", "record_count"],
        [{"country": key, "record_count": value} for key, value in sorted(country_counts.items())],
    )
    write_csv(
        TABLES_DIR / "vp35_record_counts_by_year.csv",
        ["year", "record_count"],
        [{"year": key, "record_count": value} for key, value in sorted(year_counts.items())],
    )
    write_csv(
        TABLES_DIR / "vp35_record_counts_by_host.csv",
        ["host", "record_count"],
        [{"host": key, "record_count": value} for key, value in sorted(host_counts.items())],
    )

    report = f"""# VP35 Step 1 Data Audit Report

## Inputs

- Raw protein FASTA: `{RAW_PROTEIN_FASTA.relative_to(ROOT)}`
- Raw metadata CSV: `{RAW_METADATA_CSV.relative_to(ROOT)}`
- Raw coding-region FASTA: `{RAW_CODING_FASTA.relative_to(ROOT) if RAW_CODING_FASTA.exists() else 'not found'}`
- UniProt reference FASTA: `{REFERENCE_FASTA.relative_to(ROOT)}`

## Protein FASTA filtering

- Total raw protein FASTA records: {len(protein_records)}
- Complete VP35 protein records retained: {len(kept_records)}
- Removed raw protein records: {len(removed_rows)}
- Unique complete VP35 amino-acid sequences: {len(seq_to_records)}
- Duplicate complete VP35 protein records: {len(kept_records) - len(seq_to_records)}
- Expected complete VP35 length used for quality filtering: {EXPECTED_LENGTH} aa
- Accepted length window: {MIN_COMPLETE_LENGTH}-{MAX_COMPLETE_LENGTH} aa

### Removal reasons

| Reason | Records |
|---|---:|
"""
    for reason, count in sorted(reason_counts.items()):
        if reason != "kept":
            report += f"| {reason} | {count} |\n"

    report += f"""
## Metadata matching

- Metadata rows: {len(metadata_rows)}
- Metadata rows labeled as complete VP35: {len(vp35_metadata_rows)}
- VP35 metadata rows matching retained FASTA records: {len(filtered_metadata_rows)}
- VP35 metadata rows without matching retained FASTA sequence: {len(metadata_only_rows)}
- Retained FASTA records without matching metadata row: {len(protein_only_accessions)}

Metadata-only VP35 accessions: {", ".join(row.get("Accession", "") for row in metadata_only_rows) or "none"}

Protein-only retained FASTA accessions: {", ".join(sorted(protein_only_accessions)) or "none"}

### Metadata-only VP35 records

`S32584` is present in the metadata as a 340 aa `structural protein VP35` record, but the accession is not present in the raw protein FASTA or the raw coding-region FASTA provided with this project. It was therefore not removed by the VP35 filter; it is excluded from sequence-based analyses because no FASTA sequence is available in the raw inputs. The record is preserved in `data/audit/metadata_only_vp35_records.csv` so the metadata/sequence mismatch remains documented.

## Retained VP35 metadata summary

- Countries represented: {len(country_counts)}
- Collection years represented, including missing: {len(year_counts)}
- Host labels represented: {len(host_counts)}
- Species labels represented: {len(species_counts)}
- Records with missing collection date: {year_counts.get("missing", 0)}
- Records with missing country/location: {country_counts.get("missing", 0)}

## Coding-region VP35 records

- VP35 coding-region records retained: {len(coding_vp35_records)}
- Unique VP35 coding-region sequences: {len(coding_seq_to_records)}

## Mutation summary vs UniProt reference

- Reference length: {len(reference)} aa
- Positions with at least one observed amino-acid difference: {len({pos for pos, _, _ in mutation_counter})}
- Observed mutation labels: {", ".join(row["mutation"] for row in mutation_rows) or "none"}

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
"""
    (AUDIT_DIR / "data_audit_report.md").write_text(report)

    print(f"Retained complete VP35 protein records: {len(kept_records)}")
    print(f"Unique VP35 protein haplotypes: {len(seq_to_records)}")
    print(f"Observed mutations: {', '.join(row['mutation'] for row in mutation_rows) or 'none'}")
    print(f"Audit report: {AUDIT_DIR / 'data_audit_report.md'}")


if __name__ == "__main__":
    main()
