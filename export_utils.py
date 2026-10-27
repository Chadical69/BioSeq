"""Export helpers for BioSeq."""
from __future__ import annotations

import csv
import io
import json
from typing import Any


def build_csv_report(result: dict[str, Any]) -> str:
    """Create one CSV with clearly labelled sections and metric/value columns."""
    stats = result["statistics"]
    tr = result["translation"]
    output = io.StringIO(newline="")
    writer = csv.writer(output)
    writer.writerow(["section", "metric", "value"])
    writer.writerow(["sequence", "sequence_length", stats["sequence_length"]])
    writer.writerow(["sequence", "canonical_bases", stats["canonical_bases"]])
    writer.writerow(["sequence", "ambiguous_bases", stats["ambiguous_bases"]])
    writer.writerow(["sequence", "gc_content_percent", f'{stats["gc_content"]:.4f}'])
    writer.writerow(["sequence", "at_content_percent", f'{stats["at_content"]:.4f}'])
    for base in ("A", "T", "G", "C"):
        writer.writerow(["nucleotide_count", base, stats["counts"][base]])
        writer.writerow(["nucleotide_percentage", base, f'{stats["percentages"][base]:.4f}'])
    writer.writerow(["sequence_outputs", "dna_sequence", result["sequence"]])
    writer.writerow(["sequence_outputs", "complement", result["complement"]])
    writer.writerow(["sequence_outputs", "reverse_complement", result["reverse_complement"]])
    writer.writerow(["sequence_outputs", "rna_sequence", result["rna_sequence"]])
    writer.writerow(["translation", "strand_convention", result["strand_convention"]])
    writer.writerow(["translation", "reading_frame", f'+{tr["frame"]}'])
    writer.writerow(["translation", "codon_count", tr["codon_count"]])
    writer.writerow(["translation", "trailing_bases_ignored", tr["trailing_bases_ignored"]])
    writer.writerow(["translation", "protein_sequence", tr["protein_sequence"]])
    writer.writerow(["translation", "start_codon_positions_1_based", ";".join(map(str, tr["start_codon_positions"]))])
    writer.writerow(["translation", "stop_codon_positions_1_based", ";".join(map(str, tr["stop_codon_positions"]))])
    for codon, count in tr["codon_counts"].items():
        writer.writerow(["codon_frequency", codon, count])
    return output.getvalue()


def build_text_report(result: dict[str, Any]) -> str:
    stats = result["statistics"]
    tr = result["translation"]
    lines = [
        "BIOSEQ — DNA SEQUENCE ANALYSIS REPORT",
        "=" * 42,
        "",
        "SEQUENCE STATISTICS",
        f'Sequence length: {stats["sequence_length"]} bases',
        f'Canonical bases (A/T/G/C): {stats["canonical_bases"]}',
        f'Ambiguous bases: {stats["ambiguous_bases"]}',
        "",
        "Nucleotide counts and percentages (percentages exclude ambiguous bases):",
    ]
    for base in ("A", "T", "G", "C"):
        lines.append(f'  {base}: {stats["counts"][base]} ({stats["percentages"][base]:.2f}%)')
    lines += [
        f'GC content: {stats["gc_content"]:.2f}%',
        f'AT content: {stats["at_content"]:.2f}%',
        "",
        "DNA SEQUENCE OUTPUTS",
        "Original DNA (5′→3′ as entered):",
        result["sequence"],
        "",
        "Complement (same displayed orientation):",
        result["complement"],
        "",
        "Reverse complement:",
        result["reverse_complement"],
        "",
        f'RNA (convention: {result["strand_convention"]}):',
        result["rna_sequence"],
        "",
        "TRANSLATION",
        f'Reading frame: +{tr["frame"]}',
        f'Complete codons translated: {tr["codon_count"]}',
        f'Trailing bases ignored: {tr["trailing_bases_ignored"]}',
        f'Protein sequence: {tr["protein_sequence"] or "(no complete codons)"}',
        f'Start codon positions (1-based): {", ".join(map(str, tr["start_codon_positions"])) or "None"}',
        f'Stop codon positions (1-based): {", ".join(map(str, tr["stop_codon_positions"])) or "None"}',
        "",
        "CODON FREQUENCIES",
    ]
    if tr["codon_counts"]:
        lines.extend([f"  {codon}: {count}" for codon, count in tr["codon_counts"].items()])
    else:
        lines.append("  No complete codons.")
    lines += [
        "",
        "SCIENTIFIC NOTE",
        "Translation assumes the entered sequence is a DNA coding sequence in the 5′→3′ direction.",
        "This educational output is not confirmation that the sequence encodes a functional protein.",
    ]
    return "\n".join(lines) + "\n"


def build_fasta(sequence: str, identifier: str = "BioSeq_sequence", line_width: int = 70) -> str:
    wrapped = "\n".join(sequence[i : i + line_width] for i in range(0, len(sequence), line_width))
    return f">{identifier} length={len(sequence)}\n{wrapped}\n"


def build_json_report(result: dict[str, Any]) -> str:
    """Create a JSON-serializable report with derived metrics and translation."""
    return json.dumps(result, indent=2, ensure_ascii=False)
