"""Core DNA analysis utilities for BioSeq.

These functions are deliberately separated from the Streamlit interface so they
can be tested and reused independently.
"""
from __future__ import annotations

from collections import Counter
from typing import Any

from Bio.Seq import Seq

# IUPAC nucleotide symbols accepted in DNA sequences.
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
COMPLEMENT_TABLE = str.maketrans(
    {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
    }
)
BASES = ("A", "T", "G", "C")
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"  # Standard genetic code's usual start codon for this educational tool.


def normalize_sequence(raw: str) -> str:
    """Normalize plain DNA or FASTA text, preserving unsupported symbols for validation.

    FASTA headers (lines beginning with '>') are ignored. Whitespace is removed
    and the sequence is uppercased. Invalid characters are not silently deleted.
    """
    if not raw:
        return ""
    lines = []
    for line in raw.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(">"):
            continue
        # Ignore FASTA comment lines when present.
        if stripped.startswith(";"):
            continue
        lines.append(stripped)
    return "".join("".join(lines).split()).upper()


def validate_sequence(sequence: str) -> list[str]:
    """Return validation messages; an empty list means the sequence is valid."""
    errors: list[str] = []
    if not sequence:
        return ["Enter a DNA sequence or upload a FASTA file."]
    invalid = sorted(set(sequence) - IUPAC_DNA)
    if invalid:
        errors.append(
            "Unsupported character(s): " + ", ".join(repr(ch) for ch in invalid)
            + ". Use DNA IUPAC symbols: A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N."
        )
    return errors


def sequence_statistics(sequence: str) -> dict[str, Any]:
    """Calculate length, nucleotide counts, percentages and GC/AT content.

    Percentages and GC/AT content use only unambiguous A/T/G/C bases as the
    denominator. Ambiguous IUPAC symbols are counted separately and excluded
    from those percentages.
    """
    counts_all = Counter(sequence)
    counts = {base: int(counts_all.get(base, 0)) for base in BASES}
    canonical_total = sum(counts.values())
    ambiguous_count = len(sequence) - canonical_total
    percentages = {
        base: (counts[base] / canonical_total * 100.0 if canonical_total else 0.0)
        for base in BASES
    }
    gc_count = counts["G"] + counts["C"]
    at_count = counts["A"] + counts["T"]
    gc_content = gc_count / canonical_total * 100.0 if canonical_total else 0.0
    at_content = at_count / canonical_total * 100.0 if canonical_total else 0.0
    return {
        "sequence_length": len(sequence),
        "canonical_bases": canonical_total,
        "ambiguous_bases": ambiguous_count,
        "counts": counts,
        "percentages": percentages,
        "gc_count": gc_count,
        "at_count": at_count,
        "gc_content": gc_content,
        "at_content": at_content,
    }


def complementary_strand(sequence: str) -> str:
    """Return the complementary DNA strand in the same displayed orientation."""
    return sequence.translate(COMPLEMENT_TABLE)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return complementary_strand(sequence)[::-1]


def transcribe_dna(sequence: str, strand_convention: str = "Coding strand (5′→3′)") -> str:
    """Return RNA in 5′→3′ orientation.

    For a coding strand supplied 5′→3′, T is replaced with U.
    For a template strand supplied 3′→5′, the complementary RNA is produced
    directly in 5′→3′ orientation.
    """
    if strand_convention.startswith("Template"):
        return complementary_strand(sequence).replace("T", "U")
    return sequence.replace("T", "U")


def translate_sequence(sequence: str, frame: int = 1) -> dict[str, Any]:
    """Translate a DNA coding sequence in frame +1, +2, or +3.

    Translation uses the standard genetic code and includes '*' symbols for stop
    codons. Trailing bases that do not complete a codon are reported and ignored.
    Ambiguous codons are translated as X by Biopython when not resolvable.
    """
    if frame not in (1, 2, 3):
        raise ValueError("Reading frame must be 1, 2, or 3.")
    offset = frame - 1
    framed_sequence = sequence[offset:]
    trailing_bases = len(framed_sequence) % 3
    codable_sequence = framed_sequence[: len(framed_sequence) - trailing_bases]
    codons = [codable_sequence[i : i + 3] for i in range(0, len(codable_sequence), 3)]
    protein = str(Seq(codable_sequence).translate(table="Standard", to_stop=False)) if codable_sequence else ""
    start_positions = [offset + (i * 3) + 1 for i, codon in enumerate(codons) if codon == START_CODON]
    stop_positions = [offset + (i * 3) + 1 for i, codon in enumerate(codons) if codon in STOP_CODONS]
    codon_counts = dict(sorted(Counter(codons).items(), key=lambda item: (-item[1], item[0])))
    return {
        "frame": frame,
        "offset": offset,
        "codons": codons,
        "codon_count": len(codons),
        "trailing_bases_ignored": trailing_bases,
        "protein_sequence": protein,
        "start_codon_positions": start_positions,
        "stop_codon_positions": stop_positions,
        "codon_counts": codon_counts,
    }


def analyze_sequence(raw: str, strand_convention: str = "Coding strand (5′→3′)", frame: int = 1) -> dict[str, Any]:
    """Validate and analyze raw plain-text or FASTA DNA input."""
    sequence = normalize_sequence(raw)
    errors = validate_sequence(sequence)
    if errors:
        raise ValueError(" ".join(errors))
    stats = sequence_statistics(sequence)
    translation = translate_sequence(sequence, frame)
    return {
        "sequence": sequence,
        "statistics": stats,
        "complement": complementary_strand(sequence),
        "reverse_complement": reverse_complement(sequence),
        "rna_sequence": transcribe_dna(sequence, strand_convention),
        "strand_convention": strand_convention,
        "translation": translation,
    }
