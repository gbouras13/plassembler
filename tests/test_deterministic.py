"""Deterministic (Tier-A) regression tests for pure functions.

These run without the external bioinformatics toolchain and lock the exact
output of pure helper functions against a committed fixture. They are fast and
fully reproducible, so they act as true regression guards.
"""

from pathlib import Path

from src.plassembler.utils.depth import get_contig_circularity, get_contig_lengths

GOLDEN_FASTA = Path("tests/test_data/golden/contigs.fasta")


def test_get_contig_lengths_exact():
    """Contig lengths are read exactly from the FASTA."""
    assert get_contig_lengths(GOLDEN_FASTA) == {
        "chromosome": 100,
        "plasmid00001": 60,
        "plasmid00002": 40,
    }


def test_get_contig_circularity_exact():
    """Circularity is 'circular' for the chromosome (by id) and for any contig
    whose description contains 'circular'; everything else is 'not_circular'."""
    assert get_contig_circularity(GOLDEN_FASTA) == {
        "chromosome": "circular",  # id contains "chromosome"
        "plasmid00001": "circular",  # description contains "circular"
        "plasmid00002": "not_circular",
    }
