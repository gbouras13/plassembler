"""Deterministic (Tier-A) regression tests for pure functions.

These run without the external bioinformatics toolchain and lock the exact
output of pure helper functions against a committed fixture. They are fast and
fully reproducible, so they act as true regression guards.
"""

from pathlib import Path

import numpy as np

from src.plassembler.utils.depth import (
    collate_depths,
    get_contig_circularity,
    get_contig_lengths,
)

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


def test_collate_depths_copy_number_math():
    """Copy number is mean plasmid depth divided by mean chromosome depth."""
    depths = {
        "chromosome": [10] * 100,
        "plasmid00001": [20] * 60,  # 2x the chromosome depth
    }
    contig_lengths = {"chromosome": 100, "plasmid00001": 60}

    df = collate_depths(depths, "short", contig_lengths).set_index("contig")

    assert df.loc["chromosome", "mean_depth_short"] == 10
    assert df.loc["plasmid00001", "mean_depth_short"] == 20
    # constant depth -> zero stdev
    assert df.loc["chromosome", "sd_depth_short"] == 0
    # copy numbers relative to the chromosome
    assert df.loc["chromosome", "plasmid_copy_number_short"] == 1.0
    assert df.loc["plasmid00001", "plasmid_copy_number_short"] == 2.0


def test_collate_depths_no_chromosome_is_nan_not_crash():
    """No contig named 'chromosome' -> copy number NaN, not a NameError/crash."""
    depths = {"plasmid00001": [20] * 60}
    contig_lengths = {"plasmid00001": 60}
    df = collate_depths(depths, "short", contig_lengths)
    assert df["plasmid_copy_number_short"].isna().all()
    # must stay float-convertible (downstream does .astype(float))
    df["plasmid_copy_number_short"].astype(float)


def test_collate_depths_zero_chromosome_depth_stays_float():
    """Zero chromosome depth (e.g. the fake --no_chromosome chromosome) -> inf,
    which the depth filter treats as 'keep'. Regression: a string 'NA' here broke
    the downstream .astype(float) in combine_depth_mash_tsvs (CI --no_chromosome)."""
    depths = {"chromosome": [0] * 100, "plasmid00001": [5] * 60}
    contig_lengths = {"chromosome": 100, "plasmid00001": 60}
    df = collate_depths(depths, "short", contig_lengths).set_index("contig")
    assert np.isinf(df.loc["plasmid00001", "plasmid_copy_number_short"])
    # the actual CI failure: this must not raise
    df["plasmid_copy_number_short"].astype(float)
