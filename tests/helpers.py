"""Helpers for golden-file / deterministic regression tests.

Golden files live under ``tests/test_data/golden``. To (re)generate them after an
intentional output change, run the suite with ``PLASSEMBLER_WRITE_GOLDEN=1`` set,
which makes the ``assert_*`` helpers overwrite the golden file instead of
comparing against it.
"""

import os
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO

GOLDEN_DIR = Path("tests/test_data/golden")

# Set PLASSEMBLER_WRITE_GOLDEN=1 to regenerate golden files instead of asserting.
WRITE_GOLDEN = os.environ.get("PLASSEMBLER_WRITE_GOLDEN") == "1"


def read_fasta(path):
    """Return ``{record.id: (uppercase_sequence, description)}`` for a FASTA."""
    return {
        rec.id: (str(rec.seq).upper(), rec.description)
        for rec in SeqIO.parse(str(path), "fasta")
    }


def assert_fasta_matches(actual_path, golden_path):
    """Assert a produced FASTA matches a golden FASTA by id and sequence.

    Headers/line-wrapping are normalised away; only contig ids and (uppercased)
    sequences are compared.
    """
    if WRITE_GOLDEN:
        Path(golden_path).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(actual_path, golden_path)
        return

    actual = read_fasta(actual_path)
    golden = read_fasta(golden_path)
    assert set(actual) == set(golden), (
        f"contig ids differ: {sorted(actual)} != {sorted(golden)}"
    )
    for contig_id in golden:
        assert actual[contig_id][0] == golden[contig_id][0], (
            f"sequence for {contig_id!r} differs from golden"
        )


def assert_tsv_matches(actual_path, golden_path, float_tol=0.01, sort_col=None):
    """Assert a produced TSV matches a golden TSV, tolerant of float rounding.

    :param float_tol: absolute tolerance for numeric columns.
    :param sort_col: column to sort both frames by before comparing (use when
        row order is not guaranteed).
    """
    actual = pd.read_csv(actual_path, sep="\t")

    if WRITE_GOLDEN:
        Path(golden_path).parent.mkdir(parents=True, exist_ok=True)
        actual.to_csv(golden_path, sep="\t", index=False)
        return

    golden = pd.read_csv(golden_path, sep="\t")

    assert list(actual.columns) == list(golden.columns), (
        f"columns differ: {list(actual.columns)} != {list(golden.columns)}"
    )

    if sort_col:
        actual = actual.sort_values(sort_col).reset_index(drop=True)
        golden = golden.sort_values(sort_col).reset_index(drop=True)

    assert len(actual) == len(golden), (
        f"row count differs: {len(actual)} != {len(golden)}"
    )

    for col in golden.columns:
        a, g = actual[col], golden[col]
        if pd.api.types.is_numeric_dtype(g) and pd.api.types.is_numeric_dtype(a):
            close = (a - g).abs() <= float_tol
            both_nan = a.isna() & g.isna()
            assert bool((close | both_nan).all()), (
                f"numeric column {col!r} differs beyond tolerance {float_tol}"
            )
        else:
            assert (a.astype(str) == g.astype(str)).all(), (
                f"column {col!r} differs from golden"
            )
