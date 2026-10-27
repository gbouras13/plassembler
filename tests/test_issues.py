"""Regression tests for reported GitHub issues."""

import pytest

from src.plassembler.utils.correct_reads import correction_read_type
from src.plassembler.utils.run_unicycler import build_extra_unicycler_options

# ---------------------------------------------------------------------------
# issue #83: Unicycler crashes with --spades_options but no --unicycler_options
# (a None unicycler_options was interpolated into the command as "None")
# ---------------------------------------------------------------------------


def test_unicycler_options_both_none():
    assert build_extra_unicycler_options(None, None) == ""


def test_unicycler_options_spades_only_has_no_none():
    out = build_extra_unicycler_options(None, "--memory 2000")
    assert out == '--spades_options "--memory 2000"'
    assert "None" not in out


def test_unicycler_options_unicycler_only():
    assert build_extra_unicycler_options("--mode bold", None) == "--mode bold"


def test_unicycler_options_both():
    out = build_extra_unicycler_options("--mode bold", "--memory 2000")
    assert out == '--mode bold --spades_options "--memory 2000"'


# ---------------------------------------------------------------------------
# issue #84: `plassembler long --pacbio_model pacbio-hifi` tried to correct
# already-corrected HiFi reads (with canu, which errored on them)
# ---------------------------------------------------------------------------


def test_read_type_hifi_skips_correction():
    preset, skip = correction_read_type("--pacbio-hifi")
    assert preset == "pb"
    assert skip is True


@pytest.mark.parametrize("model", ["--pacbio-raw", "--pacbio-corr"])
def test_read_type_pacbio(model):
    preset, skip = correction_read_type(model)
    assert preset == "pb"
    assert skip is False


def test_read_type_nanopore():
    preset, skip = correction_read_type("nothing")
    assert preset == "ont"
    assert skip is False
