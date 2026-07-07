"""Regression tests for reported GitHub issues."""

import pytest

from src.plassembler.utils.run_canu import canu_read_type_and_error_rate
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
# issue #84: `plassembler long --pacbio_model pacbio-hifi` passed the wrong
# read type to canu (-pacbio instead of -pacbio-hifi) and tried to correct
# already-corrected HiFi reads
# ---------------------------------------------------------------------------


def test_canu_read_type_hifi_skips_correction():
    read_type, err, skip = canu_read_type_and_error_rate("--pacbio-hifi", 0.12)
    assert read_type == "pacbio-hifi"
    assert err == 0.005
    assert skip is True


@pytest.mark.parametrize("model", ["--pacbio-raw", "--pacbio-corr"])
def test_canu_read_type_pacbio(model):
    read_type, err, skip = canu_read_type_and_error_rate(model, 0.12)
    assert read_type == "pacbio"
    assert err == 0.045
    assert skip is False


def test_canu_read_type_nanopore_passthrough():
    read_type, err, skip = canu_read_type_and_error_rate("nothing", 0.12)
    assert read_type == "nanopore"
    assert err == 0.12  # keeps the CLI --corrected_error_rate
    assert skip is False
