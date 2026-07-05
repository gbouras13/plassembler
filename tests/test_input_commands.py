"""Error-path and good-path unit tests for input validation.

All of these are pure Python (no external tools). Fatal branches route through
``logger.error``, which the shared conftest wires to ``SystemExit``, so error
paths are asserted via ``pytest.raises(SystemExit)``.
"""

from pathlib import Path

import pytest

from src.plassembler.utils.input_commands import (
    validate_fastqs_assembled_mode,
    validate_flye_assembly_info,
    validate_flye_directory,
    validate_pacbio_model,
)

TEST_DATA = Path("tests/test_data")
VAL = TEST_DATA / "validation"
FLYE_DIR = TEST_DATA / "end_to_end" / "test_flye_dir"
PLASS_CLASS = TEST_DATA / "plass_class"


# ---------------------------------------------------------------------------
# validate_pacbio_model
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "model,expected",
    [
        ("pacbio-raw", "--pacbio-raw"),
        ("pacbio-corr", "--pacbio-corr"),
        ("pacbio-hifi", "--pacbio-hifi"),
    ],
)
def test_validate_pacbio_model_good(model, expected):
    assert validate_pacbio_model(model) == expected


def test_validate_pacbio_model_bad_exits():
    with pytest.raises(SystemExit):
        validate_pacbio_model("not_a_model")


# ---------------------------------------------------------------------------
# validate_flye_assembly_info
# ---------------------------------------------------------------------------


def test_flye_assembly_info_both_present_skips():
    assembly = PLASS_CLASS / "assembly.fasta"
    info = PLASS_CLASS / "assembly_info.txt"
    assert validate_flye_assembly_info(str(assembly), str(info)) is True


def test_flye_assembly_info_assembly_without_info_no_skip():
    assembly = PLASS_CLASS / "assembly.fasta"
    assert validate_flye_assembly_info(str(assembly), "nothing") is False


def test_flye_assembly_info_info_without_assembly_no_skip():
    info = PLASS_CLASS / "assembly_info.txt"
    assert validate_flye_assembly_info("nothing", str(info)) is False


def test_flye_assembly_info_neither_no_skip():
    assert validate_flye_assembly_info("nothing", "nothing") is False


# ---------------------------------------------------------------------------
# validate_flye_directory
# ---------------------------------------------------------------------------


def test_flye_directory_valid_skips():
    assert validate_flye_directory(str(FLYE_DIR)) is True


def test_flye_directory_missing_no_skip():
    assert validate_flye_directory(str(TEST_DATA / "does_not_exist_dir")) is False


# ---------------------------------------------------------------------------
# validate_fastqs_assembled_mode
# ---------------------------------------------------------------------------


def test_fastqs_assembled_mode_good_tuple():
    longreads = VAL / "test.fastq"  # not gzipped
    s1 = VAL / "test_2.fastq.gz"
    s2 = VAL / "test_2.fastq.gz"
    result = validate_fastqs_assembled_mode(str(longreads), str(s1), str(s2))
    # (short_flag, long_flag, long_gzipped)
    assert result == (True, True, False)


def test_fastqs_assembled_mode_inconsistent_compression_exits():
    s1 = VAL / "test_2.fastq.gz"  # gzipped
    s2 = VAL / "test.fastq"  # not gzipped
    with pytest.raises(SystemExit):
        validate_fastqs_assembled_mode("nothing", str(s1), str(s2))


def test_fastqs_assembled_mode_single_short_file_exits():
    longreads = VAL / "test.fastq"
    s1 = VAL / "test.fastq"
    with pytest.raises(SystemExit):
        validate_fastqs_assembled_mode(str(longreads), str(s1), "nothing")


def test_fastqs_assembled_mode_no_inputs_exits():
    with pytest.raises(SystemExit):
        validate_fastqs_assembled_mode("nothing", "nothing", "nothing")
