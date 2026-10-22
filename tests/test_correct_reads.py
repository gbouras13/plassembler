"""Tests for long read correction in plassembler.utils.correct_reads."""

import os
import shutil
import sys
from pathlib import Path

import pytest
from loguru import logger

from src.plassembler.utils.correct_reads import (
    correction_target_bases,
    run_lrcorrect,
)

TEST_DATA = Path("tests/test_data")
READS = TEST_DATA / "plasmid_long_good.fastq"


@pytest.fixture
def captured_logs():
    """Collect WARNING and ERROR messages, without the exiting sink conftest installs.

    That sink raises SystemExit, which would stop any later sink - including this
    one - from ever seeing the message. The lrcorrect error path exits on its own,
    so dropping the sinks for the duration still exercises the real control flow.
    """
    logger.remove()
    messages = []
    logger.add(lambda m: messages.append((m.record["level"].name, str(m))))
    try:
        yield messages
    finally:
        logger.remove()
        logger.add(sys.stderr)
        logger.add(lambda _: sys.exit(1), level="ERROR")


def fake_lrcorrect(tmp_path, monkeypatch, script):
    """Put a stub `lrcorrect` first on PATH.

    The stub gets the same arguments as the real tool, so `$2`... can be used to
    find the output path; it is always the argument after `-o`.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "lrcorrect"
    stub.write_text(script)
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")


# the output path follows -o; the stubs below find it without depending on order
FIND_OUTPUT = 'while [ "$1" != "-o" ]; do shift; done; out="$2"\n'


def messages_at(captured, level):
    return "\n".join(text for lvl, text in captured if lvl == level)


def test_target_bases_keeps_canus_coverage():
    """50x of Flye's plasmid contigs, or 250x of 10 kbp when there are none:
    the maxInputCoverage/genomeSize pair canu was given."""
    assert correction_target_bases(145059) == 145059 * 50
    assert correction_target_bases(0) == 10000 * 250


def test_corrected_reads_are_used(tmp_path, monkeypatch, captured_logs):
    fake_lrcorrect(
        tmp_path,
        monkeypatch,
        "#!/bin/sh\n" + FIND_OUTPUT + 'printf "@r\\nACGT\\n+\\n####\\n" > "$out"\n',
    )
    out = tmp_path / "corrected.fastq"
    assert run_lrcorrect(2, tmp_path / "logs", READS, out, "ont", 1000) is True
    assert out.read_text().startswith("@r")
    assert messages_at(captured_logs, "WARNING") == ""


def test_lrcorrect_gets_plassemblers_arguments(tmp_path, monkeypatch, captured_logs):
    args_file = tmp_path / "args.txt"
    fake_lrcorrect(
        tmp_path,
        monkeypatch,
        f'#!/bin/sh\necho "$@" > {args_file}\n'
        + FIND_OUTPUT
        + 'printf "@r\\nACGT\\n+\\n####\\n" > "$out"\n',
    )
    out = tmp_path / "corrected.fastq"
    run_lrcorrect(3, tmp_path / "logs", READS, out, "pb", 7252950)
    args = args_file.read_text().split()
    for flag, value in (
        ("--threads", "3"),
        ("--preset", "pb"),
        ("--target-bases", "7252950"),
        ("--qual-char", "#"),
        ("-o", str(out)),
    ):
        assert args[args.index(flag) + 1] == value
    assert args[-1] == str(READS)


def test_too_little_data_carries_on(tmp_path, monkeypatch, captured_logs):
    """No reads out, with exit status 0, is lrcorrect's "too little data": a
    sample with no plasmids. The run must carry on with uncorrected reads."""
    fake_lrcorrect(
        tmp_path,
        monkeypatch,
        "#!/bin/sh\n"
        + FIND_OUTPUT
        + ': > "$out"\n'
        + 'echo "lrcorrect: no reads could be corrected: too few reads overlap" >&2\n',
    )
    out = tmp_path / "corrected.fastq"
    assert run_lrcorrect(2, tmp_path / "logs", READS, out, "ont", 1000) is False
    warnings = messages_at(captured_logs, "WARNING")
    assert "too few" in warnings
    assert "no reads could be corrected" in warnings
    assert messages_at(captured_logs, "ERROR") == ""


def test_lrcorrect_failure_is_fatal(tmp_path, monkeypatch, captured_logs):
    """A corrector that fails must stop the run and show its own message.

    The regression this replaces (#92): canu failing to run was reported as
    "don't be too concerned" and the run went on with uncorrected reads, which
    gave Unicycler a graph full of spurious contigs.
    """
    fake_lrcorrect(
        tmp_path,
        monkeypatch,
        "#!/bin/sh\necho 'lrcorrect: error: cannot parse reads.fq: bad header' >&2\nexit 1\n",
    )
    logdir = tmp_path / "logs"
    with pytest.raises(SystemExit):
        run_lrcorrect(2, logdir, READS, tmp_path / "corrected.fastq", "ont", 1000)
    errors = messages_at(captured_logs, "ERROR")
    assert "lrcorrect (return code 1)" in errors
    assert "cannot parse reads.fq: bad header" in errors
    assert str(logdir) in errors


def test_lrcorrect_killed_is_fatal(tmp_path, monkeypatch, captured_logs):
    """Killed by a signal (e.g. out of memory) is a failure, not "no plasmids"."""
    fake_lrcorrect(tmp_path, monkeypatch, "#!/bin/sh\nkill -9 $$\n")
    with pytest.raises(SystemExit):
        run_lrcorrect(
            2, tmp_path / "logs", READS, tmp_path / "corrected.fastq", "ont", 1000
        )
    assert "return code -9" in messages_at(captured_logs, "ERROR")


def test_lrcorrect_missing_is_fatal(tmp_path, monkeypatch, captured_logs):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    with pytest.raises(SystemExit):
        run_lrcorrect(
            2, tmp_path / "logs", READS, tmp_path / "corrected.fastq", "ont", 1000
        )
    assert "Error with lrcorrect" in messages_at(captured_logs, "ERROR")


@pytest.mark.requires_tool
@pytest.mark.skipif(shutil.which("lrcorrect") is None, reason="needs lrcorrect")
def test_real_lrcorrect(tmp_path, captured_logs):
    out = tmp_path / "corrected.fastq"
    assert run_lrcorrect(2, tmp_path / "logs", READS, out, "ont", 2000000) is True
    with open(out) as fh:
        header, seq, plus, qual = (fh.readline().strip() for _ in range(4))
    assert header.startswith("@") and plus == "+"
    assert len(seq) == len(qual) and set(qual) == {"#"}


@pytest.mark.requires_tool
@pytest.mark.skipif(shutil.which("lrcorrect") is None, reason="needs lrcorrect")
def test_real_lrcorrect_empty_input(tmp_path, captured_logs):
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    out = tmp_path / "corrected.fastq"
    assert run_lrcorrect(1, tmp_path / "logs", empty, out, "ont", 2500000) is False
