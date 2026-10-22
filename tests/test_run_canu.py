"""Tests for how plassembler.utils.run_canu reports canu stopping."""

import os
import sys
from pathlib import Path

import pytest
from loguru import logger

from src.plassembler.utils.run_canu import run_canu, run_canu_correct

# canu 2.3's own stderr, captured from real runs: a job that could not run
# (meryl-count under an expr that rejects `--`, as in #92), too few reads, and
# no reads at all. canu 2.2 words all three the same way.
CANU_DATA = Path("tests/test_data/canu")

# run_canu_correct and run_canu share the failure handling, so the cases that
# apply to both run against both
CANU_RUNS = [
    pytest.param(run_canu_correct, "correct reads", id="correct"),
    pytest.param(run_canu, "assemble", id="assemble"),
]


@pytest.fixture
def captured_logs():
    """Collect WARNING and ERROR records, without the exiting sink conftest installs.

    That sink raises SystemExit on the first ERROR, which would stop this one from
    ever seeing the message. The fatal canu path exits on its own, so dropping the
    sinks for the duration still exercises the real control flow.
    """
    logger.remove()
    records = []
    logger.add(lambda message: records.append(message.record), level="WARNING")
    try:
        yield records
    finally:
        logger.remove()
        logger.add(sys.stderr)
        logger.add(lambda _: sys.exit(1), level="ERROR")


def messages(records, level):
    return [record["message"] for record in records if record["level"].name == level]


def fake_canu(tmp_path, monkeypatch, stderr_file=None, returncode=1):
    """Put a stub `canu` first on PATH that replays canu's stderr and exits.

    Lets each way canu stops be driven without canu installed - or a broken expr
    to make it fail.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "canu"
    replay = f"cat '{Path(stderr_file).resolve()}' >&2\n" if stderr_file else ""
    stub.write_text(f"#!/bin/sh\n{replay}exit {returncode}\n")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")


def call(run, tmp_path):
    run(
        1,
        tmp_path / "logs",
        tmp_path / "reads.fastq",
        tmp_path / "canu",
        "nanopore",
        145059,
        0.12,
        50,
    )


@pytest.mark.parametrize("run, step", CANU_RUNS)
def test_canu_job_failure_is_fatal(tmp_path, monkeypatch, captured_logs, run, step):
    """A canu job that cannot run stops plassembler, naming the stage (#92).

    Advancing with uncorrected reads instead left 110 spurious contigs for the
    depth filter in #92, while telling the user not to be concerned.
    """
    fake_canu(tmp_path, monkeypatch, CANU_DATA / "job_failure.err")
    with pytest.raises(SystemExit):
        call(run, tmp_path)
    [error] = messages(captured_logs, "ERROR")
    assert f"could not {step}" in error
    assert "Kmer counting (meryl-count) jobs failed, tried 2 times, giving up." in error
    # the failed job's own log, which holds the real cause ("expr: syntax error"),
    # is named relative to the kept canu directory
    assert str(tmp_path / "canu") in error
    assert "./meryl-count.sh 1 > ./meryl-count.000001.out 2>&1" in error
    # ExternalTool's header lands at the end of the log, and is not canu's output
    assert "Command line:" not in error
    assert messages(captured_logs, "WARNING") == []


def test_canu_job_failure_gets_past_the_fallback_to_uncorrected_reads(
    tmp_path, monkeypatch
):
    """plassembler long catches any Exception from read correction and carries on
    with uncorrected reads. The fatal path has to get through that, here under
    the same exiting ERROR sink the CLI installs.
    """
    fake_canu(tmp_path, monkeypatch, CANU_DATA / "job_failure.err")
    with pytest.raises(SystemExit):
        try:
            call(run_canu_correct, tmp_path)
        except Exception:
            pytest.fail("a broken canu was handled as an ordinary exception")


def test_canu_crash_is_fatal(tmp_path, monkeypatch, captured_logs):
    """canu's banner for an error it calls abnormal also stops the run."""
    crash = tmp_path / "crash.err"
    crash.write_text(
        "CRASH:\n"
        "CRASH: canu 2.3\n"
        "CRASH: Please panic, this is abnormal.\n"
        "CRASH:\n"
        "CRASH:   failed to open 'correction/0-mercounts/meryl-count.sh' for writing.\n"
        "CRASH:\n"
        "CRASH: No log file supplied.\n"
        "CRASH:\n"
    )
    fake_canu(tmp_path, monkeypatch, crash)
    with pytest.raises(SystemExit):
        call(run_canu_correct, tmp_path)
    [error] = messages(captured_logs, "ERROR")
    assert "canu reported an internal error" in error
    assert "failed to open 'correction/0-mercounts/meryl-count.sh'" in error


@pytest.mark.parametrize("fixture", ["low_coverage.err", "no_reads.err"])
@pytest.mark.parametrize("run, step", CANU_RUNS)
def test_canu_low_coverage_is_not_fatal(
    tmp_path, monkeypatch, captured_logs, run, step, fixture
):
    """Too few reads is what a sample with no plasmids gives canu, so it carries on.

    canu ends these with the same ABORT banner as a real failure, which is why the
    banner cannot be what decides.
    """
    assert "ABORT: Don't panic" in (CANU_DATA / fixture).read_text()
    fake_canu(tmp_path, monkeypatch, CANU_DATA / fixture)
    call(run, tmp_path)
    assert messages(captured_logs, "ERROR") == []
    warnings = messages(captured_logs, "WARNING")
    assert f"canu stopped before it could {step}, as read coverage was too low." in (
        warnings
    )


@pytest.mark.parametrize("run, step", CANU_RUNS)
def test_canu_unrecognised_stop_shows_canus_reason(
    tmp_path, monkeypatch, captured_logs, run, step
):
    """A stop plassembler does not classify carries on, but with canu's own reason.

    canu has close to 200 ways of stopping; one this does not know must not turn
    a plasmid-free sample into a failed run - but nor should the user be told it
    was nothing to worry about.
    """
    unknown = tmp_path / "unknown.err"
    unknown.write_text(
        "ABORT:\n"
        "ABORT: canu 2.3\n"
        "ABORT: Don't panic, but a mostly harmless error occurred and Canu stopped.\n"
        "ABORT: Try restarting.  If that doesn't work, ask for help.\n"
        "ABORT:\n"
        "ABORT:   not enough memory for correction; increase corMemory to at least 12.\n"
        "ABORT:\n"
    )
    fake_canu(tmp_path, monkeypatch, unknown, returncode=1)
    call(run, tmp_path)
    assert messages(captured_logs, "ERROR") == []
    warnings = "\n".join(messages(captured_logs, "WARNING"))
    assert (
        f"canu exited with return code 1 before it could {step}, for a reason "
        "plassembler does not recognise." in warnings
    )
    assert "not enough memory for correction" in warnings


@pytest.mark.parametrize("run, step", CANU_RUNS)
def test_canu_success_reports_nothing(tmp_path, monkeypatch, captured_logs, run, step):
    fake_canu(tmp_path, monkeypatch, returncode=0)
    call(run, tmp_path)
    assert captured_logs == []
