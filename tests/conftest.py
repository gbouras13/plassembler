"""Shared pytest fixtures and configuration for the plassembler test suite."""

import sys
from pathlib import Path

import pytest
from loguru import logger

# plassembler routes fatal conditions through `logger.error`, which in the real
# CLI is wired to exit the process. Mirror that here so that error paths are
# assertable via `pytest.raises(SystemExit)`.
logger.add(lambda _: sys.exit(1), level="ERROR")

TEST_DATA = Path("tests/test_data")


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")
