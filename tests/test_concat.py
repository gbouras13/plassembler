"""Tests for the streaming concatenation helpers in plassembler.utils.concat."""

import gzip

import pytest
from Bio import SeqIO

from src.plassembler.utils.concat import (
    concatenate_single_fasta,
    concatenate_single_fastq,
)

READS_A = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n"
READS_B = "@r3\nGGGG\n+\nKKKK\n"


def records(path):
    return [
        (r.id, str(r.seq), tuple(r.letter_annotations["phred_quality"]))
        for r in SeqIO.parse(path, "fastq")
    ]


def test_concatenate_fastq_keeps_every_record(tmp_path):
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a.write_text(READS_A)
    b.write_text(READS_B)
    out = tmp_path / "out.fastq"
    concatenate_single_fastq(a, b, out)
    assert [r[0] for r in records(out)] == ["r1", "r2", "r3"]
    assert out.read_text() == READS_A + READS_B


def test_concatenate_fastq_handles_gzipped_input(tmp_path):
    """.gz inputs are decompressed, as the SeqIO version did."""
    a, b = tmp_path / "a.fastq.gz", tmp_path / "b.fastq"
    with gzip.open(a, "wt") as fh:
        fh.write(READS_A)
    b.write_text(READS_B)
    out = tmp_path / "out.fastq"
    concatenate_single_fastq(a, b, out)
    assert [r[0] for r in records(out)] == ["r1", "r2", "r3"]


def test_concatenate_fastq_with_empty_first_file(tmp_path):
    """An empty input is normal (e.g. no unmapped reads) and contributes nothing."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a.write_text("")
    b.write_text(READS_B)
    out = tmp_path / "out.fastq"
    concatenate_single_fastq(a, b, out)
    assert out.read_text() == READS_B


def test_concatenate_inserts_missing_newline_between_files(tmp_path):
    """A source with no trailing newline must not run into the next file's header."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a.write_text(READS_A.rstrip("\n"))
    b.write_text(READS_B)
    out = tmp_path / "out.fastq"
    concatenate_single_fastq(a, b, out)
    assert [r[0] for r in records(out)] == ["r1", "r2", "r3"]


def test_concatenate_fastq_rejects_a_fasta(tmp_path):
    """Wrong-format input is still caught, as BioPython's parser used to."""
    a, b = tmp_path / "a.fasta", tmp_path / "b.fastq"
    a.write_text(">contig1\nACGT\n")
    b.write_text(READS_B)
    with pytest.raises(ValueError):
        concatenate_single_fastq(a, b, tmp_path / "out.fastq")


def test_concatenate_fasta_keeps_every_contig(tmp_path):
    a, b = tmp_path / "a.fasta", tmp_path / "b.fasta"
    a.write_text(">chromosome\nACGTACGT\n")
    b.write_text(">1 circular=True\nGGGG\n>2\nTTTT\n")
    out = tmp_path / "out.fasta"
    concatenate_single_fasta(a, b, out)
    parsed = list(SeqIO.parse(out, "fasta"))
    assert [r.id for r in parsed] == ["chromosome", "1", "2"]
    # descriptions must survive: get_contig_circularity looks for "circular"
    assert "circular" in parsed[1].description


def test_concatenate_fasta_rejects_a_fastq(tmp_path):
    a, b = tmp_path / "a.fastq", tmp_path / "b.fasta"
    a.write_text(READS_A)
    b.write_text(">1\nACGT\n")
    with pytest.raises(ValueError):
        concatenate_single_fasta(a, b, tmp_path / "out.fasta")


def test_concatenate_leaves_no_output_when_second_file_is_wrong_format(tmp_path):
    """A failure must not leave a half-written file where the run expects a whole one.

    The block copy opens the output before it has seen the second input, so
    without the .tmp-and-rename the first file's reads would be left behind as a
    complete-looking FASTQ.
    """
    a, b = tmp_path / "a.fastq", tmp_path / "b.fasta"
    a.write_text(READS_A)
    b.write_text(">contig1\nACGT\n")
    out = tmp_path / "out.fastq"
    with pytest.raises(ValueError):
        concatenate_single_fastq(a, b, out)
    assert not out.exists()
    assert list(tmp_path.glob("*.tmp")) == []


def test_concatenate_leaves_an_existing_output_untouched_on_failure(tmp_path):
    """The previous output survives a failed re-run rather than being truncated."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fasta"
    a.write_text(READS_A)
    b.write_text(">contig1\nACGT\n")
    out = tmp_path / "out.fastq"
    out.write_text(READS_B)
    with pytest.raises(ValueError):
        concatenate_single_fastq(a, b, out)
    assert out.read_text() == READS_B


def test_concatenate_leaves_no_output_when_gzip_is_truncated(tmp_path):
    """A truncated .gz raises from the decompressor mid-copy, so nothing is kept."""
    a, b = tmp_path / "a.fastq.gz", tmp_path / "b.fastq"
    with gzip.open(a, "wb") as fh:
        fh.write((READS_A * 5000).encode())
    raw = a.read_bytes()
    a.write_bytes(raw[: len(raw) // 2])
    b.write_text(READS_B)
    out = tmp_path / "out.fastq"
    with pytest.raises(EOFError):
        concatenate_single_fastq(a, b, out)
    assert not out.exists()


def test_concatenate_accepts_a_string_output_path(tmp_path):
    """tests/test_plassembler.py passes os.path.join(...) rather than a Path."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    a.write_text(READS_A)
    b.write_text(READS_B)
    out = str(tmp_path / "out.fastq")
    concatenate_single_fastq(a, b, out)
    assert [r[0] for r in records(out)] == ["r1", "r2", "r3"]
