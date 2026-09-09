import gzip
import os
from pathlib import Path

from loguru import logger

# copied between files in 1 MiB blocks
COPY_BLOCK = 1024 * 1024


def concatenate_short_fastqs(out_dir):
    """moves and copies files
    :param out_dir:  Output Directory
    :param logger: logger
    :return:
    """
    # list all the inputs for concatenation
    unmapped_fastq_one_short: Path = Path(out_dir) / "unmapped_R1.fastq"
    unmapped_fastq_two_short: Path = Path(out_dir) / "unmapped_R2.fastq"
    non_chrom_fastq_one_short: Path = Path(out_dir) / "mapped_non_chromosome_R1.fastq"
    non_chrom_fastq_two_short: Path = Path(out_dir) / "mapped_non_chromosome_R2.fastq"

    # final outputs
    short_one_file: Path = Path(out_dir) / "short_read_concat_R1.fastq"
    short_two_file: Path = Path(out_dir) / "short_read_concat_R2.fastq"

    try:
        concatenate_single_fastq(
            unmapped_fastq_one_short, non_chrom_fastq_one_short, short_one_file
        )
        concatenate_single_fastq(
            unmapped_fastq_two_short, non_chrom_fastq_two_short, short_two_file
        )
    except Exception:
        logger.error("Error with concatenate_fastqs\n")


def _append_file(source: Path, out_handle, record_marker: bytes):
    """Append source to an open binary handle, decompressing it if gzipped.

    Copies in fixed blocks, so memory does not depend on file size, and ensures
    the appended block ends with a newline: without that, a source whose final
    line is unterminated would run into the next file's first header.

    A non-empty file must start with record_marker (b"@" for fastq, b">" for
    fasta). Parsing every record through BioPython used to catch a wrong-format
    input; this keeps that check at O(1) instead of O(file).

    :param source: file to append
    :param out_handle: destination, opened in binary mode
    :param record_marker: first byte every record of this format starts with
    :raises ValueError: if source is non-empty and does not start with the marker
    """
    opener = gzip.open if Path(source).suffix == ".gz" else open
    last_byte = b""
    first = True
    with opener(source, "rb") as in_handle:
        while True:
            block = in_handle.read(COPY_BLOCK)
            if not block:
                break
            if first:
                if block[:1] != record_marker:
                    raise ValueError(
                        f"{source} does not look like a "
                        f"{record_marker.decode()}-delimited file: it starts with "
                        f"{block[:1]!r}"
                    )
                first = False
            out_handle.write(block)
            last_byte = block[-1:]
    # an empty file contributes nothing, which is normal (e.g. no unmapped reads)
    if last_byte and last_byte != b"\n":
        out_handle.write(b"\n")


def _concatenate(sources, out_path, record_marker: bytes):
    """Block-copy sources into out_path, leaving no output behind on failure.

    The copy goes to a sibling .tmp and is renamed into place only once every
    source has been read, so a wrong-format second input, a truncated gzip or a
    full disk cannot leave a half-written FASTQ/FASTA where the run expects a
    complete one. The SeqIO version got this for free by parsing everything
    before it opened the output.

    :param sources: files to copy, in order
    :param out_path: destination path
    :param record_marker: first byte every record of this format starts with
    """
    out_path = Path(out_path)
    tmp_path = out_path.with_name(f"{out_path.name}.tmp")
    try:
        with open(tmp_path, "wb") as out_handle:
            for source in sources:
                _append_file(source, out_handle, record_marker)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    # same directory, so this is an atomic replace rather than a copy
    os.replace(tmp_path, out_path)


def concatenate_single_fastq(fastq_in1: Path, fastq_in2: Path, fastq_out: Path):
    """concatenates 2 fastq files

    Concatenating reads needs no parsing. Round-tripping them through
    SeqIO.parse into a list of SeqRecords cost roughly 5-10x the file size in
    RAM. Only the non-chromosomal short reads pass through here, so that is a
    transient spike of a GB or two on a typical isolate rather than the whole
    run - but it landed immediately before Unicycler, which needs the memory
    itself. A block copy is constant-memory and far faster.

    :param fastq_in1:  fastq_in1 input fastq 1
    :param fastq_in2: fastq_in1 input fastq 2
    :param fastq_out: fastq_out output fastq 2
    :return:
    """
    _concatenate([fastq_in1, fastq_in2], fastq_out, b"@")


def concatenate_single_fasta(file1: Path, file2: Path, output_file: Path):
    """concatenates 2 fasta files
    :param file1: input fasta 1
    :param file2: input fasta 2
    :param output_file: output fasta
    :return:
    """
    _concatenate([file1, file2], output_file, b">")
