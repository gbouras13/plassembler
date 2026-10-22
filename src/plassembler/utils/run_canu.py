import math
import re
import subprocess
import sys
from itertools import product
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from plassembler.utils.external_tools import ExternalTool

# canu ends every failed run with the same "ABORT: Don't panic" banner - a
# sample with too few reads gets it too - so the banner says nothing about what
# went wrong. These are the lines that do, and they read the same in 2.2 and 2.3.

# A compute job that never produced its output, which canu words the same way
# at every stage ("Kmer counting (meryl-count) jobs failed, tried 2 times, giving
# up."). Jobs only start once the reads have passed canu's coverage checks, so
# this is canu's environment failing - a job killed for memory, a missing binary,
# or a shell tool its job scripts cannot use, such as BusyBox's expr (#92).
CANU_JOB_FAILURE = re.compile(
    r"^-- (.*failed, tried \d+ times, giving up\.)$", re.MULTILINE
)
# canu's banner for an error it calls abnormal ("Please panic, this is abnormal")
CANU_CRASH = re.compile(r"^CRASH:", re.MULTILINE)
# canu's stop for thin reads, either on loading them (minInputCoverage) or
# between stages (stopOnLowCoverage). Once the chromosome is removed there are
# few reads left in a sample with no plasmids, so this is an expected outcome.
CANU_LOW_COVERAGE = re.compile(r"Read coverage \(.*\) lower than allowed")


def canu_read_type_and_error_rate(pacbio_model, default_error_rate):
    """Map a validated ``--pacbio_model`` value to canu's read-type flag, the
    corrected error rate, and whether canu read correction should be skipped.

    ``pacbio_model`` is the output of ``validate_pacbio_model`` (so one of
    ``"--pacbio-raw"``, ``"--pacbio-corr"``, ``"--pacbio-hifi"``) or
    ``"nothing"``.

    PacBio HiFi reads are already high-accuracy: canu's read type is
    ``-pacbio-hifi`` and its correction step is skipped (running ``canu -correct``
    on HiFi reads errors with "Cannot correct already corrected reads" - issue
    #84). Previously the check compared against ``"pacbio-hifi"`` (without the
    ``--`` prefix), so HiFi reads were wrongly assembled as regular ``-pacbio``.

    :return: ``(canu_read_type, corrected_error_rate, skip_correction)``
    """
    if pacbio_model == "--pacbio-hifi":
        return "pacbio-hifi", 0.005, True
    if pacbio_model != "nothing":
        return "pacbio", 0.045, False
    return "nanopore", default_error_rate, False


def _canu_stderr_tail(log_text, max_lines=25):
    """Last few lines of canu's stderr, ready to embed in a message.

    canu writes everything it reports to stderr, and the lines before its ABORT
    banner say which stage failed and name each job's own log. ExternalTool's
    "Command line:" header is dropped: it is flushed when the log closes, so it
    lands after canu's output.

    :param log_text: contents of the canu .err logfile
    :param max_lines: how many trailing lines to keep
    :return: the trailing lines, or "" if there are none
    """
    lines = [
        line
        for line in log_text.strip().splitlines()
        if not line.startswith("Command line:")
    ]
    return "\n".join(lines[-max_lines:])


def _report_canu_failure(canu, canu_output_dir, step, returncode):
    """Reports why canu exited non-zero, and stops the run if canu is broken.

    Only the failures canu itself marks as a job failing or as a crash stop the
    run. Everything else - including the coverage stop a sample with no plasmids
    is expected to hit - carries on as before: canu has close to 200 ways of
    stopping, and one this misses must not turn a plasmid-free sample into a
    failed run.

    :param canu: the ExternalTool that ran canu
    :param canu_output_dir: canu's -d directory, which holds each job's own log
    :param step: what canu was asked to do, e.g. "correct reads"
    :param returncode: canu's exit status
    """
    try:
        log_text = Path(canu.err_log).read_text(errors="replace")
    except OSError:
        log_text = ""
    stderr_tail = _canu_stderr_tail(log_text)
    job_failure = CANU_JOB_FAILURE.search(log_text)

    if job_failure or CANU_CRASH.search(log_text):
        reason = (
            job_failure.group(1) if job_failure else "canu reported an internal error"
        )
        # everything goes into a single logger.error: under the CLI that sink
        # exits the process, so a second call would never be reached
        message = (
            f"Error with canu: it could not {step}. {reason}\n"
            "This is canu failing to run, not a sign of too few reads - for "
            "example a job killed for running out of memory, or a tool its job "
            f"scripts cannot use. canu's working directory has been kept at "
            f"{canu_output_dir}, and the log of each job named below is inside "
            f"it. Please check {canu.err_log}"
        )
        if stderr_tail:
            message = f"{message}\ncanu stderr:\n{stderr_tail}"
        logger.error(message)
        # the CLI exits inside logger.error above; this covers run_canu being
        # driven as a library, where no ERROR sink is installed
        sys.exit(1)

    if CANU_LOW_COVERAGE.search(log_text):
        logger.warning(
            f"canu stopped before it could {step}, as read coverage was too low."
        )
        logger.warning(
            "This usually means few reads are left once the chromosome is "
            "removed, most often because the sample has no plasmids."
        )
    else:
        logger.warning(
            f"canu exited with return code {returncode} before it could {step}, "
            "for a reason plassembler does not recognise."
        )
        if stderr_tail:
            logger.warning(f"canu stderr:\n{stderr_tail}")
    logger.warning(
        f"Please check stdout log file: {canu.out_log} and stderr log file: "
        f"{canu.err_log}"
    )


def run_canu_correct(
    threads: Path,
    logdir: Path,
    longreads: Path,
    canu_output_dir: Path,
    canu_nano_or_pacbio: str,
    total_flye_plasmid_length: int,
    corrected_error_rate: float,
    coverage: int,
) -> None:
    """runs canu correct
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """
    # for the assembly param need to divide by a million
    total_flye_plasmid_length = round(total_flye_plasmid_length / 1000000, 5)
    canu = ExternalTool(
        tool="canu -correct",
        input="",
        output="",
        params=f" -p canu -d {canu_output_dir} genomeSize={total_flye_plasmid_length}m maxInputCoverage={coverage} stopOnLowCoverage=1 maxThreads={threads} useGrid=false -{canu_nano_or_pacbio} correctedErrorRate={corrected_error_rate} {longreads}",
        logdir=logdir,
        outfile="",
    )
    try:
        canu.run()
    except subprocess.CalledProcessError as error:
        _report_canu_failure(canu, canu_output_dir, "correct reads", error.returncode)


def run_canu(
    threads: Path,
    logdir: Path,
    longreads: Path,
    canu_output_dir: Path,
    canu_nano_or_pacbio: str,
    total_flye_plasmid_length: int,
    corrected_error_rate: float,
    coverage: int,
) -> None:
    """runs canu
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """
    # for the assembly param need to divide by a million
    total_flye_plasmid_length = round(total_flye_plasmid_length / 1000000, 5)
    canu = ExternalTool(
        tool="canu",
        input="",
        output="",
        params=f" -p canu -d {canu_output_dir} genomeSize={total_flye_plasmid_length}m maxInputCoverage={coverage} stopOnLowCoverage=1 maxThreads={threads} useGrid=false -{canu_nano_or_pacbio} correctedErrorRate={corrected_error_rate} {longreads}",
        logdir=logdir,
        outfile="",
    )
    try:
        canu.run()
    except subprocess.CalledProcessError as error:
        _report_canu_failure(canu, canu_output_dir, "assemble", error.returncode)


"""
5mer entropy
"""


def shannon_entropy_5mers(dna_sequence):
    """gets 5mer shannon entropy
    :param dna_sequence - str of DNA sequence
    :return: float entropy
    """

    # nucleotides and kmers
    nucleotides = "ACGT"
    k = 4

    # Generate all k-mer combinations
    five_mers_all = ["".join(p) for p in product(nucleotides, repeat=k)]

    kmer_counts = {}
    # total number of k-mers
    total_kmers = 0

    for kmer in five_mers_all:
        count = dna_sequence.count(kmer)
        kmer_counts[kmer] = count
        total_kmers += count

    # Calculate the probability of each k-mer
    probabilities = [
        int(count) / int(total_kmers) for count in kmer_counts.values() if count > 0
    ]

    # Calculate the Shannon entropy using the formula
    entropy = -sum(p * math.log2(p) for p in probabilities)

    # max entropy is log2(4^k). So in this case it is log2(4^4) = log2(256) = 8

    return entropy


"""
determine entropy
"""


def filter_entropy(canu_fasta, outdir):
    """runs canu
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    filtered_records = []

    for record in SeqIO.parse(canu_fasta, "fasta"):
        entropy = shannon_entropy_5mers(str(record.seq))
        if (
            entropy > 5
        ):  # reject low entropy contigs - any real sequence should have entropy more than this (close to 8)
            filtered_records.append(record)

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "canu_filtered_contigs.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fasta")

    return output_filename


def filter_entropy_fastqs(fastq: Path, output_filename: Path) -> None:
    """Filter FASTQ records based on entropy and write the filtered records to a new FASTQ file.

    :param canu_fastq: Input FASTQ file containing long read sequences.
    :param outdir: Output directory for the filtered FASTQ file.
    :return: Path to the filtered FASTQ file.
    """

    filtered_records = []

    for record in SeqIO.parse(fastq, "fastq"):
        entropy = shannon_entropy_5mers(str(record.seq))
        if entropy > 5:
            # Reject low entropy sequences - adjust the threshold as needed.
            filtered_records.append(record)

    # Write the filtered records to a new FASTQ file.
    with open(output_filename, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fastq")


"""
trim contigs
also in trycycler fwiw
https://github.com/rrwick/Trycycler/blob/main/scripts/canu_trim.py
"""


def trim_contigs(canu_filtered_fasta, outdir):
    """trims contigs based on canu trim designation
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    trimmed_records = []

    i = 1
    for record in SeqIO.parse(canu_filtered_fasta, "fasta"):
        # Extract the header
        header = record.description

        # instantiate the start and end coords
        start_coord = 0
        end_coord = len(record.seq)

        # Check if 'trim=' is present in the header
        if "trim=" in header and "suggestCircular=yes" in header:
            # Extract the coordinates from the header
            trim_info = header.split("trim=")[1]
            trim_coords = trim_info.split("-")

            # Parse the start and end coordinates
            start_coord = int(trim_coords[0])
            end_coord = int(trim_coords[1])

        circular_info = header.split("suggestCircular=")[1]
        circular_info = circular_info.split(" ")[0]

        if circular_info == "yes":
            desc = "circular=True"
        else:
            desc = ""

        # Extract the sequence
        subsequence = record.seq[start_coord:end_coord]
        record_trim = SeqRecord(seq=subsequence, id=str(i), description=desc)
        trimmed_records.append(record_trim)
        i += 1

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "canu_filtered_trimmed_contigs.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(trimmed_records, output_handle, "fasta")

    return output_filename
