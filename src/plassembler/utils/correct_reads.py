import math
import subprocess
import sys
from itertools import product
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from plassembler.utils.external_tools import ExternalTool

# corrected reads are subsampled to this many times the summed length of Flye's
# non-chromosomal contigs, so as not to give Unicycler too much read depth
CORRECTION_COVERAGE = 50
# when Flye found no plasmid contigs there is no length to aim for: aim for 250x
# of 10 kbp instead, enough for a small plasmid Flye missed
NO_PLASMID_GENOME_SIZE = 10000
NO_PLASMID_COVERAGE = 250


def correction_read_type(pacbio_model):
    """Map a validated ``--pacbio_model`` value to lrcorrect's ``--preset`` and
    whether read correction should be skipped.

    ``pacbio_model`` is the output of ``validate_pacbio_model`` (so one of
    ``"--pacbio-raw"``, ``"--pacbio-corr"``, ``"--pacbio-hifi"``) or
    ``"nothing"``.

    PacBio HiFi reads are already high-accuracy, so they are not corrected
    (issue #84). Previously the check compared against ``"pacbio-hifi"`` (without
    the ``--`` prefix), so HiFi reads were wrongly corrected as regular PacBio.

    :return: ``(lrcorrect_preset, skip_correction)``
    """
    if pacbio_model == "--pacbio-hifi":
        return "pb", True
    if pacbio_model != "nothing":
        return "pb", False
    return "ont", False


def correction_target_bases(total_flye_plasmid_length):
    """How many bases of reads to correct and hand to Unicycler.

    :param total_flye_plasmid_length: summed length of Flye's non-chromosomal contigs
    :return: the target number of bases for lrcorrect's ``--target-bases``
    """
    if total_flye_plasmid_length == 0:
        return NO_PLASMID_GENOME_SIZE * NO_PLASMID_COVERAGE
    return total_flye_plasmid_length * CORRECTION_COVERAGE


def _lrcorrect_stderr_tail(err_log, max_lines=20):
    """Last few lines of lrcorrect's stderr, ready to embed in a message.

    ExternalTool's "Command line:" header is dropped: it is flushed when the log
    closes, so it can land after lrcorrect's own output.

    :param err_log: path of the lrcorrect .err logfile
    :param max_lines: how many trailing lines to keep
    :return: the trailing lines, or "" if the log is missing or empty
    """
    try:
        log_text = Path(err_log).read_text(errors="replace")
    except OSError:
        return ""
    lines = [
        line
        for line in log_text.strip().splitlines()
        if not line.startswith("Command line:")
    ]
    return "\n".join(lines[-max_lines:])


def _fail_lrcorrect(reason, err_log):
    """Reports a fatal lrcorrect problem and stops the run.

    Everything goes into a single `logger.error`: under the CLI that sink exits
    the process, so a second call would never be reached.

    :param reason: what went wrong, e.g. "lrcorrect (return code 1)"
    :param err_log: path of the lrcorrect .err logfile
    """
    message = f"Error with lrcorrect: {reason}. Please check {err_log}"
    stderr_tail = _lrcorrect_stderr_tail(err_log)
    if stderr_tail:
        message = f"{message}\nlrcorrect stderr:\n{stderr_tail}"
    logger.error(message)
    # the CLI exits inside logger.error above; this covers being driven as a
    # library, where no ERROR sink is installed
    sys.exit(1)


def run_lrcorrect(
    threads, logdir, longreads, corrected_fastq, preset, target_bases
) -> bool:
    """Corrects long reads with lrcorrect, subsampled to ``target_bases``.

    lrcorrect exits 0 and writes no reads when there is too little data to
    correct anything - few reads are left once the chromosome is removed, most
    often because the sample has no plasmids. That is reported and the caller
    carries on with uncorrected reads. Any other failure is fatal: before, a
    broken corrector looked exactly like a sample with no plasmids (#92).

    :param threads: threads
    :param logdir: logdir
    :param longreads: entropy-filtered non-chromosomal long reads
    :param corrected_fastq: where to write the corrected reads
    :param preset: "ont" or "pb", from correction_read_type
    :param target_bases: from correction_target_bases
    :return: True if corrected reads were written, False if there were none
    """
    lrcorrect = ExternalTool(
        tool="lrcorrect",
        input=f"{longreads}",
        output=f"-o {corrected_fastq}",
        params=f"--threads {threads} --preset {preset} --target-bases {target_bases} --qual-char #",
        logdir=logdir,
        outfile="",
    )
    try:
        lrcorrect.run()
    except subprocess.CalledProcessError as error:
        _fail_lrcorrect(
            f"lrcorrect (return code {error.returncode})", lrcorrect.err_log
        )
    except OSError as error:
        _fail_lrcorrect(str(error), lrcorrect.err_log)

    if Path(corrected_fastq).stat().st_size > 0:
        return True
    logger.warning("lrcorrect could not correct any reads, as there were too few.")
    logger.warning(
        "This usually means few reads are left once the chromosome is "
        "removed, most often because the sample has no plasmids."
    )
    stderr_tail = _lrcorrect_stderr_tail(lrcorrect.err_log, max_lines=3)
    if stderr_tail:
        logger.warning(f"lrcorrect stderr:\n{stderr_tail}")
    return False


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


def filter_entropy_fastqs(fastq: Path, output_filename: Path) -> None:
    """Filter FASTQ records based on entropy and write the filtered records to a new FASTQ file.

    :param fastq: Input FASTQ file containing long read sequences.
    :param output_filename: Path of the filtered FASTQ file to write.
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
