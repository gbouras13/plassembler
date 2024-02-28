import gzip
from pathlib import Path

from Bio import SeqIO

from plassembler.utils.external_tools import ExternalTool


def run_unicycler(
    threads: int,
    logdir: Path,
    short_one: Path,
    short_two: Path,
    longreads: Path,
    unicycler_output_dir: Path,
    unicycler_options: str,
    spades_options: str,
):
    """runs Unicycler
    :param short_one: R1 short read fastq
    :param short_two: R2 short read fastq
    :param long: long read fastq
    :param unicycler_output_dir: unicycler Output Directory
    :param threads: threads
    :param logdir: logdir
    :param unicycler_options: extra unicycler options
    :param spades_options: extra spades options (for unicycler)
    :return:
    """

    if unicycler_options is None and spades_options is None:
        unicycler = ExternalTool(
            tool="unicycler",
            input="",
            output="",
            params=f" -1 {short_one} -2 {short_two} -l {longreads} -t {threads} -o {unicycler_output_dir}",
            logdir=logdir,
            outfile="",
        )
    else:
        if spades_options is None:
            unicycler = ExternalTool(
                tool="unicycler",
                input="",
                output="",
                params=f" -1 {short_one} -2 {short_two} -l {longreads} -t {threads} -o {unicycler_output_dir} {unicycler_options}",
                logdir=logdir,
                outfile="",
            )
        else:
            unicycler = ExternalTool(
                tool="unicycler",
                input="",
                output="",
                params=f' -1 {short_one} -2 {short_two} -l {longreads} -t {threads} -o {unicycler_output_dir} {unicycler_options} --spades_options "{spades_options}" ',
                logdir=logdir,
                outfile="",
            )

    ExternalTool.run_tool(unicycler, to_stdout=False)


def run_unicycler_long(
    threads: int,
    logdir: Path,
    corrected_longreads: Path,
    entropy_filtered_longreads: Path,
    unicycler_output_dir: Path,
    unicycler_options: str,
    spades_options: str,
) -> None:
    """runs Unicycler on long reads with -s -l
    :param corrected_longreads: long read fastq (subsmapled and corrected with canu)
    :param entropy_filtered_longreads: long read fastq pre correction - for scaffolding so not subsampled
    :param unicycler_output_dir: unicycler Output Directory
    :param threads: threads
    :param logdir: logdir
    :param unicycler_options: extra unicycler options
    :param spades_options: extra spades options (for unicycler)
    :return:
    """

    if unicycler_options is None and spades_options is None:
        unicycler_long = ExternalTool(
            tool="unicycler",
            input="",
            output="",
            params=f" -s {corrected_longreads} -l {entropy_filtered_longreads} -t {threads} -o {unicycler_output_dir}",
            logdir=logdir,
            outfile="",
        )
    else:
        if spades_options is None:
            unicycler_long = ExternalTool(
                tool="unicycler",
                input="",
                output="",
                params=f" -s {corrected_longreads} -l {entropy_filtered_longreads} -t {threads} -o {unicycler_output_dir} {unicycler_options}",
                logdir=logdir,
                outfile="",
            )
        else:
            unicycler_long = ExternalTool(
                tool="unicycler",
                input="",
                output="",
                params=f' -s {corrected_longreads} -l {entropy_filtered_longreads} -t {threads} -o {unicycler_output_dir} {unicycler_options} --spades_options "{spades_options}"',
                logdir=logdir,
                outfile="",
            )

    ExternalTool.run_tool(unicycler_long, to_stdout=False)


def corrected_fasta_to_fastq(input_fasta_gz: Path, output_fastq: Path) -> None:
    """
    Convert a gzipped FASTA file to FASTQ format with '#' as the quality scores and save it to a FASTQ file.

    Args:
        input_fasta_gz (Path): Path to the gzipped FASTA input file.
        output_fastq (Path): Path to the output FASTQ file.

    Returns:
        None
    """
    with gzip.open(input_fasta_gz, "rt") as fasta_file:
        with open(output_fastq, "w") as fastq_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = record.id
                sequence = str(record.seq)
                quality = "#" * len(sequence)
                fastq_record = f"@{header}\n{sequence}\n+\n{quality}\n"
                fastq_file.write(fastq_record)
