from pathlib import Path

from plassembler.utils.external_tools import ExternalTool


def build_extra_unicycler_options(unicycler_options, spades_options):
    """Build the extra-options suffix for a unicycler command.

    Each of ``unicycler_options`` / ``spades_options`` may be ``None``
    independently. Regression (issue #83): when only ``spades_options`` was set,
    a ``None`` ``unicycler_options`` was interpolated into the command as the
    literal string ``"None"``, which Unicycler rejected.
    """
    parts = []
    if unicycler_options is not None:
        parts.append(str(unicycler_options))
    if spades_options is not None:
        parts.append(f'--spades_options "{spades_options}"')
    return " ".join(parts)


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

    extra = build_extra_unicycler_options(unicycler_options, spades_options)
    unicycler = ExternalTool(
        tool="unicycler",
        input="",
        output="",
        params=(
            f" -1 {short_one} -2 {short_two} -l {longreads} -t {threads} "
            f"-o {unicycler_output_dir} {extra}"
        ),
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
    :param corrected_longreads: long read fastq (subsampled and corrected with lrcorrect)
    :param entropy_filtered_longreads: long read fastq pre correction - for scaffolding so not subsampled
    :param unicycler_output_dir: unicycler Output Directory
    :param threads: threads
    :param logdir: logdir
    :param unicycler_options: extra unicycler options
    :param spades_options: extra spades options (for unicycler)
    :return:
    """

    extra = build_extra_unicycler_options(unicycler_options, spades_options)
    unicycler_long = ExternalTool(
        tool="unicycler",
        input="",
        output="",
        params=(
            f" -s {corrected_longreads} -l {entropy_filtered_longreads} "
            f"-t {threads} -o {unicycler_output_dir} {extra}"
        ),
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(unicycler_long, to_stdout=False)
