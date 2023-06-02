from plassembler.utils.external_tools import ExternalTool


def run_unicycler(
    threads, logdir, short_one, short_two, longreads, unicycler_output_dir
):
    """runs Unicycler
    :param short_one: R1 short read fastq
    :param short_two: R2 short read fastq
    :param long: long read fastq
    :param unicycler_output_dir: unicycler Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    unicycler = ExternalTool(
        tool="unicycler",
        input="",
        output="",
        params=f" -1 {short_one} -2 {short_two} -l {longreads} -t {threads} -o {unicycler_output_dir}",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(unicycler, to_stdout=False)
