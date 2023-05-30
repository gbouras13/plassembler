from pathlib import Path

from plassembler.utils.external_tools import ExternalTool


def sam_to_bam(sam, bam, threads, logdir):
    """converts sam to bam with samtools
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" view -h -@ {threads} -b {sam}",
        logdir=logdir,
        outfile=bam,
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=True)


def sam_to_sorted_bam(sam, sorted_bam, threads, logdir):
    """converts sam to sorted bam with samtools
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" sort -@ {threads} {sam} -o {sorted_bam}",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=False)


def split_bams(outdir, threads, logdir):
    """
    ensemble function
    """
    non_chrom_bam(outdir, threads, logdir)
    unmapped_bam(outdir, threads, logdir)
    chrom_bam(outdir, threads, logdir)


def non_chrom_bam(outdir, threads, logdir):
    """gets non chrom mapped bam and bed
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    input_bam: Path = Path(outdir) / "short_read.bam"
    non_chrom_bed: Path = Path(outdir) / "non_chromosome.bed"
    non_chrom_bam: Path = Path(outdir) / "non_chromosome.bam"

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" view -b -h -@ {threads} -L {non_chrom_bed} {input_bam}",
        logdir=logdir,
        outfile=non_chrom_bam,
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=True)


def unmapped_bam(outdir, threads, logdir):
    """gets unmapped bam
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    input_bam: Path = Path(outdir) / "short_read.bam"
    unmapped_bam: Path = Path(outdir) / "unmapped_bam_file.bam"

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" view -b -h -f 4 -@ {threads} {input_bam}",
        logdir=logdir,
        outfile=unmapped_bam,
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=True)


def chrom_bam(outdir, threads, logdir):
    """gets chromosome bam and bed
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    input_bam: Path = Path(outdir) / "short_read.bam"
    chrom_bed: Path = Path(outdir) / "chromosome.bed"
    chrom_bam: Path = Path(outdir) / "chromosome.bam"

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" view -b -h -@ {threads} -L {chrom_bed} {input_bam}",
        logdir=logdir,
        outfile=chrom_bam,
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=True)


def bam_to_fastq_short(outdir, threads, logdir):
    """
    ensemble function
    """
    bam_to_fastq_unmapped(outdir, threads, logdir)
    bam_to_fastq_non_chrom(outdir, threads, logdir)


def bam_to_fastq_unmapped(outdir, threads, logdir):
    """gets fastq from unmapped bam
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    unmapped_bam: Path = Path(outdir) / "unmapped_bam_file.bam"
    unmap_fastq_one: Path = Path(outdir) / "unmapped_R1.fastq"
    unmap_fastq_two: Path = Path(outdir) / "unmapped_R2.fastq"

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" fastq -@ {threads} {unmapped_bam} -1 {unmap_fastq_one} -2 {unmap_fastq_two} -0 /dev/null -s /dev/null -n",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=False)


def bam_to_fastq_non_chrom(outdir, threads, logdir):
    """gets fastq from non chrom bam
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    non_chrom_bam: Path = Path(outdir) / "non_chromosome.bam"
    non_chrom_fastq_one: Path = Path(outdir) / "mapped_non_chromosome_R1.fastq"
    non_chrom_fastq_two: Path = Path(outdir) / "mapped_non_chromosome_R2.fastq"

    samtools = ExternalTool(
        tool="samtools",
        input="",
        output="",
        params=f" fastq -@ {threads} {non_chrom_bam} -1 {non_chrom_fastq_one} -2 {non_chrom_fastq_two} -0 /dev/null -s /dev/null -n",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(samtools, to_stdout=False)
