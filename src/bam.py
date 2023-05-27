from src.external_tools import ExternalTool
from pathlib import Path

# sam to bam


def sam_to_bam_short(outdir, threads, logdir):
    """converts sam to bam with samtools
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    sam: Path =  outdir/ f"short_read.sam"
    bam: Path =  outdir/ f"short_read.bam"

    minimap2 = ExternalTool(
        tool="samtools",
        input=f"",
        output=f"",
        params=f" view -h -@ {threads} -b {sam}",
        logdir=logdir,
        outfile = bam
    )

    # need to write to stdout
    ExternalTool.run_tool(minimap2, to_stdout = True)



def split_bams(outdir, threads, logger):
    """gets all relevant fastqs from bams
    :param long_short: "long" or "short"
        :param outdir: output directory path
    :param threads: threads
    :param logger: logger
    :return:
    """

    non_chrom_bed = os.path.join(outdir, "non_chromosome.bed")
    chrom_bed = os.path.join(outdir, "chromosome.bed")
    bam = os.path.join(outdir, "short_read.bam")

    unmapped_bam_file = os.path.join(outdir, "unmapped_bam_file.bam")
    non_chrom_bam_file = os.path.join(outdir, "non_chromosome.bam")
    chrom_bam_file = os.path.join(outdir, "chromosome.bam")

    non_chrom_bam = open(non_chrom_bam_file, "w")
    chrom_bam = open(chrom_bam_file, "w")
    unmapped_bam = open(unmapped_bam_file, "w")
    try:
        sam_to_bam = sp.Popen(
            ["samtools", "view", "-b", "-h", "-f", "4", "-@", threads, bam],
            stdout=unmapped_bam,
            stderr=sp.PIPE,
        )
        log.write_to_log(sam_to_bam.stderr, logger)
        sam_to_bam_non_chrom = sp.Popen(
            ["samtools", "view", "-b", "-h", "-@", threads, "-L", non_chrom_bed, bam],
            stdout=non_chrom_bam,
            stderr=sp.PIPE,
        )
        log.write_to_log(sam_to_bam_non_chrom.stderr, logger)
        sam_to_bam_chrom = sp.Popen(
            ["samtools", "view", "-b", "-h", "-@", threads, "-L", chrom_bed, bam],
            stdout=chrom_bam,
            stderr=sp.PIPE,
        )
        log.write_to_log(sam_to_bam_chrom.stderr, logger)
    except:
        sys.exit("Error with samtools view.\n")


def bam_to_fastq_short(outdir, threads, logger):
    """gets all relevant fastqs from bams
    :param long_short: "long" or "short"
        :param outdir: output directory path
    :param threads: threads
    :param logger: logger
    :return:
    """

    # reads that don't map to chromosome
    try:
        unmapped_bam_file = os.path.join(outdir, "unmapped_bam_file.bam")
        non_chrom_bam_file = os.path.join(outdir, "non_chromosome.bam")

        unmap_fastq_one = os.path.join(outdir, "unmapped_R1.fastq")
        unmap_fastq_two = os.path.join(outdir, "unmapped_R2.fastq")

        non_chrom_fastq_one = os.path.join(outdir, "mapped_non_chromosome_R1.fastq")
        non_chrom_fastq_two = os.path.join(outdir, "mapped_non_chromosome_R2.fastq")

        # chrom_bam_file = os.path.join(outdir, "chromosome.bam")
        # chrom_fastq_one = os.path.join(outdir, "mapped_chromosome_R1.fastq")
        # chrom_fastq_two = os.path.join(outdir, "mapped_chromosome_R2.fastq")

        extract_unmap_short_fastq = sp.Popen(
            [
                "samtools",
                "fastq",
                "-@",
                threads,
                unmapped_bam_file,
                "-1",
                unmap_fastq_one,
                "-2",
                unmap_fastq_two,
                "-0",
                "/dev/null",
                "-s",
                "/dev/null",
                "-n",
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
        )
        log.write_to_log(extract_unmap_short_fastq.stdout, logger)

        extract_non_chrom_short_fastq = sp.Popen(
            [
                "samtools",
                "fastq",
                "-@",
                threads,
                "-F",
                "4",
                non_chrom_bam_file,
                "-1",
                non_chrom_fastq_one,
                "-2",
                non_chrom_fastq_two,
                "-0",
                "/dev/null",
                "-s",
                "/dev/null",
                "-n",
            ],
            stdout=sp.PIPE,
            stderr=sp.PIPE,
        )
        log.write_to_log(extract_non_chrom_short_fastq.stdout, logger)

        # extract_chrom_short_fastq = sp.Popen(["samtools", "fastq", "-@", threads, '-F', "4", '-f', '1', chrom_bam_file, "-1", chrom_fastq_one, "-2", chrom_fastq_two, "-0", "/dev/null", "-s", "/dev/null", "-n"], stdout=sp.PIPE, stderr=sp.PIPE)
        # log.write_to_log(extract_chrom_short_fastq.stdout, logger)

    except:
        sys.exit("Error with samtools fastq.\n")
