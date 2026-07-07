import subprocess as sp
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

import pysam


def extract_long_fastqs_slow_keep_fastqs(out_dir, samname, plasmidname):
    #################################################
    # Get the single and multiple map reads as sets
    #################################################

    # get list of all read names
    read_names = []
    with pysam.AlignmentFile(samname, "r") as samfile:
        for read in samfile.fetch():
            read_names.append(read.query_name)

    # count occurrences of each read name
    count_dict = defaultdict(int)
    for item in read_names:
        count_dict[item] += 1

    # sets give O(1) membership (the per-read lookups below run once per
    # alignment, so lists here would make the whole function quadratic)
    single_read_names = {name for name, count in count_dict.items() if count == 1}
    multi_read_names = {name for name, count in count_dict.items() if count != 1}

    # ExitStack guarantees all output handles are closed even if a write or a
    # pysam call raises partway through
    with ExitStack() as stack:
        # reads mapping to plasmids, or not mapping to any contigs
        plasmidfile = stack.enter_context(open(plasmidname, "w"))
        # reads mapping to multiple contigs
        multimap_plasmid_chromosome_fastqfile = stack.enter_context(
            open(Path(out_dir) / "multimap_plasmid_chromosome_long.fastq", "w")
        )
        # chromosome fastqs
        chrom_fastqfile = stack.enter_context(
            open(Path(out_dir) / "chromosome_mapped_long.fastq", "w")
        )

        #################################################
        # process all single reads and count plasmid vs chromosome multimaps
        #################################################

        plasmid_mm_dict = defaultdict(int)
        chromosome_mm_dict = defaultdict(int)

        with pysam.AlignmentFile(samname, "r") as samfile:
            for read in samfile.fetch():
                read_name = read.query_name
                sequence = read.query_sequence
                quality = read.query_qualities
                # get contig name for the read
                contig_name = samfile.get_reference_name(read.reference_id)

                # single reads - easy :)
                if read_name in single_read_names:
                    # plasmid-mapped reads and all unmapped reads
                    if (contig_name and "plasmid" in contig_name) or read.is_unmapped:
                        plasmidfile.write(f"@{read_name}\n")
                        plasmidfile.write(f"{sequence}\n")
                        plasmidfile.write(f"+{read_name}\n")
                        plasmidfile.write("".join(chr(q + 33) for q in quality) + "\n")
                    elif contig_name and "chromosome" in contig_name:
                        chrom_fastqfile.write(f"@{read_name}\n")
                        chrom_fastqfile.write(f"{sequence}\n")
                        chrom_fastqfile.write(f"+{read_name}\n")
                        chrom_fastqfile.write(
                            "".join(chr(q + 33) for q in quality) + "\n"
                        )
                # build count dictionaries for the multimap reads (next step)
                else:
                    if contig_name and "plasmid" in contig_name:
                        plasmid_mm_dict[read_name] += 1
                    elif contig_name and "chromosome" in contig_name:
                        chromosome_mm_dict[read_name] += 1

        #################################################
        # process all multimap reads
        #################################################

        with pysam.AlignmentFile(samname, "r") as samfile:
            for read in samfile.fetch():
                read_name = read.query_name
                sequence = read.query_sequence
                quality = read.query_qualities
                flag = read.flag

                if read_name in multi_read_names:
                    # multimap to both plasmid and chromosome
                    if (
                        plasmid_mm_dict[read_name] > 0
                        and chromosome_mm_dict[read_name] > 0
                    ):
                        if quality is not None and (flag == 0 or flag == 16):
                            # get only the primary
                            multimap_plasmid_chromosome_fastqfile.write(
                                f"@{read_name}\n"
                            )
                            multimap_plasmid_chromosome_fastqfile.write(f"{sequence}\n")
                            multimap_plasmid_chromosome_fastqfile.write(
                                f"+{read_name}\n"
                            )
                            multimap_plasmid_chromosome_fastqfile.write(
                                "".join(chr(q + 33) for q in quality) + "\n"
                            )
                    # multimap to plasmid only -> plasmid file
                    elif plasmid_mm_dict[read_name] > 0:
                        if quality is not None and (flag == 0 or flag == 16):
                            plasmidfile.write(f"@{read_name}\n")
                            plasmidfile.write(f"{sequence}\n")
                            plasmidfile.write(f"+{read_name}\n")
                            plasmidfile.write(
                                "".join(chr(q + 33) for q in quality) + "\n"
                            )
                    # multimap to chromosome only -> chromosome file
                    elif chromosome_mm_dict[read_name] > 0:
                        if quality is not None and (flag == 0 or flag == 16):
                            chrom_fastqfile.write(f"@{read_name}\n")
                            chrom_fastqfile.write(f"{sequence}\n")
                            chrom_fastqfile.write(f"+{read_name}\n")
                            chrom_fastqfile.write(
                                "".join(chr(q + 33) for q in quality) + "\n"
                            )


"""
Thanks to @fanvanf
"""


def extract_long_fastqs_fast(sam_name, plasmidfile, threads):
    cmd = f'samtools  view -@ {threads} {sam_name} | awk \'{{if((($3 ~ /plas/)&& ($2 == "0"|| $2 == "16"))||($2 == "4")) print "@"$1"\\n"$10"\\n+"$1"\\n"$11}}\' > {plasmidfile}'
    # shell=True is required for the samtools | awk pipeline; check=True surfaces
    # a failing samtools instead of silently leaving an empty plasmid FASTQ.
    sp.run(cmd, shell=True, check=True)
