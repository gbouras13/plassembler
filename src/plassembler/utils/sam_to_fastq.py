import os
import subprocess as sp
from collections import defaultdict

import pysam


def extract_long_fastqs_slow_keep_fastqs(out_dir, samname, plasmidname):
    #################################################
    # Define file paths
    #################################################

    # sam_name = os.path.join(out_dir, "long_read.sam")

    # reads mapping to plasmids, or not mapping to any contigs
    plasmidfile = open(plasmidname, "w")

    # Open a FASTQ file for writing reads mapping to multiple contigs
    multimap_plasmid_chromosome_fastqfile = open(
        os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"), "w"
    )

    # chromosome fastqs
    chrom_fastqfile = open(os.path.join(out_dir, "chromosome_mapped_long.fastq"), "w")

    #################################################
    # Get the single and multiple map reads as lists
    #################################################

    # get names, single and multi as lists
    read_names = []
    single_read_names = []
    multi_read_names = []

    # open samfile
    samfile = pysam.AlignmentFile(samname, "r")

    # get list of all read names
    for read in samfile.fetch():
        read_names.append(read.query_name)

    # Create a defaultdict with int as the default factory
    count_dict = defaultdict(int)

    # Loop through the list and count occurrences
    for item in read_names:
        count_dict[item] += 1

    # Get the counts
    for key, value in count_dict.items():
        if value == 1:
            single_read_names.append(key)
        else:
            multi_read_names.append(key)

    samfile.close()

    #################################################
    # process all single reads and then get counts of plasmid vs chromosome ####
    #################################################

    samfile = pysam.AlignmentFile(samname, "r")

    # Create a defaultdict with int as the default  for the multimap reads
    plasmid_mm_dict = defaultdict(int)
    chromosome_mm_dict = defaultdict(int)

    for read in samfile.fetch():
        # Access the read's name, sequence, quality scores, etc.
        read_name = read.query_name
        sequence = read.query_sequence
        quality = read.query_qualities
        flag = read.flag
        # get contig name for the read
        contig_name = samfile.get_reference_name(read.reference_id)

        # print(read_name)
        # print(read.reference_id)
        # print(read.next_reference_id)
        # print(contig_name)
        # print(flag)

        # single reads - easy :)
        if read_name in single_read_names:
            # gets all reads that plasmid mapped reads and all unmapped reads
            if (contig_name and "plasmid" in contig_name) or read.is_unmapped:
                # Write the read to the FASTQ file
                plasmidfile.write(f"@{read_name}\n")  # Write read name
                plasmidfile.write(f"{sequence}\n")  # Write sequence
                plasmidfile.write("+{0}\n".format(read_name))  # Write quality header
                plasmidfile.write(
                    "".join(chr(q + 33) for q in quality) + "\n"
                )  # Write quality scores in ASCII format
            elif contig_name and "chromosome" in contig_name:
                # Write the read to the unmapped reads FASTQ file
                chrom_fastqfile.write(f"@{read_name}\n")  # Write read name
                chrom_fastqfile.write(f"{sequence}\n")  # Write sequence
                chrom_fastqfile.write(
                    "+{0}\n".format(read_name)
                )  # Write quality header
                chrom_fastqfile.write(
                    "".join(chr(q + 33) for q in quality) + "\n"
                )  # Write quality scores in ASCII format
        # create count dictionaries for multimap reads next step
        else:
            if contig_name and "plasmid" in contig_name:
                plasmid_mm_dict[read_name] += 1
            elif contig_name and "chromosome" in contig_name:
                chromosome_mm_dict[read_name] += 1

    samfile.close()

    #################################################
    # process all multimap reads
    #################################################

    samfile = pysam.AlignmentFile(samname, "r")

    for read in samfile.fetch():
        read_name = read.query_name
        sequence = read.query_sequence
        quality = read.query_qualities
        flag = read.flag

        # multireads
        if read_name in multi_read_names:
            if (
                plasmid_mm_dict[read_name] > 0 and chromosome_mm_dict[read_name] > 0
            ):  # multimap both plasmid and chromosome
                if quality is not None and (
                    flag == 0 or flag == 16
                ):  # get only the primary
                    multimap_plasmid_chromosome_fastqfile.write(
                        f"@{read_name}\n"
                    )  # Write read name
                    multimap_plasmid_chromosome_fastqfile.write(
                        f"{sequence}\n"
                    )  # Write sequence
                    multimap_plasmid_chromosome_fastqfile.write(
                        "+{0}\n".format(read_name)
                    )  # Write quality header
                    multimap_plasmid_chromosome_fastqfile.write(
                        "".join(chr(q + 33) for q in quality) + "\n"
                    )  # Write quality scores in ASCII format
            # write all that map to plasmid to the plasmid file
            elif plasmid_mm_dict[read_name] > 0:  # multimap plasmid
                if quality is not None and (
                    flag == 0 or flag == 16
                ):  # get only the primary
                    plasmidfile.write(f"@{read_name}\n")  # Write read name
                    plasmidfile.write(f"{sequence}\n")  # Write sequence
                    plasmidfile.write(
                        "+{0}\n".format(read_name)
                    )  # Write quality header
                    plasmidfile.write(
                        "".join(chr(q + 33) for q in quality) + "\n"
                    )  # Write quality scores in ASCII format
                    # write all that map to chromosome to the plasmid file
            elif chromosome_mm_dict[read_name] > 0:  # multimap chromosome
                if quality is not None and (
                    flag == 0 or flag == 16
                ):  # get only the primary
                    chrom_fastqfile.write(f"@{read_name}\n")  # Write read name
                    chrom_fastqfile.write(f"{sequence}\n")  # Write sequence
                    chrom_fastqfile.write(
                        "+{0}\n".format(read_name)
                    )  # Write quality header
                    chrom_fastqfile.write(
                        "".join(chr(q + 33) for q in quality) + "\n"
                    )  # Write quality scores in ASCII format

    # Close the FASTQ files
    plasmidfile.close()
    chrom_fastqfile.close()
    multimap_plasmid_chromosome_fastqfile.close()

    # Close the SAM file
    samfile.close()


"""
Thanks to @fanvanf
"""


def extract_long_fastqs_fast(sam_name, plasmidfile, threads):
    # sam_name = os.path.join(out_dir, "long_read.sam")
    # plasmidfile = os.path.join(out_dir, "plasmid_long.fastq")

    cmd = f'samtools  view -@ {threads} {sam_name} | awk \'{{if((($3 ~ /plas/)&& ($2 == "0"|| $2 == "16"))||($2 == "4")) print "@"$1"\\n"$10"\\n+"$1"\\n"$11}}\' > {plasmidfile}'
    sp.run(cmd, shell=True)
