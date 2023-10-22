import glob
import os
import shutil

####################################################
# cleanup
##########################################################


def remove_intermediate_files(out_dir, keep_chromosome, assembled_mode, long_only):
    """removes intermediate files
    :param out_dir:  Output Directory
    :return:
    """

    # find all files with the suffix "fastq"
    # find all files with the specified suffixes
    suffixes = ["fastq", "bam", "sa", "sam", "json", "bed", "msh"]
    files = []
    for suffix in suffixes:
        files.extend(glob.glob(os.path.join(out_dir, "*." + suffix)))

    # loop through the files and remove them
    for file in files:
        remove_file(file)

    remove_directory(os.path.join(out_dir, "00-assembly"))
    remove_directory(os.path.join(out_dir, "10-consensus"))
    remove_directory(os.path.join(out_dir, "20-repeat"))
    remove_directory(os.path.join(out_dir, "30-contigger"))
    remove_directory(os.path.join(out_dir, "40-polishing"))

    if assembled_mode is True:
        remove_directory(os.path.join(out_dir, "flye_output"))
        remove_file(os.path.join(out_dir, "plassembler_plasmids.fasta"))
        remove_file(os.path.join(out_dir, "plassembler_plasmids.gfa"))

    # delete intermediate mash file
    remove_file(os.path.join(out_dir, "mash.tsv"))

    # delete intermediate fasta assemble files
    remove_file(os.path.join(out_dir, "combined.fasta"))
    remove_file(os.path.join(out_dir, "flye_renamed.fasta"))
    remove_file(os.path.join(out_dir, "plasmids.fasta"))

    # delete fastq intermediate files
    remove_file(os.path.join(out_dir, "chopper_long_reads.fastq.gz"))
    remove_file(os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"))

    # long
    remove_file(os.path.join(out_dir, "plasmids_canu.fasta"))
    remove_file(os.path.join(out_dir, "plasmid_long.fastq"))

    # chromosome
    if keep_chromosome is False:
        remove_file(os.path.join(out_dir, "chromosome.fasta"))

    # canu
    # leave the canu directory for long
    # remove_directory(os.path.join(out_dir, "canu"))
    remove_file(os.path.join(out_dir, "long_combined.fasta"))
    remove_file(os.path.join(out_dir, "canu_filtered_contigs.fasta"))
    remove_file(os.path.join(out_dir, "canu_filtered_trimmed_contigs.fasta"))

    # remove dnaapler for long
    remove_directory(os.path.join(out_dir, "dnaapler"))


def move_and_copy_files(
    out_dir,
    prefix,
    unicycler_success_flag,
    keep_fastqs,
    assembled_mode,
    long_only,
    use_raven,
    skip_assembly,
    canu_flag,
):
    """moves and copies files
    :param out_dir:  Output Directory
    :param prefix: prefix
    :param unicycler_success_flag: whether or not unicycler worked
    :param assembled_mode: whether or not unicycler worked
    :param long_only: whether or not unicycler worked
    :param use_raven: whether or not unicycler worked
    :param skip_assembly: --flye_directory specified
    :param canu_flag: --canu_flag is True
    :return:
    """

    # move the flye outputs
    if assembled_mode is False:
        if use_raven is False:
            # make flye dir
            flye_dir = os.path.join(out_dir, "flye_output")
            if not os.path.exists(flye_dir):
                os.mkdir(flye_dir)
            shutil.move(os.path.join(out_dir, "assembly.fasta"), flye_dir)
            shutil.move(os.path.join(out_dir, "assembly_info.txt"), flye_dir)
            if skip_assembly is False:
                shutil.move(os.path.join(out_dir, "flye.log"), flye_dir)
                shutil.move(os.path.join(out_dir, "assembly_graph.gfa"), flye_dir)
                shutil.move(os.path.join(out_dir, "assembly_graph.gv"), flye_dir)
        else:
            # make raven dir
            raven_dir = os.path.join(out_dir, "raven_output")
            if not os.path.exists(raven_dir):
                os.mkdir(raven_dir)
            # move gfa and
            shutil.move(os.path.join(out_dir, "assembly.fasta"), raven_dir)
            shutil.move(os.path.join(out_dir, "assembly_graph.gfa"), raven_dir)

    if unicycler_success_flag is True:
        if canu_flag is False:
            # move unicycler graph output to main directory
            shutil.copy2(
                os.path.join(out_dir, "unicycler_output", "assembly.gfa"),
                os.path.join(out_dir, prefix + "_plasmids.gfa"),
            )
    else:
        # to touch empty versions of the output files if no plasmids
        touch_output_fail_files(out_dir, prefix)

    # put kept fastqs into separate directory
    # make fastqs dir
    if keep_fastqs is True:
        fastqs_dir = os.path.join(out_dir, "plasmid_fastqs")
        if not os.path.exists(fastqs_dir):
            os.mkdir(fastqs_dir)

        # always move plasmid_long
        shutil.move(
            os.path.join(out_dir, "plasmid_long.fastq"),
            os.path.join(fastqs_dir, "plasmids_long.fastq"),
        )

        if long_only is False:  # for the hybrid only
            # move fastqs
            shutil.move(
                os.path.join(out_dir, "short_read_concat_R1.fastq"),
                os.path.join(fastqs_dir, "plasmids_R1.fastq"),
            )
            shutil.move(
                os.path.join(out_dir, "short_read_concat_R2.fastq"),
                os.path.join(fastqs_dir, "plasmids_R2.fastq"),
            )

            shutil.move(
                os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"),
                os.path.join(fastqs_dir, "multimap_long.fastq"),
            )


# function to touch create a file
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


# to create empty plasmids fasta and gfa files
def touch_output_fail_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_plasmids.fasta"))
    touch_file(os.path.join(out_dir, prefix + "_plasmids.gfa"))
    touch_file(os.path.join(out_dir, prefix + "_summary.tsv"))


def touch_output_fail_files_long(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_plasmids.fasta"))
    touch_file(os.path.join(out_dir, prefix + "_summary.tsv"))


def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
