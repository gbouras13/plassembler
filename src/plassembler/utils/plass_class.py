import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from plassembler.utils.bam import sam_to_sorted_bam
from plassembler.utils.depth import (
    collate_depths,
    combine_depth_dfs,
    concatenate_chrom_plasmids,
    depth_df_single,
    get_contig_circularity,
    get_contig_lengths,
    get_depths_from_bam,
)
from plassembler.utils.mapping import minimap_long_reads, minimap_short_reads
from plassembler.utils.run_mash import get_contig_count, is_file_empty


class Plass:
    """Plassembler Output Class"""

    def __init__(
        self,
        outdir: str = "output/",
        contig_count: int = 1,
        no_plasmids_flag: bool = False,
        chromosome_flag: bool = True,
        threads: int = 1,
        depth_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        mash_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        combined_depth_mash_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        long_only: bool = False,
        unicycler_success: bool = True,
    ) -> None:
        """
        Parameters
        --------
        outdir: int, required
            output directory
        contig_count: int, required
            the number of contigs assembled by flye assembly
        no_plasmids_flag: bool, required
            flag determining whether there are plasmids at some point in the pipeline. If False means there are plasmids
        chromosome_flag: bool, required
            flag saying whether or not there was an identified chromosome in the Flye output.
        depth_df: Pandas df, required
            dataframe with depth output
        mash_df: Pandas df, required
            dataframe with mash output
        threads: int, required
            integer giving args.threads - defaults to 1.
        long_only: bool, required
            whether plassembler is in kmer mode
        unicycler_success: bool, required
            whether unicycler succeeded
        """
        self.outdir = outdir
        self.contig_count = contig_count
        self.no_plasmids_flag = no_plasmids_flag
        self.chromosome_flag = chromosome_flag
        self.threads = threads
        self.depth_df = depth_df
        self.mash_df = mash_df
        self.combined_depth_mash_df = combined_depth_mash_df
        self.long_only = long_only
        self.unicycler_success = unicycler_success

    def get_contig_count(self):
        """Counts the number of contigs assembled
        :return:
        """
        outdir = self.outdir
        fasta_file = os.path.join(outdir, "assembly.fasta")
        contig_count = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig_count += 1
        logger.info(f"Assembled {contig_count} contigs.")
        self.contig_count = contig_count

    def identify_chromosome_process_raven(self, chromosome_len):
        """Identified chromosome and processes Raven output - renames chromosome contig and the others as plasmid_1, plasmid_2 etc
        Also makes the chromosome bed file for downstream samtools mapping
        :param outdir: output directory
        :param chromosome_len: lower bound on length of chromosome from input command
        :param keep_chromosome: whether the user wants to keep the chromosome as chromosome.fasta
        :return chromosome_flag: bool whether chromosome assembles
        """
        outdir = self.outdir
        long_only = self.long_only

        fasta_file = os.path.join(outdir, "assembly.fasta")

        # get max contig length
        max_length = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_length = len(record.seq)
            if sequence_length > max_length:
                max_length = sequence_length

        # to say that the chromosome has been correctly identified
        chromosome_flag = True
        if max_length < int(
            chromosome_len
        ):  # no chromosome identified -> don't bother with the renaming
            chromosome_flag = False
        # assuming chromosome identified
        else:
            # make bed file with plasmid contigs to extract mapping reads
            bed_file = open(os.path.join(outdir, "non_chromosome.bed"), "w")
            bed_chrom_file = open(os.path.join(outdir, "chromosome.bed"), "w")
            with open(os.path.join(outdir, "flye_renamed.fasta"), "w") as rename_fa:
                # for chromosome fasta too
                with open(os.path.join(outdir, "chromosome.fasta"), "w") as chrom_fa:
                    # for plasmid numbering
                    i = 1
                    # for chromosome numbering (chromids, multiple plasmid contigs)
                    c = 1
                    for dna_record in SeqIO.parse(
                        os.path.join(outdir, "assembly.fasta"), "fasta"
                    ):
                        # if the length is over chromosome length
                        contig_len = len(dna_record.seq)
                        if contig_len > int(chromosome_len):
                            # if the first chromosome (most cases), then just chromosome. Otherwise continue.
                            if c == 1:
                                dna_header = "chromosome"
                            else:
                                if c == "2":
                                    message = "Multiple contigs above the specified chromosome length -c have been detected. \nIf you are hoping for plasmids from haploid bacteria, please check your value for -c."
                                    logger.info(message)
                                dna_header = "chromosome_" + str(c)
                            dna_description = ""
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            SeqIO.write(dna_record, chrom_fa, "fasta")
                            bed_chrom_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # chromosome
                            c += 1
                        # plasmids
                        else:
                            dna_header = "plasmid_" + str(i)
                            dna_description = ""
                            # write the updated record
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            # get length for bed file
                            # make bed file contig leng
                            bed_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # Write read name
                            i += 1
            # if lony only is true, create new unicycler_output file (fake)
            if long_only is True:
                message = "Extracting possible plasmids from Raven assembly."
                logger.info(message)
                # remove bed if exists
                if os.path.exists(os.path.join(outdir, "non_chromosome.bed")):
                    os.remove(os.path.join(outdir, "non_chromosome.bed"))
                with open(
                    os.path.join(outdir, "plasmids_initial.fasta"), "w"
                ) as rename_fa:
                    # for plasmid numbering
                    i = 1
                    for dna_record in SeqIO.parse(
                        os.path.join(outdir, "assembly.fasta"), "fasta"
                    ):
                        # if the length is over chromosome length
                        contig_len = len(dna_record.seq)
                        if contig_len < int(chromosome_len):
                            dna_header = str(i)
                            # write the updated record
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            # get length for bed file
                            # make bed file
                            bed_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # Write read name
                            i += 1
        # add to object
        self.chromosome_flag = chromosome_flag

    def identify_chromosome_process_flye(self, chromosome_len):
        """Identified chromosome and processes Flye output - renames chromosome contig and the others as plasmid_1, plasmid_2 etc
        Also makes the chromosome bed file for downstream samtools mapping
        :param outdir: output directory
        :param chromosome_len: lower bound on length of chromosome from input command
        :param keep_chromosome: whether the user wants to keep the chromosome as chromosome.fasta
        :return chromosome_flag: bool whether chromosome assembles
        """
        outdir = self.outdir
        info_file = os.path.join(outdir, "assembly_info.txt")
        col_list = [
            "seq_name",
            "length",
            "cov",
            "circ",
            "repeat",
            "mult",
            "alt_group",
            "graph_path",
        ]
        info_df = pd.read_csv(
            info_file, delimiter="\t", index_col=False, names=col_list, skiprows=1
        )
        max_length = max(info_df["length"])

        # to say that the chromosome has been correctly identified
        chromosome_flag = True
        if max_length < int(
            chromosome_len
        ):  # no chromosome identified -> don't bother with the renaming
            chromosome_flag = False
        # assuming chromosome identified
        else:
            # make bed file with plasmid contigs to extract mapping reads
            bed_file = open(os.path.join(outdir, "non_chromosome.bed"), "w")
            bed_chrom_file = open(os.path.join(outdir, "chromosome.bed"), "w")
            with open(os.path.join(outdir, "flye_renamed.fasta"), "w") as rename_fa:
                # for chromosome fasta too
                with open(os.path.join(outdir, "chromosome.fasta"), "w") as chrom_fa:
                    # for plasmid numbering
                    i = 1
                    # for chromosome numbering (chromids, multiple plasmid contigs)
                    c = 1
                    for dna_record in SeqIO.parse(
                        os.path.join(outdir, "assembly.fasta"), "fasta"
                    ):
                        # if the length is over chromosome length
                        contig_len = len(dna_record.seq)
                        if contig_len > int(chromosome_len):
                            # if the first chromosome (most cases), then just chromosome. Otherwise continue.
                            if c == 1:
                                dna_header = "chromosome"
                            else:
                                if c == "2":
                                    message = "Multiple contigs above the specified chromosome length -c have been detected. \nIf you are hoping for plasmids from haploid bacteria, please check your value for -c."
                                    logger.info(message)
                                dna_header = "chromosome_" + str(c)
                            dna_description = ""
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            SeqIO.write(dna_record, chrom_fa, "fasta")
                            bed_chrom_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # chromosome
                            c += 1
                        # plasmids
                        else:
                            dna_header = "plasmid_" + str(i)
                            dna_description = ""

                            # write the updated record
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            # get length for bed file
                            # make bed file
                            bed_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # Write read name
                            i += 1

        # add to object
        self.chromosome_flag = chromosome_flag

    def identify_chromosome_process_flye_long(self, chromosome_len):
        """Identified chromosome and processes Flye output long only - renames chromosome contig and the others as plasmid_1, plasmid_2 etc
        Also makes the chromosome bed file for downstream samtools mapping
        :param outdir: output directory
        :param chromosome_len: lower bound on length of chromosome from input command
        :param keep_chromosome: whether the user wants to keep the chromosome as chromosome.fasta
        :return chromosome_flag: bool whether chromosome assembles
        """
        outdir = self.outdir
        info_file = os.path.join(outdir, "assembly_info.txt")
        col_list = [
            "seq_name",
            "length",
            "cov",
            "circ",
            "repeat",
            "mult",
            "alt_group",
            "graph_path",
        ]
        info_df = pd.read_csv(
            info_file, delimiter="\t", index_col=False, names=col_list, skiprows=1
        )
        max_length = max(info_df["length"])

        # to say that the chromosome has been correctly identified
        chromosome_flag = True
        if max_length < int(
            chromosome_len
        ):  # no chromosome identified -> don't bother with the renaming
            chromosome_flag = False
        # assuming chromosome identified
        else:
            # make bed file with plasmid contigs to extract mapping reads
            bed_file = open(os.path.join(outdir, "non_chromosome.bed"), "w")
            bed_chrom_file = open(os.path.join(outdir, "chromosome.bed"), "w")
            with open(os.path.join(outdir, "flye_renamed.fasta"), "w") as rename_fa:
                # for chromosome fasta too
                with open(os.path.join(outdir, "chromosome.fasta"), "w") as chrom_fa:
                    # for plasmid numbering
                    i = 1
                    # for chromosome numbering (chromids, multiple plasmid contigs)
                    c = 1
                    for dna_record in SeqIO.parse(
                        os.path.join(outdir, "assembly.fasta"), "fasta"
                    ):
                        # if the length is over chromosome length
                        contig_len = len(dna_record.seq)
                        if contig_len > int(chromosome_len):
                            # if the first chromosome (most cases), then just chromosome. Otherwise continue.
                            if c == 1:
                                dna_header = "chromosome"
                            else:
                                if c == "2":
                                    message = "Multiple contigs above the specified chromosome length -c have been detected. \nIf you are hoping for plasmids from haploid bacteria, please check your value for -c."
                                    logger.info(message)
                                dna_header = "chromosome_" + str(c)
                            dna_description = ""
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            SeqIO.write(dna_record, chrom_fa, "fasta")
                            bed_chrom_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # chromosome
                            c += 1
                        # plasmids
                        else:
                            dna_header = "plasmid_" + str(i)
                            dna_description = ""
                            # write the updated record
                            dna_record = SeqRecord(
                                dna_record.seq,
                                id=dna_header,
                                description=dna_description,
                            )
                            SeqIO.write(dna_record, rename_fa, "fasta")
                            # get length for bed file
                            # make bed file
                            bed_file.write(
                                f"{dna_header}\t1\t{contig_len}\n"
                            )  # Write read name
                            i += 1
        # add to object
        self.chromosome_flag = chromosome_flag

    def check_unicycler_success(self, unicycler_dir):
        unicycler_file: Path = Path(unicycler_dir, "assembly.fasta")
        # check if unicycler succeded according to the output (it won't if no plasmids)
        unicycler_success = os.path.isfile(unicycler_file)
        if (
            os.path.isfile(unicycler_file) is True
            and os.stat(unicycler_file).st_size == 0
        ):
            unicycler_success = False
        self.unicycler_success = unicycler_success

    def get_depth(self, logdir, pacbio_model, threads):
        """wrapper function to get depth of each plasmid
        :param pacbio_model:  pacbio_model
        :param threads: threads
        :param logdir: logdir
        :return:
        """
        outdir = self.outdir
        concatenate_chrom_plasmids(outdir)

        input_long_reads: Path = Path(outdir) / "chopper_long_reads.fastq.gz"
        fasta: Path = Path(outdir) / "combined.fasta"
        sam_file: Path = Path(outdir) / "combined_long.sam"
        sorted_bam: Path = Path(outdir) / "combined_sorted_long.bam"

        # map
        minimap_long_reads(
            input_long_reads, fasta, sam_file, threads, pacbio_model, logdir
        )
        # sort
        sam_to_sorted_bam(sam_file, sorted_bam, threads, logdir)

        # short reads
        r1: Path = Path(outdir) / "trimmed_R1.fastq"
        r2: Path = Path(outdir) / "trimmed_R2.fastq"
        fasta: Path = Path(outdir) / "combined.fasta"
        sam_file: Path = Path(outdir) / "combined_short.sam"
        sorted_bam: Path = Path(outdir) / "combined_sorted_short.bam"

        # map
        minimap_short_reads(r1, r2, fasta, sam_file, threads, logdir)
        # sort
        sam_to_sorted_bam(sam_file, sorted_bam, threads, logdir)

        # get contig lengths

        fasta: Path = Path(outdir) / "combined.fasta"
        contig_lengths = get_contig_lengths(fasta)

        # depths
        short_bam_file: Path = Path(outdir) / "combined_sorted_short.bam"
        long_bam_file: Path = Path(outdir) / "combined_sorted_long.bam"
        depthsShort = get_depths_from_bam(short_bam_file, contig_lengths)
        depthsLong = get_depths_from_bam(long_bam_file, contig_lengths)

        # circular status
        circular_status = get_contig_circularity(fasta)

        summary_depth_df_short = collate_depths(depthsShort, "short", contig_lengths)
        summary_depth_df_long = collate_depths(depthsLong, "long", contig_lengths)

        # save the depth df in the class
        self.depth_df = combine_depth_dfs(
            summary_depth_df_short, summary_depth_df_long, circular_status
        )

    def get_depth_long(self, logdir, pacbio_model, threads, plas_fasta):
        """wrapper function to get depth of each plasmid
        :param pacbio_model:  pacbio_model
        :param threads: threads
        :param logdir: logdir
        :plas_fasta: plasmids from dnaapler
        :return:
        """
        outdir = self.outdir

        input_long_reads: Path = Path(outdir) / "chopper_long_reads.fastq.gz"
        chromosome: Path = Path(outdir) / "chromosome.fasta"
        combined_fasta: Path = Path(outdir) / "long_combined.fasta"
        sam_file: Path = Path(outdir) / "combined_long.sam"
        sorted_bam: Path = Path(outdir) / "combined_sorted_long.bam"

        # # write to combined fasta

        # Create a list to hold the combined sequences
        combined_sequences = []

        # Read and append sequences from the first FASTA file
        for record in SeqIO.parse(chromosome, "fasta"):
            combined_sequences.append(record)

        # Read and append sequences from the second FASTA file
        for record in SeqIO.parse(plas_fasta, "fasta"):
            combined_sequences.append(record)

        # Write the combined sequences to the output file
        SeqIO.write(combined_sequences, combined_fasta, "fasta")

        # map
        minimap_long_reads(
            input_long_reads, combined_fasta, sam_file, threads, pacbio_model, logdir
        )
        # sort
        sam_to_sorted_bam(sam_file, sorted_bam, threads, logdir)

        # get contig lengths
        contig_lengths = get_contig_lengths(combined_fasta)

        # depths
        long_bam_file: Path = Path(outdir) / "combined_sorted_long.bam"
        depthsLong = get_depths_from_bam(long_bam_file, contig_lengths)

        # circular status
        circular_status = get_contig_circularity(combined_fasta)
        summary_depth_df_long = collate_depths(depthsLong, "long", contig_lengths)

        # save the depth df in the class
        self.depth_df = depth_df_single(summary_depth_df_long, circular_status)

    def process_mash_tsv(self, plassembler_db_dir):
        """
        Process mash output
        :param outdir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        outdir = self.outdir

        # contig count of the unicycler assembly
        if self.long_only is False:
            contig_count = get_contig_count(
                os.path.join(outdir, "unicycler_output", "assembly.fasta")
            )
        else:  # for long only, just assembly.fasta
            contig_count = get_contig_count(os.path.join(outdir, "assembly.fasta"))
        # update with final plasmid count number
        self.contig_count = contig_count
        # get mash tsv output contig
        mash_tsv = os.path.join(outdir, "mash.tsv")
        col_list = [
            "contig",
            "ACC_NUCCORE",
            "mash_distance",
            "mash_pval",
            "mash_matching_hashes",
        ]

        # check if mash tsv file is empty -> no hits
        mash_empty = is_file_empty(mash_tsv)
        # instantiate tophits list
        tophits_mash_df = []

        if mash_empty is False:
            mash_df = pd.read_csv(
                mash_tsv, delimiter="\t", index_col=False, names=col_list
            )
            # get list of contigs from unicycler: 1 -> total number of contigs
            contigs = range(1, contig_count + 1)

            # instantiate tophits list
            tophits = []

            for contig in contigs:
                hit_df = (
                    mash_df.loc[mash_df["contig"] == contig]
                    .sort_values("mash_distance")
                    .reset_index(drop=True)
                )
                hits = len(hit_df["mash_distance"])
                # add only if there is a hit
                if hits > 0:
                    tmp_df = (
                        mash_df.loc[mash_df["contig"] == contig]
                        .sort_values("mash_distance")
                        .reset_index(drop=True)
                        .loc[0]
                    )
                    tophits.append(
                        [
                            tmp_df.contig,
                            "Yes",
                            tmp_df.ACC_NUCCORE,
                            tmp_df.mash_distance,
                            tmp_df.mash_pval,
                            tmp_df.mash_matching_hashes,
                        ]
                    )
                else:  # no hits append no it
                    tophits.append([contig, "", "", "", "", ""])
                # create tophits df
            tophits_mash_df = pd.DataFrame(
                tophits,
                columns=[
                    "contig",
                    "PLSDB_hit",
                    "ACC_NUCCORE",
                    "mash_distance",
                    "mash_pval",
                    "mash_matching_hashes",
                ],
            )

        else:  # empty mash file
            contigs = range(1, contig_count + 1)
            # create empty df
            tophits_mash_df = pd.DataFrame(
                columns=[
                    "contig",
                    "PLSDB_hit",
                    "ACC_NUCCORE",
                    "mash_distance",
                    "mash_pval",
                    "mash_matching_hashes",
                ]
            )
            for contig in contigs:
                tophits_mash_df.loc[contig - 1] = [contig, "", "", "", "", ""]

        # read in the plasdb tsv to get the description
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        cols = [
            "UID_NUCCORE",
            "ACC_NUCCORE",
            "Description_NUCCORE",
            "CreateDate_NUCCORE",
            "Topology_NUCCORE",
            "Completeness_NUCCORE",
            "TaxonID_NUCCORE",
            "Genome_NUCCORE",
            "Length_NUCCORE",
            "Source_NUCCORE",
            "UID_ASSEMBLY",
            "Status_ASSEMBLY",
            "SeqReleaseDate_ASSEMBLY",
            "SubmissionDate_ASSEMBLY",
            "Latest_ASSEMBLY",
            "UID_BIOSAMPLE",
            "ACC_BIOSAMPLE",
            "Location_BIOSAMPLE",
            "Coordinates_BIOSAMPLE",
            "IsolationSource_BIOSAMPLE",
            "Host_BIOSAMPLE",
            "CollectionDate_BIOSAMPLE",
            "Host_DISEASE",
            "SamplType_BIOSAMPLE",
            "taxon_name",
            "taxon_rank",
            "lineage",
            "taxon_species_id",
            "taxon_species_name",
            "taxon_genus_id",
            "taxon_genus_name",
            "taxon_family_id",
            "taxon_family_name",
            "taxon_order_id",
            "taxon_order_name",
            "taxon_class_id",
            "taxon_class_name",
            "taxon_phylum_id",
            "taxon_phylum_name",
            "taxon_superkingdom_id",
            "taxon_superkingdom_name",
            "loc_lat",
            "loc_lng",
            "loc_parsed",
            "GC_NUCCORE",
            "Identical",
            "OldVersion",
            "hits_rMLST",
            "hitscount_rMLST",
            "inclusions",
            "Host_BIOSAMPLE_processed",
            "Host_DISEASE_processed",
            "D1",
            "D2",
            "plasmidfinder",
            "pmlst",
            "relaxase_type(s)",
            "mpf_type",
        ]
        plsdb_tsv = pd.read_csv(
            plsdb_tsv_file,
            delimiter="\t",
            index_col=False,
            names=cols,
            skiprows=1,
            low_memory=False,
        )
        combined_mash_df = tophits_mash_df.merge(
            plsdb_tsv, on="ACC_NUCCORE", how="left"
        )

        self.mash_df = combined_mash_df

    def combine_depth_mash_tsvs(self, prefix):
        """
        Combine depth and mash dataframes
        :param outdir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        outdir = self.outdir
        self.depth_df["contig"] = self.depth_df["contig"].astype("str")
        self.mash_df["contig"] = self.mash_df["contig"].astype("str")
        combined_depth_mash_df = self.depth_df.merge(
            self.mash_df, on="contig", how="left"
        )
        # no hit for chromosome
        combined_depth_mash_df.loc[
            combined_depth_mash_df["contig"].str.contains("chromosome"), "PLSDB_hit"
        ] = ""
        combined_depth_mash_df.to_csv(
            os.path.join(outdir, prefix + "_summary.tsv"), sep="\t", index=False
        )
        self.combined_depth_mash_df = combined_depth_mash_df

    def finalise_contigs(self, prefix):
        """
        Renames the contigs of unicycler with the new plasmid copy numbers and outputs finalised file
        """
        outdir = self.outdir

        combined_depth_mash_df = self.combined_depth_mash_df
        combined_depth_mash_df = combined_depth_mash_df.loc[
            combined_depth_mash_df["contig"] != "chromosome"
        ].reset_index(drop=True)
        # get contigs only
        plasmid_fasta = os.path.join(outdir, "unicycler_output", "assembly.fasta")
        i = 0
        with open(os.path.join(outdir, prefix + "_plasmids.fasta"), "w") as dna_fa:
            for dna_record in SeqIO.parse(plasmid_fasta, "fasta"):
                if "circular" in dna_record.description:  # circular contigs
                    id_updated = (
                        dna_record.description.split(" ")[0]
                        + " "
                        + dna_record.description.split(" ")[1]
                        + " plasmid_copy_number_short="
                        + str(combined_depth_mash_df.plasmid_copy_number_short[i])
                        + "x plasmid_copy_number_long="
                        + str(combined_depth_mash_df.plasmid_copy_number_long[i])
                        + "x "
                        + "circular=true"
                    )
                else:  # non circular contigs
                    id_updated = (
                        dna_record.description.split(" ")[0]
                        + " "
                        + dna_record.description.split(" ")[1]
                        + " plasmid_copy_number_short="
                        + str(combined_depth_mash_df.plasmid_copy_number_short[i])
                        + "x plasmid_copy_number_long="
                        + str(combined_depth_mash_df.plasmid_copy_number_long[i])
                        + "x "
                    )
                i += 1
                record = SeqRecord(dna_record.seq, id=id_updated, description="")
                SeqIO.write(record, dna_fa, "fasta")

    def finalise_contigs_long(self, prefix):
        """
        Renames the contigs of assembly with new ones
        """
        outdir = self.outdir

        combined_depth_mash_df = self.combined_depth_mash_df
        combined_depth_mash_df = combined_depth_mash_df.loc[
            combined_depth_mash_df["contig"] != "chromosome"
        ].reset_index(drop=True)
        # get contigs only
        plasmid_fasta = os.path.join(outdir, "plasmids.fasta")
        i = 0
        with open(os.path.join(outdir, prefix + "_plasmids.fasta"), "w") as dna_fa:
            for dna_record in SeqIO.parse(plasmid_fasta, "fasta"):
                id = dna_record.id
                length = len(dna_record.seq)
                copy_number = combined_depth_mash_df.plasmid_copy_number_long[i]
                if "circular" in dna_record.description:  # circular contigs from canu
                    desc = f"len={length} plasmid_copy_number_long={copy_number}x circular=True"
                else:
                    desc = f"len={length} plasmid_copy_number_long={copy_number}x"
                i += 1
                record = SeqRecord(dna_record.seq, id=id, description=desc)
                SeqIO.write(record, dna_fa, "fasta")


class Assembly:
    """Plassembler Assembly Mode Output Class"""

    def __init__(
        self,
        outdir: str = "output/",
        long_flag: bool = True,
        short_flag: bool = True,
        chromosome_name: str = "chromosome",
        contig_count: int = 1,
        plasmid_names: list = ["1"],
        threads: int = 1,
        depth_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        mash_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        combined_depth_mash_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
    ) -> None:
        """
        Parameters
        --------
        outdir: int, required
            output directory
        contig_count: int, required
            the number of contigs assembled by flye assembly
        depth_df: Pandas df, required
            dataframe with depth output
        mash_df: Pandas df, required
            dataframe with mash output
        threads: int, required
            integer giving args.threads - defaults to 1.
        long_flag: bool, required
            whether there's long read FASTQs
        short_flag: bool, required
            whether there's short read FASTQs
        """
        self.outdir = outdir
        self.contig_count = contig_count
        self.threads = threads
        self.depth_df = depth_df
        self.mash_df = mash_df
        self.combined_depth_mash_df = combined_depth_mash_df
        self.long_flag = long_flag
        self.short_flag = short_flag

    def combine_input_fastas(self, chromosome_fasta: Path, plasmids_fasta: Path):
        """wrapper function to get depth of each plasmid
        :param prefix: prefix (default plassembler)
        :param outdir:  Output Directory
        :param threads: threads
        :param logger: logger
        :return:
        """
        # combined input fasta
        combined_fasta = Path(self.outdir) / "combined.fasta"
        chromosome_name = ""

        # rename the first contig as chromosome
        with open(chromosome_fasta, "r") as f_in, open(combined_fasta, "w") as f_out:
            # Parse the input FASTA file
            records = list(SeqIO.parse(f_in, "fasta"))
            # keep chromosome name
            chromosome_name = records[0].id
            # Rename the first record
            records[0].id = "chromosome"
            records[0].description = ""
            # Write the modified records to the output FASTA file
            SeqIO.write(records, f_out, "fasta")

        plasmid_names = []

        with open(plasmids_fasta, "r") as f_in, open(combined_fasta, "a") as f_out:
            records = list(SeqIO.parse(f_in, "fasta"))
            i = 0
            for record in records:
                i += 1
                # keep plasmid name
                plasmid_names.append(record.id)
                record.id = str(i)
                records[0].description = ""
                # Write the  records to the output FASTA file
                SeqIO.write(record, f_out, "fasta")

        self.chromosome_name = chromosome_name
        self.plasmid_names = plasmid_names

    def get_depth(self, logdir, threads, pacbio_model):
        """wrapper function to get depth of each plasmid
        :param threads: threads
        :param logdir: logdir
        :return:
        """
        outdir = self.outdir

        input_long_reads: Path = Path(outdir) / "chopper_long_reads.fastq.gz"
        fasta: Path = Path(outdir) / "combined.fasta"
        sam_file: Path = Path(outdir) / "combined_long.sam"
        sorted_bam: Path = Path(outdir) / "combined_sorted_long.bam"

        # map
        if self.long_flag is True:
            minimap_long_reads(
                input_long_reads, fasta, sam_file, threads, pacbio_model, logdir
            )
            # sort
            sam_to_sorted_bam(sam_file, sorted_bam, threads, logdir)

        # short reads
        r1: Path = Path(outdir) / "trimmed_R1.fastq"
        r2: Path = Path(outdir) / "trimmed_R2.fastq"
        fasta: Path = Path(outdir) / "combined.fasta"
        sam_file: Path = Path(outdir) / "combined_short.sam"
        sorted_bam: Path = Path(outdir) / "combined_sorted_short.bam"

        # map
        if self.short_flag is True:
            minimap_short_reads(r1, r2, fasta, sam_file, threads, logdir)
            # sort
            sam_to_sorted_bam(sam_file, sorted_bam, threads, logdir)

        # get contig lengths

        fasta: Path = Path(outdir) / "combined.fasta"
        contig_lengths = get_contig_lengths(fasta)

        # depths
        short_bam_file: Path = Path(outdir) / "combined_sorted_short.bam"
        long_bam_file: Path = Path(outdir) / "combined_sorted_long.bam"
        if self.short_flag is True:
            depthsShort = get_depths_from_bam(short_bam_file, contig_lengths)
        if self.long_flag is True:
            depthsLong = get_depths_from_bam(long_bam_file, contig_lengths)

        # circular status
        circular_status = get_contig_circularity(fasta)

        if self.short_flag is True:
            summary_depth_df_short = collate_depths(
                depthsShort, "short", contig_lengths
            )

        if self.long_flag is True:
            summary_depth_df_long = collate_depths(depthsLong, "long", contig_lengths)

        # save the depth df in the class
        if self.long_flag is True and self.short_flag is True:
            self.depth_df = combine_depth_dfs(
                summary_depth_df_short, summary_depth_df_long, circular_status
            )
        elif self.long_flag is True and self.short_flag is False:  # only long
            self.depth_df = depth_df_single(summary_depth_df_long, circular_status)
        elif self.long_flag is False and self.short_flag is True:  # only short
            self.depth_df = depth_df_single(summary_depth_df_short, circular_status)

    def process_mash_tsv(self, plassembler_db_dir, plasmid_fasta):
        """
        Process mash output
        :param outdir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        outdir = self.outdir
        # contig count of the unicycler assembly
        contig_count = get_contig_count(plasmid_fasta)
        # update with final plasmid count number
        self.contig_count = contig_count

        # get mash tsv output contig
        mash_tsv = os.path.join(outdir, "mash.tsv")
        col_list = [
            "contig",
            "ACC_NUCCORE",
            "mash_distance",
            "mash_pval",
            "mash_matching_hashes",
        ]

        # check if mash tsv file is empty -> no hits
        mash_empty = is_file_empty(mash_tsv)
        # instantiate tophits list
        tophits_mash_df = []

        if mash_empty is False:
            mash_df = pd.read_csv(
                mash_tsv, delimiter="\t", index_col=False, names=col_list
            )
            # get list of contigs from unicycler: 1 -> total number of contigs
            contigs = range(1, contig_count + 1)

            # instantiate tophits list
            tophits = []

            for contig in contigs:
                hit_df = (
                    mash_df.loc[mash_df["contig"] == contig]
                    .sort_values("mash_distance")
                    .reset_index(drop=True)
                )
                hits = len(hit_df["mash_distance"])
                # add only if there is a hit
                if hits > 0:
                    tmp_df = (
                        mash_df.loc[mash_df["contig"] == contig]
                        .sort_values("mash_distance")
                        .reset_index(drop=True)
                        .loc[0]
                    )
                    tophits.append(
                        [
                            tmp_df.contig,
                            "Yes",
                            tmp_df.ACC_NUCCORE,
                            tmp_df.mash_distance,
                            tmp_df.mash_pval,
                            tmp_df.mash_matching_hashes,
                        ]
                    )
                else:  # no hits append no it
                    tophits.append([contig, "", "", "", "", ""])
                # create tophits df
            tophits_mash_df = pd.DataFrame(
                tophits,
                columns=[
                    "contig",
                    "PLSDB_hit",
                    "ACC_NUCCORE",
                    "mash_distance",
                    "mash_pval",
                    "mash_matching_hashes",
                ],
            )

        else:  # empty mash file
            contigs = range(1, contig_count + 1)
            # create empty df
            tophits_mash_df = pd.DataFrame(
                columns=[
                    "contig",
                    "PLSDB_hit",
                    "ACC_NUCCORE",
                    "mash_distance",
                    "mash_pval",
                    "mash_matching_hashes",
                ]
            )
            for contig in contigs:
                tophits_mash_df.loc[contig - 1] = [contig, "", "", "", "", ""]

        # read in the plasdb tsv to get the description
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        cols = [
            "UID_NUCCORE",
            "ACC_NUCCORE",
            "Description_NUCCORE",
            "CreateDate_NUCCORE",
            "Topology_NUCCORE",
            "Completeness_NUCCORE",
            "TaxonID_NUCCORE",
            "Genome_NUCCORE",
            "Length_NUCCORE",
            "Source_NUCCORE",
            "UID_ASSEMBLY",
            "Status_ASSEMBLY",
            "SeqReleaseDate_ASSEMBLY",
            "SubmissionDate_ASSEMBLY",
            "Latest_ASSEMBLY",
            "UID_BIOSAMPLE",
            "ACC_BIOSAMPLE",
            "Location_BIOSAMPLE",
            "Coordinates_BIOSAMPLE",
            "IsolationSource_BIOSAMPLE",
            "Host_BIOSAMPLE",
            "CollectionDate_BIOSAMPLE",
            "Host_DISEASE",
            "SamplType_BIOSAMPLE",
            "taxon_name",
            "taxon_rank",
            "lineage",
            "taxon_species_id",
            "taxon_species_name",
            "taxon_genus_id",
            "taxon_genus_name",
            "taxon_family_id",
            "taxon_family_name",
            "taxon_order_id",
            "taxon_order_name",
            "taxon_class_id",
            "taxon_class_name",
            "taxon_phylum_id",
            "taxon_phylum_name",
            "taxon_superkingdom_id",
            "taxon_superkingdom_name",
            "loc_lat",
            "loc_lng",
            "loc_parsed",
            "GC_NUCCORE",
            "Identical",
            "OldVersion",
            "hits_rMLST",
            "hitscount_rMLST",
            "inclusions",
            "Host_BIOSAMPLE_processed",
            "Host_DISEASE_processed",
            "D1",
            "D2",
            "plasmidfinder",
            "pmlst",
            "relaxase_type(s)",
            "mpf_type",
        ]
        plsdb_tsv = pd.read_csv(
            plsdb_tsv_file,
            delimiter="\t",
            index_col=False,
            names=cols,
            skiprows=1,
            low_memory=False,
        )
        combined_mash_df = tophits_mash_df.merge(
            plsdb_tsv, on="ACC_NUCCORE", how="left"
        )

        self.mash_df = combined_mash_df

    def combine_depth_mash_tsvs(self, prefix, no_copy_numbers):
        """
        Combine depth and mash dataframes - if not --no_copy_numbers
        :param outdir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        outdir = self.outdir
        self.mash_df["contig"] = self.mash_df["contig"].astype("str")

        if no_copy_numbers is False:
            self.depth_df["contig"] = self.depth_df["contig"].astype("str")
            combined_depth_mash_df = self.depth_df.merge(
                self.mash_df, on="contig", how="left"
            )
            # no hit for chromosome
            combined_depth_mash_df.loc[
                combined_depth_mash_df["contig"].str.contains("chromosome"), "PLSDB_hit"
            ] = ""
        else:
            combined_depth_mash_df = self.mash_df

        combined_depth_mash_df.to_csv(
            os.path.join(outdir, prefix + "_summary.tsv"), sep="\t", index=False
        )
        self.combined_depth_mash_df = combined_depth_mash_df
