import os 
import pandas as pd
import depth
import mapping
import cleanup
import run_mash
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Plass:
    """Plassembler Output Class"""

    def __init__(
        self,
        out_dir: str = 'output/',
        contig_count: int = 1,
        no_plasmids_flag: bool = False,
        chromosome_flag: bool = True,
        threads: int = 1,
        depth_df:  pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]}),
        mash_df:  pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]}),
        combined_depth_mash_df: pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]}),
        kmer_mode: bool = False,
        unicycler_success: bool = True
        ) -> None:
        """
        Parameters
        --------
        out_dir: int, required
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
        kmer_mode: bool, required
            whether plassembler is in kmer mode
        mash_empty: bool, required
            whether mash output is empty (True) or not
        unicycler_success: bool, required
            whether unicycler succeeded
        """
        self.out_dir = out_dir
        self.contig_count = contig_count
        self.no_plasmids_flag = no_plasmids_flag
        self.chromosome_flag = chromosome_flag
        self.threads = threads
        self.depth_df = depth_df
        self.mash_df = mash_df
        self.combined_depth_mash_df = combined_depth_mash_df
        self.kmer_mode = kmer_mode
        self.unicycler_success = unicycler_success

    def get_contig_count(self,  logger):
        """ Counts the number of contigs assembled by flye and prints to log
        :param out_dir: output directory
        :return:
        """
        out_dir = self.out_dir
        info_file =  os.path.join(out_dir, "assembly_info.txt")
        col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
        info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
        contig_count = len(info_df['seq_name'])
        message = "Flye assembled " + str(contig_count) + " contigs."
        print(message)
        logger.info(message)
        self.contig_count = contig_count

    def make_chrom_bed(self):
        """ Creates simple bed file with chromosome, 1 and length
        :param out_dir: output directory
        :return:
        """
        out_dir = self.out_dir
        info_file =  os.path.join(out_dir, "assembly_info.txt")
        col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
        info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
        
        contig_count = len(info_df['seq_name'])


    def identify_chromosome_process_flye(self,  chromosome_len, keep_chromosome):
        """Identified chromosome and processes Flye output - renames chromosome contig and the others as plasmid_1, plasmid_2 etc
        Also makes the chromosome bed file for downstream samtools mapping
        :param out_dir: output directory
        :param chromosome_len: lower bound on length of chromosome from input command
        :param keep_chromosome: whether the user wants to keep the chromosome as chromosome.fasta
        :return chromosome_flag: bool whether chromosome assembles
        """
        out_dir = self.out_dir
        info_file =  os.path.join(out_dir, "assembly_info.txt")
        col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
        info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
        max_length = max(info_df['length'])
        # get putative chromosome contig
        chrom_contig = info_df[info_df['length'] == max_length].iloc[0]['seq_name']
        # comment out circ here
        # chrom_circ = info_df[info_df['length'] == max_length].iloc[0]['circ']

        # of chromosome is at least 90% of inputted chromosome length, flag 
        # to say that the chromosome has been correctly identified 


        chromosome_flag = True
        if max_length < int(chromosome_len)*0.9:  # no chromosome identified -> don't bother with the renaming
            chromosome_flag = False
        ### assuming chromosome identified
        else:
            # make bed file with plasmid contigs to extract mapping reads
            bed_file =  open(os.path.join(out_dir, "non_chromosome.bed"), 'w')
            bed_chrom_file =  open(os.path.join(out_dir, "chromosome.bed"), 'w')
            with open(os.path.join(out_dir, "flye_renamed.fasta"), 'w') as rename_fa:
                for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
                    i = 1
                    # chromosome
                    if dna_record.id == chrom_contig:
                        dna_header = "chromosome"
                        dna_description = ""
                        dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                        SeqIO.write(dna_record, rename_fa, 'fasta')
                        bed_chrom_file.write(f'chromosome\t1\t{max_length}\n')  # chromosome
                    # plasmids
                    else:
                        dna_header = "plasmid_" + str(i)
                        dna_description = ""
                        # get length for bed file
                        l = info_df.length.loc[info_df['seq_name'] == dna_record.id]
                        plas_len =  int(l.iloc[0])
                        # write the updated record
                        dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                        SeqIO.write(dna_record, rename_fa, 'fasta')
                        # get length for bed file
                        # make bed file
                        bed_file.write(f'{dna_header}\t1\t{plas_len}\n')  # Write read name
                        i += 1
            # keep the chromosome if specified
            if keep_chromosome == True:
                with open(os.path.join(out_dir, "chromosome.fasta"), 'w') as chrom_fa:
                    for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
                        # chromosome
                        if dna_record.id == chrom_contig:
                            dna_header = "chromosome"
                            dna_description = ""
                            dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                            SeqIO.write(dna_record, chrom_fa, 'fasta')

        # add to object
        self.chromosome_flag = chromosome_flag


    def check_unicycler_success(self):
        out_dir = self.out_dir
        unicycler_file = os.path.join(out_dir, "unicycler_output", "assembly.fasta")
        # check if unicycler succeded according to the output (it won't if no plasmids)
        unicycler_success = os.path.isfile(unicycler_file)
        if os.path.isfile(unicycler_file) == True and os.stat(unicycler_file).st_size == 0:
            unicycler_success = False
        self.unicycler_success = unicycler_success


    def get_depth(self, logger, threads, prefix):
        """ wrapper function to get depth of each plasmid in kmer mode
        :param prefix: prefix (default plassembler)
        :param out_dir:  Output Directory
        :param threads: threads
        :param logger: logger
        :return: 
        """
        out_dir = self.out_dir
        depth.concatenate_chrom_plasmids(out_dir, logger)
        depth.minimap_depth_sort_long(out_dir, threads)
        if self.kmer_mode == False:
            depth.minimap_depth_sort_short(out_dir, threads)

        contig_lengths = depth.get_contig_lengths(out_dir)
        if self.kmer_mode == False:
            depthsShort = depth.get_depths_from_bam(out_dir, "short", contig_lengths)
        depthsLong = depth.get_depths_from_bam(out_dir, "long", contig_lengths)
        circular_status = depth.get_contig_circularity(out_dir)
        if self.kmer_mode == False:
            summary_depth_df_short = depth.collate_depths(depthsShort,"short",contig_lengths)
        summary_depth_df_long = depth.collate_depths(depthsLong,"long",contig_lengths)
        # save the depth df in the class
        if self.kmer_mode == False:
            self.depth_df = depth.combine_depth_dfs(out_dir, summary_depth_df_short, summary_depth_df_long, prefix, circular_status)
        else:
            self.depth_df = depth.kmer_final_output(out_dir, summary_depth_df_long, prefix, circular_status)



    def process_mash_tsv(self,  plassembler_db_dir):
        """
        Process mash output
        :param out_dir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        out_dir = self.out_dir
        # contig count of the unicycler assembly
        contig_count = run_mash.get_contig_count( os.path.join(out_dir,"unicycler_output", "assembly.fasta"))
        # update with final plasmid count number 
        self.contig_count = contig_count
        # get mash tsv output contig 
        mash_tsv = os.path.join(out_dir,"mash.tsv")
        col_list = ["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"]

        # check if mash tsv file is empty -> no hits
        mash_empty = run_mash.is_file_empty(mash_tsv)
        # instantiate tophits list
        tophits_mash_df = []

        if mash_empty == False:   
            mash_df = pd.read_csv(mash_tsv, delimiter= '\t', index_col=False, names=col_list ) 
            # get list of contigs from unicycler: 1 -> total number of contigs
            contigs = range(1, contig_count+1)

            # instantiate tophits list
            tophits = []

            for contig in contigs:
                hit_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True)
                hits = len(hit_df['mash_distance'])
                # add only if there is a hit
                if hits > 0:
                    tmp_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True).loc[0]
                    tophits.append([tmp_df.contig, 'Yes', tmp_df.ACC_NUCCORE, tmp_df.mash_distance, tmp_df.mash_pval, tmp_df.mash_matching_hashes])
                else: # no hits append no it
                    tophits.append([contig, 'No_hit', 'No_hit', 'No_hit', 'No_hit', 'No_hit'])
                # create tophits df
            tophits_mash_df = pd.DataFrame(tophits, columns=["contig", "PLSDB_hit", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])

        else: # empty mash file
            contigs = range(1, contig_count+1)
            # create empty df
            tophits_mash_df = pd.DataFrame(columns=["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])
            for contig in contigs:
                tophits_mash_df.loc[contig-1] = [contig, 'No_hit', 'No_hit', 'No_hit', 'No_hit']

        # read in the plasdb tsv to get the description
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        cols = ["UID_NUCCORE","ACC_NUCCORE","Description_NUCCORE","CreateDate_NUCCORE","Topology_NUCCORE","Completeness_NUCCORE","TaxonID_NUCCORE","Genome_NUCCORE","Length_NUCCORE","Source_NUCCORE","UID_ASSEMBLY","Status_ASSEMBLY","SeqReleaseDate_ASSEMBLY","SubmissionDate_ASSEMBLY","Latest_ASSEMBLY","UID_BIOSAMPLE","ACC_BIOSAMPLE","Location_BIOSAMPLE","Coordinates_BIOSAMPLE","IsolationSource_BIOSAMPLE","Host_BIOSAMPLE","CollectionDate_BIOSAMPLE","Host_DISEASE","SamplType_BIOSAMPLE","taxon_name","taxon_rank","lineage","taxon_species_id","taxon_species_name","taxon_genus_id","taxon_genus_name","taxon_family_id","taxon_family_name","taxon_order_id","taxon_order_name","taxon_class_id","taxon_class_name","taxon_phylum_id","taxon_phylum_name","taxon_superkingdom_id","taxon_superkingdom_name","loc_lat","loc_lng","loc_parsed","GC_NUCCORE","Identical","OldVersion","hits_rMLST","hitscount_rMLST","inclusions","Host_BIOSAMPLE_processed","Host_DISEASE_processed","D1","D2","plasmidfinder","pmlst","relaxase_type(s)","mpf_type"]
        plsdb_tsv = pd.read_csv(plsdb_tsv_file, delimiter= '\t', index_col=False, names=cols, skiprows=1,low_memory=False) 
        combined_mash_df = tophits_mash_df.merge(plsdb_tsv, on='ACC_NUCCORE', how='left')

        self.mash_df = combined_mash_df


    def combine_depth_mash_tsvs(self,  prefix):
        """
        Combine depth and mash dataframes
        :param out_dir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        out_dir = self.out_dir
        self.depth_df['contig'] = self.depth_df['contig'].astype('str')
        self.mash_df['contig'] = self.mash_df['contig'].astype('str')
        combined_depth_mash_df = self.depth_df.merge(self.mash_df, on='contig', how='left')
        combined_depth_mash_df.to_csv(os.path.join(out_dir, prefix + "_summary.tsv"), sep="\t", index=False) 
        self.combined_depth_mash_df = combined_depth_mash_df


    def finalise_contigs(self, out_dir, prefix ):
        """
        Renames the contigs of unicycler with the new plasmid copy numbers and outputs finalised file
        :param out_dir: output directory
        :return: 
        """

        combined_depth_mash_df = self.combined_depth_mash_df
        combined_depth_mash_df = combined_depth_mash_df.loc[combined_depth_mash_df['contig'] != 'chromosome'].reset_index(drop=True)
        # get contigs only
        plasmid_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
        i = 0
        with open(os.path.join(out_dir, prefix + "_plasmids.fasta"), 'w') as dna_fa:
            for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
                if "circular" in dna_record.description: # circular contigs
                    if self.kmer_mode == False: 
                        id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_short=" + str(combined_depth_mash_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(combined_depth_mash_df.plasmid_copy_number_long[i]) + "x " + dna_record.description.split(' ')[3]
                    else: # kmer mode
                        id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_long=" + str(combined_depth_mash_df.plasmid_copy_number_long[i]) + "x " + dna_record.description.split(' ')[3]
                else: # non circular contigs
                    if self.kmer_mode == False:
                        id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_short=" + str(combined_depth_mash_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(combined_depth_mash_df.plasmid_copy_number_long[i]) + "x " 
                    else: #kmer mode
                        id_updated = dna_record.description.split(' ')[0] + " " + dna_record.description.split(' ')[1] + " plasmid_copy_number_long=" + str(combined_depth_mash_df.plasmid_copy_number_long[i]) + "x "
                i += 1
                record = SeqRecord(dna_record.seq, id=id_updated, description = "" )
                SeqIO.write(record, dna_fa, 'fasta')

class Assembly:
    """Plassembler Assembly Mode Output Class"""

    def __init__(
        self,
        out_dir: str = 'output/',
        long_flag: bool = True,
        short_flag: bool = True,
        chromosome_name: str = "chromosome",
        contig_count: int = 1,
        plasmid_names: list = ['1'],
        threads: int = 1,
        depth_df:  pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]}),
        mash_df:  pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]}),
        combined_depth_mash_df: pd.DataFrame() =  pd.DataFrame({'col1': [1, 2, 3], 'col2': [4, 5, 6]})

        ) -> None:
        """
        Parameters
        --------
        out_dir: int, required
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
        kmer_mode: bool, required
            whether plassembler is in kmer mode
        long_flag: bool, required
            whether there's long read FASTQs
        short_flag: bool, required
            whether there's short read FASTQs
        """
        self.out_dir = out_dir
        self.contig_count = contig_count
        self.no_plasmids_flag = no_plasmids_flag
        self.chromosome_flag = chromosome_flag
        self.threads = threads
        self.depth_df = depth_df
        self.mash_df = mash_df
        self.combined_depth_mash_df = combined_depth_mash_df
        self.kmer_mode = kmer_mode
        self.long_flag = long_flag
        self.short_flag = short_flag


    def combine_input_fastas(self, chromosome_fasta, plasmids_fasta):
        """ wrapper function to get depth of each plasmid 
        :param prefix: prefix (default plassembler)
        :param out_dir:  Output Directory
        :param threads: threads
        :param logger: logger
        :return: 
        """
        # combined input fasta
        combined_fasta = os.path.join(self.out_dir, "combined.fasta")

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
                i +=1
                # keep plasmid name
                plasmid_names.append(record.id)
                record.id = i
                records[0].description = ""
            # Write the  records to the output FASTA file
                SeqIO.write(record, f_out, "fasta")

        self.chromosome_name = chromosome_name
        self.plasmid_names = plasmid_names



    def get_depth(self, threads, prefix):

        out_dir = self.out_dir

        if self.long_flag == True:
            depth.minimap_depth_sort_long(out_dir, threads)
        if self.short_flag == True:
            depth.minimap_depth_sort_short(out_dir, threads)

        contig_lengths = depth.get_contig_lengths(out_dir)

        if self.long_flag == True:
            depthsLong = depth.get_depths_from_bam(out_dir, "long", contig_lengths)
        if self.short_flag == True:
            depthsShort = depth.get_depths_from_bam(out_dir, "short", contig_lengths)

        # circular status
        circular_status = depth.get_contig_circularity(out_dir)

        if self.long_flag == True:
            summary_depth_df_long = depth.collate_depths(depthsLong,"long",contig_lengths)
        if self.short_flag == True:
            summary_depth_df_short = depth.collate_depths(depthsShort,"short",contig_lengths)

        # save the depth df in the class
        if self.long_flag == True and self.short_flag == True:
            self.depth_df = depth.combine_depth_dfs(out_dir, summary_depth_df_short, summary_depth_df_long, prefix, circular_status)
        elif self.long_flag == True and self.short_flag == False:
            self.depth_df = depth.kmer_final_output(out_dir, summary_depth_df_long, prefix, circular_status)
        elif self.long_flag == False and self.short_flag == True:
            self.depth_df = depth.kmer_final_output(out_dir, summary_depth_df_short, prefix, circular_status)


    
    def process_mash_tsv(self,  plassembler_db_dir, plasmid_fasta):
        """
        Process mash output
        :param out_dir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        out_dir = self.out_dir
        # contig count of the unicycler assembly
        contig_count = run_mash.get_contig_count( plasmid_fasta)
        # update with final plasmid count number 
        self.contig_count = contig_count

        # get mash tsv output contig 
        mash_tsv = os.path.join(out_dir,"mash.tsv")
        col_list = ["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"]

        # check if mash tsv file is empty -> no hits
        mash_empty = run_mash.is_file_empty(mash_tsv)
        # instantiate tophits list
        tophits_mash_df = []

        if mash_empty == False:   
            mash_df = pd.read_csv(mash_tsv, delimiter= '\t', index_col=False, names=col_list ) 
            # get list of contigs from unicycler: 1 -> total number of contigs
            contigs = range(1, contig_count+1)

            # instantiate tophits list
            tophits = []

            for contig in contigs:
                hit_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True)
                hits = len(hit_df['mash_distance'])
                # add only if there is a hit
                if hits > 0:
                    tmp_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True).loc[0]
                    tophits.append([tmp_df.contig, 'Yes', tmp_df.ACC_NUCCORE, tmp_df.mash_distance, tmp_df.mash_pval, tmp_df.mash_matching_hashes])
                else: # no hits append no it
                    tophits.append([contig, 'No_hit', 'No_hit', 'No_hit', 'No_hit', 'No_hit'])
                # create tophits df
            tophits_mash_df = pd.DataFrame(tophits, columns=["contig", "PLSDB_hit", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])

        else: # empty mash file
            contigs = range(1, contig_count+1)
            # create empty df
            tophits_mash_df = pd.DataFrame(columns=["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])
            for contig in contigs:
                tophits_mash_df.loc[contig-1] = [contig, 'No_hit', 'No_hit', 'No_hit', 'No_hit']

        # read in the plasdb tsv to get the description
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        cols = ["UID_NUCCORE","ACC_NUCCORE","Description_NUCCORE","CreateDate_NUCCORE","Topology_NUCCORE","Completeness_NUCCORE","TaxonID_NUCCORE","Genome_NUCCORE","Length_NUCCORE","Source_NUCCORE","UID_ASSEMBLY","Status_ASSEMBLY","SeqReleaseDate_ASSEMBLY","SubmissionDate_ASSEMBLY","Latest_ASSEMBLY","UID_BIOSAMPLE","ACC_BIOSAMPLE","Location_BIOSAMPLE","Coordinates_BIOSAMPLE","IsolationSource_BIOSAMPLE","Host_BIOSAMPLE","CollectionDate_BIOSAMPLE","Host_DISEASE","SamplType_BIOSAMPLE","taxon_name","taxon_rank","lineage","taxon_species_id","taxon_species_name","taxon_genus_id","taxon_genus_name","taxon_family_id","taxon_family_name","taxon_order_id","taxon_order_name","taxon_class_id","taxon_class_name","taxon_phylum_id","taxon_phylum_name","taxon_superkingdom_id","taxon_superkingdom_name","loc_lat","loc_lng","loc_parsed","GC_NUCCORE","Identical","OldVersion","hits_rMLST","hitscount_rMLST","inclusions","Host_BIOSAMPLE_processed","Host_DISEASE_processed","D1","D2","plasmidfinder","pmlst","relaxase_type(s)","mpf_type"]
        plsdb_tsv = pd.read_csv(plsdb_tsv_file, delimiter= '\t', index_col=False, names=cols, skiprows=1,low_memory=False) 
        combined_mash_df = tophits_mash_df.merge(plsdb_tsv, on='ACC_NUCCORE', how='left')

        self.mash_df = combined_mash_df


    def combine_depth_mash_tsvs(self,  prefix):
        """
        Combine depth and mash dataframes
        :param out_dir: output directory
        :return: mash_empty: boolean whether there was a mash hit
        """
        out_dir = self.out_dir
        self.depth_df['contig'] = self.depth_df['contig'].astype('str')
        self.mash_df['contig'] = self.mash_df['contig'].astype('str')
        combined_depth_mash_df = self.depth_df.merge(self.mash_df, on='contig', how='left')
        combined_depth_mash_df.to_csv(os.path.join(out_dir, prefix + "_summary.tsv"), sep="\t", index=False) 
        self.combined_depth_mash_df = combined_depth_mash_df


