import os 
import pandas as pd
import depth
import mapping
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Plass:
    """Plassembler Output Class"""

    def __init__(
        self,
        contig_count: int = 1,
        successful_unicycler_recovery: bool = True,
        no_plasmids_flag: bool = False,
        chromosome_flag: bool = True

        ) -> None:
        """
        Parameters
        --------
        contig_count: int, required
            the number of contigs assembled by flye assembly
        successful_unicycler_recovery: bool, required
            flag determining whether unicycler finished. If False, suggests no plasmid from unicycler output
        no_plasmids_flag: bool, required
            flag determining whether there are plasmids at some point in the pipeline. If False means there are plasmids
        chromosome_flag: bool, required
            flag saying whether or not there was an identified chromosome in the Flye output. 
        case: int, required
            integer giving what case the plassembler scenario is - needs to be 1,2,3 or 4. 
        """
        self.contig_count = contig_count
        self.successful_unicycler_recovery = successful_unicycler_recovery
        self.no_plasmids_flag = no_plasmids_flag
        self.chromosome_flag = chromosome_flag

    def get_contig_count(self, out_dir, logger):
        """ Counts the number of contigs assembled by flye and prints to log
        :param out_dir: output directory
        :return:
        """
        info_file =  os.path.join(out_dir, "assembly_info.txt")
        col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
        info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
        contig_count = len(info_df['seq_name'])
        message = "Flye assembled " + str(contig_count) + " contigs."
        print(message)
        logger.info(message)
        self.contig_count = contig_count
   
    def identify_chromosome_process_flye(self, out_dir, chromosome_len):
        """Identified chromosome and processes Flye output - renames chromosome contig and the others as plasmid_1, plasmid_2 etc

        :param out_dir: output directory
        :param chromosome_len: lower bound on length of chromosome from input command
        :return: chromosome_flag = a flag saying whether or not there was an identified chromosome
        """
        
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
        else:
            with open(os.path.join(out_dir, "flye_renamed.fasta"), 'w') as rename_fa:
                for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
                    i = 1
                    # chromosome
                    if dna_record.id == chrom_contig:
                        dna_header = "chromosome"
                        dna_description = ""
                        dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                        SeqIO.write(dna_record, rename_fa, 'fasta')
                    # plasmids
                    else:
                        dna_header = "plasmid_" + str(i)
                        dna_description = ""
                        dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                        SeqIO.write(dna_record, rename_fa, 'fasta')
                        i += 1

        # add to object
        self.chromosome_flag = chromosome_flag

    def get_depth(out_dir, logger,  threads, prefix):
        """ wrapper function to get depth of each plasmid in kmer mode
        :param prefix: prefix (default plassembler)
        :param out_dir:  Output Directory
        :param threads: threads
        :param logger: logger
        :return: 
        """
        depth.concatenate_chrom_plasmids(out_dir, logger)
        mapping.index_fasta(os.path.join(out_dir, "combined.fasta"),  logger)
        depth.bwa_map_depth_sort(out_dir, threads)
        depth.minimap_depth_sort(out_dir, threads)
        contig_lengths = depth.get_contig_lengths(out_dir)
        depthsShort = depth.get_depths_from_bam(out_dir, "short", contig_lengths)
        depthsLong = depth.get_depths_from_bam(out_dir, "long", contig_lengths)
        circular_status = depth.get_contig_circularity(out_dir)
        summary_depth_df_short = depth.collate_depths(depthsShort,"short",contig_lengths)
        summary_depth_df_long = depth.collate_depths(depthsLong,"long",contig_lengths)
        depth.combine_depth_dfs(out_dir, summary_depth_df_short, summary_depth_df_long, prefix, circular_status)