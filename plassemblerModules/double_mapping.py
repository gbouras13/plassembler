import os
import sys
import subprocess as sp
import plassemblerModules


def double_mapping_analysis(out_dir, threads, logger):

    print('Extracting Reads mapping to Plasmids and Chromosome.')
    logger.info('Extracting Reads mapping to Plasmids and Chromosome.')
    extract_reads_mapping_to_plasmid_and_chromosome(out_dir, logger)
    print('Assembling Double Mapping Reads.')
    logger.info('Assembling Double Mapping Reads.')
    # assemble 
    plassemblerModules.run_unicycler(True, threads, logger,os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R1.fastq.gz"),os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R2.fastq.gz"),os.path.join(out_dir, "long_reads_mapping_to_plasmid_and_chromosome.fastq"), os.path.join(out_dir, "unicycler_plasmid_chromosome_map_output"))


###############################
# sub functions 
############################

##### Processing Sam to Bam #######

### concatenating and deduplicating



def deduplicate_fastqs(out_dir, threads, logger):
    concat_long = os.path.join(out_dir, "long_read_concat.fastq")
    dedup_long = os.path.join(out_dir, "long_read_dedup.fastq")
    concat_short_one = os.path.join(out_dir, "short_read_concat_R1.fastq.gz")
    dedup_short_one = os.path.join(out_dir, "short_read_dedup_R1.fastq.gz")
    concat_short_two = os.path.join(out_dir, "short_read_concat_R2.fastq.gz")
    dedup_short_two = os.path.join(out_dir, "short_read_dedup_R2.fastq.gz")
    try:
        dedup_long= sp.Popen(["seqkit", "rmdup", concat_long, "-j", threads, "-n", "-o",dedup_long], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(dedup_long.stdout, logger)
        dedup_s1= sp.Popen(["seqkit", "rmdup", concat_short_one, "-j", threads, "-n", "-o",dedup_short_one], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(dedup_s1.stdout, logger)
        dedup_s2= sp.Popen(["seqkit", "rmdup", concat_short_two, "-j", threads, "-n", "-o",dedup_short_two], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(dedup_s2.stdout, logger)
    except:
        sys.exit("Error with seqkit\n")  





##########################################################
#### reads mapping to both plasmids and chromosome
##########################################################

def extract_reads_mapping_to_plasmid_and_chromosome(out_dir, logger):
    long_fastq_chrom = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    long_fastq_non_chrom = os.path.join(out_dir, "long_read_non_chromosome_mapped.fastq")
    chrom_fastq_one_short = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
    chrom_fastq_two_short = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
    non_chrom_fastq_one_short = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq.gz")
    non_chrom_fastq_two_short = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq.gz")
    double_map_long = os.path.join(out_dir, "long_reads_mapping_to_plasmid_and_chromosome.fastq")
    double_map_one_short = os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R1.fastq.gz")
    double_map_two_short = os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R2.fastq.gz")

    try:
        bbmap_long= sp.Popen(["filterbyname.sh", "in="+ long_fastq_non_chrom, "names="+long_fastq_chrom, "out="+double_map_long], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(bbmap_long.stdout, logger)
        bbmap_s1= sp.Popen(["filterbyname.sh", "in="+ non_chrom_fastq_one_short, "names="+chrom_fastq_one_short, "out="+double_map_one_short], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(bbmap_s1.stdout, logger)
        bbmap_s2= sp.Popen(["filterbyname.sh", "in="+ non_chrom_fastq_two_short, "names="+chrom_fastq_two_short, "out="+double_map_two_short], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(bbmap_s2.stdout, logger)
    except:
        sys.exit("Error with bbmap\n")  



