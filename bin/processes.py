import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import logging
import glob
import gzip

def write_to_log(s, logger):
           while True:
                output = s.readline().decode()
                if output:
                    logger.log(logging.INFO, output)
                else:
                    break

def run_flye(out_dir, threads, logger):
    trim_long = os.path.join(out_dir, "porechop.fastq")
    try:
        flye = sp.Popen(["flye", "--nano-raw", trim_long, "--out-dir", out_dir, "--threads", threads ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(flye.stdout, logger)
    except:
        sys.exit("Error with Flye\n")  

def contig_count(out_dir):
    info_file =  os.path.join(out_dir, "assembly_info.txt")
    col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
    info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
    contig_count = len(info_df['seq_name'])
    print("Flye assembled " + str(contig_count) + " contigs.")
    return contig_count

def extract_chromosome(out_dir, chromosome_len):
    info_file =  os.path.join(out_dir, "assembly_info.txt")
    col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
    info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
    max_length = max(info_df['length'])
    chrom_contig = info_df[info_df['length'] == max_length].iloc[0]['seq_name']
    chrom_circ = info_df[info_df['length'] == max_length].iloc[0]['circ']
    # chromosome isn't circular or at least 80% of inputted chromosome length
    correct_chromosome = True
    if chrom_circ != 'Y' or max_length < int(chromosome_len)*0.8:
        correct_chromosome = False
    else:
        extract_chromosome_fasta(out_dir, chrom_contig)
        extract_plasmid_fastas(out_dir, chrom_contig)
    return correct_chromosome

def extract_chromosome_fasta(out_dir, contig_name):
    with open(os.path.join(out_dir, "chromosome.fasta"), 'w') as chrom_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
            # has to be contig of the chromosome
            if dna_record.id == contig_name:
                dna_header = "chromosome"
                dna_description = ""
                dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                SeqIO.write(dna_record, chrom_fa, 'fasta')


def extract_plasmid_fastas(out_dir, contig_name):
    with open(os.path.join(out_dir, "non_chromosome.fasta"), 'w') as non_chrom_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
            # has to be contig of the chromosome
            if dna_record.id != contig_name:
                dna_header = dna_record.id
                dna_description = ""
                dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                SeqIO.write(dna_record, non_chrom_fa, 'fasta')


def filtlong(input_long_reads, out_dir,min_length, min_quality,  logger):
    out_trim_long = os.path.join(out_dir, "filtlong.fastq")
    f = open(out_trim_long, "w")
    try:
        filtlong = sp.Popen(["filtlong", "--min_length", min_length, "--min_mean_q", min_quality,  input_long_reads, ], stdout=f, stderr=sp.PIPE) 
        write_to_log(filtlong.stderr, logger)
    except:
        sys.exit("Error with filtlong\n")  

def porechop(out_dir,threads, logger):
    filtlong_reads = os.path.join(out_dir, "filtlong.fastq")
    porechop_reads = os.path.join(out_dir, "porechop.fastq")
    try:
        porechop = sp.Popen(["porechop", "-i", filtlong_reads, "-o", porechop_reads, "-t", threads ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(porechop.stderr, logger)
    except:
        sys.exit("Error with filtlong\n")  


def trim_short_read(short_one, short_two, out_dir,  logger):
    out_one = os.path.join(out_dir, "trimmed_R1.fastq")
    out_two = os.path.join(out_dir, "trimmed_R2.fastq")
    try:
        fastp = sp.Popen(["fastp", "--in1", short_one, "--in2", short_two, "--out1", out_one, "--out2", out_two ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(fastp.stderr, logger)
    except:
        sys.exit("Error with Fastp\n")  


 ################################################
# number 3 - hybrid map to plasmids
################################################
def plasmid_assembly(out_dir, threads, logger):
    #### indexing contigs
    print('Indexing Plasmid Contigs.')
    logger.info("Indexing Plasmid Contigs.")
    index_fasta( os.path.join(out_dir, "non_chromosome.fasta"),  logger)
    print('Indexing Chromosome.')
    logger.info("Indexing Chromosome.")
    index_fasta( os.path.join(out_dir, "chromosome.fasta"),  logger)

    ##### Mapping #######

    #### long reads mapping
    print('Mapping Long Reads to Plasmid Contigs.')
    logger.info('Mapping Long Reads to Plasmid Contigs.')
    minimap_long_reads(False, out_dir, threads, logger)
    print('Mapping Long Reads to Chromosome.')
    logger.info('Mapping Long Reads to Chromosome.')
    minimap_long_reads(True, out_dir, threads, logger)

    #### short reads mapping
    print('Mapping Short Reads to Plasmid Contigs')
    logger.info('Mapping Short Reads to Plasmid Contigs')
    bwa_map_short_reads( out_dir, False, threads,  logger)

    print('Mapping Short Reads to Chromosome Contig')
    logger.info('Mapping Short Reads to Chromosome Contig')
    bwa_map_short_reads( out_dir, True, threads,  logger)

    #### Processing bams ######
    print('Processing Bams.')
    logger.info('Processing Bams.')

    sam_to_bam( out_dir, "chromosome_long", threads,  logger)
    sam_to_bam( out_dir, "non_chromosome_long", threads,  logger)
    sam_to_bam( out_dir, "chromosome_short", threads,  logger)
    sam_to_bam( out_dir, "non_chromosome_short", threads,  logger)

    ### extractng mapped_unmapped bams ###
    bam_to_mapped_or_unmapped(out_dir, "chromosome_long", threads, logger)
    bam_to_mapped_or_unmapped(out_dir, "non_chromosome_long", threads, logger)
    bam_to_mapped_or_unmapped(out_dir, "chromosome_short", threads, logger)
    bam_to_mapped_or_unmapped(out_dir, "non_chromosome_short", threads, logger)

    ### extracting fastqs
    print('Extracting Fastqs.')
    logger.info('Extracting Fastqs.')

    extract_long_fastq(out_dir, "chromosome", logger)
    extract_long_fastq(out_dir, "non_chromosome", logger)
    extract_short_fastq( out_dir, "chromosome",  threads,  logger)   
    extract_short_fastq( out_dir, "non_chromosome",  threads,  logger)   

    # concatenating and deduplicating fastqs

    print('Concatenating and Deduplicating Fastqs.')
    logger.info('Concatenating and Deduplicating Fastqs.')

    concatenate_fastqs(out_dir,logger)
    deduplicate_fastqs(out_dir, threads, logger)

    # running unicycler
    print('Running Unicycler')
    logger.info('Running Unicycler')
    # short only flag
    short_r1 = os.path.join(out_dir, "short_read_dedup_R1.fastq.gz")
    short_r2 = os.path.join(out_dir, "short_read_dedup_R2.fastq.gz")
    long_reads = os.path.join(out_dir, "long_read_dedup.fastq")
    unicycler(False, threads, logger, short_r1, short_r2, long_reads, os.path.join(out_dir, "unicycler_output"))


def double_mapping_analysis(out_dir, threads, logger):

    print('Extracting Reads mapping to Plasmids and Chromosome.')
    logger.info('Extracting Reads mapping to Plasmids and Chromosome.')
    extract_reads_mapping_to_plasmid_and_chromosome(out_dir, logger)

    print('Assembling Double Mapping Reads.')
    logger.info('Assembling Double Mapping Reads.')
    # assemble 
    unicycler(True, threads, logger,os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R1.fastq.gz"),os.path.join(out_dir, "short_reads_mapping_to_plasmid_and_chromosome_R2.fastq.gz"),os.path.join(out_dir, "long_reads_mapping_to_plasmid_and_chromosome.fastq"), os.path.join(out_dir, "unicycler_plasmid_chromosome_map_output"))


###############################
# sub functions 
############################

#### indexing #####

def index_fasta(fasta,  logger):
    try:
        bwa_index = sp.Popen(["bwa", "index", fasta ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(bwa_index.stdout, logger)
    except:
        sys.exit("Error with bwa index\n")  

 ##### Mapping #######

# minimap for long reads
def minimap_long_reads(flag, out_dir, threads, logger):
    input_long_reads = os.path.join(out_dir, "porechop.fastq")
    # chromosome is a flag for mapping chromosome or not
    if flag == True:
        fasta = os.path.join(out_dir, "chromosome.fasta")
        sam = os.path.join(out_dir, "long_read_chromosome.sam")
    else:
        fasta = os.path.join(out_dir, "non_chromosome.fasta")
        sam = os.path.join(out_dir, "long_read_non_chromosome.sam")
    f = open(sam, "w")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, fasta, input_long_reads ], stdout=f, stderr=sp.PIPE) 
        write_to_log(minimap.stderr, logger)
    except:
        sys.exit("Error with minimap2\n")  

# bwa for short reads 

def bwa_map_short_reads(out_dir,flag, threads, logger):
    # chromosome is a flag for mapping chromosome or not
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    if flag == False:
        fasta = os.path.join(out_dir, "non_chromosome.fasta")
        sam = os.path.join(out_dir, "short_read_non_chromosome.sam")
    else:
        fasta = os.path.join(out_dir, "chromosome.fasta")
        sam = os.path.join(out_dir, "short_read_chromosome.sam") 
    f = open(sam, "w")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, fasta, trim_one, trim_two ], stdout=f, stderr=sp.PIPE) 
        write_to_log(bwa_map.stderr, logger)
    except:
        sys.exit("Error with bwa mem\n")  

 ##### Processing Sam to Bam #######

# sam to bam
def sam_to_bam(out_dir,flag,threads, logger):
    #  flag for what bams to index
    if flag == "chromosome_long":
        sam = os.path.join(out_dir, "long_read_chromosome.sam")
        bam = os.path.join(out_dir, "long_read_chromosome.bam")
    elif flag == "non_chromosome_long":
        sam = os.path.join(out_dir, "long_read_non_chromosome.sam")
        bam = os.path.join(out_dir, "long_read_non_chromosome.bam")
    elif flag == "non_chromosome_short":
        sam = os.path.join(out_dir, "short_read_non_chromosome.sam")
        bam = os.path.join(out_dir, "short_read_non_chromosome.bam")
    else: # chromosome_short
        sam = os.path.join(out_dir, "short_read_chromosome.sam")
        bam = os.path.join(out_dir, "short_read_chromosome.bam")
    f = open(bam, "w")
    try:
        sam_to_bam = sp.Popen(["samtools", "view", "-h", "-@", threads, "-b", sam], stdout=f, stderr=sp.PIPE) 
        write_to_log(sam_to_bam.stderr, logger)
    except:
        sys.exit("Error with samtools sam_to_bam.\n")  


#### extract mapped or unmapped bams ######

def bam_to_mapped_or_unmapped(out_dir, flag, threads, logger):
        #  flag for how to reduce size of bams
    if flag == "chromosome_long":
        bam = os.path.join(out_dir, "long_read_chromosome.bam")
        map_bam = os.path.join(out_dir, "long_read_chromosome_unmapped.bam")
        map_flag = "-f"
    elif flag == "chromosome_short":
        bam = os.path.join(out_dir, "short_read_chromosome.bam")
        map_bam = os.path.join(out_dir, "short_read_chromosome_unmapped.bam")
        map_flag = "-f"
    elif flag == "non_chromosome_short":
        bam = os.path.join(out_dir, "short_read_non_chromosome.bam")
        map_bam = os.path.join(out_dir, "short_read_non_chromosome_mapped.bam")
        map_flag = "-F"
    elif flag == "non_chromosome_long":
        bam = os.path.join(out_dir, "long_read_non_chromosome.bam")
        map_bam = os.path.join(out_dir, "long_read_non_chromosome_mapped.bam")
        map_flag = "-F"
    f = open(map_bam, "w")
    try:
        bam_to_map = sp.Popen(["samtools", "view", "-h", "-@", threads, map_flag, "4", bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(bam_to_map.stderr, logger)
    except:
        sys.exit("Error with samtools bam_to_unmap. \n")  

#### extract fastqs ######
# long
def extract_long_fastq(out_dir, flag, logger):
    # want unmapped for chromosome
    if flag == "chromosome":
        bam = os.path.join(out_dir, "long_read_chromosome_unmapped.bam")
        long_fastq = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    else:
        bam = os.path.join(out_dir, "long_read_non_chromosome_mapped.bam")
        long_fastq = os.path.join(out_dir, "long_read_non_chromosome_mapped.fastq")
    f = open(long_fastq, "w")
    try:
        extract_map_fastq_long= sp.Popen(["samtools", "bam2fq", bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(extract_map_fastq_long.stderr, logger)
    except:
        sys.exit("Error with samtools bam2fastq\n")  

# short

def extract_short_fastq(out_dir, flag, threads, logger):
    if flag == "chromosome":
        bam = os.path.join(out_dir, "short_read_chromosome_unmapped.bam")
        fastq_one = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
        fastq_two = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
        map_flag = "-f"
    else: # non chromosome
        bam = os.path.join(out_dir, "short_read_non_chromosome_mapped.bam")
        fastq_one = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq.gz")
        fastq_two = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq.gz")
        map_flag = "-F"
    try:
        extract_map_fastq= sp.Popen(["samtools", "fastq", "-@", threads, map_flag, "4", bam, "-1", fastq_one, "-2", fastq_two, "-0", "/dev/null", "-s", "/dev/null", "-n"], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(extract_map_fastq.stdout, logger)
    except:
        sys.exit("Error with samtools extract_map_fastq\n")  

### concatenating and deduplicating

def concatenate_fastqs(out_dir,logger):
    long_fastq_chrom = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    long_fastq_non_chrom = os.path.join(out_dir, "long_read_non_chromosome_mapped.fastq")
    chrom_fastq_one_short = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
    chrom_fastq_two_short = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
    non_chrom_fastq_one_short = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq.gz")
    non_chrom_fastq_two_short = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq.gz")
    long_file = open(os.path.join(out_dir, "long_read_concat.fastq"), "w")
    short_one_file = open(os.path.join(out_dir, "short_read_concat_R1.fastq.gz"), "w")
    short_two_file = open(os.path.join(out_dir, "short_read_concat_R2.fastq.gz"), "w")
    try:
        concatenate_single(long_fastq_chrom, long_fastq_non_chrom, long_file, logger)
        concatenate_single(chrom_fastq_one_short, non_chrom_fastq_one_short, short_one_file, logger)
        concatenate_single(chrom_fastq_two_short, non_chrom_fastq_two_short, short_two_file, logger)
    except:
        sys.exit("Error with concatenate_fastqs\n")  
        
def concatenate_single(fastq_in_1, fastq_in_2, fastq_out, logger):
        concat_fastq = sp.Popen(["cat", fastq_in_1, fastq_in_2 ], stdout=fastq_out, stderr=sp.PIPE)
        write_to_log(concat_fastq.stderr, logger) 



def deduplicate_fastqs(out_dir, threads, logger):
    concat_long = os.path.join(out_dir, "long_read_concat.fastq")
    dedup_long = os.path.join(out_dir, "long_read_dedup.fastq")
    concat_short_one = os.path.join(out_dir, "short_read_concat_R1.fastq.gz")
    dedup_short_one = os.path.join(out_dir, "short_read_dedup_R1.fastq.gz")
    concat_short_two = os.path.join(out_dir, "short_read_concat_R2.fastq.gz")
    dedup_short_two = os.path.join(out_dir, "short_read_dedup_R2.fastq.gz")
    try:
        dedup_long= sp.Popen(["seqkit", "rmdup", concat_long, "-j", threads, "-n", "-o",dedup_long], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(dedup_long.stdout, logger)
        dedup_s1= sp.Popen(["seqkit", "rmdup", concat_short_one, "-j", threads, "-n", "-o",dedup_short_one], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(dedup_s1.stdout, logger)
        dedup_s2= sp.Popen(["seqkit", "rmdup", concat_short_two, "-j", threads, "-n", "-o",dedup_short_two], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(dedup_s2.stdout, logger)
    except:
        sys.exit("Error with seqkit\n")  


### unicycler and deduplicating

def unicycler(short_only, threads, logger, short_one, short_two, long, unicycler_output_dir):
    try:
        if short_only == False:
            unicycler= sp.Popen(["unicycler", "-1", short_one, "-2", short_two, "-l", long, "-t", threads, "-o",  unicycler_output_dir ], stdout=sp.PIPE, stderr=sp.PIPE) 
        else: 
            unicycler= sp.Popen(["unicycler", "-1", short_one, "-2", short_two, "-t", threads, "-o",  unicycler_output_dir ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(unicycler.stdout, logger)
    except:
        sys.exit("Error with Unicycler\n")  






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
        write_to_log(bbmap_long.stdout, logger)
        bbmap_s1= sp.Popen(["filterbyname.sh", "in="+ non_chrom_fastq_one_short, "names="+chrom_fastq_one_short, "out="+double_map_one_short], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(bbmap_s1.stdout, logger)
        bbmap_s2= sp.Popen(["filterbyname.sh", "in="+ non_chrom_fastq_two_short, "names="+chrom_fastq_two_short, "out="+double_map_two_short], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(bbmap_s2.stdout, logger)
    except:
        sys.exit("Error with bbmap\n")  



####################################################
# cleanup
##########################################################

def remove_intermediate_files(out_dir):

    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq.gz") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.bam") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.sa") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.sam") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.amb") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.ann") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.pac") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.bwt") ], shell=True)
    sp.run(["rm", "-rf", os.path.join(out_dir,"00-assembly") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"10-consensus") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"20-repeat") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"30-contigger") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"40-polishing") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"params.json") ])
    # delete flye assemble files
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"non_chromosome.fasta") ])


def move_and_copy_files(out_dir, prefix, fail):
    # move flye output into dir
    sp.run(["mkdir", "-p", os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"assembly.fasta"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"assembly_info.txt"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv",  os.path.join(out_dir,"flye.log"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv", os.path.join(out_dir,"assembly_graph.gfa"), os.path.join(out_dir,"flye_output") ])
    sp.run(["mv", os.path.join(out_dir,"assembly_graph.gv"), os.path.join(out_dir,"flye_output") ])
    if fail == False:
         # move unicycler output to main directory
        sp.run(["cp", os.path.join(out_dir,"unicycler_output", "assembly.fasta"), os.path.join(out_dir, prefix + "_plasmids.fasta") ])
        sp.run(["cp", os.path.join(out_dir,"unicycler_output", "assembly.gfa"), os.path.join(out_dir, prefix + "_plasmids.gfa") ])
    else:
        # to touch empty versions of the output files if no plasmids 
        touch_output_fail_files(out_dir, prefix)


# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, 'a'):
        os.utime(path, None)

# to create empty plasmids fasta and gfa files
def touch_output_fail_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_plasmids.fasta"))
    touch_file(os.path.join(out_dir, prefix + "_plasmids.gfa"))

