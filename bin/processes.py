import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import logging

def write_to_log(s, logger):
           while True:
                output = s.readline().decode()
                if output:
                    logger.log(logging.INFO, output)
                else:
                    break

def run_flye(input_long_reads, out_dir, threads, logger):
    try:
        flye = sp.Popen(["flye", "--nano-raw", input_long_reads, "--out-dir", out_dir, "--threads", threads ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(flye.stderr, logger)
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
    if chrom_circ != 'Y' or max_length < chromosome_len*0.8:
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

def trim_short_read(short_one, short_two, out_dir,  logger):
    out_one = os.path.join(out_dir, "trimmed_R1.fastq")
    out_two = os.path.join(out_dir, "trimmed_R2.fastq")
    try:
        fastp = sp.Popen(["fastp", "--in1", short_one, "--in2", short_two, "--out1", out_one, "--out2", out_two ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(fastp.stderr, logger)
    except:
        sys.exit("Error with Fastp\n")  

def index_chromosome(out_dir,  logger):
    chrom_fasta = os.path.join(out_dir, "chromosome.fasta")
    try:
        bwa_index = sp.Popen(["bwa", "index", chrom_fasta ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(bwa_index.stdout, logger)
    except:
        sys.exit("Error with bwa index\n")  

def bwa_map_chromosome(out_dir,threads, logger):
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    chrom_fasta = os.path.join(out_dir, "chromosome.fasta")
    sam = os.path.join(out_dir, "short_read.sam")
    f = open(sam, "w")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, chrom_fasta, trim_one, trim_two ], stdout=f, stderr=sp.PIPE) 
        write_to_log(bwa_map.stderr, logger)
    except:
        sys.exit("Error with bwa mem\n")  

def minimap_chromosome(input_long_reads, out_dir, threads, logger):
    chrom_fasta = os.path.join(out_dir, "chromosome.fasta")
    sam = os.path.join(out_dir, "long_read.sam")
    f = open(sam, "w")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, chrom_fasta, input_long_reads ], stdout=f, stderr=sp.PIPE) 
        write_to_log(minimap.stderr, logger)
    except:
        sys.exit("Error with minimap2\n")  

def sam_to_bam_long(out_dir, threads, logger):
    sam = os.path.join(out_dir, "long_read.sam")
    bam = os.path.join(out_dir, "long_read.bam")
    f = open(bam, "w")
    try:
        sam_to_bam = sp.Popen(["samtools", "view", "-h", "-@", threads, "-b", sam], stdout=f, stderr=sp.PIPE) 
        write_to_log(sam_to_bam.stderr, logger)
    except:
        sys.exit("Error with samtools sam_to_bam.\n")  

def sam_to_bam(out_dir, threads, logger):
    sam = os.path.join(out_dir, "short_read.sam")
    bam = os.path.join(out_dir, "short_read.bam")
    f = open(bam, "w")
    try:
        sam_to_bam = sp.Popen(["samtools", "view", "-h", "-@", threads, "-b", sam], stdout=f, stderr=sp.PIPE) 
        write_to_log(sam_to_bam.stderr, logger)
    except:
        sys.exit("Error with samtools sam_to_bam.\n")  

def bam_to_unmap(out_dir, threads, logger):
    bam = os.path.join(out_dir, "short_read.bam")
    unmap_bam = os.path.join(out_dir, "unmap.bam")
    f = open(unmap_bam, "w")
    try:
        bam_to_unmap = sp.Popen(["samtools", "view", "-h", "-@", threads, "-f", "4", bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(bam_to_unmap.stderr, logger)
    except:
        sys.exit("Error with samtools bam_to_unmap. \n")  


def bam_to_unmap_long(out_dir, threads, logger):
    bam = os.path.join(out_dir, "long_read.bam")
    unmap_bam = os.path.join(out_dir, "long_unmap.bam")
    f = open(unmap_bam, "w")
    try:
        bam_to_unmap = sp.Popen(["samtools", "view", "-h", "-@", threads, "-f", "4", bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(bam_to_unmap.stderr, logger)
    except:
        sys.exit("Error with samtools bam_to_unmap. \n")  


def extract_unmap_fastq(out_dir, threads, logger):
    unmap_bam = os.path.join(out_dir, "unmap.bam")
    unmap_one = os.path.join(out_dir, "unmapped_R1.fastq.gz")
    unmap_two = os.path.join(out_dir, "unmapped_R2.fastq.gz")
    try:
        extract_unmap_fastq= sp.Popen(["samtools", "fastq", "-@", threads, "-f", "4", unmap_bam, "-1", unmap_one, "-2", unmap_two, "-0", "/dev/null", "-s", "/dev/null", "-n"], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(extract_unmap_fastq.stderr, logger)
    except:
        sys.exit("Error with samtools extract_unmap_fastq\n")  

def extract_unmap_fastq_long(out_dir, logger):
    unmap_bam = os.path.join(out_dir, "long_unmap.bam")
    unmap_long = os.path.join(out_dir, "unmapped_long.fastq")
    f = open(unmap_long, "w")
    try:
        extract_unmap_fastq_long= sp.Popen(["samtools", "bam2fq", unmap_bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(extract_unmap_fastq_long.stderr, logger)
    except:
        sys.exit("Error with samtools extract_unmap_fastq_long\n")  

def unicycler(out_dir, threads, logger):
    unmap_one = os.path.join(out_dir, "unmapped_R1.fastq.gz")
    unmap_two = os.path.join(out_dir, "unmapped_R2.fastq.gz")
    unmap_long = os.path.join(out_dir, "unmapped_long.fastq")
    try:
        unicycler= sp.Popen(["unicycler", "-1", unmap_one, "-2", unmap_two, "-l", unmap_long, "-t", threads, "-o", os.path.join(out_dir, "unicycler") ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(unicycler.stdout, logger)
    except:
        sys.exit("Error with Unicycler\n")  


def remove_intermediate_files(out_dir):
    sp.run(["rm", "-rf", os.path.join(out_dir,"trimmed_R1.fastq") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"trimmed_R2.fastq") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"unmapped.bam") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"short_read.bam") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"short_read.sam") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"unmap.bam") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta.sa") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta.amb") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta.ann") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta.pac") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"chromosome.fasta.bwt") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"00-assembly") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"10-consensus") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"20-repeat") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"30-contigger") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"40-polishing") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"params.json") ])

#### 

def extract_plasmid_fastas(out_dir, contig_name):
    with open(os.path.join(out_dir, "non_chromosome.fasta"), 'w') as non_chrom_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
            # has to be contig of the chromosome
            if dna_record.id != contig_name:
                dna_header = dna_record.id
                dna_description = ""
                dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                SeqIO.write(dna_record, non_chrom_fa, 'fasta')

def index_non_chrom(out_dir,  logger):
    chrom_fasta = os.path.join(out_dir, "non_chromosome.fasta")
    try:
        bwa_index = sp.Popen(["bwa", "index", chrom_fasta ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(bwa_index.stdout, logger)
    except:
        sys.exit("Error with bwa index\n")  

def bwa_map_non_chrom(out_dir,threads, logger):
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    chrom_fasta = os.path.join(out_dir, "non_chromosome.fasta")
    sam = os.path.join(out_dir, "short_read_non_chrom.sam")
    f = open(sam, "w")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, chrom_fasta, trim_one, trim_two ], stdout=f, stderr=sp.PIPE) 
        write_to_log(bwa_map.stderr, logger)
    except:
        sys.exit("Error with bwa mem\n")  

def sam_to_bam_non_chrom(out_dir, threads, logger):
    sam = os.path.join(out_dir, "short_read_non_chrom.sam")
    bam = os.path.join(out_dir, "short_read_non_chrom.bam")
    f = open(bam, "w")
    try:
        sam_to_bam = sp.Popen(["samtools", "view", "-h", "-@", threads, "-b", sam], stdout=f, stderr=sp.PIPE) 
        write_to_log(sam_to_bam.stderr, logger)
    except:
        sys.exit("Error with samtools sam_to_bam.\n")  

def bam_to_map_non_chrom(out_dir, threads, logger):
    bam = os.path.join(out_dir, "short_read_non_chrom.bam")
    map_bam = os.path.join(out_dir, "map_non_chrom.bam")
    f = open(map_bam, "w")
    try:
        bam_to_map = sp.Popen(["samtools", "view", "-h", "-@", threads, "-F", "4", bam], stdout=f, stderr=sp.PIPE) 
        write_to_log(bam_to_map.stderr, logger)
    except:
        sys.exit("Error with samtools bam_to_unmap. \n")  

def extract_map_non_chrom_fastq(out_dir, threads, logger):
    map_bam = os.path.join(out_dir, "map_non_chrom.bam")
    map_one = os.path.join(out_dir, "mapped_non_chrom_R1.fastq.gz")
    map_two = os.path.join(out_dir, "mapped_non_chrom_R2.fastq.gz")
    try:
        extract_map_fastq= sp.Popen(["samtools", "fastq", "-@", threads, "-f", "4", map_bam, "-1", map_one, "-2", map_two, "-0", "/dev/null", "-s", "/dev/null", "-n"], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(extract_map_fastq.stdout, logger)
    except:
        sys.exit("Error with samtools extract_map_fastq\n")  

def unicycler_map(out_dir, threads, logger):
    map_one = os.path.join(out_dir, "unmapped_R1.fastq.gz")
    map_two = os.path.join(out_dir, "unmapped_R2.fastq.gz")
    try:
        unicycler= sp.Popen(["unicycler", "-1", map_one, "-2", map_two, "-t", threads, "-o", os.path.join(out_dir, "unicycler_map") ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(unicycler.stdout, logger)
    except:
        sys.exit("Error with Unicycler\n")  


