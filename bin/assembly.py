
#########################################
##### some code taken from & modified
##### https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/scripts/get_depths.py
#########################################
from Bio import SeqIO
import os
import sys
import subprocess as sp
import statistics
import numpy as np
import pandas as pd
import mapping
import concat
import depth
from Bio.SeqRecord import SeqRecord


def get_depth_assembly(input_chromosome,input_plasmids, out_dir, logger,  threads, prefix):
    """ wrapper function to get depth of each plasmid in assembly mode - short and long reads
    :param prefix: prefix (default plassembler)
    :param out_dir:  Output Directory
    :param threads: threads
    :param logger: logger
    :return: 
    """
    # combine input fasta
    f = open(os.path.join(out_dir,"input_tmp.fasta"), "w")

    sp.run(["cat", input_chromosome, input_plasmids ], stdout=f, stderr=sp.PIPE)

    # set input_fasta on the rest
    input_fasta = os.path.join(out_dir,"input.fasta")

    # rename the first chromosome 
    rename_first_contig_chromosome(os.path.join(out_dir,"input_tmp.fasta"), input_fasta)

    mapping.index_fasta(input_fasta,  logger)
    minimap_short_depth_sort_assembly(input_fasta, out_dir, threads)
    minimap_long_depth_sort_assembly(input_fasta, out_dir, threads)
    contig_lengths = get_contig_lengths_assembly(input_fasta)
    depthsShort = depth.get_depths_from_bam(out_dir, "short", contig_lengths)
    depthsLong = depth.get_depths_from_bam(out_dir, "long", contig_lengths)
    circular_status = get_contig_circularity_assembly(input_fasta)
    summary_depth_df_short = depth.collate_depths(depthsShort,"short",contig_lengths)
    summary_depth_df_long = depth.collate_depths(depthsLong,"long",contig_lengths)
    depth.combine_depth_dfs(out_dir, summary_depth_df_short, summary_depth_df_long, prefix, circular_status)

def get_depth_assembly_long_only(input_chromosome, input_plasmids, out_dir, logger,  threads, prefix):
    """ wrapper function to get depth of each plasmid - long only
    :param prefix: prefix (default plassembler)
    :param out_dir:  Output Directory
    :param threads: threads
    :param logger: logger
    :return: 
    """
    # combine input fasta
    f = open(os.path.join(out_dir,"input_tmp.fasta"), "w")

    sp.run(["cat", input_chromosome, input_plasmids ], stdout=f, stderr=sp.PIPE)

    # set input_fasta on the rest
    input_fasta = os.path.join(out_dir,"input.fasta")

    # rename the first chromosome 
    rename_first_contig_chromosome(os.path.join(out_dir,"input_tmp.fasta"), input_fasta)

    mapping.index_fasta(input_fasta,  logger)
    minimap_long_depth_sort_assembly(input_fasta, out_dir, threads)
    contig_lengths = get_contig_lengths_assembly(input_fasta)
    depthsLong = depth.get_depths_from_bam(out_dir, "long", contig_lengths)
    circular_status = get_contig_circularity_assembly(input_fasta)
    summary_depth_df_long = depth.collate_depths(depthsLong,"long",contig_lengths)
    kmer_final_output(out_dir, summary_depth_df_long, prefix, circular_status)


def concatenate_chrom_plasmids(out_dir, logger):
    """ concatenates chromosome and plasmids
    :param out_dir:  Output Directory
    :param logger: logger
    :return: 
    """
    chrom_fasta =os.path.join(out_dir,"chromosome.fasta")
    plas_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    concat_fasta = open(os.path.join(out_dir, "combined.fasta"), "w")
    try:
        concat.concatenate_single(chrom_fasta, plas_fasta, concat_fasta, logger)
    except:
        sys.exit("Error with concatenate_fastas\n")  


# get lengths of contigs
def get_contig_lengths_assembly(input_fasta):
    """ gets contig lengths of combined chrom and plasmids fastas in input fasta
    :param out_dir:  Output Directory
    :return: contig_lengths: dictionary of headers and lengths
    """
    contig_lengths = {}
    for dna_record in SeqIO.parse(input_fasta, 'fasta'):
        plas_len = len(dna_record.seq)
        dna_header = dna_record.id
        contig_lengths[dna_header] = plas_len
    return contig_lengths

def rename_first_contig_chromosome(input_fasta, renamed_fasta):

    # rename the first contig as chromosome 
    with open(input_fasta, "r") as f_in, open(renamed_fasta, "w") as f_out:
        # Parse the input FASTA file
        records = list(SeqIO.parse(f_in, "fasta"))
        # Rename the first record
        records[0].id = "chromosome"
        records[0].description = ""
        # Write the modified records to the output FASTA file
        SeqIO.write(records, f_out, "fasta")



# get circular status of contigs
def get_contig_circularity_assembly(input_fasta):
    """ gets circularity of contigs
    :param out_dir:  Output Directory
    :return: circular_status: dictionary of contig header and circular status
    """
    circular_status = {}
    # add circularity
    for dna_record in SeqIO.parse(input_fasta, 'fasta'):
        dna_header = dna_record.id
        # check if circular is in unicycler output description
        if "circular=true" in dna_record.description:
            circular_status[dna_header] = "circular"
        # circular chromsome
        elif "chromosome" in dna_record.id:
            circular_status[dna_header] = "circular"
        else:
            circular_status[dna_header] = "not_circular"
    return circular_status


def minimap_short_depth_sort_assembly(input_fasta, out_dir, threads):
    """ maps short reads using bwa to combined fasta and sorts bam in assembly mode
    :param out_dir:  Output Directory
    :return: threads: threads
    """
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    bam = os.path.join(out_dir, "combined_sorted.bam")
    try:
        minimap2_map = sp.Popen(["minimap2", "-ax", "sr", "-t", threads, input_fasta, trim_one, trim_two ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=minimap2_map.stdout, stderr=sp.DEVNULL ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with minim2p or samtools sort.\n")  


def minimap_long_depth_sort_assembly(input_fasta, out_dir, threads):
    """ maps long reads using minimap2 to combined fasta and sorts bam in assembly mode
    :param out_dir:  out_dir
    :param: threads: threads
    """
    input_long_reads = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    bam = os.path.join(out_dir, "combined_sorted_long.bam")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, input_fasta, input_long_reads ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        samtools_sort = sp.Popen(["samtools", "sort", "-@", threads, "-o", bam, "-" ], stdin=minimap.stdout, stderr=sp.DEVNULL ) 
        samtools_sort.communicate()[0]
    except:
        sys.exit("Error with mapping and sorting\n")  


def get_depths_from_bam(out_dir, shortFlag, contig_lengths):
    """ maps runs samtools depth on bam
    :param out_dir:  out_dir
    :param: shortFlag: string either "short" or "long"
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: depths: dictionary of contigs and depths
    """
    depths = {}
    if shortFlag == "short":
        filename = os.path.join(out_dir, "combined_sorted.bam")
    else: # long
        filename = os.path.join(out_dir, "combined_sorted_long.bam")
    for repName, repLength in contig_lengths.items():
        depths[repName] = [0] * repLength
    depthCommand = ['samtools', 'depth', filename]
    with open(os.devnull, 'wb') as devNull:
        depthOutput = sp.check_output(depthCommand, stderr=devNull).decode()
    for line in depthOutput.splitlines(): # parse output
        parts = line.strip().split('\t')
        repName = parts[0]
        depths[repName][int(parts[1])-1] = int(parts[2])
    return depths


def collate_depths(depths, shortFlag, contig_lengths):
    """ calculates summary statistics for all depths
    :param depths:  dictionary of contigs and depths from get_depths_from_bam
    :param: shortFlag: string either "short" or "long"
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: summary_df: pandas df of depth summary statistics
    """
    # define the columns of dataframe
    contig_names = []
    contig_length = []    
    mean_depth_col = []
    sd_depth_col = []
    q25_depth = []
    q75_depth = []
    # iterate over the conitgs
    for replicon_name, base_depths in depths.items():
        replicon_length = contig_lengths[replicon_name]
        unmaskedDepths = []
        for depth in enumerate(base_depths):
            unmaskedDepths.append(depth)
        try:
            mean_depth = round(statistics.mean(base_depths),2)
            depth_stdev = round(statistics.stdev(base_depths),2)
            q25, q75 = np.percentile(base_depths, [25 ,75])
            q25, q75 = int(q25), int(q75)
            # save the chromosome depth 
            if replicon_name == "chromosome":
                chromosome_depth = mean_depth
        except statistics.StatisticsError: # if can't calculate
            mean_depth, depth_stdev, q25, q75 = 'NA', 'NA', 'NA', 'NA'
        # append to list
        contig_names.append(replicon_name)
        contig_length.append(replicon_length)
        mean_depth_col.append(mean_depth)
        sd_depth_col.append(depth_stdev)
        q25_depth.append(q25)
        q75_depth.append(q75)
    # make summary df    
    if shortFlag == "short":
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'length': contig_length,
        'mean_depth_short': mean_depth_col,
        'sd_depth_short': sd_depth_col, 
        'q25_depth_short': q25_depth,
        'q75_depth_short': q75_depth
        })
        summary_df['plasmid_copy_number_short'] = round(summary_df['mean_depth_short'] / chromosome_depth,2)
    else: # long
        summary_df = pd.DataFrame(
        {'contig': contig_names,
        'mean_depth_long': mean_depth_col,
        'sd_depth_long': sd_depth_col, 
        'q25_depth_long': q25_depth,
        'q75_depth_long': q75_depth
        })
        summary_df['plasmid_copy_number_long'] = round(summary_df['mean_depth_long'] / chromosome_depth,2)
    # return df         
    return(summary_df)


def kmer_final_output(out_dir, df_long, prefix, circular_status):
    """ final output for kmer mode
    :param out_dir:  output directory
    :param df_long: long depth summary df
    :param: prefix: prefix - default plassembler
    :param: circular_status: dictionary of contig header and circular status
    """

     # add in circularity info 

    df_long['circularity'] = df_long['contig'].map(circular_status)
    out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
    with open(out_file, 'w') as f:
        df_long.to_csv(f, sep="\t", index=False, header=True)
    



from Bio import SeqIO
import pandas as pd
import os
import subprocess as sp
import log
import sys
import cleanup


def mash_sketch_assembly( out_dir, input_fasta,  logger):
    """
    Runs mash on input fastas in assembly mode
    :param input_fasta in
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    try:
        # copy fasta first
        mash_sketch = sp.Popen(["mash", "sketch",  input_fasta, "-i" ], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stdout, logger)
    except:
        sys.exit("Error with mash sketch.\n")  


def run_mash(out_dir, plasmid_sketch, plassembler_db_dir, logger):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param plassembler_db_dir: plassembler db directory
    :param logger: logger
    :return:
    """

    plsdb_sketch = os.path.join(plassembler_db_dir, "plsdb.msh")

    mash_tsv = os.path.join(out_dir,"mash.tsv")

    outFile = open(mash_tsv, "w")
    try:
        mash_sketch = sp.Popen(["mash", "dist",  plasmid_sketch, plsdb_sketch, "-v", "0.1", "-d", "0.1", "-i" ], stdout=outFile, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stderr, logger)
    except:
        sys.exit("Error with mash dist.\n")  

def get_contig_count( fasta):
    """
    Process mash output
    :param out_dir: output directory
    :return: i: int contig_count
    """
    i = 0
    for dna_record in SeqIO.parse(fasta, 'fasta'): 
            i += 1
    return i


def process_mash_tsv(out_dir, plassembler_db_dir, prefix, input_fasta):
    """
    Process mash output
    :param out_dir: output directory
    :return: mash_empty: boolean whether there was a mash hit
    """

    contig_count = get_contig_count(input_fasta)

    mash_tsv = os.path.join(out_dir,"mash.tsv")

    col_list = ["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"]

    mash_empty = is_file_empty(mash_tsv)
    # instantiate tophits list
    tophits_mash_df = []

    if mash_empty == False:   
        mash_df = pd.read_csv(mash_tsv, delimiter= '\t', index_col=False, names=col_list ) 
        # get list of genes
        contigs = range(1, contig_count+1)

        # instantiate tophits list
        tophits = []

        for contig in contigs:
            hit_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True)
            hits = len(hit_df['mash_distance'])
            # add only if there is a hit
            if hits > 0:
                tmp_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True).loc[0]
                tophits.append([tmp_df.contig, tmp_df.ACC_NUCCORE, tmp_df.mash_distance, tmp_df.mash_pval, tmp_df.mash_matching_hashes])
            # create tophits df
        tophits_mash_df = pd.DataFrame(tophits, columns=["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])

        # read in the plasdb tsv 
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        # with open(plsdb_tsv_file, 'rb') as f:
        #     result = chardet.detect(f.readline())
        #     print(result)
        cols = ["UID_NUCCORE","ACC_NUCCORE","Description_NUCCORE","CreateDate_NUCCORE","Topology_NUCCORE","Completeness_NUCCORE","TaxonID_NUCCORE","Genome_NUCCORE","Length_NUCCORE","Source_NUCCORE","UID_ASSEMBLY","Status_ASSEMBLY","SeqReleaseDate_ASSEMBLY","SubmissionDate_ASSEMBLY","Latest_ASSEMBLY","UID_BIOSAMPLE","ACC_BIOSAMPLE","Location_BIOSAMPLE","Coordinates_BIOSAMPLE","IsolationSource_BIOSAMPLE","Host_BIOSAMPLE","CollectionDate_BIOSAMPLE","Host_DISEASE","SamplType_BIOSAMPLE","taxon_name","taxon_rank","lineage","taxon_species_id","taxon_species_name","taxon_genus_id","taxon_genus_name","taxon_family_id","taxon_family_name","taxon_order_id","taxon_order_name","taxon_class_id","taxon_class_name","taxon_phylum_id","taxon_phylum_name","taxon_superkingdom_id","taxon_superkingdom_name","loc_lat","loc_lng","loc_parsed","GC_NUCCORE","Identical","OldVersion","hits_rMLST","hitscount_rMLST","inclusions","Host_BIOSAMPLE_processed","Host_DISEASE_processed","D1","D2","plasmidfinder","pmlst","relaxase_type(s)","mpf_type"]
        plsdb_tsv = pd.read_csv(plsdb_tsv_file, delimiter= '\t', index_col=False, names=cols, skiprows=1,low_memory=False) 

        combined_df = tophits_mash_df.merge(plsdb_tsv, on='ACC_NUCCORE', how='left')

        combined_df.to_csv(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"), sep="\t", index=False)        

    else:
        cleanup.touch_file(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"))


    return mash_empty
        


# check if the trna file has more than 1 line (not empty)
def is_file_empty(file):
    """
    Determines if file is empty
    :param file: file path
    :return: empty Boolean
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty


def process_mash_tsv_assembly(out_dir, plassembler_db_dir, prefix, input_fasta):
    """
    Process mash output
    :param out_dir: output directory
    :return: mash_empty: boolean whether there was a mash hit
    """

    contig_count = get_contig_count(input_fasta)

    mash_tsv = os.path.join(out_dir,"mash.tsv")

    col_list = ["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"]

    mash_empty = is_file_empty(mash_tsv)

    # instantiate tophits list
    tophits_mash_df = []

    if mash_empty == False:   
        mash_df = pd.read_csv(mash_tsv, delimiter= '\t', index_col=False, names=col_list ) 
        # get list of contigs from unique entries in 
        contigs = mash_df['contig'].unique().tolist()

        contigs = []
        for dna_record in SeqIO.parse(input_fasta, 'fasta'): 
            contigs.append(dna_record.id)

        # instantiate tophits list
        tophits = []

        for contig in contigs:
            # only if a plasmid - not chromosome
            if contig != 'chromosome':
                hit_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True)
                hits = len(hit_df['mash_distance'])
                # add only if there is a hit
                if hits > 0:
                    tmp_df = mash_df.loc[mash_df['contig'] == contig].sort_values('mash_distance').reset_index(drop=True).loc[0]
                    tophits.append([tmp_df.contig, tmp_df.ACC_NUCCORE, tmp_df.mash_distance, tmp_df.mash_pval, tmp_df.mash_matching_hashes])
                # create tophits df
        tophits_mash_df = pd.DataFrame(tophits, columns=["contig", "ACC_NUCCORE", "mash_distance", "mash_pval", "mash_matching_hashes"])

        # read in the plasdb tsv 
        plsdb_tsv_file = os.path.join(plassembler_db_dir, "plsdb.tsv")
        # with open(plsdb_tsv_file, 'rb') as f:
        #     result = chardet.detect(f.readline())
        #     print(result)
        cols = ["UID_NUCCORE","ACC_NUCCORE","Description_NUCCORE","CreateDate_NUCCORE","Topology_NUCCORE","Completeness_NUCCORE","TaxonID_NUCCORE","Genome_NUCCORE","Length_NUCCORE","Source_NUCCORE","UID_ASSEMBLY","Status_ASSEMBLY","SeqReleaseDate_ASSEMBLY","SubmissionDate_ASSEMBLY","Latest_ASSEMBLY","UID_BIOSAMPLE","ACC_BIOSAMPLE","Location_BIOSAMPLE","Coordinates_BIOSAMPLE","IsolationSource_BIOSAMPLE","Host_BIOSAMPLE","CollectionDate_BIOSAMPLE","Host_DISEASE","SamplType_BIOSAMPLE","taxon_name","taxon_rank","lineage","taxon_species_id","taxon_species_name","taxon_genus_id","taxon_genus_name","taxon_family_id","taxon_family_name","taxon_order_id","taxon_order_name","taxon_class_id","taxon_class_name","taxon_phylum_id","taxon_phylum_name","taxon_superkingdom_id","taxon_superkingdom_name","loc_lat","loc_lng","loc_parsed","GC_NUCCORE","Identical","OldVersion","hits_rMLST","hitscount_rMLST","inclusions","Host_BIOSAMPLE_processed","Host_DISEASE_processed","D1","D2","plasmidfinder","pmlst","relaxase_type(s)","mpf_type"]
        plsdb_tsv = pd.read_csv(plsdb_tsv_file, delimiter= '\t', index_col=False, names=cols, skiprows=1,low_memory=False) 

        combined_df = tophits_mash_df.merge(plsdb_tsv, on='ACC_NUCCORE', how='left')

        combined_df.to_csv(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"), sep="\t", index=False)        

    else:
        cleanup.touch_file(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"))


    return mash_empty
        

def rename_contigs_assembly(out_dir, input_fasta, prefix, short_reads):
    """
    Renames the contigs with the new plasmid copy numbers
    :param out_dir: output directory
    :return: 
    """

    depth_df = pd.read_csv(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"), delimiter= '\t', index_col=False, header=0 ) 
    depth_df = depth_df.loc[depth_df['contig'] != 'chromosome'].reset_index(drop=True)
    # get contigs only
    plasmid_fasta = input_fasta
    i = 0
    if short_reads == True:
        with open(os.path.join(out_dir, prefix + "_plasmids.fasta"), 'w') as dna_fa:
            for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
                if "chromosome" not in dna_record.description:
                    if "circular" in dna_record.description:
                        id_updated = dna_record.id + " plasmid_copy_number_short=" + str(depth_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " + "circular=True"
                    else:
                        id_updated = dna_record.id + " plasmid_copy_number_short=" + str(depth_df.plasmid_copy_number_short[i]) + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x "
                    i += 1
                    record = SeqRecord(dna_record.seq, id=id_updated, description = "" )
                    SeqIO.write(record, dna_fa, 'fasta')
    else:
        with open(os.path.join(out_dir, prefix + "_plasmids.fasta"), 'w') as dna_fa:
            for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
                if "chromosome" not in dna_record.description:
                    if "circular" in dna_record.description:
                        id_updated = dna_record.id + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x " + "circular=True"
                    else:
                        id_updated = dna_record.id + "x plasmid_copy_number_long=" + str(depth_df.plasmid_copy_number_long[i]) + "x "
                    i += 1
                    record = SeqRecord(dna_record.seq, id=id_updated, description = "" )
                    SeqIO.write(record, dna_fa, 'fasta')


#######################################################
# add PLSDB hit to copy_number_summary
#####################################################

def update_copy_number_summary_plsdb(out_dir, prefix, mash_empty):
    """
    Updates copy number summary
    :param out_dir: output directory
    :return: 
    """
    depth_df = pd.read_csv(os.path.join(out_dir, prefix + "_copy_number_summary.tsv"), delimiter= '\t', index_col=False, header=0 ) 
    
    if mash_empty == False:

        mash_df = pd.read_csv(os.path.join(out_dir, prefix + "_top_hits_mash_plsdb.tsv"), delimiter= '\t', index_col=False, header=0 )  

        mash_df_reduced =  mash_df[['contig']].copy()
        mash_df_reduced['plsdb_hit'] = 'Yes'

        mash_df_reduced['contig']=mash_df_reduced['contig'].astype(str)
        depth_df['contig']=depth_df['contig'].astype(str)

        combined_df = depth_df.merge(mash_df_reduced, on='contig', how='left')
        combined_df['plsdb_hit'] = combined_df['plsdb_hit'].fillna("No")

        # overwrite the file
        out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
        with open(out_file, 'w') as f:
            combined_df.to_csv(f, sep="\t", index=False, header=True)
    # empty mash - update with 
    else:
        depth_df['plsdb_hit'] = 'No'
            # overwrite the file
        out_file = os.path.join(out_dir, prefix + "_copy_number_summary.tsv")
        with open(out_file, 'w') as f:
            depth_df.to_csv(f, sep="\t", index=False, header=True)



def remove_intermediate_files(out_dir):
    """ removes intermediate files
    :param out_dir:  Output Directory
    :return: 
    """
    sp.run(["rm -rf "+ os.path.join(out_dir,"input.fasta*") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.fastq.gz") ], shell=True)
    sp.run(["rm -rf "+ os.path.join(out_dir,"*.bam") ], shell=True)
    sp.run(["rm", "-rf", os.path.join(out_dir,"params.json") ])
    sp.run(["rm", "-rf", os.path.join(out_dir,"input_tmp.fasta") ])

    # delete mash
    sp.run(["rm", "-rf", os.path.join(out_dir,"mash.tsv") ])
