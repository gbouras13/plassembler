from Bio import SeqIO
import pandas as pd
import os
import subprocess as sp
import log
import sys
import cleanup


def mash_sketch(out_dir, logger):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    plasmid_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")

    try:
        mash_sketch = sp.Popen(["mash", "sketch",  plasmid_fasta, "-i" ], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stdout, logger)
    except:
        sys.exit("Error with mash sketch.\n")  


def run_mash(out_dir, plassembler_db_dir, logger):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param plassembler_db_dir: plassembler db directory
    :param logger: logger
    :return:
    """

    plsdb_sketch = os.path.join(plassembler_db_dir, "plsdb.msh")

    plasmid_sketch = os.path.join(out_dir,"unicycler_output", "assembly.fasta.msh")

    mash_tsv = os.path.join(out_dir,"mash.tsv")

    outFile = open(mash_tsv, "w")
    try:
        mash_sketch = sp.Popen(["mash", "dist",  plasmid_sketch, plsdb_sketch, "-v", "0.1", "-d", "0.1", "-i" ], stdout=outFile, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stderr, logger)
    except:
        sys.exit("Error with mash dist.\n")  

#mash dist test_mash/unicycler_output/assembly.fasta.msh ../plsdb_110222_plassembler_v0.1.4_databases/plsdb.msh -v 0.1 -d 0.1 -i

def get_contig_count(out_dir):
    """
    Process mash output
    :param out_dir: output directory
    :return: i: int contig_count
    """
    plasmid_fasta = os.path.join(out_dir,"unicycler_output", "assembly.fasta")
    i = 0
    for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
            i += 1
    return i


def process_mash_tsv(out_dir, plassembler_db_dir, prefix):
    """
    Process mash output
    :param out_dir: output directory
    :return: mash_empty: boolean whether there was a mash hit
    """

    contig_count = get_contig_count(out_dir)

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