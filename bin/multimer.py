import os
import sys
import subprocess as sp
import log

def minimap2_unicycler_vs_flye_plasmids(out_dir, prefix, logger):

    # flye plasmids

    flye_plasmids = os.path.join(out_dir,"flye_renamed.fasta")

    # will be the same as output
    unicycler_plasmids = os.path.join(out_dir, prefix + "_plasmids.fasta")

    paf = os.path.join(out_dir, "mapping.paf")

    f = open(paf, "w")

    try:
        minimap = sp.Popen(["minimap2",   unicycler_plasmids, flye_plasmids ], stdout=f, stderr=sp.PIPE) 
        log.write_to_log(minimap.stderr, logger)
    except:
        sys.exit("Error with minimap2\n")  




