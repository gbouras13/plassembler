Outputs
-------
plassembler will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in FASTA format, and a `_plasmids.gfa` file, which will contain the assembly graph from Unicycler that can be visualised in [Bandage](https://github.com/rrwick/Bandage). 

plassembler also outputs a `copy_number_summary.tsv` file, which gives the estimated copy number for each plasmid, for both short reads and long reads (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) for more details about plasmid copy numbers) and also a `_top_hits_mash_plsdb.tsv` file, which gives each contig's top hit by mash distance in the PLSDB (if there is a hit), along with all its supporting information. 

If there is no PLSDB hit for a contig, there will be no entry in the `_top_hits_mash_plsdb.tsv` file.

If plassembler fails to find any plasmids, these files will still exist, but will be empty (to ensure plassembler can be easily integrated into workflow managers like Snakemake).

plassembler will also output a log file, a `flye_output` directory, which contains the output from Flye (it may be useful to decide whether you need more sequencing reads, or some strange assembly artifact occured) and a `unicycler_output` directory containing the output from Unicycler.

**Other Benefits**
If plassembler results in many non-circular contigs (particularly those that, with the help of something like BLAST, map to bacterial chromosomes), it is likely because your read sets do not come from the same isolate! Accordingly, plassembler can also be used for QC purposes in this way.

**What happens if plassembler fails to find a plasmid?**

* There are two reasons why plassembler will fail to find a plasmid:

1. Where Flye assembles a complete chromosome, but fails to find any plasmids. Most of the time, this simply means that there are no plasmids in the your bacterial isolate. 
2. Where there is insufficient coverage for Flye to assemble a complete circular chromosome. In these cases, it is recommended that you either do more long read sequencing so that you can assemble a complete circular chromosome, or check that your -c parameter is correcting set.