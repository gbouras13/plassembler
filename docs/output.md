Outputs
-------
`plassembler` will output a `_plasmids.fasta` file, which will contain the assembled plasmid sequence(s) in FASTA format (including long and short read copy numbers in the header), and a `_plasmids.gfa` file, which will contain the assembly graph from Unicycler that can be visualised in [Bandage](https://github.com/rrwick/Bandage). 

`plassembler` also outputs a `_summary.tsv` file, which gives the estimated copy number for each plasmid, for both short reads and long reads (see this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2) for more details about plasmid copy numbers) and also gives each contig's top hit by mash distance in the PLSDB (if there is a hit), along with all its supporting information. 

If `plassembler` fails to assemble any plasmids at all in `_plasmids.fasta`, all these files will still exist, but will be empty (to ensure plassembler can be easily integrated into workflow managers like Snakemake).

`plassembler` will also output a log file, a `flye_output` directory, which contains the output from Flye (it may be useful to decide whether you need more sequencing reads, or some strange assembly artifact occured) and a `unicycler_output` directory containing the output from Unicycler. If `--use_raven` is specified, a `raven_output` director will exist instead of `flye_output`.



Other Outputs
------------

If you use the `--keep_fastqs` flag, `plassembler` will keep FASTQ files containing all reads that went into the Unicycler assembly. These will be kept in the `plasmid_fastqs` directory. It will also keep all long reads that map to more than one contig as `multimap_long.fastq`. 

If you use `--keep_chromosome`, Plassembler will keep the unpolished chromosome(s) as `chromosome.fasta`.
	
**What happens if `plassembler` fails to find a plasmid?**

* There are two reasons why `plassembler` will fail to find a plasmid:

1. Where Flye assembles a complete chromosome, and then Plassembler fails to find any plasmids using the short reads. Most of the time, this simply means that there are no plasmids in the your bacterial isolate. Plassembler will still make output files (to ensure you can easily use it with workflow managers like Snakemake), but these will be empty.
2. Where there is insufficient coverage for Flye to assemble a complete circular chromosome. In these cases, it is recommended that you either do more long read sequencing so that you can assemble a complete circular chromosome, or check that your -c parameter is correcting set.