plassembler creates a number of output files

The main outputs are `_plasmids.fasta` and `_plasmids.gfa`, which contain the Unicycler assembly outputs in fasta and graph format. The `.gfa` file can be visualised is programs such as Bandage (https://github.com/rrwick/Bandage/).

plassmbler also outputs `copy_number_summary.tsv`, which gives plasmid copy number statistics for both short and long read sets.

**Note**

* If plassembler fails to find any plasmids, these files will exist, but will be empty (to ensure plassembler can be easily integrated into workflow managers like Snakemake).

plassembler will contain all flye output in the `flye_output` directory and all unicycler output in the `unicycler_output` directory.

Additionally, as a work in progress, there will also be a directory `unicycler_plasmid_chromosome_map_output` that contains the short read only assemblies from all short reads that map to both the chromosome and plasmids (https://www.pnas.org/doi/10.1073/pnas.2008731118).
