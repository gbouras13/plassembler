plassembler is a pipeline to assemble accurate extra-chromosomal plasmids from hybrid ONT/short read bacterial sequencing data.

plassembler uses Flye (Kolmogorov et al 2019) to create an initial long-read only assembly. Following this, minimap2 (Li 2018) and BWA are used to map long and short read sets to the Flye assembly.

Reads that map to the non-chromosome contigs (i.e. the putative plasmids), and reads that do not map to the chromosome (to recover possible plasmids Flye has missed) are then extracted, combined and de-duplicated.

These reads are then assembled using Unicycler (Wick et al 2017). The resulting output should ideally include accurate circularised plasmid contigs.

plassembler also calculated plasmid copy number for both long and short read sets.
