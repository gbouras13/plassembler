plassembler is tool designed for automated assembly of plasmids in haploid bacterial genomes that have been hybrid sequenced with Oxford Nanopore Technologies long read (Ideally R9.4.1, also should work with R 10.4.1) & paired-end short read sequencing.

plassembler uses Flye (Kolmogorov et al 2019) to create an initial long-read only assembly. Following this, minimap2 (Li 2018) and BWA are used to map long and short read sets to the Flye assembly.

Reads that map to the non-chromosome contigs (i.e. the putative plasmids), and reads that do not map to the chromosome (to recover possible plasmids Flye has missed) are extracted, combined and de-duplicated.

These reads are then assembled using Unicycler (Wick et al 2017). The resulting output should ideally include accurate circularised plasmid contigs.

plassembler also calculates plasmid copy number for both long and short read sets.
