Plassembler is designed for automated assembly of plasmids in haploid bacterial genomes that have been hybrid sequenced with  long read (Nanopore, but also with Pacbio if `--pacbio_model` is specified) & paired-end short read sequencing.

Plassembler uses Flye (Kolmogorov et al 2019) or optionally Raven (Vaser and Šikić 2021) to create an initial long-read only assembly. Following this, minimap2 (Li 2018) is used to map long and short read sets.

Reads that map to the non-chromosome contigs (i.e. the putative plasmids), and reads that are unmapped  (to recover possible plasmids Flye has missed) are extracted and combined.

These reads are then assembled using Unicycler (Wick et al 2017). The resulting output should ideally include accurate circularised plasmid contigs plassembler calculates plasmid copy number for both long and short read sets.

Finally, Plassembler uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 
