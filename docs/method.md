Method
-------

1. Long reads are filtered using [nanofilt](https://github.com/wdecoster/nanofilt) .
2. Long-read only assembly is conducted with [Flye](https://github.com/fenderglass/Flye).
3. If the resulting assembly is checked. If the largest contig is over 90% of the length of the provided chromosome size -c, then it is identified as the chromosome and extracted. Any other contigs are extracted as putative plasmid contigs, if Flye assembled any. Otherwise, plassembler will exit - you probably need to get some more long reads to complete your assembly (or check -c).
4. Short reads are filtered using [fastp](https://github.com/OpenGene/fastp).
5. Long reads are mapped to any putative plasmid contigs using [minimap2](https://github.com/lh3/minimap2#uguide), and short reads are mapped using [bwa](https://github.com/lh3/bwa).
6. Long reads are mapped to the chromosome using minimap2 and short reads are mapped using bwa. This is done to identify reads that do not map to the chromosome (for any plasmids that Flye may have missed assembling).
7. All reads that map to the putative plasmid contigs and all reads that do not map the chromosome are extracted, combined and de-duplicated using [seqkit](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163962).
8. The de-deduplicated read sets are assembled using the hybrid assembler [Unicycler](https://github.com/rrwick/Unicycler) to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculated using a modified version of code found [here](https://github.com/rrwick/Small-plasmid-Nanopore). See also this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).
10. Plasmid copy number is calculated by dividing the plasmid read depth by the chromosome read depth.



Other Features (Work in Progress)

1. All reads that map to both the chromosome and plasmid are extracted and assembled (short read only assembly).
