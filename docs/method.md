Method
-------

1. Long reads are filtered using nanofilt (https://github.com/wdecoster/nanofilt) .
2. Long-read assembly is conducted with Flye (https://github.com/fenderglass/Flye).
3. If the resulting assembly has more than 1 contig, the largest contig is checked. If it is over 90% of the length of the provided chromosome size, or is circular, then it is identified as the chromosome and extracted. All other contigs are extracted as putative plasmid extra-chromosomal contigs.
4. Short reads are filtered using fastp (https://github.com/OpenGene/fastp).
5. Long reads are mapped to the extra-chromosomal contigs using minimap2 (https://github.com/lh3/minimap2#uguide), and short reads are mapped using bwa (https://github.com/lh3/bwa).
6. Long reads are mapped to the chromosome using minimap2 and short reads are mapped using bwa. This is done to identify reads that do not map to the chromosome (for any plasmids that Flye may have missed assembling).
7. All reads that map to the extra-chromosomal contigs and all reads that do not map the chromosome are extracted, combined and de-duplicated.
8. These are assembled using the hybrid assembler Unicycler to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculating using a modified version of code found in (https://github.com/rrwick/Small-plasmid-Nanopore, https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2). 
10. Plasmid copy number is calculated by dividing the plasmid depth by the chromosome depth.


Other Features (Work in Progress)

1. All reads that map to both the chromosome and plasmid are extracted and assembled (short read only assembly).
