Method
-------

1. Long reads are filtered using [chopper](https://github.com/wdecoster/chopper) .
2. Long reads are subsampled to 30x to reduce runtime (subsampling depth be changed with `-s` or skipped with `--no_subsample`).
3. Long-read only assembly is conducted with [Flye](https://github.com/fenderglass/Flye).
4. If the resulting assembly is checked. Contigs bigger than the provided chromosome size `-c`, are identified as chromosomal and extracted. Any other contigs are extracted as putative plasmid contigs, if Flye assembled any. If no contigs were larger than `-c`, plassembler will exit - you probably need to get some more long reads to complete your assembly (or check `-c` wasn't too big).
5. Short reads are filtered using [fastp](https://github.com/OpenGene/fastp).
6. Long and short reads are mapped to a reference containing the chromosomal contigs plus putative plasmid contigs using [minimap2](https://github.com/lh3/minimap2#uguide).
7. All reads that map to the putative plasmid contigs and all reads that are unmapped are extracted and combined.
8. These reads are assembled using the hybrid assembler [Unicycler](https://github.com/rrwick/Unicycler) to generate final plasmid contigs.
9. Average read coverage depth for each plasmid is calculated using a modified version of code found [here](https://github.com/rrwick/Small-plasmid-Nanopore). See also this [paper](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2).
10. Plasmid copy number is calculated by dividing the plasmid read depth by the chromosome read depth.
11. All plasmid contigs are compared against [PLSDB](https://doi.org/10.1093/nar/gkab1111) using [mash](https://github.com/marbl/Mash) with a cutoff maximum mash distance of 0.1.

