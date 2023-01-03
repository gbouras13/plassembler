Why Use plassembler?
---------

In long-read first assembled bacterial genomes, small plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see https://f1000research.com/articles/8-2138 https://github.com/rrwick/Trycycler/wiki/Clustering-contigs).

plassembler was created as an automated tool to ensure plasmids assemble correctly without duplicated regions for high-throughput uses - and to provide some useful statistics as well (like copy number). Additionally, it will likely recover small plasmids that long read assemblers like Flye simply miss.

plassembler is primarily intended for high-throughput users with hybrid ONT/short read bacterial sequencing data who want accurate plasmid assemblies.

I developed plassembler because, if conducting Flye assembly alone (and using short reads to polish the assemblies), plasmid assemblies would often be duplicated or triplicated, which would bias downstream analyses involving measures such as plasmid copy numbers.

I also noticed that short read assemblies of plasmids alone would often fail to assemble some plasmids into one circularised contig, due to repetitive elements such as transposons. 

Why Not Use Unicycler?
----

Unicycler is awesome and still probably the best way to assemble plasmids from hybrid sequencing - plassembler uses it! But there are a few reasons to use plassembler instead:

1. Time. plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so is much faster than simply running Unicycler. 
2. Plassembler will output only the plasmids, which may be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://preprints.scielo.org/index.php/scielo/preprint/view/5053) so plassembler can get you only what is necessary from Unicycler.
3. Plassembler will give you summary coverage stats for both long and short reads.
4. Plassembler can be used as fast-ish quality control to check if your short and long reads come from the same sample - if plassembler results in many non-circular contigs (particularly those that, with the help of something like BLAST, map to bacterial chromosomes), it is likely because your read sets do not come from the same isolate! 
5. As of v 0.1.4, you will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. Additionally, if there are many contigs with no PLSDB hits, this can help you determine that your long and short reads do not come from the same isolate, or that there may be some other biological phenomenon in your "isolate" (for example, perhaps there is some other non-plasmid mobile genetic element present in low abundance). For less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.