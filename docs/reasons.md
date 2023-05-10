Why Does Plassembler Exist?
---------

In long-read first assembled bacterial genomes, small plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see [this](https://doi.org/10.1371/journal.pcbi.1010905) and [this](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs) ).

Plassembler was therefore created as an automated tool to ensure plasmids assemble correctly without duplicated regions for high-throughput uses - and to provide some useful statistics as well (such as estimate plasmid copy numbers for both long and short read sets). Plassembler will likely also recover small plasmids that long read assemblers like Flye simply miss.

As it turns out (though this wasn't a motivation for making it!), Plassembler also assembles plasmids more accurately than Unicycler. I think this is because subsampling short reads sets can improve assembly quality (see [this](https://doi.org/10.1093/bioinformatics/btv311) and [this](https://doi.org/10.1371/journal.pone.0060204) and [this](https://academic.oup.com/bioinformatics/article/27/4/479/198367) - thanks Michael Hall for the references!), so throwing away chromosomal reads probably has a similar effect.

Plassembler also uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 

Additionally, due to its mapping approach, Plassembler can also be used as a quality control tool for checking whether your long and short read sets come from the same isolate. This may be particularly useful if your read sets come from different extractions, or you have multiplexed many samples (& want to avoid mislabelling).  

Why Not Just Use Unicycler?
---------

Unicycler is awesome and still a good way to assemble plasmids from hybrid sequencing - plassembler uses it! But there are a few reasons to use plassembler instead:

1. Time. Plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so it is much faster (4-20x, and will be higher if you have lots of long reads). 
2. Accuracy. Benchmarking has shown Plassembler is more accurate than Unicycler in terms of indels and mismatches. I honestly wasn't even aiming for this when I wrote Plassembler, just a big speed-up, but it is a nice result of course!
3. Plassembler will output only the likely plasmids, and can more easily be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://doi.org/10.1371/journal.pcbi.1010905) so plassembler can get you only what is necessary from Unicycler.
4. Plassembler will give you summary depth and copy number stats for both long and short reads.
5. Plassembler can be used as a quality control to check if your short and long reads come from the same sample - if plassembler results in many non-circular contigs (particularly those that have no hits in PLSDB), it is likely because your read sets do not come from the same isolate! 
6. You will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. 
* Note: Especially for less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.
