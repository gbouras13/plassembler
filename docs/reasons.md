# Why Does Plassembler Exist?

In long-read first assembled bacterial genomes, small plasmids are often difficult to assemble correctly with long read assemblers such as Flye. They often have circularisation issues and can be duplicated or missed (see [this](https://doi.org/10.1371/journal.pcbi.1010905), [this](https://f1000research.com/articles/8-2138) and [this](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs) ).

Plassembler was therefore created as an automated tool to ensure plasmids assemble correctly without duplicated regions for high-throughput uses - and to provide some useful statistics as well (such as estimate plasmid copy numbers for both long and short read sets). Plassembler will likely also recover small plasmids that long read assemblers like Flye  miss.

As it turns out (though this wasn't a motivation for making it!), Plassembler also recovers more small plasmids than the existing gold standard tool Unicycler. I think this is because it throws away chromosomal reads, similar to subsampling short reads sets which can improve assembly quality (see [this](https://doi.org/10.1093/bioinformatics/btv311) and [this](https://doi.org/10.1371/journal.pone.0060204) and [this](https://academic.oup.com/bioinformatics/article/27/4/479/198367) - thanks Michael Hall for the references!). In other words, there are more plasmid reads a proportion of the overall read set, so a higher chance of recovering smaller plasmids.

You can see this increase in accuracy and speed in the benchmarking results for [simulated](docs/benchmarking_results_sim.md) and [real](docs/benchmarking_results_real.md) datasets.

Plassembler also uses [mash](https://github.com/marbl/Mash) as a quick way to determine whether each assembled contig has any similar hits in [PLSDB](https://doi.org/10.1093/nar/gkab1111). 

Additionally, due to its mapping approach, Plassembler can also be used as a quality control tool for checking whether your long and short read sets come from the same isolate. This may be particularly useful if your read sets come from different extractions, or you have multiplexed many samples (& want to avoid mislabelling).  

# Why Not Just Use Unicycler?

Unicycler is awesome and still a good way to assemble plasmids from hybrid sequencing - plassembler uses it! But there are a few reasons to use plassembler instead:

1. Time. Plassember throws away all the chromosomal reads (i.e. most of them) before running Unicycler, so it is much faster (wall clock 3-10x faster generally). 
2. Accuracy. Benchmarking has shown Plassembler is better than Unicycler in terms of recovering small plasmids.
3. Plassembler will output only the likely plasmids, and can more easily be integrated into pipelines. You shouldn't be assembling the chromosome using Unicycler [anymore](https://doi.org/10.1371/journal.pcbi.1010905) so Plassembler can get you only what is necessary from Unicycler.
4. Plassembler will give you summary depth and copy number stats for both long and short reads.
5. Plassembler can be used as a quality control to check if your short and long reads come from the same sample - if plassembler results in many non-circular contigs (particularly those that have no hits in PLSDB), it is likely because your read sets do not come from the same isolate! See [Quality Control](#quality-control).
6. You will get information whether each assembled contig has a similar entry in [PLSDB](https://doi.org/10.1093/nar/gkab1111). Especially for common pathogen species that are well represented in databases, this will likely tell you specifically what plasmid you have in your sample. 
* Note: Especially for less commonly sequenced species, I would not suggest that that absence of a PLSDB hit is necessary meaningful, especially for circular contigs - those would likely be novel plasmids uncaptured by PLSDB.

