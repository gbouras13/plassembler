# Plassembler & Assembler Non Determinism 

All following data was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS. 

In the course of bechmarking Plassembler, I ran into some interesting problems regarding non-determinism.

It was flagged in review that Plassembler was non-deterministic, which I didn't actually know about (all the testing I'd done had been with the same amount of threads - a good point for anyone testing their tools in future!)

Previously (perhaps naively), I had assumed assemblers ( [SPAdes](https://github.com/ablab/spades), [Flye](https://github.com/fenderglass/Flye) and [Unicycler](https://github.com/rrwick/Unicycler) (which uses SPAdes)) were deterministic. But I quickly found out they were not.

It is know by the developers that [SPAdes](https://github.com/ablab/spades/issues/111) and [Flye](https://github.com/fenderglass/Flye/pull/546) were non-deterministc - and Flye even has a `--deterministic` option.

Interestingly, I found out that even Flye using `--deterministic` isn't completely deterministic between threads, so there is some deeper stochastisity at play.

The following Tables show the results of Flye assemblies of real read sets for Wick et al [here](https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/method.md), Houtak et al (C222), [here](https://doi.org/10.1101/2023.03.28.534496) and [here](https://github.com/gbouras13/CRS_Saureus_Evolutionary_Landscape) and Mathers et al (CAV1217) [here](https://doi.org/10.1128/AAC.01823-16).

I didn't do SPAdes because, as you can see in the main benchmarking analysis [here](), there is no difference in Unicycler between threads, so I didn't bother investigate further - indicating SPAdes is pretty close to deterministic at least.

You can find the full pipeline used to generate these results [here](https://github.com/gbouras13/plassembler_simulation_benchmarking).

First I have included the Flye assemblies using a version of Plassembler (v 1.0.0 as per the release and bioconda) with Flye v2.9.2 without `--deterministic`. The exact same readset for all 3 runs (1, 8 and 16 threads) are used as input (as chopper and rasusa are deterministic).

The exact command used within Plassembler is:

```
flye --nano-hq <input fastq> --out-dir <output directory> --threads <thread count>
```

Flye without `--deterministic`
=================================
