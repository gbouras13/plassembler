# Plassembler Real Reads Benchmarking Results 

All benchmarking was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS. 

The following Tables show the summary statistics of `plassembler` (with Flye and Raven) and Unicycler run from real read sets simulated from 6 samples from Wick, Judd, Wyres et al 2021. These were compared against the independent Trycyler derived and manually curated ground truth you can find [here](https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/method.md). 

You can find the full pipeline used to generate these results [here](https://github.com/gbouras13/plassembler_simulation_benchmarking).

You can find the full details on how the real read sets were obtained and subsampled [here](fastqs.md).

Time and Memory Usage
===============

As can be seen below, `plassembler` was consistently faster than Unicycler in terms of wall clock time, with the largest increase coming single-threaded. `plassembler` with Raven was also faster than with Flye, but less so than in the [simulated benchmarking](benchmarking_results_sim.md).


| Threads | Program                | Median Wall Clock Time (s) | Min Wall Clock Time (s) | Max Wall Clock Time (s) | Median Max Memory (MB) | Min Max Memory (MB) | Max Max Memory (MB) |
| ------- | ---------------------- | -------------------------- | ----------------------- | ----------------------- | ---------------------- | ------------------- | ------------------- |
| 1       | Plassembler with Flye  | 7063                       | 4559                    | 7865                    | 9921                   | 9116                | 10024               |
| 1       | Plassembler with Raven | 5918                       | 3430                    | 7028                    | 2623                   | 2214                | 3039                |
| 1       | Unicycler              | 48325                      | 37282                   | 58823                   | 4671                   | 3583                | 5832                |
| 8       | Plassembler with Flye  | 1493                       | 1250                    | 1688                    | 10172                  | 9710                | 10897               |
| 8       | Plassembler with Raven | 1126                       | 709                     | 1356                    | 6626                   | 2872                | 8374                |
| 8       | Unicycler              | 7500                       | 4509                    | 9659                    | 7535                   | 7003                | 8128                |
| 16      | Plassembler with Flye  | 975                        | 779                     | 1163                    | 11803                  | 10344               | 13467               |
| 16      | Plassembler with Raven | 699                        | 397                     | 840                     | 6574                   | 3343                | 14784               |
| 16      | Unicycler              | 3876                       | 2944                    | 5036                    | 14041                  | 13647               | 14062               |

Assembly Accuracy
==================

All results below were taken from each program run with 8 threads.

`plassembler` and Unicycler had identical genome fractions and low indel and mismatch rates. Similarly to the [simulated benchmarking](benchmarking_results_sim.md), `plassembler` recovered 2 additional small plasmids missed by Unicycler of lengths 1934 bp  _(_K. variicola_ INF345) and 10697 bp (_K. oxytoca_ MSB1 2C). 

It turned out that the 10697 bp plasmid recovered in Klebsiella oxytoca MSB1 2C was not recovered using the Trycycler based long-read first assembly method in any sample by Wick, Judd, Wyers, et al., 2021, and annotation with Bakta v1.7.0 (Schwengers et al., 2021) revealed that this plasmid contains a Type III toxin-antitoxin system (Supplementary Table 9) along with other plasmid replication genes (repA etc) - it is possibly a real plasmid missed in that study!

Brackets indicate (minimum and maximum). Indels and Mismatches are expressed as median, and per 100 kbp from QUAST.

| Program                | Additional Plasmids Recovered | Total Missed Plasmids | Total Fragmented Plasmids | Total Misassemblies | Genome Fraction                           | Indels (per 100kbp)    | Mismatches (per 100kbp) |
| ---------------------- | ----------------------------- | --------------------- | ------------------------- | ------------------- | ----------------------------------------- | ---------------------- | ----------------------- |
| Plassembler with Flye  | 2                             | 2                     | 1                         | 0                   | 97.02 (mean), 100 (median),  (82.14, 100) | 0.46 (0 min, 8.6 max)  | 0.2 (0 min, 135.23 max) |
| Plassembler with Raven | 2                             | 2                     | 1                         | 0                   | 97.02 (mean), 100 (median),  (82.14, 100) | 0.64 (0 min, 1.67 max) | 0.58 (0 min, 7.94 max)  |
| Unicycler              | 0                             | 2                     | 1                         | 0                   | 97.02 (mean), 100 (median),  (82.14, 100) | 0.18 (0 min, 1.67 max) | 0.18 (0 min, 2.65 max)  |
