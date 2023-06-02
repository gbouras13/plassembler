# Plassembler Simulated Reads Benchmarking Results 

All benchmarking was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS. 

The following Tables show the summary statistics of `plassembler` (with Flye and Raven) and Unicycler run from read sets simulated from 20 assemblies from Wick, Judd, Wyres et al 2021 [here](https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/method.md), Houtak et al (C222), [here](https://doi.org/10.1101/2023.03.28.534496) and [here](https://github.com/gbouras13/CRS_Saureus_Evolutionary_Landscape), Mathers et al (CAV1217) [here](https://doi.org/10.1128/AAC.01823-16), and for De Maio et al [here](https://doi.org/10.1099/mgen.0.000294).

Accuracy was assessed using QUAST comparing plasmid assemblies against the ground truth.

You can find the full pipeline used to generate these results [here](https://github.com/gbouras13/plassembler_simulation_benchmarking).

Time and Memory Usage
===============

As can be seen below, Plassembler was consistently faster than Unicycler by 3-10x in terms of wall clock time, with the largest increase coming single-threaded. Plassembler with Raven was also faster than with Flye.

| Threads | Program                | Median Wall Clock Time (s) | Min Wall Clock Time (s) | Max Wall Clock Time (s) | Median Max Memory (MB) | Min Max Memory (MB) | Max Max Memory (MB) |
| ------- | ---------------------- | -------------------------- | ----------------------- | ----------------------- | ---------------------- | ------------------- | ------------------- |
| 1       | Plassembler with Flye  | 7012                       | 1926                    | 28103                   | 3039                   | 2442                | 5275                |
| 1       | Plassembler with Raven | 3893                       | 815                     | 22106                   | 1764                   | 1307                | 2464                |
| 1       | Unicycler              | 32619                      | 13183                   | 66880                   | 2804                   | 1404                | 3784                |
| 8       | Plassembler with Flye  | 1563                       | 531                     | 4852                    | 4842                   | 2274                | 7712                |
| 8       | Plassembler with Raven | 700                        | 132                     | 3051                    | 2579                   | 1587                | 7891                |
| 8       | Unicycler              | 4411                       | 2139                    | 10003                   | 6568                   | 5610                | 6826                |
| 16      | Plassembler with Flye  | 1019                       | 497                     | 2675                    | 5892                   | 3832                | 15086               |
| 16      | Plassembler with Raven | 430                        | 114                     | 1749                    | 3517                   | 2360                | 15195               |
| 16      | Unicycler              | 2554                       | 1347                    | 5098                    | 12967                  | 6509                | 13549               |

Assembly Accuracy
==================

All results below were taken from each program run with 8 threads.

As can be seen below, `plassembler` missed only 1 plasmid (2370 bp plasmid in _Enterobacter kobei_ MSB1 1B) while Unicycler missed 7, of which were small (under 10kbp). These were _S. aureus_ C222 (2473 bp), _C. koseri_ MINF 9D  (9294 bp), _K. oxytoca_ MSB1 2C (4574 bp), _K. variicola_ INF345  (5783 bp), _E. cloacae_ RBHSTW-00059 (2495 bp) and _K. pneumoniae_ RBHSTW-00128 (3980 bp).

`plassembler` also had more fragmented assemblies than Unicycler (3 and 5 v 1). All assemblers fragmented the linear plasmid in K. variicola_ INF345, while Plassembler with Flye also fragmented a 3764 bp plasmid in _C. gillenii_ RBHSTW-00142, and a 4841 bp plasmid in _K. oxytoca_ RBHSTW-00167.

Brackets indicate (minimum and maximum). Indels and Mismatches are expressed as median, and per 100 kbp from QUAST.


| Program                | Missed Plasmids | Incomplete/Fragmented Plasmids | Misassemblies | Genome Fraction                             | Indels      | Mismatches      |
| ---------------------- | --------------- | ------------------------------ | ------------- | ------------------------------------------- | ----------- | --------------- |
| Plassembler with Flye  | 1               | 3                              | 0             | 99.78 (mean), 99.97 (median),  (98.16, 100) | 0 (0, 1.37) | 0.91 (0, 11.12) |
| Plassembler with Raven | 1               | 5                              | 0             | 99.04 (mean), 99.96 (median),  (87.86, 100) | 0 (0, 1.37) | 1.04 (0, 11.12) |
| Unicycler              | 7               | 1                              | 1             | 93.81 (mean), 99.88 (median),  (0, 100)     | 0 (0, 1.54) | 0.88 (0, 7.28)  |