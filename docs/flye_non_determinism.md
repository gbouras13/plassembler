# Plassembler & Assembler Non Determinism 

All following data was conducted on a Intel® Core™ i7-10700K CPU @ 3.80GHz on a machine running Ubuntu 20.04.6 LTS. 

In the course of bechmarking `plassembler`, I ran into some interesting problems regarding non-determinism.

It was flagged in manuscript review that `plassembler` was non-deterministic, which I didn't actually know about. All my testing I'd done had been with the same amount of threads - so this will be something I take into future benchmarking analyses.

Previously (perhaps naively), I had assumed assemblers ( [SPAdes](https://github.com/ablab/spades), [Flye](https://github.com/fenderglass/Flye) and [Unicycler](https://github.com/rrwick/Unicycler) (which uses SPAdes)) were deterministic. But after some digging, I quickly found out they were not.

It turns out essentially every assembler is not deterministic between threads due to multi-threading:

* [Canu](https://github.com/marbl/canu/issues/889)
* [megahit](https://github.com/voutcn/megahit/issues/170)
* [Shasta](https://github.com/chanzuckerberg/shasta/issues/296)
* [Flye](https://github.com/fenderglass/Flye/issues/509 ) and [here](https://github.com/fenderglass/Flye/issues/546)
* [SPAdes](https://github.com/ablab/spades/issues/111)

In practice from my benchmarking however, I found that [Unicycler](https://github.com/rrwick/Unicycler) (which uses SPAdes), is "practically" deterministic between thread counts (i.e. it always gave the same output). So I did not investigate any further.

In Flye, the developers recently added a `--deterministic` option.

Interestingly, I found out that even Flye using `--deterministic` isn't deterministic between threads either, so there is some deeper stochastisity at play in Flye. So I decided to investigate this.

It turns out that `--deterministic` doesn't really make a difference to Flye's determinism - so I would probably not use it as it slows Flye down!

The following Tables show the results of Flye assemblies of real read sets for Wick et al [here](https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/method.md) and Houtak et al (C222, subsampled to 30x and 60x), [here](https://doi.org/10.1101/2023.03.28.534496).

I have included the Flye v2.9.2 assembly lengths with and without without `--deterministic` run using Plassembler. The exact same readset for all 3 runs (1, 8 and 16 threads) are used as input (as chopper is deterministic).

Note that this is with Plassembler v1.0.0, so if you dig into the details on Zenodo the log files will look different.

The exact command used within `plassembler` is:

```
flye --nano-hq <input fastq> --out-dir <output directory> --threads <thread count>
```

Flye without `--deterministic`
=================================


| Isolate                          | Threads                    | 1                                                      | 8 | 16 |
| -------------------------------- | -------------------------- | ------------------------------------------------------ | ------------------------------------------------------ | ------------------------------------------------------ |
| A. baumannii J9                  | Flye Contig Lengths (bp)   | 3797823; 145027; 12145                                 | 3797806; 145027; 12086 | 3797803; 145026; 12144 |
|                                  | Plassembler Contig Lengths | 145071; 6078                                           | 145069; 6078 | 145071; 6078 |
|                                  |                            |                                                        |  |  |
| S. aureus C222 30x               | Flye Contig Lengths (bp)   | 2788497; 2472                                          | 2788497; 2472 | 2788499; 2472 |
|                                  | Plassembler Contig Lengths | 2473                                                   | 2473 | 2473 |
|                                  |                            |                                                        |  |  |
| S. aureus C222 60x               | Flye Contig Lengths (bp)   | 2788510; 2469                                          | 2788510; 2471 | 2788509; 1895 |
|                                  | Plassembler Contig Lengths | 2473                                                   | 2473 | 2473 |
|                                  |                            |                                                        |  |  |
| C. koseri MINF 9D                | Flye Contig Lengths (bp)   | 4757420; 64915; 18570                                  | 4757424; 64916; 9283 | 4757426; 64913; 18593 |
|                                  | Plassembler Contig Lengths | 64962; 9294                                            | 64962; 9294 | 64962; 9294 |
|                                  |                            |                                                        |  |  |
| E. kobei MSB1 1B                 | Flye Contig Lengths (bp)   | 4836290; 132437; 108265; 14925; 7173                   | 4836272; 244830; 4650 | 4836272; 136413; 108366; 4650 |
|                                  | Plassembler Contig Lengths | 136482; 108411; 4665; 3715; 2370                       | 136482; 108411; 4665; 3715; 2370 | 136482; 108411; 4665; 3715; 2370 |
|                                  |                            |                                                        |  |  |
| Haemophilus sp002998595 M1C132 1 | Flye Contig Lengths (bp)   | 2051487; 39372; 21394; 19915                           | 2051476; 39371; 27888; 21412; 13546; 3160 | 2051483; 39373; 21404; 19929 |
|                                  | Plassembler Contig Lengths | 39345; 10719; 9975                                     | 39345; 10719; 9975 | 39345; 10719; 9975 |
|                                  |                            |                                                        |  |  |
| K. oxytoca MSB1 2C               | Flye Contig Lengths (bp)   | 5802199; 118088; 58430; 9135                           | 5802339; 118091; 58430; 9141 | 5802182; 118088; 58430; 4402 |
|                                  | Plassembler Contig Lengths | 118161; 58382; 10697; 9975; 4574                       | 118161; 58382; 10697; 9975; 4574 | 118161; 58382; 10697; 9975; 4574 |
|                                  |                            |                                                        |  |  |
| K. variicola INF345              | Flye Contig Lengths (bp)   | 5415051; 250848; 243518; 31785; 11548                  | 5415038; 250870; 243809; 31785; 11547; 1947 | 5415038; 250848; 243515; 30857; 11550 |
|                                  | Plassembler Contig Lengths | 250272; 243534; 30412 + 663 (linear); 5738; 3514; 1934 | 250914; 243534; 30412 + 663 (linear); 5783; 3514; 1934 | 250914; 243534; 30412 + 663 (linear); 5783; 3514; 1934 |

Flye with `--deterministic`
=================================


| Isolate                          | Threads                    | 1                                                      | 8                                                      | 16                                                     |
| -------------------------------- | -------------------------- | ------------------------------------------------------ | ------------------------------------------------------ | ------------------------------------------------------ |
| A. baumannii J9                  | Flye Contig Lengths (bp)   | 3797822; 144252; 11505                                 | 3797808; 145027; 18231                                 | 3797822; 144252; 18053                                 |
|                                  | Plassembler Contig Lengths | 145069; 6078                                           | 145069; 6078                                           | 145071; 6078                                           |
|                                  |                            |                                                        |                                                        |                                                        |
| S. aureus C222 30x               | Flye Contig Lengths (bp)   | 2788497; 2472                                          | 2788496; 2471                                          | 2788497; 2471                                          |
|                                  | Plassembler Contig Lengths | 2473                                                   | 2473                                                   | 2473                                                   |
|                                  |                            |                                                        |                                                        |                                                        |
| S. aureus C222 60x               | Flye Contig Lengths (bp)   | 2788514; 2471                                          | 2788510; 2471                                          | 2788509; 2469                                          |
|                                  | Plassembler Contig Lengths | 2473                                                   | 2473                                                   | 2473                                                   |
|                                  |                            |                                                        |                                                        |                                                        |
| C. koseri MINF 9D                | Flye Contig Lengths (bp)   | 4757416; 64916; 18569                                  | 4757419; 64916; 17890                                  | 4757436; 64915; 18569                                  |
|                                  | Plassembler Contig Lengths | 64962; 9294                                            | 64962; 9294                                            | 64962; 9294                                            |
|                                  |                            |                                                        |                                                        |                                                        |
| E. kobei MSB1 1B                 | Flye Contig Lengths (bp)   | 4836273; 128941; 107534; 9421; 7516; 811               | 4836272; 244870; 8916                                  | 4836292; 136412; 108364; 9304                          |
|                                  | Plassembler Contig Lengths | 136482; 108411; 4665; 3715; 2370                       | 136482; 108411; 4665; 3715; 2370                       | 136482; 108411; 4665; 3715; 2370                       |
|                                  |                            |                                                        |                                                        |                                                        |
| Haemophilus sp002998595 M1C132 1 | Flye Contig Lengths (bp)   | 2080062; 48768; 39374; 21401; 19934; 5348              | 2051484; 39371; 21404; 19932                           | 2051487; 39372; 30363; 19931                           |
|                                  | Plassembler Contig Lengths | 48620 (likely prophage); 39398; 10719; 9975            | 39345; 10719; 9975                                     | 39345; 10719; 9975                                     |
|                                  |                            |                                                        |                                                        |                                                        |
| K. oxytoca MSB1 2C               | Flye Contig Lengths (bp)   | 5802201; 118088; 58433; 9133                           | 5802179; 118081; 58428; 9097                           | 5802238; 118089; 58430; 9133                           |
|                                  | Plassembler Contig Lengths | 118161; 58381; 10697; 9975; 4574                       | 118161; 58382; 10697; 9975; 4574                       | 118161; 58382; 10697; 9975; 4574                       |
|                                  |                            |                                                        |                                                        |                                                        |
| K. variicola INF345              | Flye Contig Lengths (bp)   | 5415037; 250849; 243507; 31779; 5772; 1499             | 5415055; 250850; 243514; 31778; 11570                  | 5415037; 250843; 243613; 31172; 17176; 11547; 6130     |
|                                  | Plassembler Contig Lengths | 250272; 243534; 30412 + 663 (linear); 5783; 3514; 1934 | 250272; 243534; 30412 + 663 (linear); 5783; 3514; 1934 | 250272; 243534; 30412 + 663 (linear); 5783; 3514; 1934 |