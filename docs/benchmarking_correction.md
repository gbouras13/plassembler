# Long read correction benchmarking (v1.9.0)

From v1.9.0, `plassembler long` corrects its non-chromosomal reads with [lrcorrect](https://github.com/gbouras13/lrcorrect) instead of `canu -correct`. This page records the benchmarking behind that change. All inputs, outputs and scripts are on Zenodo: [doi:10.5281/zenodo.23656913](https://doi.org/10.5281/zenodo.23656913).

## Summary

* On 16 real read sets with curated reference genomes (68 plasmids), Unicycler recovered as many complete circular plasmids from lrcorrect's reads as from canu's: **41 vs 40**. Without correction it recovered 23.
* lrcorrect's corrected reads were more accurate than canu's on 15 of the 16 real read sets and tied on the other.
* lrcorrect was faster on 14 of the 16 (by 2.6-8.5x on 8 of them; about as fast as canu on the high-copy-plasmid read sets), and needed 0.6-3.6 GB of RAM against canu's 1.7-4.2 GB.
* lrcorrect's output is identical across reruns and thread counts. canu's was not: 4 runs gave 4 different outputs on one read set.
* lrcorrect exits cleanly with no reads when there is too little data. canu exits with an error, which `plassembler` could not tell apart from canu failing to run ([#92](https://github.com/gbouras13/plassembler/issues/92)).

## Why replace canu

canu needs perl and java. It took ~2 minutes even on a few Mbp of reads, it needs a genome size, and it broke in the Biocontainers image ([#92](https://github.com/gbouras13/plassembler/issues/92)): its job scripts call `expr --`, which BusyBox rejects, and its perl needs `PERL5LIB`, which `singularity exec` does not set. When it failed, `plassembler` carried on with uncorrected reads, which give Unicycler a much worse graph (see below).

## What was compared

Each corrector was given the same reads: those `plassembler long` would hand to `canu -correct`. That is, after chopper, Flye, mapping to Flye's chromosome and the entropy filter. They were run with 8 threads on one node of the University of Adelaide's Phoenix HPC (Intel Ice Lake).

| Corrector | Command |
|---|---|
| canu 2.3 | `canu -correct -p canu -d DIR genomeSize=G maxInputCoverage=C stopOnLowCoverage=1 maxThreads=8 useGrid=false -nanopore correctedErrorRate=0.12 READS` (exactly as `plassembler` v1.8.5; `-pacbio correctedErrorRate=0.045` for CLR) |
| lrcorrect 0.1.0 | `lrcorrect --threads 8 --preset ont --target-bases T --qual-char '#' -o OUT READS` (exactly as `plassembler` v1.9.0) |
| racon 1.5.0 | random subsample to T, then `minimap2 -x ava-ont --dual=yes` and `racon -f` |
| hifiasm 0.25.0 | random subsample to T, then `hifiasm --ont --write-ec` (ONT read sets only) |
| No correction | the reads as they are, which is what `plassembler` fell back to when canu failed |

G is the summed length of Flye's non-chromosomal contigs (10 kbp if there were none), C is 50 (250 if there were none), and T = G x C.

Each corrector's output was then assembled exactly as `plassembler long` does: `unicycler -s CORRECTED -l UNCORRECTED -t 8` (Unicycler 0.5.1, SPAdes 4.3.0).

Measurements:

* **Read identity**: corrected reads mapped to the reference with `minimap2 -x map-ont -c --eqx`. This is matching bases over all alignment columns, pooled across reads.
* **Complete circular plasmid**: a single Unicycler contig marked circular that covers at least 99% of a reference plasmid and is within 1% of its length. Errors are mismatches and indel bases against the reference.
* **Time and RAM**: wall time and peak RSS from GNU `time -v`.

## Real read sets

| Read set | Chemistry | Raw read identity | Source | Reference |
|---|---|---|---|---|
| *A. baumannii* J9, *C. koseri* MINF_9D, *E. kobei* MSB1_1B, *Haemophilus* M1C132_1, *K. oxytoca* MSB1_2C, *K. variicola* INF345 | R9.4.1, rapid (tech rep 1) and ligation (tech rep 2) libraries | 93.8-95.0% | [Wick et al. 2021](https://doi.org/10.1099/mgen.0.000631), [Monash Bridges](https://bridges.monash.edu/articles/dataset/Small_plasmid_Nanopore_data/13543754) | Wick et al. 2021 curated assemblies |
| *K. variicola* INF345 | not recorded (raw identity suggests R9.4.1, or R10.4.1 fast/HAC) | 94.5% | ERR15284527 (PRJEB93892, 2025) | Wick et al. 2021 |
| *K. pneumoniae* MGH78578 | R9.4.1 SUP | 95.0% | ERR8958826 (PRJEB51164, [Sanderson et al.](https://doi.org/10.1099/mgen.0.000910)) | GenBank CP000647-CP000652 |
| *K. pneumoniae* MGH78578 | R10.4 SUP | 98.0% | ERR8958710 (PRJEB51164) | GenBank CP000647-CP000652 |
| *K. pneumoniae* MGH78578 | R10.4.1 | 97.3% | ERR12342517 (PRJEB51164) | GenBank CP000647-CP000652 |

Each read set was subsampled to 100x of its reference with `rasusa reads -b <100 x genome size> -s 1`, then run through `plassembler long` (v1.8.5 code) up to read correction with `-c 1000000`. The MGH78578 reference was sequenced in 2007 and differs a little from the strain sequenced in 2022, so errors on its largest plasmid partly reflect real differences.

### Complete circular plasmids

| Read set | Plasmids | canu | lrcorrect | racon | No correction |
|---|---|---|---|---|---|
| J9, R9.4.1 rapid | 2 | 2 | 2 | 1 | 1 |
| J9, R9.4.1 ligation | 2 | 2 | 2 | 2 | 0 |
| MINF_9D, R9.4.1 rapid | 2 | 2 | 2 | 2 | 2 |
| MINF_9D, R9.4.1 ligation | 2 | 2 | 2 | 2 | 1 |
| MSB1_1B, R9.4.1 rapid | 6 | 3 | **5** | 3 | 3 |
| MSB1_1B, R9.4.1 ligation | 6 | 3 | 3 | 0 | 0 |
| M1C132_1, R9.4.1 rapid | 5 | 2 | 2 | 2 | 3 |
| M1C132_1, R9.4.1 ligation | 5 | 2 | 2 | 2 | 1 |
| MSB1_2C, R9.4.1 rapid | 4 | 1 | 1 | 3 | 1 |
| MSB1_2C, R9.4.1 ligation | 4 | 2 | **3** | 1 | 1 |
| INF345, R9.4.1 rapid | 5 | 3 | 3 | 3 | 1 |
| INF345, R9.4.1 ligation | 5 | **4** | 2 | 2 | 2 |
| INF345, 2025 | 5 | 2 | 2 | 2 | 2 |
| MGH78578, R9.4.1 SUP | 5 | **5** | 3 | 2 | 2 |
| MGH78578, R10.4 SUP | 5 | 2 | **3** | 0 | 3 |
| MGH78578, R10.4.1 | 5 | 3 | **4** | 3 | 0 (Unicycler crashed) |
| **Total** | **68** | **40** | **41** | **30** | **23** |

### Corrected reads and run time

| Read set | Input reads (Mbp) | canu identity | lrcorrect identity | canu reads kept | lrcorrect reads kept | canu time (s) / RAM (GB) | lrcorrect time (s) / RAM (GB) |
|---|---|---|---|---|---|---|---|
| J9, rapid | 8.7 | 99.87% | 99.91% | 565 | 913 | 110 / 3.6 | 70 / 1.4 |
| J9, ligation | 10.4 | 99.77% | 99.89% | 234 | 296 | 81 / 3.6 | 12 / 0.9 |
| MINF_9D, rapid | 29.3 | 99.76% | 99.88% | 197 | 339 | 86 / 1.7 | 28 / 1.0 |
| MINF_9D, ligation | 12.5 | 99.72% | 99.85% | 269 | 432 | 35 / 3.6 | 4 / 0.7 |
| MSB1_1B, rapid | 13.6 | 99.64% | 99.84% | 1841 | 2370 | 496 / 3.8 | 116 / 1.4 |
| MSB1_1B, ligation | 11.0 | 99.57% | 99.87% | 594 | 801 | 334 / 3.6 | 56 / 1.1 |
| M1C132_1, rapid | 46.0 | 99.73% | 99.86% | 345 | 686 | 75 / 2.1 | 42 / 1.0 |
| M1C132_1, ligation | 27.9 | 99.61% | 99.80% | 181 | 280 | 29 / 2.0 | 4 / 0.6 |
| MSB1_2C, rapid | 29.8 | 99.29% | 99.71% | 3339 | 4503 | 375 / 4.2 | 146 / 2.2 |
| MSB1_2C, ligation | 10.0 | 99.65% | 99.83% | 482 | 565 | 62 / 3.5 | 16 / 1.0 |
| INF345, rapid | 49.7 | 99.70% | 99.88% | 2183 | 3929 | 474 / 3.9 | 461 / 3.6 |
| INF345, ligation | 26.1 | 99.66% | 99.86% | 1414 | 4209 | 385 / 3.6 | 264 / 2.2 |
| INF345, 2025 | 52.0 | 99.77% | 99.94% | 1650 | 4038 | 364 / 3.9 | 274 / 2.4 |
| MGH78578, R9.4.1 SUP | 78.5 | 99.78% | 99.88% | 1473 | 3049 | 404 / 3.8 | 419 / 3.2 |
| MGH78578, R10.4 SUP | 46.6 | 99.95% | 99.95% | 1216 | 2705 | 302 / 3.8 | 295 / 2.3 |
| MGH78578, R10.4.1 | 37.5 | 99.84% | 99.95% | 1991 | 3783 | 392 / 3.9 | 437 / 3.0 |

* canu keeps only the longest ~40x of its input (its `corOutCoverage` default), not the 50x `plassembler` asks for. lrcorrect keeps a random 50x subsample whole, which matters for small plasmids, whose reads are short.
* lrcorrect is slowest on the MGH78578 and INF345 read sets. These are dominated by small plasmids at 1,000-2,400x depth: the number of overlaps grows with the square of depth.
* On MSB1_2C rapid, Flye's chromosome was fragmented, so 2.7 Mbp of chromosomal contigs counted as "plasmid" and no subsampling happened. Every corrector gave a very fragmented assembly (48-143 contigs).

### Small plasmids that every corrector misses

Many small plasmids are missed whatever the corrector: M1C132_1's 7.4 and 5.7 kbp, MSB1_1B's 3.7 kbp, and MSB1_2C's 10 kbp in most runs. In these cases the corrected reads cover the plasmid, but Unicycler does not resolve it. MGH78578's 4.3 and 3.5 kbp plasmids share an 820 bp region at 96% identity, and most runs, with any corrector, recover only one of them.

## Simulated read sets

Simulated with [Badread](https://github.com/rrwick/Badread) v0.4.2:

* **ab95**: the 145 kbp *A. baumannii* J9 plasmid in `tests/test_data`, ~95% identity, 7 Mbp.
* **tp95**: the 44 and 9.3 kbp plasmids in `tests/test_data`, ~95% identity, 2.1 Mbp.
* **R9-like, R10-like, CLR**: J9's 145 kbp plasmid plus the 44 and 9.3 kbp plasmids at relative depths 1:3:10, 14.8 Mbp each, with `--small_plasmid_bias`:
  * R9-like: `nanopore2020` model, identity 92,98,4.
  * R10-like: `nanopore2023` model, identity 99,99.9,0.5.
  * CLR: `pacbio2016` model, identity 87,95,3.

| Read set | Corrector | Read identity | Reads kept | Time (s) | RAM (GB) | Unicycler (errors per plasmid) |
|---|---|---|---|---|---|---|
| ab95 | canu | 99.94% | 263/480 | 126 | 3.6 | 1 circular (0) |
| ab95 | lrcorrect | 100.00% | 478/480 | 52 | 1.1 | 1 circular (0) |
| ab95 | none | 94.64% | 480/480 | - | - | 4.3 kbp fragment |
| tp95 | canu | 99.91% | 104/172 | 50 | 1.7 | 2 circular (0, 0) |
| tp95 | lrcorrect | 100.00% | 168/172 | 11 | 0.7 | 2 circular (0, 0) |
| tp95 | none | 94.77% | 172/172 | - | - | SPAdes failed |
| R9-like | canu | 99.69% | 367/934 | 242 | 3.6 | Unicycler crashed |
| R9-like | lrcorrect | 99.97% | 646/934 | 77 | 1.5 | 3 circular (41, 5, 9) |
| R9-like | hifiasm | 93.35% | 647/934 | 10 | 16.9 | 2 circular, 1 partial |
| R10-like | canu | 99.9995% | 391/967 | 182 | 3.6 | 3 circular (0, 0, 0) |
| R10-like | lrcorrect | 100.00% | 656/967 | 86 | 1.5 | 3 circular (0, 0, 0) |
| R10-like | none | 98.72% | 967/967 | - | - | 3 circular (1, 0, 0); Unicycler 1.6x slower |
| CLR | canu | 98.55% | 412/943 | 123 | 3.6 | 3 circular (154, 11, 24) |
| CLR | lrcorrect | 99.83% | 629/943 | 131 | 1.4 | 3 circular (58, 5, 19) |
| CLR | racon | 99.45% | 645/943 | 157 | 0.2 | 3 circular (31, 0, 5) |

hifiasm's ONT correction did little on these reads, and by default it reserves ~16 GB for its bloom filter whatever the input size (`-f0` avoids this).

### Too little data

| Input | canu | lrcorrect |
|---|---|---|
| empty, 1 read, 5 reads, a plasmid-free sample (2 reads) | exits 1 | exits 0, no reads, says why |
| 3x of the 145 kbp plasmid | exits 1 (below its 10x minimum) | 27 reads; Unicycler gives a partial plasmid |
| 8x of the 145 kbp plasmid | exits 1; Unicycler fails on the raw reads `plassembler` falls back to | 79 reads; Unicycler gives one circular 145 kbp contig covering 98% of the plasmid |

### Determinism

Four runs of each corrector on the same input, with 8 threads:

* **canu**: the same output each time on ab95, but 4 different outputs on the R9-like set.
* **racon, hifiasm, lrcorrect**: the same output each time on both.
* **lrcorrect**: byte-identical output at 1, 4 and 8 threads, and so were the Unicycler assemblies built from it.

## Data and scripts

The [Zenodo record](https://doi.org/10.5281/zenodo.23656913) holds:

* the scripts;
* the exact input reads each corrector received, with each set's truth reference;
* every corrector's output, Unicycler assembly, log and metrics;
* `results_long.tsv`, one row per read set and corrector.

The original real read sets are not re-hosted. They are at the accessions above, and `fetch_reads.sh` downloads them.
