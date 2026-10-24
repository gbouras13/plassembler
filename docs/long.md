# Long v1.3.0 Updates (24 October 2023)

 Note: as of v1.9.0, reads are corrected with [lrcorrect](https://github.com/gbouras13/lrcorrect) instead of `canu -correct`, and the old canu-based method (`--canu_flag`) has been removed. See [Read correction](#read-correction-v190) below.

While I was mostly happy with with the `plassembler long` update in v 1.2.0, while testing it a lot of real life data benchmarking [hybracter](https://github.com/gbouras13/hybracter), I found some strange instances where it would give way too many contigs as output - seemingly, it would sometimes assemble the same plasmid multiple times with slight variations (probably due to various factors or too deep read sets). In any case, I wanted a better automated solution. 

Inspired by this [tweet](https://twitter.com/rrwick/status/1548926644085108738/) by Ryan Wick, I decided to experiment with treating long reads as both short reads (in the sense of creating a de Brujin graph based assembly) and long reads (for scaffolding) in Unicycler - and the results were great.

As you can see in the results, I am confident that assuming you have good quality R10 Nanopore data, `plassembler long` should now recover small plasmids.

`plassembler long` now implements the following steps (after the Flye assembly and getting the plasmid long reads as previous):

* Removes extremely low entropy repetitive reads
* Runs `lrcorrect` (before v1.9.0, `canu -correct`) to simultaneously error correct the long reads and subsample them to 50x of the length of Flye's non-chromosomal contigs (so as not to give Unicycler too much read depth - see [this](https://doi.org/10.1093/bioinformatics/btv311) and [this](https://doi.org/10.1371/journal.pone.0060204) and [this](https://academic.oup.com/bioinformatics/article/27/4/479/198367)) - major thanks to Ryan Wick for suggesting suggesting this!
* Runs Unicycler as follows: `unicycler -s {error_corrected_longreads} -l {all_longreads}`

# Read correction (v1.9.0)

From v1.9.0, `plassembler long` corrects the non-chromosomal reads with [lrcorrect](https://github.com/gbouras13/lrcorrect) rather than `canu -correct`. lrcorrect is a single Rust binary with minimap2 built in. It finds all-vs-all overlaps, aligns each read's best partners to it and calls a majority-vote consensus, over two rounds. It needs no perl, java or genome size. It is deterministic: the output is identical whatever the thread count.

The reads are subsampled at random to 50x of the length of Flye's non-chromosomal contigs, or to 2.5 Mbp when Flye found none (the targets canu was given), then corrected and given to Unicycler with `-s`. If there are too few reads to correct anything, `plassembler` carries on with the uncorrected reads, as before. If lrcorrect fails to run, `plassembler` stops and shows lrcorrect's error.

Comparison on the same non-chromosomal read sets: 8 threads on one Phoenix (Adelaide HPC) node, Unicycler v0.5.1, and identity measured by mapping to the truth with minimap2. All read sets are simulated with [Badread](https://github.com/rrwick/Badread):

* ab95: the 145 kbp *A. baumannii* plasmid in `tests/test_data`, ~95% read identity.
* tp95: the 44 kbp and 9.3 kbp plasmids in `tests/test_data`, ~95%.
* R9-like, R10-like and CLR: 145 + 44 + 9.3 kbp plasmids at relative depths 1:3:10, at ~92%, ~98.7% and ~87% identity.

| Read set | Corrector | Read identity | Reads kept | Time (s) | Peak RAM (GB) | Unicycler result (errors per plasmid) |
|---|---|---|---|---|---|---|
| ab95 | canu | 99.94% | 263/480 | 126 | 3.6 | 1 circular (0) |
| ab95 | lrcorrect | 100.00% | 478/480 | 52 | 1.1 | 1 circular (0) |
| tp95 | canu | 99.91% | 104/172 | 50 | 1.7 | 2 circular (0, 0) |
| tp95 | lrcorrect | 100.00% | 168/172 | 11 | 0.7 | 2 circular (0, 0) |
| R9-like | canu | 99.69% | 367/934 | 242 | 3.6 | Unicycler crashed |
| R9-like | lrcorrect | 99.97% | 646/934 | 77 | 1.5 | 3 circular (41, 5, 9) |
| R10-like | canu | 99.9995% | 391/967 | 182 | 3.6 | 3 circular (0, 0, 0) |
| R10-like | lrcorrect | 100.00% | 656/967 | 86 | 1.5 | 3 circular (0, 0, 0) |
| R10-like | none | 98.72% | 967/967 | - | - | 3 circular (1, 0, 0), Unicycler 1.6x slower |
| CLR | canu | 98.55% | 412/943 | 123 | 3.6 | 3 circular (154, 11, 24) |
| CLR | lrcorrect | 99.83% | 629/943 | 131 | 1.4 | 3 circular (58, 5, 19) |

On real reads, with curated reference genomes, the comparison is closer. The read sets are:

* the six isolates from [Wick et al. 2021](https://doi.org/10.1099/mgen.0.000631), each sequenced with both an R9.4.1 rapid and a ligation library (12 read sets);
* *K. pneumoniae* MGH78578 on R9.4.1 SUP, R10.4 SUP and R10.4.1 ([Sanderson et al.](https://doi.org/10.1099/mgen.0.000910));
* *K. variicola* INF345 from a 2025 ONT benchmarking run (PRJEB93892), whose raw reads are ~94.5% accurate.

Each set was subsampled to ~100x and run through `plassembler long` up to read correction. Counts are complete circular plasmids:

| Read set | Plasmids | canu | lrcorrect | No correction | canu time (s) | lrcorrect time (s) |
|---|---|---|---|---|---|---|
| *A. baumannii* J9, R9.4.1 rapid / ligation | 2 / 2 | 2 / 2 | 2 / 2 | 1 / 0 | 110 / 81 | 70 / 12 |
| *C. koseri* MINF_9D, R9.4.1 rapid / ligation | 2 / 2 | 2 / 2 | 2 / 2 | 2 / 1 | 86 / 35 | 28 / 4 |
| *E. kobei* MSB1_1B, R9.4.1 rapid / ligation | 6 / 6 | 3 / 3 | 5 / 3 | 3 / 0 | 496 / 334 | 116 / 56 |
| *Haemophilus* M1C132_1, R9.4.1 rapid / ligation | 5 / 5 | 2 / 2 | 2 / 2 | 3 / 1 | 75 / 29 | 42 / 4 |
| *K. oxytoca* MSB1_2C, R9.4.1 rapid / ligation | 4 / 4 | 1 / 2 | 1 / 3 | 1 / 1 | 375 / 62 | 146 / 16 |
| *K. variicola* INF345, R9.4.1 rapid / ligation | 5 / 5 | 3 / 4 | 3 / 2 | 1 / 2 | 474 / 385 | 461 / 264 |
| *K. variicola* INF345, 2025 run | 5 | 2 | 2 | 2 | 364 | 274 |
| *K. pneumoniae* MGH78578, R9.4.1 SUP | 5 | 5 | 3 | 2 | 404 | 419 |
| *K. pneumoniae* MGH78578, R10.4 SUP | 5 | 2 | 3 | 3 | 302 | 295 |
| *K. pneumoniae* MGH78578, R10.4.1 | 5 | 3 | 4 | 0 (Unicycler crashed) | 392 | 437 |
| **Total** | **68** | **40** | **41** | **23** | | |

* lrcorrect's reads were more accurate than canu's on 15 of the 16 real read sets and tied on the other, and it used 0.6-3.6 GB of RAM against canu's 1.7-4.2 GB.
* lrcorrect was faster on 13 of the 16 sets, usually 2-8x. On the read sets dominated by a high-copy small plasmid it took about as long as canu.
* Many small plasmids are missed whichever corrector is used. In these cases the corrected reads are there, but Unicycler does not resolve the plasmid, for example two small MGH78578 plasmids that share an 820 bp region.
* The reads in all of these sets are 94-98% accurate. None is modern R10.4.1 (Q20+) data, where the simulations suggest correction mostly saves Unicycler time.

* canu only keeps the longest ~40x of reads (its `corOutCoverage`), not the 50x `plassembler` asked for. lrcorrect keeps the whole random 50x subsample.
* Without correction, Unicycler failed or gave a fragment on the ~95%, ~92% and CLR sets. Modern R10.4.1 reads (~Q20) assemble without correction, but more slowly.
* With few reads (one or five reads, a sample with no plasmids, empty input), canu exits with an error. lrcorrect exits cleanly with no reads, which `plassembler` reports as "too few reads to correct". At 8x depth canu refuses to run (below its 10x minimum), and Unicycler fails on the raw reads it falls back to. Correcting with lrcorrect gives a single 143 kbp contig instead.
* canu's output changed between reruns on the R9-like set (4 runs, 4 different outputs). lrcorrect's did not, and was byte-identical at 1, 4 and 8 threads.

# Results

I tested `plassembler long` on 60x simulated reads generated using badread with the `nanopore2023` model on the same isolates 
from Wick et al, C222 (Houtak et al) and Cav1217 (Mathers et al) used in the Plassembler manuscript. 

Overall, it seems to work almost perfectly and quite performantly. In particular, the new approach seems to provide a speed-up vs v1.2.0  for more complicated assemblies (e.g. _Klebsiella variicola_ INF345 ) as the Unicycler step is quite fast. 

See [here](#Long-v1.2.0-(superceded)) for the old benchmarking data and times run on the same Macbook machine, although note they were done with different read sets (simulated R9 vs R10) so it's not truly a fair comparison.

The only misassembly is of the linear plasmid in _K variicola_ , which is a known issue of the Unicycler approach. If you have linear plasmids, please use a long-read assembly appraoch e.g. with Flye.

Everything was run on my Macbook Pro M1 (2020) with 8 threads. I was too lazy to run QUAST but I would recommend polishing the output with your favourite long read polisher (e.g. Medaka) anyway as implemented in [hybracter](https://github.com/gbouras13/hybracter). 

## By Isolate

(c) = Unicycler marked the plasmid as circular 

| Isolate                     | Ground Truth                                   | `plassembler long` v1.3.0                                            | Time (s) | Flye (Within `plassembler long` )                  |
| --------------------------- | ---------------------------------------------- | ------------------------------------------------------------------ | -------- | ------------------------------------------ |
| C222                        | 2473                                           | 2473 (c)                                                              | 889     | Nothing - missed 2473                      |
| _Acinetobacter baumannii_ J9  | 145059; 6078                                   | 145058 (c); 6078 (c)                                               | 1179      | 145059; 6077                              |
| CAV1217                     | 181436;  70606; 44015; 9294                    | 181435 (c);  70605 (c); 44015 (c); 9293 (c)                        | 1582     | 181433; 70609; 44015; 9294                |
| _Citrobacter koseri_ MINF 9D  | 64962; 9294                                    | 64961 (c); 9294 (c)                   | 1328     | 64962; 18088                               |
| _Enterobacter kobei_ MSB1 1B  | 136482;  108411;  4665;  3715;  2370           | 136480 (c);  108410 (c);  4665 (c);  3715 (c);  2368 (c)           | 1579     | 136481; 108410 - missed 3 small plasmids   |
| _Klebsiella oxytoca_ MSB1 2C  | 118161;  58472; 9975; 4574                     | 118160 (c);  58471 (c); 9975 (c); 4574 (c)  | 1273     | 118161;  58472; 9975 - missed 4574         |
| _Klebsiella variicola_ INF345 | 250980;  243620;  31780 (linear);  5783;  3514 | 250976 (c);  243612 (c); 30408 (linear incomplete); 5783 (c); 3514 (c)  | 1206     | 250979; 243618; 31742 - missed 5783 + 3514 |

## Overall

| Total                     | Missed Small Plasmids | Missassembled         |
| ------------------------- | --------------------- | --------------------- |
| `plassembler long` v1.3.0   | 0     | 1 _K oxtyoca_ linear plasmid   |
| Flye (Within `plassembler`) | 5                     | 0 |

## Summary

While I'd still recommend short reads if you can get them, I am now confident that if your isolate has small plasmids, `plassembler long` should find them.

```
Usage: plassembler long [OPTIONS]

  Plassembler with long reads only

Options:
  -h, --help                    Show this message and exit.
  -V, --version                 Show the version and exit.
  -d, --database PATH           Directory of PLSDB database.  [required]
  -l, --longreads PATH          FASTQ file of long reads.  [required]
  -c, --chromosome INTEGER      Approximate lower-bound chromosome length of
                                bacteria (in base pairs).  [default: 1000000]
  -o, --outdir PATH             Directory to write the output to.  [default:
                                plassembler.output/]
  -m, --min_length TEXT         minimum length for filtering long reads with
                                chopper.  [default: 500]
  -q, --min_quality TEXT        minimum quality q-score for filtering long
                                reads with chopper.  [default: 9]
  -t, --threads TEXT            Number of threads.  [default: 1]
  -f, --force                   Force overwrites the output directory.
  -p, --prefix TEXT             Prefix for output files. This is not required.
                                [default: plassembler]
  --skip_qc                     Skips qc (chopper and fastp).
  --pacbio_model TEXT           Pacbio model for Flye.  Must be one of pacbio-
                                raw, pacbio-corr or pacbio-hifi.  Use pacbio-
                                raw for PacBio regular CLR reads (<20 percent
                                error), pacbio-corr for PacBio reads that were
                                corrected with other methods (<3 percent
                                error) or pacbio-hifi for PacBio HiFi reads
                                (<1 percent error).
  -r, --raw_flag                Use --nano-raw for Flye.  Designed for Guppy
                                fast configuration reads.  By default, Flye
                                will assume SUP or HAC reads and use --nano-
                                hq.
  --keep_chromosome             If you want to keep the chromosome assembly.
  --flye_directory PATH         Directory containing Flye long read assembly.
                                Needs to contain assembly_info.txt and
                                assembly_info.fasta. Allows Plassembler to
                                Skip Flye assembly step.
  --flye_assembly PATH          Path to file containing Flye long read
                                assembly FASTA. Allows Plassembler to Skip
                                Flye assembly step in conjunction with
                                --flye_info.
  --flye_info PATH              Path to file containing Flye long read
                                assembly info text file. Allows Plassembler to
                                Skip Flye assembly step in conjunction with
                                --flye_assembly.
  --no_chromosome           Run Plassembler assuming no chromosome can be
                            assembled. Use this if your reads only contain
                            plasmids that you would like to assemble.
```