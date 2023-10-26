# Long v1.3.0 Updates (24 October 2023)

 Note: as of v1.3.0, if still you want the old canu-based `plassembler long` method implemented in v1.2.0 for some reason, you can do this with `--canu_flag`.

While I was mostly happy with with the `plassembler long` update in v 1.2.0, while testing it a lot of real life data benchmarking [hybracter](https://github.com/gbouras13/hybracter), I found some strange instances where it would give way too many contigs as output - seemingly, it would sometimes assemble the same plasmid multiple times with slight variations (probably due to various factors or too deep read sets). In any case, I wanted a better automated solution. 

Inspired by this [tweet](https://twitter.com/rrwick/status/1548926644085108738/) by Ryan Wick, I decided to experiment with treating long reads as both short reads (in the sense of creating a de Brujin graph based assembly) and long reads (for scaffolding) in Unicycler - and the results were great.

As you can see in the results, I am confident that assuming you have good quality R10 Nanopore data, `plassembler long` should now recover small plasmids.

`plassembler long` now implements the following steps (after the Flye assembly and getting the plasmid long reads as previous):

* Removes extremely low entropy repetitive reads
* Runs `canu -correct` to simultaneously error correct the long reads and subsample to 100x estimated depth (so as not to give Unicycler too much read depth - see [this](https://doi.org/10.1093/bioinformatics/btv311) and [this](https://doi.org/10.1371/journal.pone.0060204) and [this](https://academic.oup.com/bioinformatics/article/27/4/479/198367)) - major thanks to Ryan Wick for suggesting suggesting this!
* Runs Unicycler as follows: `unicycler -s {error_corrected_longreads} -l {all_longreads}`

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
| C222                        | 2473 (c)                                          | 2473                                                               | 889     | Nothing - missed 2473                      |
| _Acinetobacter baumannii_ J9  | 145059; 6078                                   | 145058 (c); 6078 (c)                                               | 1179      | 145059; 6077 (c)                              |
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
  --canu_flag                   Runs canu instead of Unicycler (aka replicates
                                v1.2.0). As of v1.3.0, Unicycler is the
                                assembler for long reads. Canu is only
                                recommended if you have low quality reads
                                (e.g. ONT R9).
  --corrected_error_rate FLOAT  Corrected error rate parameter for canu
                                -correct. For advanced users only.
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
```