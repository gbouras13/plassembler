# Long

Note: This method has been superceded by v1.3.0 but I am retaining the below out of interest. If you want the old `plassembler long` method implemented from v1.3.0 onwards for some reason, you can do this with `--canu_flag`.

After reading this [paper](https://doi.org/10.1099/mgen.0.001024) by Johnson et al, it seemed interesting that while most assemblers failed to recover small <10kbp> plasmids, Canu always did - albeit with multiplication. However, multplicated contigs in [Canu](https://github.com/marbl/canu) indicate the smallest repeating sequence unit in the header, and so can be trimmed easily (see e.g. Ryan Wick's [script](https://github.com/rrwick/Trycycler/blob/main/scripts/canu_trim.py) in Trycycler).

Therefore, from v1.2.0, `plassembler long` will follow all the same steps as `plassembler run` (i.e. hybrid mode), but instead of Unicycler, it will run Canu to recover plasmids in the unmapped reads.

Of course, another big issue with long read is size selection - if your read set doesn't have many reads small enough to be a part of your <10kbp plasmids>, then nothing will work - your best bet is still to get some short read sequencing data. 

But assuming they are in the read set, `plassembler long` should hopefully recover your small plasmids.

# Results

I tested `plassembler long` on the 60x simulated reads I generated for benchmarking `plassembler` for the manuscript and can be found [here](https://zenodo.org/record/7996690). from Wick et al, C222 (Houtak et al) and Cav1217 (Mathers et al) (find more information [here](https://plassembler.readthedocs.io/en/latest/benchmarking_results_sim) and [here](https://github.com/gbouras13/plassembler_simulation_benchmarking)).

Overall, it seems to work pretty well - not perfect (it missed a linear and small plasmid in _K variicola_ and also tends to assemble some non-circular chimeric contigs), but it seems anything circular is a real plasmid from the ground truth. 

# Other Options

You could try running Canu on the entire assembly and use Ryan Wick's [script](https://github.com/rrwick/Trycycler/blob/main/scripts/canu_trim.py) to trim the multiplicated plasmid contigs.

Ryan also [thinks](https://rrwick.github.io/2023/05/05/ont-only-accuracy-with-r10.4.1.html) Canu is a pretty good assembler generally (it is perhaps more accurate than Flye). The downside is time though, I have found that it takes 5-10x more time to run than Flye on my isolates.


# Results

Everything was run on my Macbook Pro M1 (2020) with 8 threads. I was too lazy to run QUAST but I would recommend polishing the output with your favourite long read polisher (Medaka) anyway. Below are the results (in terms of lengths).

## By Isolate

(c) = Canu marked the plasmid as circular 

| Isolate                     | Ground Truth                                   | `plassembler long` v1.2.0                                            | Time (s) | Flye (Within `plassembler long` )                  |
| --------------------------- | ---------------------------------------------- | ------------------------------------------------------------------ | -------- | ------------------------------------------ |
| C222                        | 2473                                           | 2473                                                               | 1917     | Nothing - missed 2473                      |
| _Acinetobacter baumannii_ J9  | 145059; 6078                                   | 145058 (c); 6078 (c)                                               | 967      | 145059; 10771                              |
| CAV1217                     | 181436;  70606; 44015; 9294                    | 181429 (c);  70603 (c); 44015 (c); 9294 (c)                        | 1230     | 181433;  70605; 44015; 9294                |
| _Citrobacter koseri_ MINF 9D  | 64962; 9294                                    | 93661 (chimera includes 64962 plasmid); 9294 (c)                   | 1196     | 64961; 9294                                |
| _Enterobacter kobei_ MSB1 1B  | 136482;  108411;  4665;  3715;  2370           | 136477 (c);  108402 (c);  4665 (c);  3715 (c);  2367 (c)           | 1374     | 136477; 108408 - missed 3 small plasmids   |
| _Klebsiella oxytoca_ MSB1 2C  | 118161;  58472; 9975; 4574                     | 118159 (c);  58467 (c); 9973 (c); 4573 (c); extra 18290 (chimera)  | 1646     | 118161;  58472; 9975 - missed 4574         |
| _Klebsiella variicola_ INF345 | 250980;  243620;  31780 (linear);  5783;  3514 | 250968 (c);  243616 (c);  3514 (c) - missing linear plasmid + 5783 | 3544     | 249712; 243617; 31645 - missed 5783 + 3514 |


## Overall

| Total                     | Missed Small Plasmids | Missassembled         |
| ------------------------- | --------------------- | --------------------- |
| `plassembler long` v1.2.0   | 1 + linear plasmid    | 2 (_C. koseri_ 64962) and _K oxtyoca_ chimera   |
| Flye (Within `plassembler`) | 6                     | 1 (_A. baumannii_ 6078) |

## Summary

It's pretty good! Not perfect and I'd still recommend short reads if you can get them, but not too bad. Flye still seems the best for linear plasmids as well. 

```
Usage: plassembler long [OPTIONS]
Plassembler with long reads only

Options:
  -h, --help                Show this message and exit.
  -V, --version             Show the version and exit.
  -d, --database PATH       Directory of PLSDB database.  [required]
  -l, --longreads PATH      FASTQ file of long reads.  [required]
  -c, --chromosome INTEGER  Approximate lower-bound chromosome length of
                            bacteria (in base pairs).  [default: 1000000]git
  -o, --outdir PATH         Directory to write the output to.  [default:
                            plassembler.output/]
  -m, --min_length TEXT     minimum length for filtering long reads with
                            chopper.  [default: 500]
  -q, --min_quality TEXT    minimum quality q-score for filtering long reads
                            with chopper.  [default: 9]
  -t, --threads TEXT        Number of threads.  [default: 1]
  -f, --force               Force overwrites the output directory.
  -p, --prefix TEXT         Prefix for output files. This is not required.
                            [default: plassembler]
  --skip_qc                 Skips qc (chopper and fastp).
  --pacbio_model TEXT       Pacbio model for Flye.  Must be one of pacbio-raw,
                            pacbio-corr or pacbio-hifi.  Use pacbio-raw for
                            PacBio regular CLR reads (<20 percent error),
                            pacbio-corr for PacBio reads that were corrected
                            with other methods (<3 percent error) or pacbio-
                            hifi for PacBio HiFi reads (<1 percent error).
  -r, --raw_flag            Use --nano-raw for Flye.  Designed for Guppy fast
                            configuration reads.  By default, Flye will assume
                            SUP or HAC reads and use --nano-hq.
  --keep_chromosome         If you want to keep the chromosome assembly.
  --flye_directory PATH     Directory containing Flye long read assembly.
                            Needs to contain assembly_info.txt and
                            assembly_info.fasta. Allows Plassembler to Skip
                            Flye assembly step.
```