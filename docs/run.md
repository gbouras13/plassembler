To run plassembler

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length>`

* -c will default to 2500000 (a lower bound for the estimated genome length of staphylococcus aureus).

To specify threads:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads>`

To specify a prefix for the output files:

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix>`

To specify a minimum length and minimum read quality Q-score for nanofilt :

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality>`

* -m will default to 1000 and -q will default to 9. Note that for some tiny plasmids, -m may need to be reduced (see [https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000631#tab2)).

To overwrite an existing output directory, use -f

` plassembler.py -l <long read fastq> -o <output dir> -s1 < short read R1 fastq> -s2 < short read R2 fastq>  -c <estimated chromosome length> -t <threads> -p <prefix> -m <min length> -q <min quality> -f`

plassembler defaults to 8 threads.

**What happens if plassembler fails to find a plasmid?**

* There are two reasons why plassembler will fail to find a plasmid:

1. Where Flye assembles a complete chromosome, but fails to find any plasmids. Most of the time, this simply means that there are no plasmids in the your bacterial isolate. For now, plassembler will just finish. I am working on adding a module to identify reads that do not map to the chromosome and assembling those for plasmids - to identify plasmids that Flye may miss.
2. Where there is insufficient coverage for Flye to assemble a complete circular chromosome. In these cases, it is recommended that you either do more long read sequencing so that you can assemble a complete circular chromosome, or use Unicycler directly to create a hybrid assembly.
