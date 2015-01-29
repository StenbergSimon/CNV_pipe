CNV_pipe
========

Pearl scripts involved in CNV calling from NGS data

Developed by: Simon Stenberg

Dependencies
--------
Bedtools (Developed in v2.17.0)

Samtools (Developed in 0.1.18)

Bamtools (Developed in 2.3.0)

R (developed in R version 3.1.0 (2014-04-10) -- "Spring Dance") - if -rall option is used

R library ggplot2 (if -rall options is used)

Install
-------

Copy plot_cnv.r into ~/bin/

Preceeding data handling
-------
* You need at least 1 set of reads aligned to a reference genome/contigs
* Preferably use a refrence alignment of to compare with, for example ancestral genome (reference) vs evolved genome (sample)
* Bam files must be sorted and indexed
* Preferably use some kind of GC-correction. Deeptools is one example that can do this.


Usage
--------

Call the script followed with your sample .bam alignment and the preceding options

cnv_pipe.pl sample.bam [options] > output.tsv

The output is a file with the called windows that is atleast the cutoff value below or above 0.

Options
--------

-w [int] sliding window size, bigger window decrease computanional time, decrease sensitivity and decrease false posetive rate (scales with read depth) Default: 500

-co [int] Cutoff value of log2 change compared to reference or name median to be called and placed in file

-rall prints out all window values even if cutoff is met or not, placed in separate folder (/reports/), also produces .pdf files with plots of all names

-incr [int] increments of the sliding window, controls overlap. Default: [windowsize] 

-mapq [int] mapping quality cutoff to include in the analysis. Default: [0]

-r Reference.bam reference alignment

-regex [regex] supply own regex to recognize the names in the alignment file, does not work at the moment. Default: (SN:)([^\s]+)

-h show help

--help show help

Additional Info
--------
You might find interesting scripts in variant_calling repo to complement the cnv_pipe.pl


Known Bugs
-------

Increments =! Windowsize will currently give wrongly scaled outputs in position


cnv_caller.py
-------

Things to know:

Minimal reference coverage assumed at each base for the reference = 1

Without reference log2ratio is calculated with chromosome median coverage
