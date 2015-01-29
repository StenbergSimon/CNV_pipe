cnv_caller.py
========

Python script to call copy number variations from platform independent NGS data

Usage: cnv_caller.py [options]

Options:
  -h, --help            

                        show this help message and exit
  
  -w INT, --windowSize=INT
                        
                        Windowsize (bp) to be used to calculate log2ratio
                        [Default:500]
  
  -m INT, --mappingQuality=INT
                        
                        Mapping quality cutoff for reads to be used in the
                        calculation [Default:0]
  
  -f FILE, --file=FILE  
              
                        Input bam file to be analyzed, should be sorted and
                        indexed
  
  -o PATH, --ouput=PATH
                        
                        Output path
  
  -l FILE, --name-list=FILE
                        
                        List of bam headers in order as they should be
                        plotted, [Default:/Users/Simon/git/CNV_pipe/chr.list]
  
  -a BOOLEAN, --plot=BOOLEAN
                        
                        pecify if plotting should be done using DNAcopy
                        [Default:True]
  
  -r FILE, --reference=FILE
                        
                        Bam file to be used as refernce / control

Developed by: Simon Stenberg

Dependencies
--------

pysam (install: pip install pysam)

R (developed in R version 3.1.0 (2014-04-10) -- "Spring Dance")

R library DNAcopy http://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html

Install
-------

Clone repo and go!

Preceeding data handling
-------
* You need at least 1 set of reads aligned to a reference genome/contigs
* Preferably use a refrence alignment of to compare with, for example ancestral genome (reference) vs evolved genome (sample)
* Bam files must be sorted and indexed
* Preferably use some kind of GC-correction. Deeptools is one example that can do this.

Caveats
-------

Things to know about algorithm

* Minimal reference coverage assumed at each base for the reference = 1

* Without reference log2ratio is calculated with chromosome median coverage

* Headers in both bams need to be identical

* Output path does not have to exist, it will be created if it can be

* List the bam headers (chromosome names) in a list to specify the order of chromosomes in plot. Example is included in repo (chr.list)
