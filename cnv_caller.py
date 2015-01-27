import pysam
import numpy
import os
import math

samfile = "/Users/Simon/up_91_3.Founder_003_trimmed.sorted_rmdup_gc-corrected.bam"
mapq_cutoff = 1
window_size = 100
path = "/Users/Simon/output/"

outfile = file("test.out", "w")

sam = pysam.AlignmentFile(samfile, "rb")

def covScanner(sam, mapq_cutoff, name, pos_start, pos_end):
    WINDOW_COVS = []
    for pile in sam.pileup(name, pos_start, pos_end, truncate=True):
	cov = 0
        for reads in pile.pileups:
	    if reads.alignment.mapq >= mapq_cutoff:
                cov = cov + 1
	WINDOW_COVS.append(cov)
    window_mean = float(sum(WINDOW_COVS)) / float(len(WINDOW_COVS))
    return window_mean

def windowIterator(pos_start, pos_end, window_size):
    pos_end = pos_end + window_size
    pos_start = pos_start + window_size
    return pos_start, pos_end

def getNames(sam):
    header_dict = sam.header
    NAMES = []
    LENGTHS = []
    i = 0
    while i <= len(header_dict.get('SQ')):
        NAMES.append(header_dict.get('SQ')[i].get('SN'))
	LENGTHS.append(header_dict.get('SQ')[i].get('LN'))
    	i = i + 1
	return NAMES, LENGTHS

def median(lst):
    return numpy.median(numpy.array(lst))

def getMedianCov(sam, NAMES, LENGTH, mapq_cutoff):
    medians_dict = {}
    for name, ln in zip(NAMES, LENGTH): 
       COVS = []
       for pile in sam.pileup(name, 0, ln):
           cov = 0
           for reads in pile.pileups:
              if reads.alignment.mapq >= mapq_cutoff:
                   cov = cov + 1
   	   COVS.append(cov)
       chr_median = median(COVS)
       medians_dict[name] = chr_median
    return medians_dict

def getLogratio(window_mean, name, chr_medians):
    median = chr_medians[name]
    win_logratio = window_mean / median
    if win_logratio == 0: 
        win_logratio = 0
    else:
        win_logratio = math.log(win_logratio, 2)
    return win_logratio

	
NAMES, LENGTH = getNames(sam)
chr_medians = getMedianCov(sam, NAMES, LENGTH, mapq_cutoff)

for name, ln  in zip(NAMES, LENGTH):
    	pos_start = 0
	pos_end = window_size
	while pos_end <= ln:
	   pos_start, pos_end = windowIterator(pos_start, pos_end, window_size)
           window_mean = covScanner(sam, mapq_cutoff, name, pos_start, pos_end)
           out = file(os.path.join(path, name)
	   print >>  out, pos_start, pos_end, getLogratio(window_mean, name, chr_medians)
        


"""for pile in sam.pileup("Chr4", 13000, 13100, truncate=True):
        cov = 0
	print pile.n, pile.pos, 
        for reads in pile.pileups:
            if reads.alignment.mapq >= mapq_cutoff:
                cov = cov + 1
         
        print cov
"""
