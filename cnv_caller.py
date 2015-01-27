import pysam
import numpy
import os
import math

samfile = "/Users/Simon/up_91_3.Founder_003_trimmed.sorted_rmdup_gc-corrected.bam"
mapq_cutoff = 1
window_size = 100
path = "/Users/Simon/output/"

sam = pysam.AlignmentFile(samfile, "rb")
"""
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
"""

def getNames(sam):
    header_dict = sam.header
    NAMES = []
    LENGTHS = []
    i = 0
    while i < len(header_dict.get('SQ')):
        NAMES.append(header_dict.get('SQ')[i].get('SN'))
        LENGTHS.append(header_dict.get('SQ')[i].get('LN'))
        i = i + 1
    sam.close()
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

def writeRscript(FILELIST):
    with open(os.path.join(path, "run_DNA_copy.r"), "w") as out:
        out.write("#!/usr/bin/env Rscript")
        out.write("library(DNAcopy)")
	p = 1
	for line in FILELIST:
	   out.write("".join(["c", str(p), "<- read.table(\"", line]))	

class CovScanner():
  
   def __init__(self):
       self.window_size = ()
       self.pos_start = 0
       self.pos_end = ()
       self.samfile = ()
       self.name = ()
       self.mapq_cutoff = ()
       self.MEANS = []
    
   def setWindowSize(self, win):
       self.window_size = win
       self.pos_end = win
   
   def setSamFile(self, sam):
       self.samfile = sam
   
   def setMapq(self, mapq):
       self.mapq_cutoff = mapq

   def setName(self, name):
       self.name = name
 
   def move(self, ln):
       if self.pos_end < ln:
          windowMean()
          self.pos_start = self.pos_start + self.window_size
          self.pos_end = self.pos_end + self.window_size
       else:
          return self.MEANS
  
   def windowMean(self):
       sam = pysam.AlignmentFile(samfile, "rb")
       WINDOW_COVS = []
       for pile in sam.pileup(self.name, self.pos_start, self.pos_end, truncate=True):
          cov = 0
          for reads in pile.pileups:
             if reads.alignment.mapq >= mapq_cutoff:
                cov = cov + 1
       WINDOW_COVS.append(cov)
       sam.close()
       window_mean = numpy.mean(WINDOW_COVS)
       self.MEANS.append(window_mean)  

class FilePrinter():

   def __init__(self, path):
       self._path = path
       self._fh = None

   def __enter__(self): 

       self._fh = open(self._path, "w")
       return self

   def __exit__(self, *args):

       self._fh.close()
       self._fh = None

   def printToFile(self, *data):

       if (self._fh is None):
           raise Exception("Only use within with statement blocks")

       line = '\t'.join(map(str, data)) + "\n"
       self._fh.write(line)
       

""" START """
if __name__ == "__main__":
	
    NAMES, LENGTH = getNames(sam)
    chr_medians = getMedianCov(sam, NAMES, LENGTH, mapq_cutoff)

    if not os.path.exists(path):
                    os.makedirs(path)

    for name, ln  in zip(NAMES, LENGTH):
        with FilePrinter(os.path.join(path, name)) as out:
            scan = CovScanner()
	    scan.setSamFile(samfile)
	    scan.setMapq(mapq_cutoff)
            scan.setWindowSize(window_size)
	    scan.setName(name)
	    WINDOW_MEANS = scan.move(ln)
            out.printToFile(WINDOW_MEANS)
    
"""
for name, ln  in zip(NAMES, LENGTH):
    with FilePrinter(os.path.join(path, name)) as out:
	pos_start = 0
	pos_end = window_size
	while pos_end <= ln:
	   pos_start, pos_end = windowIterator(pos_start, pos_end, window_size)
           window_mean = covScanner(sam, mapq_cutoff, name, pos_start, pos_end)
	   out.printToFile(pos_start, pos_end, getLogratio(window_mean, name, chr_medians))
        
"""

"""for pile in sam.pileup("Chr4", 13000, 13100, truncate=True):
        cov = 0
	print pile.n, pile.pos, 
        for reads in pile.pileups:
            if reads.alignment.mapq >= mapq_cutoff:
                cov = cov + 1
         
        print cov
"""
