import pysam
import numpy
import os
import math
from optparse import OptionParser as opt

samfile = "/Users/Simon/up_91_3.Founder_003_trimmed.sorted_rmdup_gc-corrected.bam"
mapq_cutoff = 1
window_size = 100
path = "/Users/Simon/output/"

"""
#Set the options that need to be set
prsr = opt()
prsr.add_option("-w", "--windowSize", dest="winsize", metavar="INT", default=500, help="Windowsize (bp) to be used to calculate log2ratio [Default:%default]")
prsr.add_option("-m", "--mappingQuality", dest="mapq", metavar="INT", default=0, help="Mapping quality cutoff for reads to be used in the calculation [Default:%default]")
prsr.add_option("-f", "--file", dest="bam", metavar="FILE", help="Input bam file to be analyzed")
prsr.add_option("-o", "--ouput", dest="path", metavar="PATH", default=os.getcwd() ,help="Output path")
prsr.add_option("-p", "--prefix", dest="prefix", metavar="DIR", default="Output_CNV_Caller", help="Specify the name of the output folder that will be created in output PATH [Default:%default]")

# Get options
(options, args) = prsr.parse_args()

def checkFile(test_file):
        if not os.path.isfile(test_file):
                quit("Could not find the file:" + test_file)

def checkValidArgs(options):
        if options.bam == None:
                quit("ERROR: No BAM file submitted")
        checkFile(options.bam)
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
    return NAMES, LENGTHS
    
def median(lst):
    return numpy.median(numpy.array(lst))

def getMedianCov(sam, NAMES, LENGTH, mapq_cutoff):
    medians_dict = {}
    for name, ln in zip(NAMES, LENGTH): 
       COVS = []
       for pile in sam.pileup(name, 0, ln):
           cov = pile.n
   	   COVS.append(cov)
       chr_median = median(COVS)
       medians_dict[name] = chr_median
    return medians_dict

def getLogratios(WINDOW_MEANS, name, chr_medians):
    median = chr_medians[name]
    LOG2RATIOS = []
    for window_mean in WINDOW_MEANS:
        win_logratio = window_mean / median
        if win_logratio == 0: 
           win_logratio = 0
        else:
           win_logratio = math.log(win_logratio, 2)
   #        print win_logratio, window_mean, median
        LOG2RATIOS.append(win_logratio)
    return LOG2RATIOS

def writeRscript(FILELIST):
    with open(os.path.join(path, "run_DNA_copy.r"), "w") as out:
        out.write("#!/usr/bin/env Rscript")
        out.write("library(DNAcopy)")
	p = 1
	for line in FILELIST:
	   out.write("".join(["c", str(p), "<- read.table(\"", line]))	

class RWriter():
   def __init__(self):
      self.data = ()
      
   def writeHeader(self):
      return None
       


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
       while self.pos_end < ln:
          self.windowMean()
          self.pos_start = self.pos_start + self.window_size
          self.pos_end = self.pos_end + self.window_size
       return self.MEANS
  
   def windowMean(self):
       sam = pysam.AlignmentFile(samfile, "rb")
       WINDOW_COVS = []
       cov = ()
       for pile in sam.pileup(self.name, self.pos_start, self.pos_end, truncate=True):
           cov = 0
	   for reads in pile.pileups:
               if reads.alignment.mapq >= mapq_cutoff:
                   cov = cov + 1
	   WINDOW_COVS.append(cov)
       window_mean = numpy.mean(WINDOW_COVS)    
       self.MEANS.append(window_mean)
       sam.close()

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

   def printListToFile(self, LIST):

       if (self._fh is None):
           raise Exception("Only use within with statement blocks")
       for line in LIST:
           self._fh.write("%s\n" % line)
       

""" START """
if __name__ == "__main__":
    
    bam = pysam.AlignmentFile(samfile, "rb")	
    NAMES, LENGTH = getNames(bam)
    chr_medians = getMedianCov(bam, NAMES, LENGTH, mapq_cutoff)

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
            WINDOW_MEANS = getLogratios(WINDOW_MEANS, name, chr_medians)
	    out.printListToFile(WINDOW_MEANS)

    bam.close()    
