import pysam
import numpy
import sys
import os
import math
from optparse import OptionParser as opt
from subprocess import call

#Set the options that need to be set
prsr = opt()
prsr.add_option("-w", "--windowSize", dest="winsize", metavar="INT", default=500, help="Windowsize (bp) to be used to calculate log2ratio [Default:%default]")
prsr.add_option("-m", "--mappingQuality", dest="mapq", metavar="INT", default=0, help="Mapping quality cutoff for reads to be used in the calculation [Default:%default]")
prsr.add_option("-f", "--file", dest="bam", metavar="FILE", help="Input bam file to be analyzed, should be sorted and indexed")
prsr.add_option("-o", "--ouput", dest="path", metavar="PATH", default=os.getcwd() ,help="Output path")
prsr.add_option("-l", "--name-list", dest="order", metavar="FILE", default=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "chr.list"),
help="List of bam headers in order as they should be plotted, [Default:%default]")
prsr.add_option("-a", "--plot", dest="plot", metavar="BOOLEAN", default=True, help="Specify if plotting should be done using DNAcopy [Default:%default]")

# Get options
(options, args) = prsr.parse_args()
DEVNULL = open(os.devnull, 'wb')

def checkFile(test_file):
        if not os.path.isfile(test_file):
                quit("Could not find the file:" + test_file)

def checkValidArgs(options):
        if options.bam == None:
                quit("ERROR: No BAM file submitted")
        checkFile(options.bam)

def setAbsPath(options):
    options.path = os.path.abspath(options.path)
    return options.path   

def getNames(bam):
    header_dict = bam.header
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

def getMedianCov(bam, NAMES, LENGTH):
    medians_dict = {}
    for name, ln in zip(NAMES, LENGTH): 
       COVS = []
       for pile in bam.pileup(name, 0, ln):
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
        LOG2RATIOS.append(win_logratio)
    return LOG2RATIOS

class RWriter():
   
   def __init__(self, name_list):
      self.data = ()
      self.fh = None
      self.name = ()
      self.head = ()
      self.mid = ()
      self.merge = ()
      self.main = ()
      self.name_list = open(name_list, "r")
      self.counter = 0
      self.path = ()
      for name in self.name_list:
          self.counter = self.counter + 1
      self.name_list.seek(0)
   
   def getName(self):
       return self.name
 
   def setPath(self, path):
       self.path = path

   def writeHeader(self):
      HEAD = ["#!/usr/bin/env Rscript", "library(DNAcopy)"]
      self.head = '\n'.join(HEAD) + '\n'

   def writeMid(self):
      MID = []
      n = 1
      for name in self.name_list:
          name = name.rstrip()
	  MID.append("c%s <- read.table(\"%s\")" % (n, os.path.join(self.path, name)))
          n = n + 1
      n = 1
      self.name_list.seek(0)
      for name in self.name_list:
          name = name.rstrip()
          MID.append("c%s$chr <- (\"%s\")" % (n, name))      
          n = n + 1
      self.mid = '\n'.join(MID) + '\n'
      self.name_list.seek(0)

   def writeMerger(self):
      MERGE = []
      MERGE.append("merged <- rbind(")
      n = 1
      for name in self.name_list:
          name = name.rstrip()
          while n <= self.counter:
             if n == 1: 
	        MERGE.append("c")
		MERGE.append(str(n))
             else:
		MERGE.append(",c")
		MERGE.append(str(n))
             n = n + 1
      MERGE.append(")")
      MERGE.append("\npdf(\"%s\")" % os.path.join(self.path, "cnv_report.pdf"))
      self.merge = ''.join(MERGE) + '\n'
      self.name_list.seek(0)     
  
   def writeMain(self):
      MAIN = []
      MAIN.append("CNA.object <- CNA(merged$V1, merged$chr, merged$V2, data.type=(\"logratio\"), presorted=TRUE)")
      MAIN.append("CNA.smooth <- smooth.CNA(CNA.object)")
      MAIN.append("CNA.segm <- segment(CNA.smooth)")
      MAIN.append("plot(CNA.segm, plot.type=\"w\")")      
      MAIN.append("plot(CNA.segm, plot.type=\"s\")")
      MAIN.append("plot(CNA.segm, plot.type=\"p\")")
      n = 1
      for name in self.name_list:
         MAIN.append("c%s.object <- CNA(c%s$V1, c%s$chr, c%s$V2, data.type=(\"logratio\"), presorted=TRUE)" % (n,n,n,n))
         MAIN.append("c%s.smooth <- smooth.CNA(c%s.object)" % (n,n))
	 MAIN.append("c%s.segm <- segment(c%s.smooth)" % (n,n))
         MAIN.append("plot(c%s.segm, plot.type=\"s\")" % n)
	 n = n + 1
      MAIN.append("dev.off()")
      self.main = '\n'.join(MAIN)

   def setName(self, name):
      self.name = name

   def assembler(self):
      self.writeHeader()
      self.writeMid()
      self.writeMerger()
      self.writeMain()
      self.writer(self.name)
 
   def writer(self, name):
      self.fh = open(name, "w")
      self.fh.write(self.head)
      self.fh.write(self.mid)
      self.fh.write(self.merge)
      self.fh.write(self.main)
      self.fh.close()

class CovScanner():
  
   def __init__(self):
       self.window_size = ()
       self.pos_start = 0
       self.pos_end = ()
       self.bamfile = ()
       self.name = ()
       self.mapq_cutoff = ()
       self.MEANS = []
       self.POS = []
    
   def setWindowSize(self, win):
       self.window_size = win
       self.pos_end = win
   
   def setSamFile(self, bam):
       self.bamfile = bam
   
   def setMapq(self, mapq):
       self.mapq_cutoff = int(mapq)

   def setName(self, name):
       self.name = name
 
   def move(self, ln):
       while self.pos_end < ln:
          self.windowMean()
          self.POS.append(self.pos_end)
          self.pos_start = self.pos_start + self.window_size
          self.pos_end = self.pos_end + self.window_size
       return self.MEANS, self.POS
  
   def windowMean(self):
       bam = pysam.AlignmentFile(self.bamfile, "rb")
       WINDOW_COVS = []
       cov = ()
       for pile in bam.pileup(self.name, self.pos_start, self.pos_end, truncate=True):
           cov = 0
	   for reads in pile.pileups:
               if reads.alignment.mapq >= self.mapq_cutoff:
                   cov = cov + 1
	   WINDOW_COVS.append(cov)
       if sum(WINDOW_COVS) > 0:
           window_mean = numpy.mean(WINDOW_COVS)    
       else:
           window_mean = 0
       self.MEANS.append(window_mean)
       bam.close()

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

       for line in LIST:
           self._fh.write("%s\n" % line)
   
   def printZipListToFile(self, LIST1, LIST2):
       
       for line1, line2 in zip(LIST1, LIST2):
           self._fh.write("%s\t%s\n" % (line1, line2))    

""" START """

if __name__ == "__main__":
    
    options.path = setAbsPath(options)
    
    if not os.path.exists(options.path):
       os.makedirs(options.path)
 
    bam = pysam.AlignmentFile(options.bam, "rb")	
    NAMES, LENGTH = getNames(bam)
    
    if bool(options.plot) == True:
       rscripter = RWriter(options.order)
       rscripter.setName(os.path.join(options.path, "run_DNAcopy.r"))
       rscripter.setPath(options.path)
       rscripter.assembler()
 
    chr_medians = getMedianCov(bam, NAMES, LENGTH)

    for name, ln  in zip(NAMES, LENGTH):
        with FilePrinter(os.path.join(options.path, name)) as out:
            scan = CovScanner()
	    scan.setSamFile(options.bam)
	    scan.setMapq(options.mapq)
            scan.setWindowSize(options.winsize)
	    scan.setName(name)
	    WINDOW_MEANS, POS = scan.move(ln)
            WINDOW_MEANS = getLogratios(WINDOW_MEANS, name, chr_medians)
	    out.printZipListToFile(WINDOW_MEANS, POS)
 
    if bool(options.plot) == True:
       cmd = ["Rscript", rscripter.getName()] 
       call(cmd, stdout=DEVNULL, stderr=DEVNULL)
       rm_cmd = ["rm", rscripter.getName()]
       call(rm_cmd)

    bam.close()
    DEVNULL.close()
