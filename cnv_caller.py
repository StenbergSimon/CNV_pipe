
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
prsr.add_option("-r", "--reference", dest="ref", metavar="FILE", help="Bam file to be used as refernce / control")
prsr.add_option("-s", "--zstart", dest="zstart", metavar="INT", help="Zoom: Start chromosomal location")
prsr.add_option("-e", "--zend", dest="zend", metavar="INT", help="Zoom: End chromosomal location")
prsr.add_option("-c", "--zchrom", dest="zchrom", metavar="STR", help="Zoom: Chromosome")
prsr.add_option("-z", "--zoom", dest="z", action="store_true", help="Runs in zoom-mode on a prerun project")

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

def getNormalizer(bam, ref, NAMES, LENGTH):
    REF = []
    BAM = []
    for name, ln in zip(NAMES,LENGTH):
        for bam_pile in bam.pileup(name, 0, ln):
            BAM.append(bam_pile.n)
        for ref_pile in ref.pileup(name, 0, ln):    
            REF.append(ref_pile.n)
    normalizer = float(sum(REF)) / float(sum(BAM))
    return normalizer

class RWriter():
   
   def __init__(self, options):
      self.options = options
      self.data = ()
      self.fh = None
      self.name = ()
      self.head = ()
      self.mid = ()
      self.merge = ()
      self.main = ()
      self.name_list = open(options.order, "r")
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
      print "writing R script..."
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
      MERGE.append("CNA.object <- CNA(merged$V1, merged$chr, merged$V2, data.type=(\"logratio\"), presorted=TRUE)")
      MERGE.append("CNA.smooth <- smooth.CNA(CNA.object)")
      MERGE.append("CNA.segm <- segment(CNA.smooth)")
      self.merge = ''.join(MERGE) + '\n'
      self.name_list.seek(0)     
  
   def writeMain(self):
      MAIN = []
      MAIN.append("plot(CNA.segm, plot.type=\"w\")")      
      MAIN.append("plot(CNA.segm, plot.type=\"s\")")
      MAIN.append("plot(CNA.segm, plot.type=\"p\")")
      MAIN.append("write.table(segments.summary(CNA.segm), file = \"%s\", sep=\"\\t\")" % os.path.join(self.path, "segments_summary.tsv"))
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
      
   def writeZoom(self):
      ZOOM = []
      ZOOM.append("\npdf(\"%s\")" % os.path.join(self.path, str(options.zchrom) + "_" + str(options.zstart) + "_" + str(options.zend) + ".pdf"))
      ZOOM.append("zoomIntoRegion(CNA.segm, chrom=%s, maploc.start=%s, maploc.end=%s, sampleid=\"Sample.1\"") % options.zchrom, options.zstart, options.zend 
      ZOOM.append("dev.off()")
      self.zoom = ''.join(ZOOM) + '\n'

   def aseembleZoom(self):
       self.writeHeader()
       self.writeMid()
       self.writeMerger()
       self.writeZoom()
   
   def zoomWriter(self, name):
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
       self.ref = None
       self.name = ()
       self.mapq_cutoff = ()
       self.RATIOS = []
       self.POS = []
       self.medians = None
       self.normalizer = ()
    
   def setWindowSize(self, win):
       self.window_size = win
       self.pos_end = win
 
   def setChrMedians(self, medians):
       self.medians = medians
   
   def setSamFile(self, bam, ref):
       self.bamfile = bam
       self.ref = ref
   
   def setNormalizer(self, normalizer):
       self.normalizer = normalizer
       
   def setMapq(self, mapq):
       self.mapq_cutoff = int(mapq)

   def setName(self, name):
       self.name = name
 
   def move(self, ln):
       while self.pos_end < ln:
          window_mean = self.windowMean(self.bamfile)
          self.POS.append(self.pos_end)
          self.RATIOS.append(self.getLogratios(window_mean))
          self.pos_start = self.pos_start + self.window_size
          self.pos_end = self.pos_end + self.window_size
       return self.RATIOS, self.POS
  
   def getLogratios(self, window_mean):
       if self.ref == None:
          logratio = window_mean / self.medians[self.name]
       else:
          win_mean_ref = self.windowMean(self.ref)
          if win_mean_ref == 0:
             win_mean_ref = 1
          logratio  = (window_mean * self.normalizer) / win_mean_ref
       if logratio == 0:
          logratio = 0
       else:
          logratio = math.log(logratio, 2)
       return logratio

   def windowMean(self, bam):
       bam = pysam.AlignmentFile(bam, "rb")
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
       bam.close()
       return window_mean


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
    normalizer = 0 
    if not os.path.exists(options.path):
       os.makedirs(options.path)
 
    bam = pysam.AlignmentFile(options.bam, "rb")	
    NAMES, LENGTH = getNames(bam)
    
    if bool(options.plot) == True:
       rscripter = RWriter(options)
       rscripter.setName(os.path.join(options.path, "run_DNAcopy.r"))
       rscripter.setPath(options.path)
       rscripter.assembler()
 
    if options.ref == None:
       chr_medians = getMedianCov(bam, NAMES, LENGTH)
    else:
       ref = pysam.AlignmentFile(options.ref, "rb")
       normalizer = getNormalizer(bam, ref, NAMES, LENGTH)     
    print "calculating log ratios..." 
    for name, ln  in zip(NAMES, LENGTH):
        with FilePrinter(os.path.join(options.path, name)) as out:
            scan = CovScanner()
            if options.ref == None:
               scan.setChrMedians(chr_medians)
	    scan.setSamFile(options.bam, options.ref)
	    scan.setMapq(int(options.mapq))
            scan.setWindowSize(int(options.winsize))
            scan.setNormalizer(normalizer)
	    scan.setName(name)
	    WINDOW_RATIOS, POS = scan.move(ln)
	    out.printZipListToFile(WINDOW_RATIOS, POS)
 
    if bool(options.plot) == True:
       print "running plot scripts..."
       cmd = ["Rscript", rscripter.getName()] 
       call(cmd, stdout=DEVNULL, stderr=DEVNULL)
       rm_cmd = ["rm", rscripter.getName()]
       call(rm_cmd)

    bam.close()
    DEVNULL.close()
