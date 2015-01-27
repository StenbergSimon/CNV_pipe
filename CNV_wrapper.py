#!/usr/bin/env python 

import os
import sys
from optparse import OptionParser as opt
import subprocess

#Set the options that need to be set
prsr = opt()
prsr.add_option("-w", "--windowSize", dest="winsize", metavar="INT", default="500", help="Windowsize (bp) to be used to calculate log2ratio [Default:%default]")
prsr.add_option("-m", "--mappingQuality", dest="mapq", metavar="INT", default="0", help="Mapping quality cutoff for reads to be used in the calculation [Default:%default]")
prsr.add_option("-c", "--cutOff", dest="co", metavar="INT", default="0.5", help="Log2Ratio cut-off to be used to mark significant hits in the plots [Default:%default]")
prsr.add_option("-f", "--file", dest="file", metavar="FILE", help="Input bam file to be analyzed")
prsr.add_option("-o", "--ouput", dest="out", metavar="FILE", default="cnv_report.tsv",help="Report file to write output to")

# Get options
(options, args) = prsr.parse_args()

if options.file == None:
	quit("ERROR: No file submitted")

def setDefaults(options):
	if options.winsize == None:
		options.winsize = 500
	if options.mapq == None:
		options.mapq = 0
	if options.co == None:
		options.co = 0.5
	if options.out == None:
		options.out = "cnv_report.tsv"
	return (options.winsize, options.mapq, options.co, options.out)

def writeRscript(options):
	print

""" PROGRAM START """
winsize, mapq, co, out= setDefaults(options)

#Get folder name

folder=os.path.dirname(options.file)


#CNV ouputfile:
cnv_f = open(out, "w")
#Call CNV_pipe with the options - save to cnv_file. Options are called from a list
cnv_cmd = ["cnv_pipe.pl", options.file, "-w", winsize, "-mapq", mapq, "-co", co, "-rall"]
p = subprocess.Popen(cnv_cmd, cwd=folder, stdout = cnv_f)
p.communicate()
 
