#!/usr/bin/env Rscript

library(DNAcopy)

options(warn=-1)
args <-commandArgs(trailingOnly = TRUE)

setwd(args[1])

c1 <- read.table("ChrI")
c2 <- read.table("ChrII")
c3 <- read.table("ChrIII")
c4 <- read.table("ChrIV")
c5 <- read.table("ChrV")
c6 <- read.table("ChrVI")
c7 <- read.table("ChrVII")
c8 <- read.table("ChrVIII")
c9 <- read.table("ChrIX")
c10 <- read.table("ChrX")
c11 <- read.table("ChrXI")
c12 <- read.table("ChrXII")
c13 <- read.table("ChrXIII")
c14 <- read.table("ChrXIV")
c15 <- read.table("ChrXV")
c16 <- read.table("ChrXVI")

c1$chr <- 1
c2$chr <- 2
c3$chr <- 3
c4$chr <- 4
c5$chr <- 5
c6$chr <- 6
c7$chr <- 7
c8$chr <- 8
c9$chr <- 9
c10$chr <- 10
c11$chr <- 11
c12$chr <- 12
c13$chr <- 13
c14$chr <- 14
c15$chr <- 15
c16$chr <- 16


merged <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16)

CNA.object <- CNA(merged$V1, merged$chr, merged$V2, data.type=("logratio"), presorted=TRUE)
CNA.smooth <- smooth.CNA(CNA.object)
CNA.segm <- segment(CNA.smooth)

pdf("cnv_report.pdf")

plot(CNA.segm, plot.type="w")
plot(CNA.segm, plot.type="s")
plot(CNA.segm, plot.type="p")

dev.off()








