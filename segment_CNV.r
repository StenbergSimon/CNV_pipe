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

pdf("cnv_report.pdf")
CNA.object <- CNA(merged$V1, merged$chr, merged$V2, data.type=("logratio"), presorted=TRUE)
CNA.smooth <- smooth.CNA(CNA.object)
CNA.segm <- segment(CNA.smooth)

plot(CNA.segm, plot.type="w")
plot(CNA.segm, plot.type="s")
plot(CNA.segm, plot.type="p")

c1.object <- CNA(c1$V1, c1$chr, c1$V2, data.type=("logratio"), presorted=TRUE)
c1.smooth <- smooth.CNA(c1.object)
c1.segm <- segment(c1.smooth)

plot(c1.segm, plot.type="s")

c2.object <- CNA(c2$V1, c2$chr, c2$V2, data.type=("logratio"), presorted=TRUE)
c2.smooth <- smooth.CNA(c2.object)
c2.segm <- segment(c2.smooth)

plot(c2.segm, plot.type="s")

c3.object <- CNA(c3$V1, c3$chr, c3$V2, data.type=("logratio"), presorted=TRUE)
c3.smooth <- smooth.CNA(c3.object)
c3.segm <- segment(c3.smooth)

plot(c3.segm, plot.type="s")

c4.object <- CNA(c4$V1, c4$chr, c4$V2, data.type=("logratio"), presorted=TRUE)
c4.smooth <- smooth.CNA(c4.object)
c4.segm <- segment(c4.smooth)

plot(c4.segm, plot.type="s")

c5.object <- CNA(c5$V1, c5$chr, c5$V2, data.type=("logratio"), presorted=TRUE)
c5.smooth <- smooth.CNA(c5.object)
c5.segm <- segment(c5.smooth)

plot(c5.segm, plot.type="s")

c6.object <- CNA(c6$V1, c6$chr, c6$V2, data.type=("logratio"), presorted=TRUE)
c6.smooth <- smooth.CNA(c6.object)
c6.segm <- segment(c6.smooth)

plot(c6.segm, plot.type="s")

c7.object <- CNA(c7$V1, c7$chr, c7$V2, data.type=("logratio"), presorted=TRUE)
c7.smooth <- smooth.CNA(c7.object)
c7.segm <- segment(c7.smooth)

plot(c7.segm, plot.type="s")

c8.object <- CNA(c8$V1, c8$chr, c8$V2, data.type=("logratio"), presorted=TRUE)
c8.smooth <- smooth.CNA(c8.object)
c8.segm <- segment(c8.smooth)

plot(c8.segm, plot.type="s")

c9.object <- CNA(c9$V1, c9$chr, c9$V2, data.type=("logratio"), presorted=TRUE)
c9.smooth <- smooth.CNA(c9.object)
c9.segm <- segment(c9.smooth)

plot(c9.segm, plot.type="s")

c10.object <- CNA(c10$V1, c10$chr, c10$V2, data.type=("logratio"), presorted=TRUE)
c10.smooth <- smooth.CNA(c10.object)
c10.segm <- segment(c10.smooth)

plot(c10.segm, plot.type="s")

c11.object <- CNA(c11$V1, c11$chr, c11$V2, data.type=("logratio"), presorted=TRUE)
c11.smooth <- smooth.CNA(c11.object)
c11.segm <- segment(c11.smooth)

plot(c11.segm, plot.type="s")

c12.object <- CNA(c12$V1, c12$chr, c12$V2, data.type=("logratio"), presorted=TRUE)
c12.smooth <- smooth.CNA(c12.object)
c12.segm <- segment(c12.smooth)

plot(c12.segm, plot.type="s")

c13.object <- CNA(c13$V1, c13$chr, c13$V2, data.type=("logratio"), presorted=TRUE)
c13.smooth <- smooth.CNA(c13.object)
c13.segm <- segment(c13.smooth)

plot(c13.segm, plot.type="s")

c14.object <- CNA(c14$V1, c14$chr, c14$V2, data.type=("logratio"), presorted=TRUE)
c14.smooth <- smooth.CNA(c14.object)
c14.segm <- segment(c14.smooth)

plot(c14.segm, plot.type="s")

c15.object <- CNA(c15$V1, c15$chr, c15$V2, data.type=("logratio"), presorted=TRUE)
c15.smooth <- smooth.CNA(c15.object)
c15.segm <- segment(c15.smooth)

plot(c15.segm, plot.type="s")

c16.object <- CNA(c16$V1, c16$chr, c16$V2, data.type=("logratio"), presorted=TRUE)
c16.smooth <- smooth.CNA(c16.object)
c16.segm <- segment(c16.smooth)

plot(c16.segm, plot.type="s")

dev.off()








