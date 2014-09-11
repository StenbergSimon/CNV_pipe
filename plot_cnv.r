#!/usr/bin/env Rscript

#Enable arguments to be passed
options(warn=-1)
args <-commandArgs(trailingOnly = TRUE)

suppressMessages(library("ggplot2"))
suppressMessages(library("plyr"))

sink(file = "cnv_pipe_temps/output", type = c("output", "message"))

#Load report file into x
x<- read.table(args[1])

#Calculate mean position of each window
x$meanpos <- rowMeans(subset(x, select = c(V2,V3)), na.rm = TRUE)

col_headings <- c('Log2_ratio', 'V2', 'V3', 'Cutoff', 'Position')
names(x) <- col_headings

wi <- (nrow(x) * 0.0154) + 8.1692

#cutoff <- data.frame( z = c(0,(max(x))), b = 0, cutoff = factor(0) )
#print (cutoff)
#Save as pdf
filename<-paste(args[1], '.pdf', sep="")
pdf(filename, height=9, width=wi)

#Plot the log2 values over the chr

plot <- ggplot(x, aes(x = Position, y = Log2_ratio, colour = Cutoff)) 
plot <- plot + ylim(-4,4)
plot <- plot + geom_point(stat="identity") 
#plot <- plot + geom_line(aes( z, b, linetype = cutoff ), cutoff)
plot <- plot + geom_segment(aes(x=0, xend=(max(x)), y = 0, yend = 0), colour = "black" )
plot <- plot + theme( legend.position = "none")
plot <- plot + scale_colour_gradientn(colours=c("red2","#132B43","red2"),limits=c(-1,1))
plot <- plot + theme(axis.title.x = element_text(colour="#132B43", size=15), axis.title.y = element_text(colour="#132B43", size=15)) + labs(y=expression(paste(Log [2]," ", "Ratio")), x="Position (bp)")
print(plot)

dev.off()
