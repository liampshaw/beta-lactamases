
library(ape)
library(ggplot2)
library(RColorBrewer)
suppressMessages(library(ggtree))
options(ignore.negative.edge=TRUE) # for ggtree
options(warn = -1)

args = commandArgs(trailingOnly=TRUE)

# Two arguments
# 1 - input alignment
# 2 - output pdf
# 3 width
# 4 height

source('treePlot.R')


pdf(file=args[2], width=args[3], height=args[4])
p.tree = treePlot(args[1])
p.tree
dev.off()
