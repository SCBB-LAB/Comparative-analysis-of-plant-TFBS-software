
args <- commandArgs(trailingOnly = TRUE)
library(seqLogo)
m <- read.table(args[1],sep=",",header=F,row.names=c(1))
p <- makePWM(m)
pdf(args[2],width=5,height=4)
seqLogo(p, ic.scale=FALSE)
dev.off()

