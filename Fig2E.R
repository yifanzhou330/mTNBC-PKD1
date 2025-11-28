rm(list = ls())
graphics.off()
load("data4analy.Rdata")
library(tidyverse)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(circlize)

clinical.ori.f = clinical.ori[,c("Tumor_Sample_Barcode","Biopsied_met4plot")]
recur.mut =  c("TP53","PIK3CA","KMT2D","PTEN","RB1","ARID1A","ERBB3","NOTCH1")

somatic.mtx.f = as.data.frame(t(somatic.mtx[recur.mut,]))

data4chord = dcast(data4sankey[,c(1,2)],meta.site~mut.gene)
rownames(data4chord) = data4chord[,1] ; data4chord = data4chord[,-1]


grid.col = NULL
grid.col[colnames(data4chord)]= rev(brewer.pal(8,name = "Set1"))
grid.col[rownames(data4chord)]= brewer.pal(6,name = "Set2")

pdf(file = "./chord_ori1.pdf",width = 10)
circos.par(gap.degree=c(rep(2,nrow(data4chord)-1),10, rep(2,ncol(data4chord)-1),10),start.degree=180)
chordDiagram(as.matrix(data4chord),
             #directional=TRUE,
             diffHeight = 0.06,
             grid.col = grid.col, 
             transparency = 0.2,
             annotationTrack = c("name","grid"),
             annotationTrackHeight = c(0.05,0.12))

legend("right",pch=20,legend=colnames(data4chord),
       col=grid.col[colnames(data4chord)],bty="n",
       cex=1,pt.cex=3,border="black")
circos.clear()
dev.off()


