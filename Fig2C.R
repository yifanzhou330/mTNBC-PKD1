rm(list = ls())
graphics.off()
load("data4analy.Rdata")
load("cnv_merge.Rdata")
library(ComplexHeatmap)

cnv = tmp
tmp = intersect(colnames(cnv), clinical.ori$Tumor_Sample_Barcode)
cnv = cnv[,tmp]

for (i in colnames(cnv)){
  mat = cnv[,i]
  mat[mat == "neutral"] = ""
  cnv[,i] = mat
}


wid = 1.2
hei = 1
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = "#CCCCCC", col = "#CCCCCC"))
  },
  del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["del"], col = col["del"]))
  },
  amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["amp"], col = col["amp"]))
  }
)
col = c("del" = "#366A9C", "amp" = "#D0342B")

gene_ord = gene_ord[cnvlist]

pdf(file = "./Fig2C.pdf",height = 4.5,width = 9.5)
oncoPrint(cnv[names(gene_ord)[1:10],] , alter_fun = alter_fun,col = col)
dev.off()




