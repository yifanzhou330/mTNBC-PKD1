rm(list = ls())
graphics.off()
library(RColorBrewer)
library(maftools)
load("data4analy.Rdata")
load("data4analy.mskcc.Rdata")

color.re3.Met.num4plot = c("#A52B29","#DBBEBE") #
names(color.re3.Met.num4plot) = c(">= 3", "<= 2")

color.Biopsied_met4plot = brewer.pal(n = 6 ,"Set2") 
names(color.Biopsied_met4plot) = c("Breast", "Chest_wall", "Liver", "Lung", "Lymph_nodes", "Others")

color.DFS4plot = c("#D1A0C7","#E0C7DC","#E6DBE4") #pink
names(color.DFS4plot) = c(">3y", "1-3y", "<1y")

colors.Lines_of_CT4plot = c("#45749E","#5489BA","#80A7D0","#B4C9E2","#C2D9F7")
names(colors.Lines_of_CT4plot) = c("> 3","3", "2", "1","0")


colors.ls = list(Lines_of_CT4plot = colors.Lines_of_CT4plot,
                 re3.Met.num4plot = color.re3.Met.num4plot,
                 DFS4plot = color.DFS4plot,
                 Biopsied_met4plot = color.Biopsied_met4plot)

vc_cols = c("#E31F28","#4B718B","#4FAF4B",
            "black",
            "#F4801F","#F5ED3A","#A55627","#F084B5")
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
                   'Multi_Hit',
                   'Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del')


pdf("Fig2A.pdf",height = 8, width = 12)
oncoplot(maf = somatic.maf, annotationColor = colors.ls, colors = vc_cols,
         sortByAnnotation = F, 
         clinicalFeatures = c("Biopsied_met4plot","Lines_of_CT4plot","re3.Met.num4plot",
                              "DFS4plot"),
         removeNonMutated = F,draw_titv = TRUE, bgCol = "#F4EFE6",borderCol = NULL,
         writeMatrix = T
         )
dev.off()


pdf("FigS2C.pdf",height = 8, width = 12)
oncoplot(maf = somatic.maf, annotationColor = colors.ls, colors = vc_cols,
         groupAnnotationBySize = F,
         sortByAnnotation = T,
         clinicalFeatures = c("Biopsied_met4plot"),
         removeNonMutated = F,draw_titv = TRUE, bgCol = "#F4EFE6",borderCol = NULL,
         writeMatrix = T
)
dev.off()

pdf("FigS3A.pdf",height = 8, width = 12)
oncoplot(maf = msk.maf, 
         #annotationColor = colors.ls, 
         colors = vc_cols,
         sortByAnnotation = T, 
         clinicalFeatures = c("Metastatic.Site"),
         removeNonMutated = F,draw_titv = TRUE, bgCol = "#F4EFE6",borderCol = NULL,
         writeMatrix = T
)
dev.off()

#############################################
#############################################
rm(list = ls())
graphics.off()
load("data4analy.Rdata")

mut.tmp = cbind(somatic.data,mut.ID = paste(somatic.data$Hugo_Symbol,somatic.data$aaChange, sep = ":") )
hotspot = as.data.frame(table(mut.tmp$mut.ID))
hotspot = hotspot[order(hotspot$Freq,decreasing = T),]
hotspot.top = as.character(hotspot$Var1[c(1,3:14)]) 

hotspot.mat = matrix(ncol = length(unique(somatic.data$Tumor_Sample_Barcode)), nrow = length(hotspot.top))
colnames(hotspot.mat) = unique(somatic.data$Tumor_Sample_Barcode)
rownames(hotspot.mat) = hotspot.top
hotspot.mat = as.data.frame(hotspot.mat)

for (i in hotspot.top) {
  tmp = mut.tmp[mut.tmp$mut.ID %in% i,]
  hotspot.mat[i,match(tmp$Tumor_Sample_Barcode, colnames(hotspot.mat))] = tmp$Variant_Classification
}


wid = 1
hei = 1
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = "#CCCCCC", col = "#CCCCCC"))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Missense_Mutation"], col = col["Missense_Mutation"]))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Nonsense_Mutation"], col = col["Nonsense_Mutation"]))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["In_Frame_Del"], col = col["In_Frame_Del"]))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["In_Frame_Ins"], col = col["In_Frame_Ins"]))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Multi_Hit"], col = col["Multi_Hit"]))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Splice_Site"], col = col["Splice_Site"]))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Frame_Shift_Del"], col = col["Frame_Shift_Del"]))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Frame_Shift_Ins"], col = col["Frame_Shift_Ins"]))
  }
)
col = c("Missense_Mutation" = "#366A9C", "Nonsense_Mutation" = "#BBDE93",
        "In_Frame_Del" = "#EE8632" ,"In_Frame_Ins" = "#D0342B",
        "Multi_Hit" = "black", "Splice_Site" = "#F3C17B", "Frame_Shift_Del" = "#AECDE1",
        "Frame_Shift_Ins" = "#ED9E9B")

hotspot.mat = as.matrix(hotspot.mat)
pdf(file = "Fig2B.pdf",height = 4.5,width = 9.5)
oncoPrint(hotspot.mat,alter_fun = alter_fun,col = col,column_order = sample.ord,)
dev.off()
