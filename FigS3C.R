rm(list = ls())
graphics.off()
pathway.gene = fread(file = "pathway.gene.csv",data.table = F)

load("data4analy.Rdata")
oncoplot(maf = somatic.maf, genes = pathway.gene$HUGO_Symbol,
         writeMatrix = T,removeNonMutated = F,
         top = 100)
onco_matrix = read.table(file = "./onco_matrix.txt", sep = "\t",
                         na.strings = c("",0))
new_onco_matrix = as.data.frame(matrix(nrow = length(unique(pathway.gene$`24 Categories`)),
                                       ncol = ncol(onco_matrix))
)
rownames(new_onco_matrix) = unique(pathway.gene$`24 Categories`)
colnames(new_onco_matrix) = colnames(onco_matrix)

for ( i in rownames(new_onco_matrix) ){
  gene = pathway.gene$HUGO_Symbol[pathway.gene$`24 Categories` == i]
  mat = onco_matrix[gene,]
  
  for ( j in 1:ncol(mat)){
    if ( all(is.na(unique(mat[,j])) )  ){ new_onco_matrix[i,j] = NA
    } else {
      
      if (any ( is.na(unique(mat[,j])) ) & length(unique(mat[,j])) == 2){
        new_onco_matrix[i,j] = unique(mat[,j])[ !is.na(  unique(mat[,j])  )  ]
      } else {
        new_onco_matrix[i,j] = "Multi_Hit"
      }}}
}

new_onco_matrix_FUSCCmeta = new_onco_matrix

col = c("Missense_Mutation" = "#366A9C", "Nonsense_Mutation" = "#BBDE93",
        "In_Frame_Del" = "#EE8632" ,"In_Frame_Ins" = "#D0342B",
        "Multi_Hit" = "black", "Splice_Site" = "#F3C17B", "Frame_Shift_Del" = "#AECDE1",
        "Frame_Shift_Ins" = "#ED9E9B",
        "Nonstop_Mutation" = "#339900")

wid = 0.89
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
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "mm"), h-unit(hei, "mm"),
              gp = gpar(fill = col["Nonstop_Mutation"], col = col["Nonstop_Mutation"]))
  }
)

tmp = pathway.gene[,c(2,3)][!duplicated(pathway.gene[,c(2,3)]$`24 Categories`),]
annot4row = data.frame(row.names = tmp$`24 Categories`, pathway.categories = tmp$`8 Categories` )
col.annot = brewer.pal(8,name = "Set3") ; names(col.annot) = unique(annot4row$pathway.categories)
col.annot = list(pathway.categories = col.annot)
rAno = rowAnnotation(df = annot4row, col =  col.annot,
                     show_annotation_name = F) 
                     
pdf(file = paste0("./FigS3C_1.pdf"),height = 4.5,width = 9.5)
oncoPrint(new_onco_matrix,alter_fun = alter_fun,col = col,
          remove_empty_rows = F,
          remove_empty_columns = F,
          left_annotation = rAno
          #heatmap_height = unit(1, "npc")
          )
dev.off()
  
tmp = oncoPrint(new_onco_matrix,alter_fun = alter_fun,col = col, 
                #remove_empty_rows = F
                remove_empty_columns = F)
path.order = rownames(new_onco_matrix)[tmp@row_order]

###############
rm(list = ls()[ls() != "path.order"])
graphics.off()
load("./data4analy.mskcc.Rdata")

oncoplot(maf = msk.maf, genes = pathway.gene$HUGO_Symbol,
         writeMatrix = T,removeNonMutated = F,
         top = 100)

onco_matrix = read.table(file = "./onco_matrix.txt", sep = "\t",
                         na.strings = c("",0))
new_onco_matrix = as.data.frame(matrix(nrow = length(unique(pathway.gene$`24 Categories`)),
                                       ncol = ncol(onco_matrix))
)
rownames(new_onco_matrix) = unique(pathway.gene$`24 Categories`)
colnames(new_onco_matrix) = colnames(onco_matrix)

for ( i in rownames(new_onco_matrix) ){
  gene = pathway.gene$HUGO_Symbol[pathway.gene$`24 Categories` == i]
  mat = onco_matrix[gene,]
  
  for ( j in 1:ncol(mat)){
    if ( all(is.na(unique(mat[,j])) )  ){ new_onco_matrix[i,j] = NA
    } else {
      
      if (any ( is.na(unique(mat[,j])) ) & length(unique(mat[,j])) == 2){
        new_onco_matrix[i,j] = unique(mat[,j])[ !is.na(  unique(mat[,j])  )  ]
      } else {
        new_onco_matrix[i,j] = "Multi_Hit"
      }}}
}


col = c("Missense_Mutation" = "#366A9C", "Nonsense_Mutation" = "#BBDE93",
        "In_Frame_Del" = "#EE8632" ,"In_Frame_Ins" = "#D0342B",
        "Multi_Hit" = "black", "Splice_Site" = "#F3C17B", "Frame_Shift_Del" = "#AECDE1",
        "Frame_Shift_Ins" = "#ED9E9B")
wid = 0.89
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

tmp = pathway.gene[,c(2,3)][!duplicated(pathway.gene[,c(2,3)]$`24 Categories`),]
annot4row = data.frame(row.names = tmp$`24 Categories`, pathway.categories = tmp$`8 Categories` )
col.annot = brewer.pal(8,name = "Set3") ; names(col.annot) = unique(annot4row$pathway.categories)
col.annot = list(pathway.categories = col.annot)
rAno = rowAnnotation(df = annot4row, col =  col.annot,
                     show_annotation_name = F) 

pdf(file = paste0("./FigS3C_2.pdf"),height =4.5,width = 7.2)
oncoPrint(new_onco_matrix,alter_fun = alter_fun,col = col,
          remove_empty_columns = F, remove_empty_rows = F,
          row_order = path.order,
          left_annotation = rAno
          #heatmap_height = unit(0.1, "points")
          )
dev.off()
