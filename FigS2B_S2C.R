rm(list = ls())
graphics.off()
load("data4analy.Rdata")
load("data4analy.FUSCC.primary.Rdata")

index = permutations(nrow(onco_matrix_FUSCCmeta),2,rownames(onco_matrix_FUSCCmeta),repeats = T)

res.fisher = as.data.frame(matrix(nrow = nrow(index),ncol = 6))
colnames(res.fisher) = c("path1","path2","p.value","adj.p","OR","logOR")

for ( i in 1:nrow(res.fisher)){                                  
  if (index[i,][1] != index[i,][2]){                             
    mat = as.data.frame(t(onco_matrix_FUSCCmeta[index[i,],]))
    colnames(mat) = c("path1","path2")
    mytab = xtabs(~path1 + path2 , data = mat)
    if ( any(as.numeric(mytab) == 0)  ) { mytab = mytab +1 } 
    res = fisher.test(mytab)
    res.fisher[i,1] = index[i,][1]
    res.fisher[i,2] = index[i,][2]
    res.fisher[i,3] = res[["p.value"]]
    res.fisher[i,5] = res[["estimate"]][["odds ratio"]]
  } else{
    res.fisher[i,1] = index[i,][1]
    res.fisher[i,2] = index[i,][2]
  }
}
res.fisher[,4] = p.adjust(res.fisher[,3] ,method = "bonferroni")
res.fisher[,6] = log(res.fisher[,5] )  

## plot
input.data = res.fisher[,c(1,2,6)]
input.data = dcast(input.data , path1~path2)
rownames(input.data) = input.data$path1
input.data = input.data[,-1]
input.data = as.matrix(input.data)

input.pvalue = res.fisher[,c(1,2,3)]
input.pvalue = dcast(input.pvalue, path1 ~ path2)
rownames(input.pvalue) = input.pvalue$path1
input.pvalue = input.pvalue[,-1]
input.pvalue = as.matrix(input.pvalue)
rownames(input.pvalue) = colnames(input.pvalue)


pdf("./FigS2B.pdf")
corrplot(input.data, 
         type = "upper", 
         order="hclust",
         col=rev(brewer.pal(n=8, name="RdBu")),
         tl.col="black", 
         tl.cex = 0.5, 
         tl.srt = 90, 
         is.corr = FALSE, 
         diag = F,  
         p.mat = input.pvalue,
         sig.level = c(.001, .05), 
         insig = c("label_sig"), 
         pch.cex = 1, 
         font = 3) 
dev.off()

######################################################
######################################################

rm(list = ls())
graphics.off()
load("data4analy.Rdata")
load("data4analy.FUSCC.primary.Rdata")

index = permutations(nrow(onco_matrix_FUSCCmeta),2,rownames(onco_matrix_FUSCCmeta),repeats = T)

res.fisher = as.data.frame(matrix(nrow = nrow(index),ncol = 6))
colnames(res.fisher) = c("path1","path2","p.value","adj.p","OR","logOR")

for ( i in 1:nrow(res.fisher)){                                  
  if (index[i,][1] != index[i,][2]){                             
    mat = as.data.frame(t(onco_matrix_FUSCCmeta[index[i,],]))
    colnames(mat) = c("path1","path2")
    mytab = xtabs(~path1 + path2 , data = mat)
    if ( any(as.numeric(mytab) == 0)  ) { mytab = mytab +1 } 
    res = fisher.test(mytab)
    res.fisher[i,1] = index[i,][1]
    res.fisher[i,2] = index[i,][2]
    res.fisher[i,3] = res[["p.value"]]
    res.fisher[i,5] = res[["estimate"]][["odds ratio"]]
  } else{
    res.fisher[i,1] = index[i,][1]
    res.fisher[i,2] = index[i,][2]
  }
}
res.fisher[,4] = p.adjust(res.fisher[,3] ,method = "bonferroni")
res.fisher[,6] = log(res.fisher[,5] )  

## plot
input.data = res.fisher[,c(1,2,6)]
input.data = dcast(input.data , path1~path2)
rownames(input.data) = input.data$path1
input.data = input.data[,-1]
input.data = as.matrix(input.data)

input.pvalue = res.fisher[,c(1,2,3)]
input.pvalue = dcast(input.pvalue, path1 ~ path2)
rownames(input.pvalue) = input.pvalue$path1
input.pvalue = input.pvalue[,-1]
input.pvalue = as.matrix(input.pvalue)
rownames(input.pvalue) = colnames(input.pvalue)

pdf("./FigS2C.pdf")
corrplot(input.data, 
         type = "upper",
         order="hclust",
         col=rev(brewer.pal(n=8, name="RdBu")),
         tl.col="black",
         tl.cex = 0.5, 
         tl.srt = 90, 
         is.corr = FALSE,
         diag = F,  
         p.mat = input.pvalue, 
         sig.level = c(.001, .05),
         insig = c("label_sig"), 
         pch.cex = 1, 
         font = 3) 
dev.off()
