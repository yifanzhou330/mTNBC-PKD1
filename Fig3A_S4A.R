rm(list = ls())
graphics.off()
load("./data4analy.FUSCC.primary.Rdata")
load("./data4analy.Rdata")

common.gene.pri = intersect( unique(fuscc.pri.mut$Hugo_Symbol) , unique(somatic.data$Hugo_Symbol) )
somatic.maf.pri.f = read.maf(somatic.data[somatic.data$Hugo_Symbol %in% common.gene.pri,])
fuscc.pri.mut.f = read.maf(fuscc.pri.mut[fuscc.pri.mut$Hugo_Symbol %in% common.gene.pri,]) 

common.gene.pri = intersect(unique(somatic.maf.pri.f@data$Hugo_Symbol),unique(fuscc.pri.mut.f@data$Hugo_Symbol))
tmp1 = somatic.maf.pri.f@data
tmp2 = fuscc.pri.mut.f@data
somatic.maf.pri.f = read.maf(tmp1[tmp1$Hugo_Symbol %in% common.gene.pri,])
fuscc.pri.mut.f = read.maf(tmp2[tmp2$Hugo_Symbol %in% common.gene.pri,]) 

somatic.maf.pri.f.mtx = mutCountMatrix(somatic.maf.pri.f, includeSyn = F, 
                                       countOnly = NULL, removeNonMutated = F) 
somatic.maf.pri.f.mtx = as.data.frame(apply(somatic.maf.pri.f.mtx, 2, function(x){ifelse(x == 0, 0, 1)} ))

fuscc.pri.mut.f.mtx = mutCountMatrix(fuscc.pri.mut.f, includeSyn = F, 
                               countOnly = NULL, removeNonMutated = F) 
fuscc.pri.mut.f.mtx = as.data.frame(apply(fuscc.pri.mut.f.mtx, 2, function(x){ifelse(x == 0, 0, 1)} ))
fuscc.pri.mut.f.mtx = fuscc.pri.mut.f.mtx[rownames(somatic.maf.pri.f.mtx),]

priVSmeta.res = as.data.frame(matrix(nrow = nrow(somatic.maf.pri.f.mtx), ncol = 4))
rownames(priVSmeta.res) = row.names(somatic.maf.pri.f.mtx)
colnames(priVSmeta.res) = c("FUSCC.pri","FUSCC.meta","p.val","adj.p")

priVSmeta.res$FUSCC.pri = apply(fuscc.pri.mut.f.mtx,1,function(x){sum(x)/length(x)})
priVSmeta.res$FUSCC.meta = apply(somatic.maf.pri.f.mtx,1,function(x){sum(x)/length(x)})

for ( i in rownames(priVSmeta.res)){
  tmp1 = t(rbind(rep("FUSCC.pri",ncol(fuscc.pri.mut.f.mtx)),fuscc.pri.mut.f.mtx[i,]))
  tmp2= t(rbind(rep("FUSCC.meta",ncol(somatic.maf.pri.f.mtx)),somatic.maf.pri.f.mtx[i,]))
  tmp = as.data.frame(rbind(tmp1,tmp2))
  colnames(tmp) = c("cohort","gene")
  tmp$cohort = as.factor(tmp$cohort)
  tmp$gene = as.factor(tmp$gene)
  p = chisq_test(cohort~gene,data= tmp, distribution=approximate())
  priVSmeta.res[i,3] <- pvalue(p)[1]
}
priVSmeta.res[,4] = p.adjust(priVSmeta.res[,3],method = "fdr")

ggplot(priVSmeta.res,aes(x=FUSCC.meta, y = FUSCC.pri)) + 
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.5, color = '#747475')+
  geom_point(alpha = 1, size = 2.5, col="#76A1C5") + 
  geom_point(data = priVSmeta.res[priVSmeta.res$p.val < 0.05,] , 
             alpha = 1, size = 2.5, col="red") +
  xlim(0,0.1) + ylim(0,0.1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  theme_bw()+
  geom_text_repel(data= priVSmeta.res[priVSmeta.res$p.val < 0.05,] ,
                  aes(label=rownames(priVSmeta.res[priVSmeta.res$p.val < 0.05,])),
                  col="black",alpha = 1,size=3) 
  ggsave(filename = "./Fig3A.pdf")

##########################################
rm(list = ls())
graphics.off()
load("./data4analy.FUSCC.primary.Rdata")
load("./data4analy.Rdata")
sample =clinical.ori$Tumor_Sample_Barcode[clinical.ori$Lines_of_CT == 0]
somatic.data = somatic.data[somatic.data$Tumor_Sample_Barcode %in% sample , ]

common.gene.pri = intersect( unique(fuscc.pri.mut$Hugo_Symbol) , unique(somatic.data$Hugo_Symbol) )
somatic.maf.pri.f = read.maf(somatic.data[somatic.data$Hugo_Symbol %in% common.gene.pri,])
fuscc.pri.mut.f = read.maf(fuscc.pri.mut[fuscc.pri.mut$Hugo_Symbol %in% common.gene.pri,]) 

common.gene.pri = intersect(unique(somatic.maf.pri.f@data$Hugo_Symbol),unique(fuscc.pri.mut.f@data$Hugo_Symbol))
tmp1 = somatic.maf.pri.f@data
tmp2 = fuscc.pri.mut.f@data
somatic.maf.pri.f = read.maf(tmp1[tmp1$Hugo_Symbol %in% common.gene.pri,])
fuscc.pri.mut.f = read.maf(tmp2[tmp2$Hugo_Symbol %in% common.gene.pri,]) 

somatic.maf.pri.f.mtx = mutCountMatrix(somatic.maf.pri.f, includeSyn = F, 
                                       countOnly = NULL, removeNonMutated = F) 
somatic.maf.pri.f.mtx = as.data.frame(apply(somatic.maf.pri.f.mtx, 2, function(x){ifelse(x == 0, 0, 1)} ))

fuscc.pri.mut.f.mtx = mutCountMatrix(fuscc.pri.mut.f, includeSyn = F, 
                                     countOnly = NULL, removeNonMutated = F) 
fuscc.pri.mut.f.mtx = as.data.frame(apply(fuscc.pri.mut.f.mtx, 2, function(x){ifelse(x == 0, 0, 1)} ))
fuscc.pri.mut.f.mtx = fuscc.pri.mut.f.mtx[rownames(somatic.maf.pri.f.mtx),]

priVSmeta.res = as.data.frame(matrix(nrow = nrow(somatic.maf.pri.f.mtx), ncol = 4))
rownames(priVSmeta.res) = row.names(somatic.maf.pri.f.mtx)
colnames(priVSmeta.res) = c("FUSCC.pri","FUSCC.meta","p.val","adj.p")

priVSmeta.res$FUSCC.pri = apply(fuscc.pri.mut.f.mtx,1,function(x){sum(x)/length(x)})
priVSmeta.res$FUSCC.meta = apply(somatic.maf.pri.f.mtx,1,function(x){sum(x)/length(x)})

for ( i in rownames(priVSmeta.res)){
  tmp1 = t(rbind(rep("FUSCC.pri",ncol(fuscc.pri.mut.f.mtx)),fuscc.pri.mut.f.mtx[i,]))
  tmp2= t(rbind(rep("FUSCC.meta",ncol(somatic.maf.pri.f.mtx)),somatic.maf.pri.f.mtx[i,]))
  tmp = as.data.frame(rbind(tmp1,tmp2))
  colnames(tmp) = c("cohort","gene")
  tmp$cohort = as.factor(tmp$cohort)
  tmp$gene = as.factor(tmp$gene)
  p = chisq_test(cohort~gene,data= tmp, distribution=approximate())
  priVSmeta.res[i,3] <- pvalue(p)[1]
}
priVSmeta.res[,4] = p.adjust(priVSmeta.res[,3],method = "fdr")

ggplot(priVSmeta.res,aes(x=FUSCC.meta, y = FUSCC.pri)) + 
  geom_abline(intercept=0,slope=1, linetype="dashed", size = 1.5, color = '#747475')+
  geom_point(alpha = 1, size = 2.5, col="#76A1C5") + 
  geom_point(data = priVSmeta.res[priVSmeta.res$p.val < 0.05,] , 
             alpha = 1, size = 2.5, col="red") +
  xlim(0,0.1) + ylim(0,0.1) + xlab("FUSCC.meta.Line0")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  theme_bw()+
  geom_text_repel(data= priVSmeta.res[priVSmeta.res$p.val < 0.05,] ,
                  aes(label=rownames(priVSmeta.res[priVSmeta.res$p.val < 0.05,])),
                  col="black",alpha = 1,size=3) 
ggsave(filename = "./FigS4A.pdf")

