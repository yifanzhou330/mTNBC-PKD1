rm(list = ls()) ; graphics.off()
library(maftools)
library(tidyverse)
library(plyr)

clinic = openxlsx::read.xlsx("clinic.xlsx",1)
somatic.maf = read.maf("/Users/ZYF/Desktop/PKD1/rawdata/配对样本/somt.maf",clinicalData = clinic)

vc_cols = c("#E31F28","#4B718B","#4FAF4B",
            "black",
            "#F4801F","#F5ED3A","#A55627","#F084B5")
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation',
                   'Multi_Hit',
                   'Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del')

somatic.mtx =  mutCountMatrix(somatic.maf, includeSyn = F, countOnly = NULL, removeNonMutated = F)
somatic.mtx = apply(somatic.mtx, 2, function(x){ifelse(x == 0, 0, 1)} )
somatic.mtx.PKD1 = as.data.frame(somatic.mtx["PKD1",])
colnames(somatic.mtx.PKD1) = "PKD1_mut"
somatic.mtx.PKD1$P_M = clinic[rownames(somatic.mtx.PKD1),"P_M"]

cluster.merge4plot = as.data.frame(table(somatic.mtx.PKD1$P_M,somatic.mtx.PKD1$PKD1_mut))
colnames(cluster.merge4plot)[1:2] = c("P_M","PKD1_mut")

cluster.merge4plot2 = ddply(cluster.merge4plot,"P_M", transform,
                            percent = Freq / sum(Freq) *100)
color =  c("0" = "#DBBEBE" ,"1" = "#A52B29"  )

tmp = fisher.test(table(somatic.mtx.PKD1$PKD1_mut,somatic.mtx.PKD1$P_M))
tmp2 = as.numeric(table(somatic.mtx.PKD1$P_M))
p = ggplot(cluster.merge4plot2, aes(x=P_M, y = percent, fill = PKD1_mut )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  ggtitle(paste0("p = ",round(tmp[["p.value"]],4)," sum=",sum(tmp2), " group=",paste(tmp2, collapse = ",")) ) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1.5), 
        axis.ticks = element_line(size = 1),
        axis.text.x=element_text(size = 12),
        axis.text.y=element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) 
p


############ kick sample
somatic.mtx.PKD1 = somatic.mtx.PKD1[ !str_detect(rownames(somatic.mtx.PKD1),"2318522"),]
somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2405985"),]
somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2204510"),]
somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2209317"),]
somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2104753"),]
#somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2312392"),]
#somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2308576"),]
#somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"1902345"),]
# somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2303809"),]
# somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2303805"),]
# somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"2110427"),]

somatic.mtx.PKD1 = somatic.mtx.PKD1[ ! str_detect(rownames(somatic.mtx.PKD1),"1902291"),]

cluster.merge4plot = as.data.frame(table(somatic.mtx.PKD1$P_M,somatic.mtx.PKD1$PKD1_mut))
colnames(cluster.merge4plot)[1:2] = c("P_M","PKD1_mut")

cluster.merge4plot2 = ddply(cluster.merge4plot,"P_M", transform,
                            percent = Freq / sum(Freq) *100)
#cluster.merge4plot2$total_res = factor(cluster.merge4plot2$total_res,levels = c("res", "non_res" )) 
cnv_color =  c("0" = "#DBBEBE" ,"1" = "#A52B29"  )

set.seed(000)
tmp = fisher.test(table(somatic.mtx.PKD1$PKD1_mut,somatic.mtx.PKD1$P_M),simulate.p.value = T)
tmp2 = as.numeric(table(somatic.mtx.PKD1$P_M))
p = ggplot(cluster.merge4plot2, aes(x=P_M, y = percent, fill = PKD1_mut )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cnv_color)+
  ggtitle(paste0("p = ",round(tmp[["p.value"]],4)," sum=",sum(tmp2), " group=",paste(tmp2, collapse = ",")) ) +
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1.5), #调整边框粗细
        axis.ticks = element_line(size = 1),#调整刻度线粗细
        axis.text.x=element_text(size = 12),
        axis.text.y=element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) 
p

ggsave(p,filename =  paste0("./Fig3B.pdf"),height = 6,width = 5)
write.csv(cluster.merge4plot2,file = "Fig3B.csv")
