library(readxl)
library(ggplot2)
library(dplyr)
library(ggprism)


tbm<-read_xlsx("df_b2m_tmb_metachr.xlsx")
ihc<-read_xlsx("ihc_data_NEW.xlsx")

clinical<-read.delim("clinical_history_cancer.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)


ihc
ihc$IHC<-""
ihc[ihc$b2m0loss1weak2strong3014125026==0,]$IHC<-"Negative"
ihc[ihc$b2m0loss1weak2strong3014125026%in%c(1,2,4),]$IHC<-"Positive"
ihc[ihc$b2m0loss1weak2strong3014125026%in%c(3,5,6),]$IHC<-"Het"


positive<-ihc[ihc$IHC=="Positive",]
tmbclinical<-merge(tbm,clinical,by.x="Tumor_Sample_Barcode",by.y="name")

clinicalpositive<-tmbclinical[tmbclinical$Tumor_Sample_Barcode%in%positive$Tumor_Sample_Barcode,]
dim(clinicalpositive)



output_plot = 'Plot/'


###### FIGURE A

col_CD8 =  "cd3p_cd8p_dens_overall_ct"    


subset = clinicalpositive[,c('Tumor_mutation_burden', col_CD8)]
subset = na.omit(subset)

# Calcolo della correlazione
cor_test <- cor.test(subset[,'Tumor_mutation_burden'], subset[,col_CD8], method = "pearson")

# Estrai r e p-value
r_value <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 3)


pTMB_CD8 <-ggplot(subset,aes(Tumor_mutation_burden,subset[,col_CD8]))+
  geom_point(shape=21,fill="orange",size=3)+
  #ggpubr::stat_cor(method = "spearman")+ylab(col)+
  #stat_cor(method = "pearson",label.x = 80)+
  #geom_smooth(method = "lm")+
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  annotate("text", x = min(subset[,1]), y = max(subset[,2])*1.2,
           label = paste0("R = ", r_value, ", p = ", p_value),
           hjust = 0, vjust = 1, size = 5)+
  theme_prism(base_size = 15)+
  xlab("TMB")+
  ylab('CD8 T cells in tumor core')




png(paste0(output_plot,  col_CD8 ,"_TMB_correlation.png"), res=300,width = 1600,height =1500)
print(pTMB_CD8)
dev.off()

pdf(paste0(output_plot, col_CD8, "_TMB_correlation.pdf"), width = 9,height =8)
print(pTMB_CD8)
dev.off()

svg(paste0(output_plot, col_CD8, "_TMB_correlation.svg"), width = 5,height =4)
print(pTMB_CD8)
dev.off()


###### FIGURE B
##### NOUS209 

subset_2 = clinicalpositive[,c('NOUS209', col_CD8)]
subset_2 = na.omit(subset_2)

# Calcolo della correlazione
cor_test2 <- cor.test(subset_2[,'NOUS209'], subset_2[,col_CD8], method = "pearson")

# Estrai r e p-value
r_value <- round(cor_test2$estimate, 2)
p_value <- signif(cor_test2$p.value, 3)


p1NOUS209_CD8 <-ggplot(subset_2,aes(NOUS209,subset_2[,col_CD8]))+
  geom_point(shape=21,fill="orange",size=3)+
  #ggpubr::stat_cor(method = "spearman")+ylab(col)+
  #stat_cor(method = "pearson",label.x = 80)+
  #geom_smooth(method = "lm")+
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  annotate("text", x = min(subset[,1]), y = max(subset[,2])*1.2,
           label = paste0("R = ", r_value, ", p = ", p_value),
           hjust = 0, vjust = 1, size = 5)+
  theme_prism(base_size = 15)+
  scale_x_continuous(breaks=seq(0, 100,by=20))+
  xlab("# NOUS209")+
  ylab('CD8 T cells in tumor core')




png(paste0(output_plot, col_CD8, "_NOUS209_correlation.png"), res=300,width = 1600,height =1500)
print(p1NOUS209_CD8)
dev.off()

pdf(paste0(output_plot, col_CD8, "_NOUS209_correlation.pdf"), width = 9,height =8)
print(p1NOUS209_CD8)
dev.off()

svg(paste0(output_plot, col_CD8, "_NOUS209_correlation.svg"), width =5,height =4)
print(p1NOUS209_CD8)
dev.off()




##### FIGURE C
##### TMB 



col_CD8 =  "cd3p_cd8p_dens_overall_im"    


subset = clinicalpositive[,c('Tumor_mutation_burden', col_CD8)]
subset = na.omit(subset)

# Calcolo della correlazione
cor_test <- cor.test(subset[,'Tumor_mutation_burden'], subset[,col_CD8], method = "pearson")

# Estrai r e p-value
r_value <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 3)


pTMB_CD8 <-ggplot(subset,aes(Tumor_mutation_burden,subset[,col_CD8]))+
  geom_point(shape=21,fill="orange",size=3)+
  #ggpubr::stat_cor(method = "spearman")+ylab(col)+
  #stat_cor(method = "pearson",label.x = 80)+
  #geom_smooth(method = "lm")+
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  annotate("text", x = min(subset[,1]), y = max(subset[,2])*1.2,
           label = paste0("R = ", r_value, ", p = ", p_value),
           hjust = 0, vjust = 1, size = 5)+
  theme_prism(base_size = 15)+
  xlab("TMB")+
  ylab('CD8 T cells at invasive margin')




png(paste0(output_plot,  col_CD8 ,"_TMB_correlation.png"), res=300,width = 1600,height =1500)
print(pTMB_CD8)
dev.off()

pdf(paste0(output_plot, col_CD8, "_TMB_correlation.pdf"), width = 9,height =8)
print(pTMB_CD8)
dev.off()

svg(paste0(output_plot, col_CD8, "_TMB_correlation.svg"), width = 5,height =4)
print(pTMB_CD8)
dev.off()




##### FIGURE D
##### NOUS209 

subset_2 = clinicalpositive[,c('NOUS209', col_CD8)]
subset_2 = na.omit(subset_2)

# Calcolo della correlazione
cor_test2 <- cor.test(subset_2[,'NOUS209'], subset_2[,col_CD8], method = "pearson")

# Estrai r e p-value
r_value <- round(cor_test2$estimate, 2)
p_value <- signif(cor_test2$p.value, 3)


p1NOUS209_CD8 <-ggplot(subset_2,aes(NOUS209,subset_2[,col_CD8]))+
  geom_point(shape=21,fill="orange",size=3)+
  #ggpubr::stat_cor(method = "spearman")+ylab(col)+
  #stat_cor(method = "pearson",label.x = 80)+
  #geom_smooth(method = "lm")+
  geom_smooth(method = "lm", se=TRUE, color = "blue") +
  annotate("text", x = min(subset[,1]), y = max(subset[,2])*1.2,
           label = paste0("R = ", r_value, ", p = ", p_value),
           hjust = 0, vjust = 1, size = 5)+
  theme_prism(base_size = 15)+
  scale_x_continuous(breaks=seq(0, 100,by=20))+
  xlab("# NOUS209")+
  ylab('CD8 T cells at invasive margin')




png(paste0(output_plot, col_CD8, "_NOUS209_correlation.png"), res=300,width = 1600,height =1500)
print(p1NOUS209_CD8)
dev.off()

pdf(paste0(output_plot, col_CD8, "_NOUS209_correlation.pdf"), width = 9,height =8)
print(p1NOUS209_CD8)
dev.off()


svg(paste0(output_plot, col_CD8, "_NOUS209_correlation.svg"), width = 5,height =4)
print(p1NOUS209_CD8)
dev.off()


