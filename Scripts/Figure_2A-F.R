library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)

############ Fig 2A

clinical<-read.delim("./Data/clinical_history_CRC.tsv")
infouro<-read.delim("./Data/UC_clinical_data.txt")

infouro<-infouro[infouro$Location%in%c("Bladder","UTUC"),]
summary_urothelial<-read.delim("./Data/NOUS209FSP_UC.txt")
summary_urothelial$Patient<-gsub("_blood","",summary_urothelial$Patient)

summaryUrothelial<-summary_urothelial[summary_urothelial$Patient%in%infouro$SampleID,]

urothelial_info<-merge(infouro,summaryUrothelial,by.x = "SampleID",by.y="Patient")
dim(urothelial_info)


clinical<-clinical[,c("N_tumor","NOUS209","MNR","tumor_sample")]
clinical$Tumor<-"CRC"
subsetUrot<-urothelial_info[,c("N_tumor","FSP_lookup","MNR_lookup","SampleID")]
colnames(subsetUrot)<-c("N_tumor","NOUS209","MNR","tumor_sample")
subsetUrot$Tumor<-"UC"
clinical<-rbind(clinical,subsetUrot)
ntumor1<-clinical[clinical$N_tumor==1,]
ntumor2<-clinical[clinical$N_tumor==2,]
ntumor3<-clinical[clinical$N_tumor==3,]





png("./Output/Fig_2A_mosaic.png",width = 1500,height = 1400,res = 380)
mosaicplot(table(clinical$N_tumor,clinical$Tumor),color=c("#D63447","orange"),cex.axis = 1,main = "")

dev.off()

pdf("./Output/Fig_2A_mosaic.pdf",width = 10,height =5)
mosaicplot(table(clinical$N_tumor,clinical$Tumor),color=c("#D63447","orange"),cex.axis = 1,main = "")

dev.off()

svg("./Output/Fig_2A_mosaic.svg",width = 6,height =5)
mat<-table(clinical$N_tumor,clinical$Tumor)
rownames(mat)<-c("Group 1","Group 2","Group 3")
mosaicplot(mat,
           color=c("#D63447","orange"),cex.axis = 1.2,main = "")


text(0.25, 0.65, "31",cex=1.5)
text(0.63, 0.65, "19",cex=1.5)
text(0.9, 0.65, "8",cex=1.5)

text(0.25, 0.03, "4",cex=1.5)
text(0.63, 0.04, "4",cex=1.5)
text(0.90, 0.06, "7",cex=1.5)

dev.off()



############# Fig 2 B
crc_clinical<-read.delim("./Data/clinical_history_CRC.tsv")
urothelial_clinical<-read_xlsx("./Data/UC_clinical_data.xlsx")
dim(urothelial_clinical)


names(urothelial_clinical)
names(crc_clinical)
crc<-crc_clinical[,c("patient_id","Age","N_tumor")]
uro<-urothelial_clinical[,c("Sample","Age at sample collection","Ntumor")]
names(uro)<-c("patient_id","Age","N_tumor")

tot<-rbind(crc,uro)
uro

svg("./Output/Fig_2B_AgeGroup.svg", width = 5, height = 4)
ggplot(tot,aes(factor(N_tumor),Age))+
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom(shape=21,size=3,aes(fill=factor(N_tumor)))+
  stat_compare_means(comparisons = list(c(1,2),
                                        c(2,3),
                                        c(1,3)))+
  xlab("")+theme_prism(base_size = 15)+
  ggsci::scale_fill_jama()+theme(legend.position = "None")+
  scale_x_discrete(labels=c("Group 1","Group 2","Group 3"))

dev.off()


ggsave("./Output/Fig_2B_AgeGroup.pdf",dpi=300,
       width =1500,height = 1300,units = "px")

ggsave("./Output/Fig_2B_AgeGroup.jpeg",dpi=300,
       width =1500,height = 1300,units = "px")


############# Fig 2 C

lines<-data.frame(N_tumor=c("Group 1","Group 2","Group 3"),
                  hline_y=c(median(ntumor1$NOUS209),
                            median(ntumor2$NOUS209),
                            median(ntumor3$NOUS209) ))

annotations <- data.frame(
  N_tumor = c("Group 1", "Group 2", "Group 3"),  # Variabile del facet
  x = c(nrow(ntumor1)-4, nrow(ntumor2)-4, nrow(ntumor3)-4),           # Posizione sull'asse x
  y = c(median(ntumor1$NOUS209)+6, median(ntumor2$NOUS209)+6, median(ntumor3$NOUS209)+6),        # Posizione sull'asse y
  label = c(paste0("Median ",median(ntumor1$NOUS209)),
            paste0("Median ",median(ntumor2$NOUS209)), 
            paste0("Median ",median(ntumor3$NOUS209)))  # Testo
)
clinical$N_tumor<-paste0("Group ",clinical$N_tumor)
svg("./Output/Fig_2C_NOUS209group.svg",width = 11,height = 4)

ggplot(clinical)+
  geom_bar(aes(x=reorder(tumor_sample,-NOUS209),y=NOUS209,fill=factor(N_tumor)),stat="identity",
           color="black",width = 0.6)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+
  facet_grid(~N_tumor,scales = "free_x",space = "free")+
  ggsci::scale_fill_jama()+
  xlab("")+ylab("# Nous-209")+#ggtitle("Nous-209 in Incident CRC and Urothelial Cancers \n(n=73)")+ylim(c(0,100))+
  geom_hline(data = lines, aes(yintercept = hline_y), color = "darkred", linetype = "dashed",size=1)+
  # coord_cartesian(xlim=c(1,nrow(ntumor3)-1))+scale_x_discrete(labels=paste0("S",1:nrow(ntumor3)-1))+
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            color = "tomato3", size = 6,)+
  scale_x_discrete(labels=rep("",58))+
  theme(legend.position = "None")
dev.off()

ggsave("./Output/Fig_2C_NOUS209group.pdf",dpi=300,
       width =4500,height = 1300,units = "px")

ggsave("./Output/Fig_2C_NOUS209group.jpeg",dpi=300,
       width =4500,height = 1300,units = "px")
ggsave("./Output/Fig_2C_NOUS209group.svg",dpi=300,
       width =4500,height = 1300,units = "px")



############# Fig 2 D

svg("./Output/Fig_2D_NOUS209grouptumor.svg",height =4,
       width=7)

ggplot(clinical,aes(factor(N_tumor),NOUS209))+
  geom_boxplot(outliers = F)+
  xlab("")+ylab("# Nous-209 FSP")+
  ggbeeswarm::geom_quasirandom(aes(fill = factor(N_tumor)),shape=21,size=3)+
  ggprism::theme_prism(base_size = 15)+ylim(c(0,150))+
  ggsci::scale_fill_jama()+stat_compare_means(comparisons = 
                                                list(c("Group 1","Group 2"),
                                                     c("Group 1","Group 3"),
                                                     c("Group 2","Group 3")))+
  facet_wrap(~Tumor)+theme(legend.position = "None")
dev.off()

ggsave("./Output/Fig_2D_NOUS209grouptumor.pdf",dpi=300,
       width =2800,height = 1300,units = "px")

ggsave("./Output/Fig_2D_NOUS209grouptumor.jpeg",dpi=300,
       width =2800,height = 1300,units = "px")
ggsave("./Output/Fig_2D_NOUS209grouptumor.svg",dpi=300,
       width =2800,height = 1300,units = "px")





################## Fig 2 E e 2F
clinicalcrc<-read.delim("./Data/clinical_history_CRC.tsv")
clinicalcrc$Stagebis<-ifelse(clinicalcrc$Stage=="I","I","II-III")

svg("./Output/Fig_2E_NOUS209StageCRC.svg",height=4,
       width = 3.5)


ggplot(clinicalcrc[!is.na(clinicalcrc$Stage),],aes(Stagebis,NOUS209))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(comparisons = list(c("I","II-III")))+
  ggbeeswarm::geom_quasirandom(shape=21,size=3,aes(fill=Stagebis))+
  ggprism::theme_prism(base_size = 15)+ggsci::scale_fill_simpsons()+
  theme(legend.position = "None")+ylab("# Nous-209")+
  xlab("Stage")
dev.off()
ggsave("./Output/Fig_2E_NOUS209StageCRC.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_2E_NOUS209StageCRC.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_2E_NOUS209StageCRC.svg",dpi=300,
       width =1400,height = 1300,units = "px")



clinicaluro<-read.delim("./Data/UC_clinical_data.txt")
clinicaluro<-clinicaluro[clinicaluro$Location%in%c("Bladder","UTUC"),]

nousuro<-read.delim("./Data/Summary_MNR_NOUS209FSP_UROTHELIAL_LOOKUP_hg38.txt")
nousuro$Patient<-gsub("_blood","",nousuro$Patient)
clinicaluro<-merge(clinicaluro,nousuro,by.x="SampleID",by.y="Patient")

clinicaluro$Stagebis<-""
clinicaluro[clinicaluro$Stage=="Tis",]$Stagebis<-"pTis"
clinicaluro[clinicaluro$Stage%in%c("Ta","T1","Tis"),]$Stagebis<-"pTis-pTa-pT1"
clinicaluro[clinicaluro$Stage%in%c("T2","T3","T4"),]$Stagebis<-"pT2-pT3-pT4"


clinicalcrc$Stagebis%>%table()
clinicaluro$Stagebis%>%table()

svg("./Output/Fig_2E_NOUS209StageUC.svg",height=4,
    width = 3.5)
ggplot(clinicaluro,aes(factor(Stagebis,
                              levels=c("pTis-pTa-pT1","pT2-pT3-pT4")),FSP_lookup))+
  geom_boxplot(outliers = F)+
  ggbeeswarm::geom_quasirandom(shape=21,size=3,aes(fill=factor(Stagebis,
                                                               levels=c("pTis-pTa-pT1","pT2-pT3-pT4"))))+
  
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("NOUS209")+
  stat_compare_means(comparison=list(
    c("pTis-pTa-pT1","pT2-pT3-pT4")))+ggsci::scale_fill_simpsons()+
  theme(legend.position = "None")+
  ylab("# Nous-209")+
  xlab("Stage")

dev.off()

ggsave("./Output/Fig_2E_NOUS209StageUC.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_2E_NOUS209StageUC.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_2E_NOUS209StageUC.svg",dpi=300,
       width =1400,height = 1300,units = "px")



################################################################################



clinicalcrc<-read.delim("./Data/clinical_history_CRC.tsv")

svg("./Output/Fig_2F_NOUS209Lymph.svg",height=4,
    width = 3.5)
ggplot(clinicalcrc[!is.na(clinicalcrc$Lymphovascular_invasion),],aes(Lymphovascular_invasion,NOUS209))+
  geom_boxplot(outliers = F)+
  ggpubr::stat_compare_means(comparisons = list(c("Yes","No")))+
  ggbeeswarm::geom_quasirandom(shape=21,size=3,aes(fill=factor(Lymphovascular_invasion,levels = c("Yes","No"))))+
  ggprism::theme_prism()+ggsci::scale_fill_simpsons()+
  theme(legend.position = "None")+xlab("")+ylab("# Nous-209")

dev.off()
ggsave("./Output/Fig_2F_NOUS209LymphoInvasion.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_2F_NOUS209LymphoInvasion.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_2F_NOUS209LymphoInvasion.svg",dpi=300,
       width =1400,height = 1300,units = "px")




clinicaluro<-read.delim("./Data/UC_clinical_data.txt")
clinicaluro<-clinicaluro[clinicaluro$Location%in%c("Bladder","UTUC"),]

nousuro<-read.delim("./Data/NOUS209FSP_UC.txt")
nousuro$Patient<-gsub("_blood","",nousuro$Patient)
clinicaluro<-merge(clinicaluro,nousuro,by.x="SampleID",by.y="Patient")

clinicaluro$Stagebis<-""
clinicaluro[clinicaluro$Stage=="Tis",]$Stagebis<-"pTis"
clinicaluro[clinicaluro$Stage%in%c("Ta","T1","Tis"),]$Stagebis<-"pTis-pTa-pT1"
clinicaluro[clinicaluro$Stage%in%c("T2","T3","T4"),]$Stagebis<-"pT2-pT3-pT4"



svg("./Output/Fig_2F_NOUS209SInvasionUC.svg",width = 3.5,height = 4)
ggplot(clinicaluro,aes(factor(Invasion,levels = c("non-MI","MI")),FSP_lookup))+
  geom_boxplot(outliers = F)+
  ggbeeswarm::geom_quasirandom(shape=21,size=3,aes(fill=factor(Invasion,levels = c("non-MI","MI"))))+
  
  ggprism::theme_prism()+xlab("")+
  stat_compare_means(comparison=list(
    c("MI","non-MI")))+ggsci::scale_fill_simpsons()+
  theme(legend.position = "None")+ylab("# Nous-209")+
  xlab("Invasion")

dev.off()

ggsave("./Output/Fig_2F_NOUS209InvasionUC.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_2F_NOUS209InvasionUC.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_2F_NOUS209SInvasionUC.svg",dpi=300,
       width =1400,height = 1300,units = "px")





















