library(dplyr)
library(readxl)
library(corrplot)
library(ggplot2)
library(ggsci)
library(gghalves)
library(ggpubr)
library(gridExtra)
library(stringr)
library(ggvenn)

########## FIG 1 A

clinical<-read.delim("./Data/clinical_history_cancer.tsv")
clinical<-clinical[order(clinical$NOUS209,decreasing = T),]
clinical$ID<-paste0("",1:nrow(clinical))
p1<-ggplot(clinical)+geom_bar(aes(x=reorder(ID,-NOUS209),y=NOUS209),
                              stat="identity",fill="#47A992",color="black",width = 0.6,
                              position = "dodge")+
  geom_hline(yintercept = median(clinical$NOUS209),color="tomato3",linetype="dashed",linewidth=1)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+
  xlab("")+ylab("# Nous-209")+ggtitle("Nous-209 in Incident CRC\n(n=58)")+
  ylim(c(0,100))+
  annotate("text",x=53,y=median(clinical$NOUS209)+10,
           label=paste0("Median: ",median(clinical$NOUS209)),
           color = "tomato3", size = 5, fontface = "italic")+
  theme(axis.text.x = element_text(size = 10))
p1


png("./Output/Fig_1A_NOUS209CRC.png",width = 4000,height = 1200,res = 380)
print(p1)
dev.off()

pdf("./Output/Fig_1A_NOUS209CRC.pdf",width = 15,height =5)
print(p1)
dev.off()

svg("./Output/Fig_1A_NOUS209CRC.svg",width = 10,height =4)
print(p1)
dev.off()

##### FIG 1 D

msisensor<-read_excel("./Data/CRC_MSI_Sensor.xlsx")
names(msisensor)

clinical$name<-gsub("_tumor","",clinical$tumor_sample)
setdiff(clinical$name,msisensor$Sample)
n_distinct(intersect(clinical$name,msisensor$Sample))
msisensor<-merge(msisensor,clinical,by.x="Sample",by.y="name")
msisensor$Color<-ifelse(msisensor$MSIscore<=20,"MSS","MSI-H")

msisensor$ID<-factor(msisensor$ID,levels = clinical[order(clinical$NOUS209,decreasing = T),]$ID)
p1<-ggplot(msisensor,aes(ID,MSIscore))+
  geom_bar(stat = "identity",color="black",aes(fill=Color),width = 0.7)+
  geom_hline(linetype="dashed",size=1,color="black",yintercept = 20)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+
  xlab("")+scale_fill_manual(values = c("MSI-H"="#DB456F",
                                        "MSS"="#FFBABA"))+
  scale_y_continuous(breaks = seq(0,100,by=20),limits = c(0,100))+theme(legend.position = "None")+
  theme(axis.text.x = element_text(size = 10))

png("./Output/Fig_1D_MsiCRC.png",width = 4000,height = 1200,res = 380)
print(p1)
dev.off()

pdf("./Output/Fig_1D_MsiCRC.pdf",width = 15,height =5)
print(p1)
dev.off()
svg("./Output/Fig_1D_MsiCRC.svg",width = 10,height =4)
print(p1)
dev.off()


########### Fig 1 B
clinical<-read.delim("./Data/Summary_MNR_NOUS209FSP_UROTHELIAL_LOOKUP_hg38.txt")
info<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
info$SampleID<-sub("^([^\\-]+-[^\\-]+)-.*", "\\1", info$Sample)
clinical$PatientID<-sub("_.*", "\\1", clinical$Patient)

clinical$Patient
info$SampleID

info<-info[info$Location %in% c("UTUC","Bladder"),]
bladder<-info[info$Location %in% c("Bladder"),]$SampleID
utuc<-info[info$Location %in% c("UTUC"),]$SampleID
infoclinical<-clinical[clinical$PatientID%in%info$SampleID,]
infoclinical$Location<-""
infoclinical[infoclinical$PatientID%in%bladder,]$Location<-"B"
infoclinical[infoclinical$PatientID%in%utuc,]$Location<-"U"


tempbladder<-infoclinical[infoclinical$Location=="B",]
temputuc<-infoclinical[infoclinical$Location=="U",]

tempbladder<-tempbladder[order(tempbladder$FSP_lookup,decreasing = T),]
tempbladder$ID<-paste0("B",1:nrow(tempbladder))

temputuc<-temputuc[order(temputuc$FSP_lookup,decreasing = T),]
temputuc$ID<-paste0("U",1:nrow(temputuc))
infoclinical<-rbind(temputuc,tempbladder)

p1<-ggplot(infoclinical)+geom_bar(aes(x=reorder(ID,-FSP_lookup),y=FSP_lookup),
                                  stat="identity",fill="#47A992",color="black",width = 0.6,
                                  position = "dodge")+
  geom_hline(yintercept = median(infoclinical$FSP_lookup),color="tomato3",linetype="dashed",linewidth=1)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+
  xlab("")+ylab("# Nous-209")+ggtitle("Nous-209 in UC\n(n=15)")+
  ylim(c(0,100))+
  annotate("text",x=13,y=median(infoclinical$FSP_lookup)+15,label=paste0("Median: ",median(infoclinical$FSP_lookup)),
           color = "tomato3", size = 5, fontface = "italic")



png("./Output/Fig_1B_NOUS209UC.png",width = 2000,height = 1200,res = 380)
print(p1)
dev.off()

pdf("./Output/Fig_1B_NOUS209UC.pdf",width = 15,height =5)
print(p1)
dev.off()

svg("./Output/Fig_1B_NOUS209UC.svg",width = 6,height =4)
print(p1)
dev.off()



######### FIG 1 E

msisensor<-read.delim("./Data/Summary_MSISensor_NEW_DATA_UROTHELIAL_Seppala_2025.tsv")

clinical<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
head(clinical)

clinical$SampleID<-sub("^([^\\-]+-[^\\-]+)-.*", "\\1", clinical$Sample)


clinical<-clinical[clinical$Location%in%c("UTUC","Bladder"),]
msisensor$PatientID<-sub("_.*", "\\1", msisensor$Patient)
msisensor[msisensor$PatientID%in%clinical$SampleID,]%>%dim()


msisensor[msisensor$MSIscore<=20&
            msisensor$Type=="Tumor-Normal",]
setdiff(clinical$SampleID,msisensor$PatientID)




plotMSI<-msisensor[msisensor$Type=="Tumor-Normal" &
                     msisensor$PatientID%in%clinical$SampleID|
                     msisensor$Patient=="LS-830",]

plotMSI<-merge(plotMSI,infoclinical,by.x="Patient",by.y="Patient")
dim(plotMSI)
setdiff(plotMSI$Patient,infoclinical$Patient)
plotMSI$Color<-ifelse(plotMSI$MSIscore<=20,"MSS","MSI-H")

fakeinfo<-infoclinical
fakeinfo[fakeinfo$ID=="UTUC1",]$FSP_lookup<-fakeinfo[fakeinfo$ID=="UTUC1",]$FSP_lookup-1
plotMSI$ID<-factor(plotMSI$ID,levels = fakeinfo[order(fakeinfo$FSP_lookup,decreasing = T),]$ID)
p1<-ggplot(plotMSI,aes(ID,MSIscore))+
  geom_bar(stat = "identity",color="black",aes(fill=Color),width = 0.7)+
  geom_hline(linetype="dashed",size=1,color="black",yintercept = 20)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+xlab("")+
  scale_fill_manual(values = c("MSI-H"="#DB456F",
                               "MSS"="#FFBABA"))+
  scale_y_continuous(breaks = seq(0,100,by=20),limits = c(0,100))+theme(legend.position = "None")
p1


png("./Output/Fig_1E_MsiUC.png",width = 4000,height = 1200,res = 380)
print(p1)
dev.off()

pdf("./Output/Fig_1E_MsiUC.pdf",width = 15,height =5)
print(p1)
dev.off()

svg("./Output/Fig_1E_MsiUC.svg",width = 7.5,height =4)
print(p1)
dev.off()








Urothelial<-read.delim("./Data/NOUS209_Urothelial_lookup_allsamples.tsv")
info<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
info_keep<-info[info$Location%in%c("UTUC","Bladder"),]
info_remove<-info[!info$Location%in%c("UTUC","Bladder"),]

info_remove$SampleID<-sub("^([^\\-]+-[^\\-]+)-.*", "\\1", info_remove$Sample)
info_remove$SampleID
Urothelial$PatientID_folder<-sub("_blood", "",Urothelial$Patient)

Urothelial<-Urothelial[!Urothelial$PatientID_folder%in%info_remove$SampleID,]


CRC<-read.delim("./Data/total_FSP_NOUS209SEPPALA.tsv")
clinical<-read.delim("./Data/clinical_history_cancer.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)
CRC<-CRC[CRC$SAMPLE%in%clinical$name,]


svg("../Paper_Definitivo/Output/Fig_1C_Venn_CRC_Urothelial.svg",height = 4,width = 7)

ggvenn::ggvenn(list("CRC"=CRC$X26,
                    "UC"=Urothelial$VECTOR),set_name_size = 11,
               ,text_size = 8,fill_color = c("tomato3","orange"),auto_scale = T)
dev.off()



















