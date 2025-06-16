FSP<-read.delim("./Data/FSP_NOUS209_CRC.tsv")
FSP<-FSP[FSP$SAMPLE%in%clinical$name,]


CRC_lenght<-data.frame()

for (sample in unique(FSP$SAMPLE)){
  
  FSPpt<-FSP[FSP$SAMPLE==sample,]$X26
  print(FSPpt%>%n_distinct())
  temp<-find_len(FSPpt,nous209len)
  temp$Patient<-sample
  CRC_lenght<-rbind(CRC_lenght,temp)
}


msisensor<-read_excel("./Data/CRC_MSI_Sensor.xlsx")
msisensor[msisensor$MSIscore<=20,]

summaryTot<-CRC_lenght%>%
  group_by(Patient)%>%
  summarise(Tot_len=sum(Len))

summaryTot$MSI<-ifelse(summaryTot$Patient=="K18795-12","MSS","MSI-H")


clinical$name<-gsub("_tumor","",clinical$tumor_sample)
clinical<-clinical[order(clinical$NOUS209,decreasing = T),]
clinical$ID<-paste0(1:nrow(clinical))
summaryTot<-merge(summaryTot,clinical,by.x="Patient",by.y="name")


summaryTot$ID<-factor(summaryTot$ID,levels = clinical[order(clinical$NOUS209,decreasing = T),]$ID)

p1<-ggplot(summaryTot,aes(ID,Tot_len))+
  geom_bar(stat="identity",width = 0.7,color="black",aes(fill=MSI))+
  ggprism::theme_prism(axis_text_angle = 45,base_size = 15)+
  geom_hline(yintercept = 400,color="darkblue",linetype="dashed",size=1)+
  geom_hline(yintercept = median(summaryTot$Tot_len),size=1,color="tomato3",linetype="dashed")+
  geom_text(aes(x=55,y=median(summaryTot$Tot_len)+100,
                label=paste0("Median:", median(summaryTot$Tot_len))),
            size=5,color="tomato3")+
  geom_text(aes(x=57.5,y=400+100),
            label="400",size=5,color="blue")+xlab("")+ylab("Cumulative FSP length")+
  scale_fill_manual(values = c("MSI-H"="#DB456F",
                               "MSS"="#FFBABA"))+
  scale_y_continuous(breaks = seq(0,2800,by=400),limits = c(0,2800))

p1


jpeg("./Output/Fig_supp1_CRC.jpeg",width = 4500,height = 1200,res = 350)
print(p1)
dev.off()
pdf("./Output/Fig_supp1_CRC.pdf",width = 4500,height = 1200)
print(p1)
dev.off()

svg("./Output/Fig_supp1_CRC.svg",width = 14,height = 4)
print(p1)
dev.off()


##############################
FSPuro<-read.delim("./Data/NOUS209_UC.tsv")
clinicaluro<-read.delim("./Data/UC_clinical_data.txt")
clinicaluro<-clinicaluro%>%filter(Location%in%c("Bladder","UTUC"))
FSPuro<-FSPuro[FSPuro$SampleID%in%clinicaluro$SampleID,]
clinicaluro$PatientID
URO_lenght<-data.frame()



for (sample in unique(FSPuro$SampleID)){
  
  FSPpt<-FSPuro[FSPuro$SampleID==sample,]$VECTOR
  print(FSPpt%>%n_distinct())
  temp<-find_len(FSPpt,nous209len)
  temp$Patient<-sample
  URO_lenght<-rbind(URO_lenght,temp)
}

summaryTotUro<-URO_lenght%>%
  group_by(Patient)%>%
  summarise(Tot_len=sum(Len))

infoclinical$SampleID<-gsub("_blood","",infoclinical$Patient)
summaryTotUro$MSI<-ifelse(summaryTotUro$Patient%in%c("LS-130","LS-459_B","LS-1130"),
                          "MSS","MSI-H")

setdiff(summaryTotUro$Patient,infoclinical$Patient)
summaryTotUro<-merge(summaryTotUro,infoclinical,by.x="Patient",by.y="SampleID")


summaryTotUro$ID<-factor(summaryTotUro$ID,levels = infoclinical[order(infoclinical$FSP_lookup,decreasing = T),]$ID)


p1<-summaryTotUro%>%
  
  ggplot(aes(ID,Tot_len))+
  geom_bar(stat="identity",width = 0.7,color="black",aes(fill=MSI))+
  ggprism::theme_prism(axis_text_angle = 45,base_size = 15)+
  geom_hline(yintercept = 400,color="darkblue",linetype="dashed",size=1)+
  geom_hline(yintercept = median(summaryTotUro$Tot_len),size=1,color="tomato3",linetype="dashed")+
  geom_text(aes(x=14.5,y=median(summaryTotUro$Tot_len)+80,
                label=paste0("Median:", median(summaryTotUro$Tot_len))),
            size=5,color="tomato3")+
  geom_text(aes(x=14.5,y=400+90),
            label="400",size=5,color="blue")+xlab("")+ylab("Cumulative FSP length")+
  ylim(c(0,1500))+scale_fill_manual(values = c("MSI-H"="#DB456F",
                                               "MSS"="#FFBABA"))+
scale_y_continuous(breaks = seq(0,1600,by=400),limits = c(0,1600))
p1
jpeg("./Output/Fig_supp2_UC.jpeg",width = 4500,height = 1200,res = 350)
print(p1)
dev.off()
pdf("./Output/Fig_supp2_UC.pdf",width = 4500,height = 1200)
print(p1)
dev.off()

svg("./Output/Fig_supp2_UC.svg",width = 10,height = 4)
print(p1)
dev.off()




###############################################################################



temp1<-summaryTot[,c("tumor_sample","Tot_len","ID")]
names(temp1)[1]<-"Patient"
temp2<-summaryTotUro[,c("Patient","Tot_len","ID")]
temptot<-rbind(temp1,temp2)

provatot<-merge(clinical,temptot,by.x = "tumor_sample",by.y = "Patient")

lines<-data.frame(N_tumor=c("Group 1","Group 2","Group 3"),
                  hline_y=c(median(provatot [provatot$N_tumor==1,]$Tot_len),
                            median(provatot[provatot$N_tumor==2,]$Tot_len),
                            median(provatot[provatot$N_tumor==3,]$Tot_len) ))

annotations <- data.frame(
  N_tumor = c("Group 1", "Group 2", "Group 3"),  # Variabile del facet
  x = c(nrow(ntumor1)-8, nrow(ntumor2)-8, nrow(ntumor3)-10),           # Posizione sull'asse x
  y = c(median(provatot [provatot$N_tumor==1,]$Tot_len)+120, 
        median(provatot [provatot$N_tumor==2,]$Tot_len)+120, 
        median(provatot [provatot$N_tumor==3,]$Tot_len)+120),        # Posizione sull'asse y
  label = c(paste0("Median ",median(provatot [provatot$N_tumor==1,]$Tot_len)),
            paste0("Median ",median(provatot [provatot$N_tumor==2,]$Tot_len)), 
            paste0("Median ",median(provatot [provatot$N_tumor==3,]$Tot_len)))  # Testo
)


annotations2 <- data.frame(
  N_tumor = c("Group 1", "Group 2", "Group 3"),  # Variabile del facet
  x = c(30, 18, 7),           # Posizione sull'asse x
  y = c(490, 
        490, 
        490),        # Posizione sull'asse y
  label = c("400","400","400")  # Testo
)

provatot$N_tumor<-paste0("Group ",provatot$N_tumor)
p1<-ggplot(provatot)+
  geom_bar(aes(x=reorder(ID.x,-Tot_len),y=Tot_len,fill=factor(N_tumor)),stat="identity",
           color="black",width = 0.6)+
  ggprism::theme_prism(base_size = 15,axis_text_angle = 45)+
  facet_grid(~N_tumor,scales = "free_x",space = "free")+
  ggsci::scale_fill_jama()+
  xlab("")+ylab("Total length")+
  geom_hline(data = lines, aes(yintercept = hline_y), 
             color = "tomato3", linetype = "dashed",size=1)+
  # coord_cartesian(xlim=c(1,nrow(ntumor3)-1))+scale_x_discrete(labels=paste0("S",1:nrow(ntumor3)-1))+
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            color = "tomato3", size = 5,)+
  geom_text(data = annotations2,aes(x=x,y=y,label = label),
            label="400",size=5,color="blue")+xlab("")+ylab("Cumulative FSP length")+
  scale_x_discrete(labels=rep("",58))+
  theme(legend.position = "None")+ylim(c(0,3000))+
  geom_hline(yintercept = 400,color="black",linetype="dashed",
                                                             size=1)+
  scale_y_continuous(breaks = seq(0,2800,by=400),limits = c(0,2800))

p1
jpeg("./Output/Fig_supp3_group.jpeg",width = 5000,height = 1200,res = 350)
print(p1)
dev.off()
pdf("./Output/Fig_supp3_group.pdf",width = 5000,height = 1200)
print(p1)
dev.off()

svg("./Output/Fig_supp3_group.svg",width = 15,height = 4)
print(p1)
dev.off()







setwd("./Daticlinici_SEPPALA/Article_Seppala/Paper_Definitivo/")

clinical<-read.delim("./Data/clinical_history_cancer.tsv")
dim(clinical)
clinical$name<-gsub("_tumor","",clinical$tumor_sample)

library(readxl)
ihc<-read_xlsx("./Data/ihc_data_NEW.xlsx")
ihc_crc<-ihc[ihc$Tumor_Sample_Barcode%in%clinical$name,]

table(ihc$b2m0loss1weak2strong3014125026)
