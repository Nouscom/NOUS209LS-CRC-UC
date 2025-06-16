library(ggplot2)
library(ggvenn)
library(ggprism)
library(ggpubr)

nous209len<-read.delim("./Data/Nous209_sequence.txt")
nous209len$Len<-nchar(nous209len$FSP.sequence)
extract_parts <- function(vec) {
  result <- c()
  
  for (elem in vec) {
    parts <- unlist(strsplit(elem, "_"))
    if (length(parts) >= 3) {
      prefix <- paste(parts[1:3], collapse = "_")
      
      if (prefix=="FSP_CEP76_1;FSP"){
        prefix="FSP_CEP76_1;FSP_CEP76_2"
        result <- c(result, prefix)
      }else if(prefix=="FSP_SMAD7_1"){
        prefix="FSP_SMAD7_1_2_a"
        result <- c(result, prefix)
        prefix="FSP_SMAD7_1_2_b"
        result <- c(result, prefix)
      }else{
        
        result <- c(result, prefix)
      }
      if (grepl("_1_2", elem)) {
        
        newelem<-gsub("_1_2", "_2", elem)
        parts2 <- unlist(strsplit(newelem, "_"))
        newresult<-paste(parts2[1:3], collapse = "_")
        result <- c(result, newresult)
      }
      if (grepl("_1_2_3", elem)) {
        
        newelem<-gsub("_1_2", "_2", elem)
        parts2 <- unlist(strsplit(newelem, "_"))
        newresult<-paste(parts2[1:3], collapse = "_")
        result <- c(result, newresult)
        
        newelem<-gsub("_1_2_3", "_3", elem)
        parts2 <- unlist(strsplit(newelem, "_"))
        newresult<-paste(parts2[1:3], collapse = "_")
        result <- c(result, newresult)
        
        
      }
      
    }
  }
  
  return(unique(result))
}

find_len <- function(vec,dfnous) {
  df<-data.frame()
  vec_extr<-extract_parts(vec)
  for (elem in vec_extr) {
    
    temp<-dfnous[grepl(elem,dfnous$FSP_ID),]
    
    if (nrow(temp)>0){
      df<-rbind(df,temp)
    }else{
      print(elem)
    }
  }
  return(unique(df))
}




FSP<-read.delim("./Data/total_FSP_NOUS209SEPPALA.tsv")

clinical<-read.delim("./Data/clinical_history_cancer.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)

HLA<-read.delim("./Data/CRC_HLA_TOTAL_COHORT.tsv")

paired<-read.delim("../Analysis/Data/PairedSample_df.tsv")
hlaint<-HLA[HLA$Patient%in%c(1086,1013,
                             paired$patient_id),]
predictions<-read.delim("./Data/Total_predictions.tsv")


LostFSP<-data.frame()
PredPatient<-data.frame()
FSPgroup<-data.frame()
for (pt in c(1013,1086,767,197,651,798,237)){
  
  
  mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  meta1<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymin,]$name
  meta2<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymax,]$name
  
  lost<-setdiff(FSP[FSP$SAMPLE==meta1,]$X26,
                FSP[FSP$SAMPLE==meta2,]$X26)
  gained<-setdiff(FSP[FSP$SAMPLE==meta2,]$X26,
                  FSP[FSP$SAMPLE==meta1,]$X26)
  
  keep<-intersect(FSP[FSP$SAMPLE==meta2,]$X26,
                  FSP[FSP$SAMPLE==meta1,]$X26)
  
  LostFSP<-rbind(LostFSP,data.frame("Pt"=rep(pt,n_distinct(lost)),
                                    "FSP"=lost))
  HLApt<-HLA[HLA$Patient==pt,]
  HLApt$HLA<-paste0("HLA-",HLApt$HLA)
  
  temp<-data.frame(Cond=c("Lost","Kept","Gained"),
                   NFSP=c(length(lost),length(keep),length(gained)),
                   Pt=c(pt,pt,pt))
  FSPgroup<-rbind(FSPgroup,temp)
  
  
  for (fs in lost){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                            grepl(fs_name,predictions$FSP),]
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Lost"
      predtemp$Pt<-pt
      PredPatient<-rbind(PredPatient,predtemp)
    }
  }
  
  
  for (fs in gained){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                            grepl(fs_name,predictions$FSP),]
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Gained"
      predtemp$Pt<-pt
      PredPatient<-rbind(PredPatient,predtemp)
    }
  }
  
  for (fs in keep){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                            grepl(fs_name,predictions$FSP),]
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Kept"
      predtemp$Pt<-pt
      
      PredPatient<-rbind(PredPatient,predtemp)
    }
  }
  
  
}


FSPURO<-read.delim("./Data/NOUS209_Urothelial_lookup_allsamples.tsv")
FSPURO

for (pt in c(1122,392)){
  
  
  if (pt ==1122){
    sample1<-clinical[clinical$patient_id==pt,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$X26,
                  sample_uro_fps$VECTOR)
    gained<-setdiff(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    temp<-data.frame(Cond=c("Lost","Kept","Gained"),
                     NFSP=c(length(lost),length(keep),length(gained)),
                     Pt=c(pt,pt,pt))
    
    FSPgroup<-rbind(FSPgroup,temp)
    
    
    LostFSP<-rbind(LostFSP,data.frame("Pt"=rep(pt,n_distinct(lost)),
                                      "FSP"=lost))
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    
    
    for (fs in lost){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Lost"
        predtemp$Pt<-pt
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    
    for (fs in gained){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Gained"
        predtemp$Pt<-pt
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    for (fs in keep){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Kept"
        predtemp$Pt<-pt
        
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    
    
    
    
  }
  
  
  if (pt ==392){
    
    mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    sample1<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymin,]$name
    sample2<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymax,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$X26,
                  sample_uro_fps$VECTOR)
    gained<-setdiff(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    temp<-data.frame(Cond=c("Lost","Kept","Gained"),
                     NFSP=c(length(lost),length(keep),length(gained)),
                     Pt=c("392_PU","392_PU","392_PU"))
    
    FSPgroup<-rbind(FSPgroup,temp)
    
    LostFSP<-rbind(LostFSP,data.frame("Pt"=rep("392_PU",n_distinct(lost)),
                                      "FSP"=lost))
    
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    
    
    for (fs in lost){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Lost"
        predtemp$Pt<-"392_PU"
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    
    for (fs in gained){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Gained"
        predtemp$Pt<-"392_PU"
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    for (fs in keep){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Kept"
        predtemp$Pt<-"392_PU"
        
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    
    
    
    lost<-setdiff(sample_uro_fps$VECTOR,
                  FSP[FSP$SAMPLE==sample2,]$X26
    )
    gained<-setdiff(FSP[FSP$SAMPLE==sample2,]$X26,
                    sample_uro_fps$VECTOR
    )
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample2,]$X26)
    temp<-data.frame(Cond=c("Lost","Kept","Gained"),
                     NFSP=c(length(lost),length(keep),length(gained)),
                     Pt=c("392_UM","392_UM","392_UM"))
    LostFSP<-rbind(LostFSP,data.frame("Pt"=rep("392_UM",n_distinct(lost)),
                                      "FSP"=lost))
    FSPgroup<-rbind(FSPgroup,temp)
    
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    
    
    for (fs in lost){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Lost"
        predtemp$Pt<-"392_UM"
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    
    for (fs in gained){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Gained"
        predtemp$Pt<-"392_UM"
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
    for (fs in keep){
      fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
      
      predtemp<-predictions[predictions$allele%in%HLApt$HLA&
                              grepl(fs_name,predictions$FSP),]
      
      if (nrow(predtemp)>0){
        predtemp$Cond<-"Kept"
        predtemp$Pt<-"392_UM"
        
        PredPatient<-rbind(PredPatient,predtemp)
      }
    }
    
  }
}



###
hla<-read.delim("./Data/CRC_HLA_TOTAL_COHORT.tsv")
hlapaired<-hla[hla$Patient%in%c(1086,1013,798,197,1122,767,651,"392",237),]
predictions<-read.delim("./Data/Total_predictions.tsv")
predictions<-predictions[predictions$length==9,]
hlapaired%>%
  group_by(HLA)%>%
  summarize(npt=n_distinct(Patient))%>%
  arrange(-npt)


Patient_A03<-hlapaired[hlapaired$HLA=="A*03:01",]$Patient%>%unique()
Patient_A02<-hlapaired[hlapaired$HLA=="A*02:01",]$Patient%>%unique()



LostFSP<-LostFSP[LostFSP$Pt%in%
                   c(1086,1013,798,197,1122,767,651,"392_PU",237),]


LostFSP_03<-LostFSP%>%
  filter(Pt%in%c(237,1122,1013))%>%
  group_by(FSP)%>%
  summarise(N_Patients=n_distinct(Pt))%>%
  arrange(-N_Patients)


pred_03<-predictions[predictions$allele=="HLA-A*03:01",]

Lost03_with_pred<-data.frame()

for (fsp in unique(pred_03$FSP)){
  
  if (nrow(LostFSP_03[grepl(fsp,LostFSP_03$FSP),])>0){
    print(fsp)
    temp<-LostFSP_03[grepl(fsp,LostFSP_03$FSP),]
    temp$N_epitopes<-sum(pred_03$FSP==fsp)
    Lost03_with_pred<-rbind(Lost03_with_pred,temp)
  }
}


Lost03_with_pred$FSP
View(Lost03_with_pred)
write.table(Lost03_with_pred[order(Lost03_with_pred$N_Patients,decreasing = T),],"FSP_lost_patient_A0301_withPred.tsv",sep="\t",row.names = F,quote = F)


LostFSP_02<-LostFSP%>%
  filter(Pt%in%c(Patient_A02,"392_PU"))%>%
  group_by(FSP)%>%
  summarise(N_Patients=n_distinct(Pt))%>%
  arrange(-N_Patients)

write.table(LostFSP_02,"FSP_lost_patient_A0201.tsv",sep="\t",row.names = F,quote = F)



pred_02<-predictions[predictions$allele=="HLA-A*02:01",]

Lost02_with_pred<-data.frame()
for (fsp in unique(pred_02$FSP)){
  
  if (nrow(LostFSP_02[grepl(fsp,LostFSP_02$FSP),])>0){
    temp<-LostFSP_02[grepl(fsp,LostFSP_02$FSP),]
    temp$N_epitopes<-sum(pred_02$FSP==fsp)
    Lost02_with_pred<-rbind(Lost02_with_pred,temp)
  }
}

Lost02_with_pred$FSP
View(Lost02_with_pred)
write.table(Lost02_with_pred[order(Lost02_with_pred$N_Patients,decreasing = T),],"FSP_lost_patient_A0201_withPred.tsv",sep="\t",row.names = F,quote = F)


PredPatient$Cond<-factor(PredPatient$Cond,
                         levels = c("Lost","Kept","Gained"))
plot1<-PredPatient%>%
  filter(length==9)%>%
  group_by(Cond,Pt)%>%
  summarise(Predictions=n_distinct(peptide))

plot1$Pt<-factor(plot1$Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))
svg("./Output/Fig_5D_barplot.svg",
       height = 4,width =8)

ggplot(plot1[plot1$Pt!="392_UM",],aes(factor(Pt),Predictions,fill=Cond))+
  geom_bar(stat="identity",color="black",position = position_dodge())+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylim(c(0,150))+ylab("# Predicted epitopes")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")
dev.off()




provaplot1<-PredPatient%>%
  filter(length==9)%>%
    group_by(Cond,Pt,FSP)%>%
  summarise(Predictions=n_distinct(peptide))%>%
  summarise(Predictions=median(Predictions))
ggplot(provaplot1[plot1$Pt!="392_UM",],aes(factor(Pt),Predictions,fill=Cond))+
  geom_bar(stat="identity",color="black",position = position_dodge())+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("# Predicted epitopes")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")



ggsave("./Output/Fig_5D_barplot.svg",dpi = 300,units = "cm",
       height = 7,width =25)

ggsave("./Output/Fig_5D_barplot.jpeg",dpi = 300,units = "cm",
       height = 7,width =25 )

ggsave("./Output/Fig_5D_barplot.pdf",dpi = 300,units = "cm",
       height = 7,width =25 )




###############################

FSPgroup$Cond<-factor(FSPgroup$Cond,levels = c("Lost","Kept","Gained"))
plot2<-FSPgroup
svg("./Output/Fig_5B_boxNous209.svg",
       height = 4,width =3.5)
ggplot(plot2[plot2$Pt!="392_UM",],aes(factor(Cond,levels=c("Lost","Kept","Gained")),NFSP))+
  geom_boxplot()+
  geom_line(aes(group = factor(Pt)),size=0.8,color="grey",linetype = "dashed")+
  geom_point(aes(fill=factor(Pt)),shape=21,size=3)+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("# Nous-209")+
  stat_compare_means(comparisons = list(c("Lost","Kept"),
                                        
                                        c("Kept","Gained"),
                                        c("Lost","Gained")),
                     paired=T)+ylim(c(0,80))+
  theme(legend.position = "None")
dev.off()


ggsave("./Output/Fig_5B_boxNous209.jpeg",dpi = 300,units = "cm",
       height = 10,width =14 )
ggsave("./Output/Fig_5B_boxNous209.pdf",units = "cm",
       height = 10,width =14 )



plot2[plot2$Pt!="392_UM",]$NFSP%>%median()





svg("./Output/Fig_5D_GoodBinder.svg",width = 3.5,height = 4)
PredPatient%>%
  filter(length==9)%>%
  filter(Pt!="392_UM")%>%
  group_by(Cond,Pt)%>%
  summarise(Predictions=n_distinct(peptide))%>%
  arrange(Pt)%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),Predictions))+
  geom_boxplot()+
  geom_line(aes(group = factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),size=0.8,color="grey",linetype = "dashed")+
  geom_point(aes(fill=factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),shape=21,size=3)+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("# Predicted epitopes")+
  stat_compare_means(comparisons = list(
    c("Lost","Kept"),
    
    c("Kept","Gained"),
    c("Lost","Gained")),
    paired=T)+ylim(c(0,200))+
  theme(legend.position = "None")
dev.off()

ggsave("./Output/Fig_5D_GoodBinder.jpeg",dpi = 300,units = "cm",
       height = 10,width =14 )
ggsave("./Output/Fig_5D_GoodBinder.pdf",units = "cm",
       height = 10,width =14 )





fsplenght<-read.delim("./Data/Nous209_sequence.txt")
fsplenght$Len<-nchar(fsplenght$FSP.sequence)


LenPatient<-data.frame()
for (pt in c(1013,1086,767,197,651,798,237)){
  
  pt
  mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  meta1<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymin,]$name
  meta2<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymax,]$name
  
  lost<-setdiff(FSP[FSP$SAMPLE==meta1,]$X26,
                FSP[FSP$SAMPLE==meta2,]$X26)
  gained<-setdiff(FSP[FSP$SAMPLE==meta2,]$X26,
                  FSP[FSP$SAMPLE==meta1,]$X26)
  
  keep<-intersect(FSP[FSP$SAMPLE==meta2,]$X26,
                  FSP[FSP$SAMPLE==meta1,]$X26)
  
  
  lost
  
  #for (fs in lost){
  # fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
  
  predtemp<- find_len(lost,nous209len)
  if (nrow(predtemp)>0){
    predtemp$Cond<-"Lost"
    predtemp$Pt<-pt
    LenPatient<-rbind(LenPatient,predtemp)
  }
  # }
  
  
  #for (fs in gained){
  fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
  
  predtemp<-find_len(gained,nous209len)
  
  if (nrow(predtemp)>0){
    predtemp$Cond<-"Gained"
    predtemp$Pt<-pt
    LenPatient<-rbind(LenPatient,predtemp)
  }
  #}
  
  #for (fs in keep){
  # fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
  
  predtemp<-find_len(keep,nous209len)
  
  if (nrow(predtemp)>0){
    predtemp$Cond<-"Kept"
    predtemp$Pt<-pt
    
    LenPatient<-rbind(LenPatient,predtemp)
  }
  # }
  
  
}


FSP$SAMPLE

############ COMPARE SAMPLES WITH UROTHELIAL

for (pt in c(1122,392)){
  
  
  if (pt ==1122){
    sample1<-clinical[clinical$patient_id==pt,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$X26,
                  sample_uro_fps$VECTOR)
    gained<-setdiff(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    
    
    
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    
    # for (fs in lost){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(lost,nous209len)
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Lost"
      predtemp$Pt<-pt
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    
    #for (fs in gained){
    # fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(gained,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Gained"
      predtemp$Pt<-pt
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    #for (fs in keep){
    # fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(keep,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Kept"
      predtemp$Pt<-pt
      
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    
  }
  
  
  if (pt ==392){
    
    mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    sample1<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymin,]$name
    sample2<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymax,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$X26,
                  sample_uro_fps$VECTOR)
    gained<-setdiff(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample1,]$X26)
    #for (fs in lost){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(lost,nous209len)
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Lost"
      predtemp$Pt<-"392_PU"
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    
    #for (fs in gained){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(gained,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Gained"
      predtemp$Pt<-"392_PU"
      LenPatient<-rbind(LenPatient,predtemp)
    }
    # }
    
    #for (fs in keep){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(keep,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Kept"
      predtemp$Pt<-"392_PU"
      
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    
    
    lost<-setdiff(sample_uro_fps$VECTOR,
                  FSP[FSP$SAMPLE==sample2,]$X26
    )
    gained<-setdiff(FSP[FSP$SAMPLE==sample2,]$X26,
                    sample_uro_fps$VECTOR
    )
    
    keep<-intersect(sample_uro_fps$VECTOR,
                    FSP[FSP$SAMPLE==sample2,]$X26)
    temp<-data.frame(Cond=c("Lost","Kept","Gained"),
                     NFSP=c(length(lost),length(keep),length(gained)),
                     Pt=c("392_UM","392_UM","392_UM"))
    
    # for (fs in lost){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(lost,nous209len)
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Lost"
      predtemp$Pt<-"392_UM"
      LenPatient<-rbind(LenPatient,predtemp)
    }
    # }
    
    
    #for (fs in gained){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(gained,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Gained"
      predtemp$Pt<-"392_UM"
      LenPatient<-rbind(LenPatient,predtemp)
    }
    #}
    
    #for (fs in keep){
    fs_name<-sapply(strsplit(fs, "_"), function(x) paste(x[1:2], collapse = "_"))
    
    predtemp<-find_len(keep,nous209len)
    
    if (nrow(predtemp)>0){
      predtemp$Cond<-"Kept"
      predtemp$Pt<-"392_UM"
      
      LenPatient<-rbind(LenPatient,predtemp)
    }
    # }
    
    
  }
}
svg("./Output/Fig_5C_CumulativeLength.svg",width = 3.5,height = 4)
LenPatient%>%
  group_by(Cond,Pt)%>%
  summarise(Cumulative=sum(Len))%>%
  filter(Pt!="392_UM")%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),Cumulative))+
  geom_boxplot()+
  geom_line(aes(group = factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),size=0.8,color="grey",linetype = "dashed")+
  geom_point(aes(fill=factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),shape=21,size=3)+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("Cumulative FSP length")+
  stat_compare_means(comparisons = list(
    c("Lost","Kept"),
    
    c("Kept","Gained"),
    c("Lost","Gained")),paired=T)+theme(legend.position = "None")
dev.off()


ggsave("./Output/Fig_5C_CumulativeLength.jpeg",dpi = 300,units = "cm",
       height = 10,width =14 )
ggsave("./Output/Fig_5C_CumulativeLength.pdf",units = "cm",
       height = 10,width =14 )




svg("./Output/Fig_5C_MedianLength.svg",width = 3.5,height = 4)
LenPatient%>%
  group_by(Cond,Pt,FSP_ID)%>%
  summarise(Cumulative=sum(Len))%>%
  summarise(Cumulative_M=median(Cumulative))%>%
  filter(Pt!="392_UM")%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),Cumulative_M))+
  geom_boxplot()+
  geom_line(aes(group = factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),size=0.8,color="grey",linetype = "dashed")+
  geom_point(aes(fill=factor(Pt,levels = c("1086",1013,798,197,1122,767,651,"392_PU","392_UM",237))),shape=21,size=3)+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("Median FSP length")+
  stat_compare_means(comparisons = list(
    c("Lost","Kept"),
    
    c("Kept","Gained"),
    c("Lost","Gained")),paired=T)+theme(legend.position = "None")

dev.off()
head(LenPatient)


LenPatient%>%
  group_by(Cond,Pt,FSP_ID)%>%
  filter(Cond=="Gained")%>%
  filter(Pt!="392_UM")%>%
  summarise(Cumulative=sum(Len))%>%
  summarise(Cumulative_M=max(Cumulative))%>%
  
  summarise(median(Cumulative_M))


######## Annotation for BOXPLOT GOOD BINDERS

annotation<-read_xlsx("./Data/Annotation_IHC_perBarplot.xlsx")
annotation$Patient<-as.character(annotation$Patient)
#annotation<-annotation[annotation$Patient!="392_UM",]

svg("./Output/B2M_annotation.svg",width = 8,height = 0.5)
ggplot(annotation)+
  geom_rect(aes(xmin  = X,xmax = X+.5,ymin = Y-0.3,ymax=Y+0.3,fill=B2M))+
  theme_void()+scale_fill_manual(values = c("Negative"="grey70","Positive"="orange"))+
  theme(legend.position = "None")
dev.off()

ggsave("./Output/B2M_annotation.jpeg",width =2000,height = 100 ,dpi=300,units = "px")
ggsave("./Output/B2M_annotation.pdf",width =2000,height = 100 ,dpi=300,units = "px")
ggsave("./Output/B2M_annotation.svg",width =2000,height = 100 ,dpi=300,units = "px")











