library(readxl)
library(ggplot2)
library(dplyr)
library(pheatmap)
FSP<-read.delim("./Data/total_FSP_NOUS209SEPPALA.tsv")

predictions<-read.delim("./Data/Total_predictions.tsv")

clinical<-read.delim("./Data/clinical_history_cancer.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)

HLA<-read.delim("./Data/CRC_HLA_TOTAL_COHORT.tsv")

FSP$FSPname<-sub(":.*", "", FSP$X26)

FSPLOST<-data.frame()
# 1086, 767 e 798 RIMOSSI
for (pt in c(1013,197,237,651)){
  
  
  mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
  meta1<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymin,]$name
  meta2<-clinical[clinical$patient_id==pt&
                    clinical$crc_ordinal_number==mymax,]$name
  
  lost<-setdiff(FSP[FSP$SAMPLE==meta1,]$FSPname,
                FSP[FSP$SAMPLE==meta2,]$FSPname)
  #gained<-setdiff(FSP[FSP$SAMPLE==meta2,]$FSP$FSPname,
  #               FSP[FSP$SAMPLE==meta1,]$FSP$FSPname)
  
  #keep<-intersect(FSP[FSP$SAMPLE==meta2,]$FSP$FSPname,
  #               FSP[FSP$SAMPLE==meta1,]$FSP$FSPname)
  
  
  HLApt<-HLA[HLA$Patient==pt,]
  HLApt$HLA<-paste0("HLA-",HLApt$HLA)
  HLApt<-unique(HLApt)
  
  #for ( ls in lost){
  
  # lssearch<-paste0(stringr::str_split(ls, "_")[[1]][1:3], collapse = "_")
  #if( nrow(predictions[predictions$allele%in%HLApt$HLA&
  #                    grepl(lssearch,predictions$FSP),])>0){
  
  #goodallele<-unique(predictions[predictions$allele%in%HLApt$HLA&
  #                                grepl(lssearch,predictions$FSP),]$allele)
  
  #temp<-data.frame(FSP=rep(ls,length(goodallele)),
  #                Pt=rep(pt,length(goodallele)),
  #               Allele=goodallele)
  #predictions[predictions$allele%in%HLApt$HLA&
  #              grepl(lssearch,predictions$FSP),]
  #temp<-merge(temp,HLApt[HLApt$HLA%in%predictions[predictions$allele%in%HLApt$HLA&
  #                                                 grepl(lssearch,predictions$FSP),]$allele,],by.x="Pt",by.y="Patient")
  #FSPLOST<-rbind(FSPLOST,temp)
  #}
  #}
  
  
  temp<-data.frame(FSP=lost,Pt=rep(pt,length(lost)))
  #temp<-merge(temp,HLApt,by.x="Pt",by.y="Patient")
  FSPLOST<-rbind(FSPLOST,temp)
  
  
}

FSPURO<-read.delim("./Data/NOUS209_Urothelial_lookup_allsamples.tsv")
FSPURO

FSPURO$FSPname<-sub(":.*", "", FSPURO$VECTOR)
for (pt in c(1122,392)){
  
  
  if (pt ==1122){
    sample1<-clinical[clinical$patient_id==pt,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$FSPname,
                  sample_uro_fps$FSPname)
    
    
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    HLApt<-unique(HLApt)
    
    for ( ls in lost){
      
      lssearch<-paste0(stringr::str_split(ls, "_")[[1]][1:3], collapse = "_")
      if( nrow(predictions[predictions$allele%in%HLApt$HLA&
                           grepl(lssearch,predictions$FSP),])>0){
        goodallele<-unique(predictions[predictions$allele%in%HLApt$HLA&
                                         grepl(lssearch,predictions$FSP),]$allele)
        
        temp<-data.frame(FSP=rep(ls,length(goodallele)),
                         Pt=rep(pt,length(goodallele)),
                         Allele=goodallele)
        predictions[predictions$allele%in%HLApt$HLA&
                      grepl(lssearch,predictions$FSP),]
        #temp<-merge(temp,HLApt[HLApt$HLA%in%predictions[predictions$allele%in%HLApt$HLA&
        #                                                 grepl(lssearch,predictions$FSP),]$allele,],by.x="Pt",by.y="Patient")
        #FSPLOST<-rbind(FSPLOST,temp)
      }
    }
    
    temp<-data.frame(FSP=lost,Pt=rep(pt,length(lost)))
    #temp<-merge(temp,HLApt,by.x="Pt",by.y="Patient")
    FSPLOST<-rbind(FSPLOST,temp)
    
  }
  
  
  if (pt ==392){
    
    mymin<-min(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    mymax<-max(clinical[clinical$patient_id==pt,]$crc_ordinal_number)
    sample1<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymin,]$name
    sample2<-clinical[clinical$patient_id==pt&
                        clinical$crc_ordinal_number==mymax,]$name
    sample_uro_fps<-FSPURO[FSPURO$SampleID==paste0("LS-",pt),]
    
    lost<-setdiff(FSP[FSP$SAMPLE==sample1,]$FSPname,
                  sample_uro_fps$FSPname)
    
    #lost<-setdiff(sample_uro_fps$FSPname,FSP[FSP$SAMPLE==sample2,]$FSPname)
    
    
    HLApt<-HLA[HLA$Patient==pt,]
    HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    HLApt<-unique(HLApt)
    
    for ( ls in lost){
      
      lssearch<-paste0(stringr::str_split(ls, "_")[[1]][1:3], collapse = "_")
      if( nrow(predictions[predictions$allele%in%HLApt$HLA&
                           grepl(lssearch,predictions$FSP),])>0){
        goodallele<-unique(predictions[predictions$allele%in%HLApt$HLA&
                                         grepl(lssearch,predictions$FSP),]$allele)
        
        temp<-data.frame(FSP=rep(ls,length(goodallele)),
                         Pt=rep(pt,length(goodallele)),
                         Allele=goodallele)
        # temp<-merge(temp,HLApt[HLApt$HLA%in%predictions[predictions$allele%in%HLApt$HLA&
        #                                                 grepl(lssearch,predictions$FSP),]$allele,],by.x="Pt",by.y="Patient")
        # FSPLOST<-rbind(FSPLOST,temp)
        #FSPLOST[FSPLOST$Pt==392,]$Pt<-"392_PU"
      }
    }
    
    temp<-data.frame(FSP=lost,Pt=rep("392",length(lost)))
    #temp<-merge(temp,HLApt,by.x="Pt",by.y="Patient")
    FSPLOST<-rbind(FSPLOST,temp)
    
    
    #lost<-setdiff(sample_uro_fps$FSPname,
    #             FSP[FSP$SAMPLE==sample2,]$FSPname
    #)
    
    #HLApt<-HLA[HLA$Patient==pt,]
    #HLApt$HLA<-paste0("HLA-",HLApt$HLA)
    #HLApt<-unique(HLApt)
    
    #for ( ls in lost){
    
    # lssearch<-paste0(stringr::str_split(ls, "_")[[1]][1:3], collapse = "_")
    #if( nrow(predictions[predictions$allele%in%HLApt$HLA&
    #                  grepl(lssearch,predictions$FSP),])>0){
    #goodallele<-unique(predictions[predictions$allele%in%HLApt$HLA&
    #                                   grepl(lssearch,predictions$FSP),]$allele)
    
    # temp<-data.frame(FSP=rep(ls,length(goodallele)),
    #                 Pt=rep(pt,length(goodallele)),
    #                Allele=goodallele)
    #temp<-merge(temp,HLApt[HLApt$HLA%in%predictions[predictions$allele%in%HLApt$HLA&
    #                                                 grepl(lssearch,predictions$FSP),]$allele,],by.x="Pt",by.y="Patient")
    #FSPLOST<-rbind(FSPLOST,temp)
    #FSPLOST[FSPLOST$Pt=="392",]$Pt<-"392_UM"
    #  }
    #}
    
    #temp<-data.frame(FSP=lost,Pt=rep("392_UM",length(lost)))
    #temp<-merge(temp,HLApt,by.x="Pt",by.y="Patient")
    #FSPLOST<-rbind(FSPLOST,temp)
    
  }
  
  
}


FSPLOST$Pt

keep<-FSPLOST%>%
  group_by(FSP)%>%
  summarize(np=n_distinct(Pt))%>%
  filter(np>1)

mat<-matrix(0,nrow=n_distinct(FSPLOST$Pt),ncol = n_distinct(FSPLOST[FSPLOST$FSP%in%keep$FSP,]$FSP))
dim(mat)
rownames(mat)<-unique(FSPLOST$Pt)
colnames(mat)<-unique(FSPLOST[FSPLOST$FSP%in%keep$FSP,]$FSP)


for (pt in rownames(mat)){
  for (col in colnames(mat)){
    
    if (nrow(FSPLOST[FSPLOST$Pt==pt&FSPLOST$FSP==col,])>0){
      mat[as.character(pt),col]<-1
      
    }
  }
}

matSum<-colSums(mat)
myorder<-order(matSum,decreasing = T)


pheatmap::pheatmap(mat[,myorder],color=c("white","#135D66"),
                   cellwidth = 9,cellheight = 14,cluster_rows = F,cluster_cols = F,legend = F)


mat2<-mat
for (row in rownames(mat)){
  for (col in colnames(mat)){
    if(mat[row,col]==0){
      
      if (row!="392"){
        # if (row!="392_PU"&row!="392_UM"){
        mymin<-min(clinical[clinical$patient_id==row,]$crc_ordinal_number)
        
        tum1<-clinical[clinical$patient_id==row&
                         clinical$crc_ordinal_number==mymin,]$name
        
        print(tum1)
        if (nrow(FSP[FSP$SAMPLE==tum1&FSP$FSPname==col,])>0){
          mat2[row,col]=-1
        }
        
        
      }
      if(row=="392"){
        mymin<-min(clinical[clinical$patient_id==392,]$crc_ordinal_number)
        
        tum1<-clinical[clinical$patient_id==392&
                         clinical$crc_ordinal_number==mymin,]$name
        
        
        if (nrow(FSP[FSP$SAMPLE==tum1&FSP$FSPname==col,])>0){
          mat2[row,col]=-1
        }
        
        
        
      }
      
      if(row=="392_UM"){
        
        
        
        
        if (nrow(FSPURO[FSPURO$SampleID=="LS-392"&FSPURO$FSPname==col,])>0){
          mat2[row,col]=-1
        }
        
      }
      
    }
  }
}

matSum2<-colSums(mat==1)
myorder2<-order(matSum2,decreasing = T)

rownames(mat2)[6]<-"392_PU"

split_names <- strsplit(colnames(mat2), "_")

# Estrai le prime tre parti e uniscile con "_"
result <- sapply(split_names, function(x) paste(x[1:min(3, length(x))], collapse = "_"))

colnames(mat2)<-result

jpeg("./Output/Fig_5F_Heatmap.jpeg",units = "px",width = 5200,height = 2300,res = 300)
pdf("./Output/Fig_5F_Heatmap.pdf",width = 5200,height = 2300)


svg("./Output/Fig_5F_Heatmap.svg",width = 13,height = 6)
print(pheatmap::pheatmap(mat2[,myorder2],color=c("#2ca02c","white","#1f77b4"),
                   cellwidth = 9,cellheight = 12,fontsize_col = 7,
                   cluster_rows = F,cluster_cols = F,legend = F))

dev.off()



cat(colnames(mat2[,myorder2]),sep="\n")



ResultsHLA<-data.frame()

for(col in colnames(mat2)){
  
  flagHLA="No"
  flagHLAneg=""
  neginposlist=""
  ptnegwithhla=""
  HLAlist=""
  
  
  neg_pat=names(which(mat2[,col]==-1))
  pos_pat=names(which(mat2[,col]==1))
  pos_pat=gsub("_PU","",gsub("_UM","",pos_pat))
  neg_pat=gsub("_PU","",gsub("_UM","",neg_pat))
  
  flagNeg=ifelse(length(neg_pat)>0,"Yes","No")
  
  HLApos<-HLA[HLA$Patient%in%pos_pat,]
  HLAneg<-HLA[HLA$Patient%in%neg_pat,]
  
  sumHLA<-HLApos%>%
    group_by(HLA)%>%
    summarise(npt=n_distinct(Patient))
  
  maxHLA=max(sumHLA$npt)
  
  if( maxHLA==n_distinct(pos_pat)){
    
    flagHLA="Yes"
    HLAlist=paste0(sumHLA[sumHLA$npt==maxHLA,]$HLA,collapse = ";")
    HLAunique=sumHLA[sumHLA$npt==maxHLA,]$HLA
    
    if (!any(HLAneg$HLA%in%HLAunique)){
      flagHLAneg="Yes"
    }else{
      neginposlist=paste0(HLAneg$HLA[HLAneg$HLA%in%HLAunique],collapse = ";")
      neginposunique=HLAneg$HLA[HLAneg$HLA%in%HLAunique]
      ptnegwithhla=paste0(HLAneg[HLAneg$HLA%in%neginposunique,]$Patient,collapse = ";")
    }
    
  }
  
  
  temp=data.frame(FSP=col,
                  N_pos_pt=n_distinct(pos_pat),
                  N_neg_pt=n_distinct(neg_pat),
                  HLA_common_in_pos=flagHLA,
                  Which_HLA_in_common=HLAlist,
                  Non_in_neg=flagHLAneg,
                  HLA_in_neg=neginposlist,
                  In_which_neg_pt=ptnegwithhla)
  ResultsHLA<-rbind(ResultsHLA,temp)
  
  
}




write.table(ResultsHLA,"./Output/Resuls_HLA_LOST_1704.tsv",sep="\t",quote=F,
            row.names = F)



filteredRes<-ResultsHLA[ResultsHLA$Non_in_neg=="Yes"&
                          ResultsHLA$HLA_common_in_pos=="Yes",]



Resprediction=data.frame()
Res_paper<-data.frame()
for (row in 1:nrow(filteredRes)){
  hlalistgoodpred=""
  hla=strsplit(filteredRes[row,]$Which_HLA_in_common,";")[[1]]
  fsp_name=sub("_chr.*","",filteredRes[row,]$FSP)
  
  if(grepl("_1",fsp_name)&grepl("_2",fsp_name)){
    
    fsp_name=c(gsub("_2","",fsp_name),gsub("_1","",fsp_name))
    
  }
  
  
  hla_pred=paste0("HLA-",hla)
  
  predmatch<-predictions[predictions$allele%in%hla_pred &
                           predictions$FSP%in%fsp_name,]
  
  if (nrow(predmatch)>0){
    Res_paper<-rbind(Res_paper,predmatch)
    hlalistgoodpred=paste0(unique(predmatch$allele),collapse = ";")
  }
  temp=data.frame(FSP=filteredRes[row,]$FSP,
                  HLAwithPredictions=hlalistgoodpred)
  Resprediction<-rbind(Resprediction,temp)
  rm(hla_pred)
  rm(predmatch)
  rm(hlalistgoodpred)
  
}

Resprediction


ResultsHLA_pred<-merge(ResultsHLA,Resprediction,by="FSP",all.x = T)


write.table(ResultsHLA_pred,"./Output/Resuls_HLA_with_PredictionsLOST_174.tsv",sep="\t",quote=F,
            row.names = F)





############

HLA_filt<-HLA[HLA$Patient%in%c(1013,197,237,651,392,1122),]
HLA_filt%>%
  group_by(HLA)%>%
  summarize(Pt=n_distinct(Patient))%>%
  arrange(-Pt)


