library(dplyr)
library(ggpubr)
library(ggplot2)
clinical<-read.delim("./Data/clinical_history_CRC.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)

primary<-clinical[clinical$patient_id==392&
                    clinical$Metachronous_status=="Primary",]$name


metachronous<-clinical[clinical$patient_id==392&
                         clinical$Metachronous_status=="Metachronous",]$name

urothelial<-"LS-392"
primary
metachronous


crc<-read.delim("./Data/FSP_NOUS209_CRC.tsv")
uro<-read.delim("./Data/NOUS209_UC.tsv")

fsp_prim<-crc[crc$SAMPLE==primary,]$X26
fsp_meta<-crc[crc$SAMPLE==metachronous,]$X26
fsp_uro<-uro[uro$PatientID==urothelial,]$VECTOR


library(ggvenn)

ggvenn(list("Primary"=fsp_prim,
            "Urothelial"=fsp_uro,
            "Metachronous"=fsp_meta
),text_size = 6)+ggsci::scale_fill_bmj()

ggsave("./Output/Fig_6A_venn.pdf",width = 8,height = 8)
ggsave("./Output/Fig_6A_venn.svg",width = 8,height = 8,dpi=300)
ggsave("./Output/Fig_6A_venn.jpeg",width = 8,height = 8,dpi=300)


nous209len<-read.delim("./Data/Nous209_sequence.txt")
nous209len$Len<-nchar(nous209len$FSP.sequence)


fsp_prim
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

only43prim<-setdiff(fsp_prim,union(fsp_meta,fsp_uro))
only22meta<-setdiff(fsp_meta,union(fsp_prim,fsp_uro))
only10uro<-setdiff(fsp_uro,union(fsp_prim,fsp_meta))

only10common<-intersect(fsp_meta,intersect(fsp_prim,fsp_uro))
only15commonPrimMeta<-setdiff(intersect(fsp_meta,fsp_prim),fsp_uro)
only10commonUroMeta<-setdiff(intersect(fsp_meta,fsp_uro),fsp_prim)
only13commonPrimUro<-setdiff(intersect(fsp_prim,fsp_uro),fsp_meta)

lost_total<-c(only43prim,only10uro,only13commonPrimUro)
gained_total<-only22meta
kept_total<-c(only10common,only15commonPrimMeta,only10commonUroMeta)


only43prim_df<-find_len(only43prim,nous209len)
only22meta_df<-find_len(only22meta,nous209len)
only10uro_df<-find_len(only10uro,nous209len)
only10common_df<-find_len(only10common,nous209len)
only15commonPrimMeta_df<-find_len(only15commonPrimMeta,nous209len)
only10commonUroMeta_df<-find_len(only10commonUroMeta,nous209len)
only13commonPrimUro_df<-find_len(only13commonPrimUro,nous209len)


only43prim_df$Cond<-"Only Prim"
only22meta_df$Cond<-"Only Meta"
only10uro_df$Cond<-"Only Uro"
only10common_df$Cond<-"All tumor"
only15commonPrimMeta_df$Cond<-"PandM"
only10commonUroMeta_df$Cond<-"UandM"
only13commonPrimUro_df$Cond<-"PandU"
only10common_df$Cond<-"Shared"
only15commonPrimMeta_df$Cond<-"Shared"
only10commonUroMeta_df$Cond<-"Shared"
only13commonPrimUro_df$Cond<-"Shared"



lost_total<-c(only43prim,only10uro,   only13commonPrimUro)
gained_total<-only22meta
kept_total<-c(only10common,only15commonPrimMeta,only10commonUroMeta)

lost_df<-find_len(lost_total,nous209len)
gained_df<-find_len(gained_total,nous209len)
kept_df<-find_len(kept_total,nous209len)


lost_df$Cond<-"Lost"
gained_df$Cond<-"Gained"
kept_df$Cond<-"Kept"




tot_392<-rbind(lost_df,gained_df,kept_df)



svg("./Output/Fig_6B_NOUS.svg",width = 4,height = 5)
tot_392%>%
  group_by(Cond)%>%
  summarise(FSP=n_distinct(FSP_ID))%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),FSP))+
 # geom_line(color="grey",size=0.8,aes(group = 1),linetype="dashed")+
  #geom_point(stat="identity",shape=21,aes(fill=Cond),size=3,color="black")+
  geom_bar(stat="identity",aes(fill=Cond),color="black")+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("# Nous-209")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")+ylim(c(0,80))
dev.off()

ggsave("./Output/Fig_6B_NOUS.pdf",width = 4,height = 4)

ggsave("./Output/Fig_6B_NOUS.jpeg",width = 4,height = 4,dpi=300)


names(tot_392)
svg("./Output/Fig_6C_Len.svg",width = 4,height = 5)
tot_392%>%
  group_by(Cond)%>%
  summarise(FSPLEN=sum(Len))%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),FSPLEN))+
  #geom_line(color="grey",size=0.8,aes(group = 1),linetype="dashed")+
  #geom_point(stat="identity",shape=21,aes(fill=Cond),size=3,color="black",position = position_dodge())+
  geom_bar(stat="identity",aes(fill=Cond),color="black")+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("Cumulative FSP length")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")+ylim(c(0,2000))
dev.off()

ggsave("./Output/Fig_6C_Len.pdf",width = 4,height = 4)

ggsave("./Output/Fig_6C_Len.jpeg",width = 4,height = 4,dpi=300)






predictions<-read.delim("./Data/Total_predictions.tsv")
hla<-read.delim("./Data/CRC_HLA_TOTAL_COHORT.tsv")

hlapt392<-paste0("HLA-",hla[hla$Patient==392,]$HLA)
hlapt392


predictions_392<-predictions[predictions$allele%in%hlapt392,]


Pred_Res<-data.frame()
for ( i in 1:nrow(tot_392)){
  
  
  temp<-predictions_392[predictions_392$FSP==tot_392[i,]$FSP_ID,]
  
  if(nrow(temp)>0){
    temp$Cond<-tot_392[i,]$Cond
    
    Pred_Res<-rbind(Pred_Res,temp)
  }
  
}


Pred_Res<-Pred_Res[Pred_Res$length==9,]
svg("./Output/Fig_6D_GoodBind.svg",width = 4,height = 5)
Pred_Res%>%
  group_by(Cond)%>%
  summarise(Good=n_distinct(peptide))%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),Good))+
  #geom_line(color="grey",size=0.8,aes(group = 1),linetype="dashed")+
  #geom_point(stat="identity",shape=21,aes(fill=Cond),size=3,color="black",position = position_dodge())+
  geom_bar(stat="identity",aes(fill=Cond),color="black")+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("# Predicted epitopes")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")+ylim(c(0,120))
dev.off()


ggsave("./Output/Fig_6D_GoodBind.pdf",width = 4,height = 4)

ggsave("./Output/Fig_6D_GoodBind.jpeg",width = 4,height = 4,dpi=300)


svg("./Output/Fig_5I_MedianLen.svg",width = 4,height = 5)
tot_392%>%
  group_by(Cond,FSP_ID)%>%
  summarise(FSPLEN=sum(Len))%>%
  summarise(FSPLEN=median(FSPLEN))%>%
  ggplot(aes(factor(Cond,levels=c("Lost","Kept","Gained")),FSPLEN))+
  #geom_line(color="grey",size=0.8,aes(group = 1),linetype="dashed")+
  #geom_point(stat="identity",shape=21,aes(fill=Cond),size=3,color="black",position = position_dodge())+
  geom_bar(stat="identity",aes(fill=Cond),color="black")+
  ggprism::theme_prism(base_size = 15)+xlab("")+ylab("Median FSP length")+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")+ylim(c(0,20))
dev.off()

ggsave("./Output/Fig_5I_Len.pdf",width = 4,height = 4)

ggsave("./Output/Fig_5I_Len.jpeg",width = 4,height = 4,dpi=300)



