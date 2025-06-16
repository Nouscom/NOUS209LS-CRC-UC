library(ggplot2)
library(dplyr)
########### Fig 4 A
pietot<-data.frame(Var1=c("B2m alone","B2M with other mutations"),
                   Freq=c(1,12))
df3 <- pietot%>% 
  mutate(
    cs = rev(cumsum(rev(Freq))),
    label=round(Freq/sum(Freq)*100,1), 
    pos = Freq/2 + lead(cs, 1),
    pos = if_else(is.na(pos), Freq/2, pos))

svg("./Output/Fig_4A_pie.svg",  width = 4, height = 4)
ggplot(data = pietot,
       aes(x = "", y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", color = "black") +

  coord_polar("y") +
  theme_void()+
  ggtitle("B2M ")+ 
  theme(legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5) 
  )+
  ggsci::scale_fill_bmj(alpha = 0.7)
dev.off()
ggsave("./Output/Fig_4A_pie.pdf", width = 8, height = 6)
ggsave("./Output/Fig_4A_pie.jpeg", width = 8, height = 6)
ggsave("./Output/Fig_4A_pie.svg", width = 8, height = 6)

######### Fig 4 B
library(plotly)
library(plotme)

as.sunburstDF <- function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF <- names(DF)
  
  if(is.data.table(DF)){
    DT <- copy(DF)
  } else {
    DT <- data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT <- names(DT)
  hierarchy_columns <- setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList <- list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns <- colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT <- unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT <- DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <- rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns <- setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}


library(readxl)
maf<-read.delim("./Data/maf_lof_or_pathogenic_CRC.tsv")
tmb<-read_xlsx("./Data/df_b2m_tmb_metachr.xlsx")
clinical<-read.delim("./Data/clinical_history_CRC.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)
mafpaired<-maf[maf$Tumor_Sample_Barcode%in%clinical$name,]
selectedgenes<-mafpaired[mafpaired$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]

mafuro<-read.delim("./Data/Urothalial_maf_lof_or_pathogenic.tsv")
selectedgenesuro<-mafuro[mafuro$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]
clinicaluro<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
clinicaluro<-clinicaluro[clinicaluro$Location%in%c("Bladder","UTUC"),]
uromapping<-read.delim("./Data/Patient_ID_sampleCegat_pair_tumor_normal.txt")
uromapping$sample_pair<-paste0("S000021_",uromapping$sample_pair)
selectedgenesuro<-merge(selectedgenesuro,uromapping,by.x="Tumor_Sample_Barcode",
                        by.y = "sample_pair")

matmatrix<-matrix(0,ncol =length(unique(mafpaired$Tumor_Sample_Barcode)),
                  nrow = length(c("MSH3","MSH6","POLD1","POLE")))
colnames(matmatrix)<-c(unique(mafpaired$Tumor_Sample_Barcode))
rownames(matmatrix)<-c("MSH3","MSH6","POLD1","POLE")

for (patient in colnames(matmatrix)){
  for (mut in rownames(matmatrix) ){
    
    
    if ( patient %in% selectedgenes$Tumor_Sample_Barcode){
      
      
      mut_present<-selectedgenes[selectedgenes$Tumor_Sample_Barcode==patient&
                                   selectedgenes$Hugo_Symbol==mut,]
      
      if (nrow(mut_present)>0){
        
        if (mut=="B2M"){
          matmatrix[mut,patient]<-paste0("B2M",";")
        }else if (mut=="POLD1"){
          matmatrix[mut,patient]<-paste0("POLD1",";")
        }else if (mut=="POLE"){
          matmatrix[mut,patient]<-paste0("POLE",";")
        }else if (mut=="MSH3"){
          matmatrix[mut,patient]<-paste0("MSH3",";")
        }
        else if (mut=="MSH6"){
          matmatrix[mut,patient]<-paste0("MSH6",";")
        }
        
        
      }
      
    }
    
  }
}

matmatrix[matmatrix=="0"]=";"

matmatrix2<-matrix(0,ncol =15,
                   nrow = length(c("MSH3","MSH6","POLD1","POLE")))
colnames(matmatrix2)<-unique(clinicaluro$SampleID)
rownames(matmatrix2)<-c("MSH3","MSH6","POLD1","POLE")



for (patient in colnames(matmatrix2)){
  for (mut in rownames(matmatrix2) ){
    
    
    if ( patient %in% selectedgenesuro$Patient_ID.y){
      
      
      mut_present<-selectedgenesuro[selectedgenesuro$Patient_ID.y==patient&
                                      selectedgenesuro$Hugo_Symbol==mut,]
      
      if (nrow(mut_present)>0){
        
        if (mut=="B2M"){
          matmatrix2[mut,patient]<-paste0("B2M",";")
        }else if (mut=="POLD1"){
          matmatrix2[mut,patient]<-paste0("POLD1",";")
        }else if (mut=="POLE"){
          matmatrix2[mut,patient]<-paste0("POLE",";")
        }else if (mut=="MSH3"){
          matmatrix2[mut,patient]<-paste0("MSH3",";")
        }
        else if (mut=="MSH6"){
          matmatrix2[mut,patient]<-paste0("MSH6",";")
        }
        
        
      }
      
    }
    
  }
}
matmatrix2[matmatrix2=="0"]=";"

matmatrixtot<- cbind(matmatrix,matmatrix2)

mutations_per_sample <- apply(matmatrixtot, 2, function(col) sum(col != ";"))

MutDF<-as.data.frame(mutations_per_sample)
MutDF$Patient<-row.names(MutDF)

MutDF<-merge(MutDF,infosorted,by.x="Patient",by.y = "name")
MutDF$MutSample<-""
MutDF[MutDF$mutations_per_sample==0,]$MutSample<-"No Mutation"
MutDF[MutDF$mutations_per_sample>=1,]$MutSample<-"One or more Mutations"


mut_b2m_crc<-mafpaired[mafpaired$Hugo_Symbol=="B2M",]$Tumor_Sample_Barcode
MutDF$B2M<-ifelse(MutDF$Patient%in%mut_b2m_crc,"B2M mut","B2M wt")



test<-MutDF[MutDF$N_tumor=="Group 1",]%>% 
  count(Tumor,  B2M,MutSample)
sunburstDF <- as.sunburstDF(test, value_column = "n", add_root = TRUE)
plot1<-plot_ly(data = sunburstDF, ids = ~ids, 
        labels= ~labels, 
        parents = ~parents,
        values= ~values, type='sunburst', 
        branchvalues = 'total',
        textinfo="label+percent parent")


save_image(plot1,"Fig_4B_piegroup1_kaleido.svg")

reticulate::py_install("kaleido")
test<-MutDF[MutDF$N_tumor=="Group 2",]%>% 
  count(Tumor,  B2M,MutSample)
sunburstDF <- as.sunburstDF(test, value_column = "n", add_root = TRUE)
svg("Fig_4B_piegroup2.svg")

plot1<-plot_ly(data = sunburstDF, ids = ~ids, 
        labels= ~labels, 
        parents = ~parents,
        values= ~values, type='sunburst', 
        branchvalues = 'total',
        textinfo="label+percent parent")

dev.off()
kaleido(plot1, "surface-plot.svg")

test<-MutDF[MutDF$N_tumor=="Group 3",]%>% 
  count(Tumor,  B2M,MutSample)
sunburstDF <- as.sunburstDF(test, value_column = "n", add_root = TRUE)
plot_ly(data = sunburstDF, ids = ~ids, 
        labels= ~labels, 
        parents = ~parents,
        values= ~values, type='sunburst', 
        branchvalues = 'total',
        textinfo="label+percent parent")



#########
maf<-read.delim("./Data/maf_lof_or_pathogenic_CRC.tsv")
maf<-maf[maf$Tumor_Sample_Barcode%in%clinical$name,]
n_distinct(maf$Tumor_Sample_Barcode)
ihc<-read_excel("./Data/IHC_data_CRC.xlsx")
ihc<-ihc[ihc$Tumor_Sample_Barcode%in%clinical$name,]
b2mmut<-unique(maf[maf$Hugo_Symbol=="B2M",]$Tumor_Sample_Barcode)

b2mmut%in%ihc$Tumor_Sample_Barcode

ihc$NGS<-ifelse(ihc$Tumor_Sample_Barcode%in%b2mmut,"Mut","Wt")
ihc$IHC<-""
ihc[ihc$b2m0loss1weak2strong3014125026==0,]$IHC<-"Negative"
ihc[ihc$b2m0loss1weak2strong3014125026%in%c(1,2,4),]$IHC<-"Positive"
ihc[ihc$b2m0loss1weak2strong3014125026%in%c(3,5,6),]$IHC<-"Het"


ihc_uro<-read_xlsx("./Data/IHC_data_UC.xlsx")
ihc_uro$IHC<-ihc_uro$B2M_results
ihc_uro[ihc_uro$B2M_results=="Heterogeneous",]$IHC<-"Het"
ihc_uro$NGS<-"Wt"

names(ihc_uro)[5]<- "Tumor_Sample_Barcode"
ihc_tot<-rbind(ihc[,c("Tumor_Sample_Barcode","IHC","NGS")],
               ihc_uro[,c("Tumor_Sample_Barcode","IHC","NGS")])

ihc_tot[!duplicated(ihc_tot$Tumor_Sample_Barcode),]%>%dim()

pdf("./Output/Fig_4C_mosaic.pdf",width = 10,height = 5)

jpeg("./Output/Fig_4C_mosaic.jpeg",width = 15,height = 10,units = "cm",res = 300)
svg("./Output/Fig_4C_mosaic.svg",width = 7.5,height = 4)
mosaicplot(table(factor(ihc_tot$IHC,
                        levels=c("Negative","Het","Positive")),ihc_tot$NGS),
           color=c("springgreen3","khaki"),main = "",cex=1.2)
dev.off()


table(factor(ihc_tot$IHC,
             levels=c("Negative","Het","Positive")),ihc_tot$NGS)


mergemut<-merge(maf[maf$Hugo_Symbol=="B2M",],ihc,by="Tumor_Sample_Barcode")
library(ggpubr)
svg("./Output/Fig_4D_allelefreq.svg",width = 3.5,height = 4)
mergemut%>%
  ggplot(aes(IHC,t_AF))+
  geom_boxplot()+
  geom_point(aes(fill=IHC),shape=21,size=3)+stat_compare_means(comparisons = list(c("Positive","Negative")))+
  ggsci::scale_fill_jama()+theme_prism(base_size = 15)+ylab("Allele Frequency")+ylim(c(0,0.8))+xlab("")+
  theme(legend.position = "None")
dev.off()

ggsave("./Output/Fig_4D_allelefreq.pdf",width = 8,height = 7)
ggsave("./Output/Fig_4D_allelefreq.svg",width = 8,height = 7)
ggsave("./Output/Fig_4D_allelefreq.jpeg",width = 8,height = 7)









#############################################
library(readxl)
library(dplyr)
clinical_crc<-read.delim("./Data/clinical_history_CRC.tsv")
ihc_crc<-read_xlsx("./Data/IHC_data_CRC.xlsx")  

clinical_crc$name<-gsub("_tumor","",clinical_crc$tumor_sample)


Crc_IHC<-merge(ihc_crc,clinical_crc,by.x="Tumor_Sample_Barcode",by.y="name")

Crc_IHC$Stage2<-ifelse(Crc_IHC$Stage%in%c("II","III"),"II-III",Crc_IHC$Stage)
Crc_IHC$IHC<-""
Crc_IHC[Crc_IHC$b2m0loss1weak2strong3014125026==0,]$IHC<-"Negative"
Crc_IHC[Crc_IHC$b2m0loss1weak2strong3014125026%in%c(1,2,4),]$IHC<-"Positive"
Crc_IHC[Crc_IHC$b2m0loss1weak2strong3014125026%in%c(3,5,6),]$IHC<-"Heterogeneous"
subCRC<-Crc_IHC[,c("Tumor_Sample_Barcode","IHC","Stage2","MMR_status")]
nous_ihc<-read.delim("./Data/Summary_MNR_NOUS209FSP_UROTHELIAL_LOOKUP_hg38.txt")
nous_ihc$Patient<-gsub("_blood","",nous_ihc$Patient)
uro<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
ihc_uro<-read_xlsx("./Data/Nouscom LS-UC_B2M_Staining results_20Maggio_NEW.xlsx")



IHC_URO<-merge(uro,ihc_uro,by.x="Sample",by.y ="ID" )

dim(IHC_URO)

nous_ihc<-read.delim("./Data/NOUS209FSP_CRC.txt")
nous_ihc$Patient<-gsub("_blood","",nous_ihc$Patient)
IHC_URO<-merge(IHC_URO,nous_ihc,by.x="SampleID",by.y="Patient")

IHC_URO$Stage2<-ifelse(IHC_URO$Stage%in%c("Tis","Ta","T1"),"Tis-Ta-T1",
                       "T2-T3-T4")
subURO<-IHC_URO[,c("PatientID","B2M result","Stage2","MMR.IHC")]
names(subURO)<-c("Tumor_Sample_Barcode","IHC","Stage2","MMR_status")

subCRC$Stage<-ifelse(subCRC$Stage2=="II-III","Advanced","Early")
subURO$Stage<-ifelse(subURO$Stage2=="T2-T3-T4","Advanced","Early")
subURO$MMR_status<-ifelse(subURO$MMR_status=="MLH1 loss","MLH1","MSH2 or MSH6")
subURO<-subURO[subURO$IHC!="TO DO",]
subURO$IHC<-ifelse(subURO$IHC=="negative/loss","Negative","Positive")


subCRC$Cancer<-"CRC"
subURO$Cancer<-"UC"
ihc_stage<-rbind(subCRC,subURO)

table(ihc_stage[ihc_stage$IHC=="Negative",]$Stage,
      ihc_stage[ihc_stage$IHC=="Negative",]$Cancer)


pietot3<-table(ihc_stage[ihc_stage$IHC=="Negative",]$Stage)%>%as.data.frame()
df2 <- pietot3%>% 
  mutate(
    cs = rev(cumsum(rev(Freq))),
    label=round(Freq/sum(Freq)*100,1), 
    pos = Freq/2 + lead(cs, 1),
    pos = if_else(is.na(pos), Freq/2, pos))
P1<-ggplot(data = pietot3,
           aes(x = "", y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", color = "black") +
  #ggrepel::geom_label_repel(aes(y = pos, label = paste0(Var1,"\n",label, "%")),
   #                         data = df2, size=6, show.legend = F, nudge_x = 1) +
  coord_polar("y") +
  theme_void()+
  theme(legend.title = element_blank())+
  ggtitle("Negative B2M IHC")+ # Titolo del grafico
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5,size = 18) # Impostazioni per il titolo
  )+scale_fill_manual(values = c("violet","orange"))




svg("./Output/Fig_4F_IHC_Stage.svg",width = 4,height = 5)
print(P1)
dev.off()


ggsave("../../Paper_Definitivo/Output/Fig_IHC_Stage.pdf", plot = P1, width = 8, height = 6)
ggsave("../../Paper_Definitivo/Output/Fig_IHC_Stage.svg", plot = P1, width = 8, height = 6)
ggsave("../../Paper_Definitivo/Output/Fig_IHC_Stage.jpeg",dpi = 300, plot = P1, width = 8, height = 6)





