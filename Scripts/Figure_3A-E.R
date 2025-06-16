library(ComplexHeatmap)
library(dplyr)
library(readxl)
library(ggplot2)

maf<-read.delim("./Data/maf_lof_or_pathogenic_CRC.tsv")
tmb<-read_xlsx("./Data/df_b2m_tmb_metachr.xlsx")
clinical<-read.delim("./Data/clinical_history_CRC.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)
mafpaired<-maf[maf$Tumor_Sample_Barcode%in%clinical$name,]
selectedgenes<-mafpaired[mafpaired$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]

mafuro<-read.delim("./Data/UC_maf_lof_or_pathogenic.tsv")
selectedgenesuro<-mafuro[mafuro$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]
clinicaluro<-read.delim("./Data/UC_clinical_data.txt")
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


col_fun <- c(
  # ";" = "white",           # Nessuna mutazione
  "B2M" = "skyblue3",            # Mutazione tipo 1
  "POLE"= "red",             # Mutazione tipo 2
  "POLD1"  = "green4",
  "MSH3"="khaki",# Mutazione tipo 3
  "MSH6"="lightcoral"
)



clinical_crc<-clinical[,c("patient_id","Metachronous_status","N_tumor","name")]
clinicaluro_sub<-clinicaluro[,c("PatientID","N_tumor","SampleID")]
clinicaluro_sub$Metachronous_status<-""

names(clinicaluro_sub)<-c("patient_id","N_tumor","name","Metachronous_status")
clinicaltot<-rbind(clinical_crc,clinicaluro_sub)


clinicaltot$Tumor<-"CRC"
clinicaltot[grepl("LS-",clinicaltot$name),]$Tumor<-"UC"
dim(clinicaltot)




matmatrixtot<- cbind(matmatrix,matmatrix2)


tmburo<-read.delim("./Data/TMB_Urothelial_cohort.tsv")
tmburo<-merge(tmburo,uromapping,by.x="Sample",by.y = "sample_pair")
names(tmburo)[c(6,5)]<-c("Tumor_Sample_Barcode","Tumor_mutation_burden")

tmbtot<-rbind(tmburo[,c(6,5)],
              data.frame("Tumor_Sample_Barcode"="LS-830",
                         "Tumor_mutation_burden"=21.17),
              tmb[,c(1,3)])


infosorted<-merge(clinicaltot,tmbtot,by.x="name",by.y = "Tumor_Sample_Barcode",sort = F)


urosummary<-read.delim("./Data/NOUS209FSP_UC.txt")
urosummary<-urosummary[,c("Patient","MNR_lookup","FSP_lookup")]
urosummary$Patient<-gsub("_blood","",urosummary$Patient)

crcsummary<-clinical[,c("name","MNR","NOUS209")]
names(urosummary)<-c("name","MNR","NOUS209")

addinfo<-rbind(crcsummary,urosummary)

infosorted<-merge(infosorted,addinfo,by="name",all.x=T)




# Conta il numero di valori diversi da ";"
count_non_semicolon <- colSums(matmatrixtot!= ";")

# Creiamo un dataframe con il nome della colonna e il numero di valori diversi da ";"
df <- data.frame(Colonna = colnames(matmatrixtot), Valori_Diversi = count_non_semicolon)


infosorted<-merge(infosorted,df,by.x="name",by.y="Colonna")
dim(infosorted)

infosorted<-infosorted%>%
  arrange(N_tumor,Valori_Diversi)
infosorted$MSI<-"MSI-H"
infosorted[infosorted$name%in%c("K18795-12","LS-130","LS-459_B","LS-1130"),]$MSI<-"Not confirmed MSI-H"

infosorted$color_TMB<-ifelse(infosorted$MSI=="MSI-H","#bd5140","grey50")
infosorted$color_NOUS209<-ifelse(infosorted$MSI=="MSI-H","#47A992","grey50")
infosorted$N_tumor<-paste0("Group ",infosorted$N_tumor)
top_annotation = HeatmapAnnotation(
  
  TMB = anno_barplot(infosorted$Tumor_mutation_burden,
                     gp = gpar(fill = infosorted$color_TMB),
                     
                     height = unit(3, "cm")),
  NOUS209 = anno_barplot(infosorted$NOUS209,gp = gpar(fill = infosorted$color_NOUS209),
                         height = unit(3, "cm")) ,
  
  gap = unit(5, "mm"),
  Group= infosorted$N_tumor,
  Tumor= infosorted$Tumor,
  col = list(
    Group = c("Group 1" = "#374E55", "Group 2" = "#DF8F44","Group 3"="#00A1D5"),
    Tumor = c("CRC" = "tomato3", "UC" = "khaki"),
    TMB=c("darkred"),
    NOUS209="#47A992"
  ),annotation_label = c("TMB","Nous-209","Group","Tumor"),
  annotation_name_side = "left"
)



mat_ordinata<-matmatrixtot[,infosorted$name]








jpeg("./Output/Fig_3A_Oncoplot.jpeg",res = 300,units = "px",width = 4800,height = 3000)
pdf("./Output/Fig_3A_Oncoplot.pdf",width = 48,height = 30)

svg("./Output/Fig_3A_Oncoplot.svg",width = 15,height = 15)


oncoPrint(mat_ordinata,
          alter_fun = list(
            background = function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill = "grey95", col = NA))
            },
            "B2M" =  function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill =col_fun ["B2M"], col = NA))
            },
            "POLE" =  function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill =col_fun ["POLE"], col = NA))
            },
            "POLD1" =  function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill =col_fun ["POLD1"], col = NA))
            },
            "MSH3" =  function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill =col_fun ["MSH3"], col = NA))
            },
            "MSH6" =  function(x, y, w, h) {
              grid.rect(x, y, w-unit(2, "pt"), h-unit(3, "pt"), 
                        gp = gpar(fill =col_fun ["MSH6"], col = NA))
            }),alter_fun_is_vectorized=F,
          width = unit(27, "cm"), height = unit(6, "cm"),
          col = col_fun,
          column_title = "",show_row_names = T,
          row_names_side = "left",column_order = infosorted$name,
          pct_side = "right",column_gap = unit(4, "mm"),
          
          bottom_annotation = top_annotation,
          column_split = infosorted$N_tumor, 
          row_names_gp = gpar(fontsize = 12),
          top_annotation = HeatmapAnnotation("column_barplot"=anno_oncoprint_barplot( border = F, # only MUT
                                                                                   height = unit(3, "cm"))
                                                                                
                                                                                   
                                          
          )
          )


dev.off()

################## Fig 3 B e C
mutations_per_sample <- apply(matmatrixtot, 2, function(col) sum(col != ";"))


MutDF<-as.data.frame(mutations_per_sample)
MutDF$Patient<-row.names(MutDF)

MutDF<-merge(MutDF,infosorted,by.x="Patient",by.y = "name")


geneCRC<-mafpaired[mafpaired$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]$Tumor_Sample_Barcode
geneURO<-mafuro[mafuro$Hugo_Symbol%in%c("MSH3","MSH6","POLD1","POLE"),]$Patient_ID

MutDF$MutSample<-""
MutDF[MutDF$mutations_per_sample==0,]$MutSample<-"No Mutation"
MutDF[MutDF$mutations_per_sample>=1,]$MutSample<-"One or\nmore Mutations"

table(MutDF$MutSample)
53/73
library(ggprism)

svg("./Output/Fig_3B_TMB_mutations.svg",height = 4,width = 4.5)
ggplot(MutDF,aes(MutSample,Tumor_mutation_burden))+
  geom_boxplot(outliers = F)+
  ggbeeswarm::geom_quasirandom(aes(fill=MutSample),shape=21,size=3)+theme_prism()+
  ggpubr::stat_compare_means(comparisons = list(c("No Mutation","One or\nmore Mutations")))+
  xlab("")+ggsci::scale_fill_locuszoom()+ylab("TMB")+
  theme(legend.position = "None")+ylim(c(0,150))
dev.off()


ggsave("./Output/Fig_3B_TMB_mutations.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_3B_TMB_mutations.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_3B_TMB_mutations.svg",dpi=300,
       width =1400,height = 1300,units = "px")



svg("./Output/Fig_3C_NOUS209_mutations.svg",height = 4,width = 4.5)
ggplot(MutDF,aes(MutSample,NOUS209))+
  geom_boxplot(outliers = F)+
  ggbeeswarm::geom_quasirandom(aes(fill=MutSample),shape=21,size=3)+theme_prism()+
  ggpubr::stat_compare_means(comparisons = list(c("No Mutation","One or\nmore Mutations")))+
  xlab("")+ggsci::scale_fill_locuszoom()+ylab("# Nous-209")+
  theme(legend.position = "None")

dev.off()

ggsave("./Output/Fig_3C_NOUS209_mutations.pdf",dpi=300,
       width =1400,height = 1300,units = "px")

ggsave("./Output/Fig_3C_NOUS209_mutations.jpeg",dpi=300,
       width =1400,height = 1300,units = "px")
ggsave("./Output/Fig_3C_NOUS209_mutations.svg",dpi=300,
       width =1400,height = 1300,units = "px")



##### Fig 3 D E 

MutDF_CRC<-MutDF[MutDF$Tumor=="CRC",]
MutDF_UC<-MutDF[MutDF$Tumor=="UC",]


MutDF_CRC<-merge(MutDF_CRC,clinical[,c("Stage","name")],by.x="Patient",by.y="name")
MutDF_CRC$Stage2<-ifelse(MutDF_CRC$Stage%in%c("II","III"),"II-III",MutDF_CRC$Stage)
MutDF_CRC$MutSample2<-ifelse(MutDF_CRC$MutSample%in%c("No Mutation"),"Absent","Present")

jpeg("./Output/Fig_3D_mosaicCRC.jpeg",res = 300,units = "px",width = 1500,
     height =1400 )
pdf("./Output/Fig_3D_mosaicCRC.pdf")
svg("./Output/Fig_3D_mosaicCRC.svg",height = 4,width =6 )

table(MutDF_CRC$Stage2,MutDF_CRC$MutSample2)%>%mosaicplot(cex=1.5,col=c("coral","#CA5FA6"),
                                                         main = "CRC")
dev.off()


MutDF_UC<-merge(MutDF_UC,clinicaluro[,c("Stage","SampleID","Invasion")],by.x="Patient",by.y="SampleID")
MutDF_UC$Stage2<-ifelse(MutDF_UC$Stage%in%c("T3","T4","T2"),"pT2-pT3-pT4","pTis-pTa-pT1")

MutDF_UC$Stage2<-factor(MutDF_UC$Stage2,levels = c("pTis-pTa-pT1","pT2-pT3-pT4"))
MutDF_UC$MutSample2<-ifelse(MutDF_UC$MutSample%in%c("No Mutation"),"Absent","Present")
MutDF_UC$Invasion<-factor(MutDF_UC$Invasion,
                          levels = c("non-MI","MI"))
jpeg("./Output/Fig_3E_mosaicUC.jpeg",res = 300,units = "px",width = 1500,
     height =1400 )
pdf("./Output/Fig_3E_mosaicUC.pdf")
svg("./Output/Fig_3E_mosaicUC.svg",height = 4,width =6 )
table(MutDF_UC$Invasion,MutDF_UC$MutSample2)%>%mosaicplot(cex=1.5,col=c("coral","#CA5FA6"),main = "UC")
dev.off()














