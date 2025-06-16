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















library(tidyverse)
library(ggforce)

table(ihc_stage$MMR_status,ihc_stage$IHC)




data <- matrix(c(15, 39, 0, 15), ncol = 2, byrow = TRUE)
colnames(data) <- c("Negative", "Positive")
rownames(data) <- c("MLH1", "MSH2/MSH6")

df <- as.data.frame(data) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = Negative:Positive, names_to = "Status", values_to = "Count") %>%
  mutate(Ring = ifelse(Status == "Negative", "Inner", "Outer"))


library(webr)


PieDonutC<-function (data, mapping, start = getOption("PieDonut.start", 
                                                      0), 
                     addPieLabel = TRUE, addDonutLabel = TRUE, 
                     showRatioDonut = TRUE, 
                     showRatioPie = TRUE, ratioByGroup = TRUE, 
                     showRatioThreshold = getOption("PieDonut.showRatioThreshold", 
                                                    0.02), labelposition = getOption("PieDonut.labelposition", 
                                                                                     2),
                     labelpositionThreshold = 0.1,
                     r0 = getOption("PieDonut.r0",                                                                                                                                                            0.3), 
                     r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2",                                                                                                                                                                                                                  1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                     color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                     showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                     pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                     explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                     family = getOption("PieDonut.family", "")) 
{
  (cols = colnames(data))
  if (use.labels) 
    data = addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                             "\n(", scales::percent(df$ratio), ")"), as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  mainCol = c("#d62728","#1f77b4")
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$ratio, 
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
                                    r = as.character(r1), start = "start1", end = "end1", 
                                    fill = pies), alpha = pieAlpha, color = color, data = df) + 
    transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"),color="black", data = df) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", vjust = "vjust"), 
                size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), color="black",data = df[df$ratio < 
                                                                                                      labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                                                                                                                                                        y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                                        vjust = "vjust"), size = pieLabelSize, data = df[df$ratio < 
                                                                                                                                                                                                           labelpositionThreshold, ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                                                                                                                                                                                                              y = "labely", label = "label"), size = pieLabelSize, 
                                                                                                                                                                                                                                                                   data = df[df$ratio >= labelpositionThreshold, ], 
                                                                                                                                                                                                                                                                   family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no", 
                                        explode = "focus"), alpha = donutAlpha, color = color, 
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no"), 
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3) + 
        geom_text(aes_string(x = "segxend", y = "segyend", 
                             label = "label", hjust = "hjust", vjust = "vjust"), 
                  size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely", 
                                      label = "label"), size = donutLabelSize, data = df3, 
                           family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                                                                                           labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                                                                                                                                             y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                             vjust = "vjust"), size = donutLabelSize, data = df3[df3$ratio1 < 
                                                                                                                                                                                                   labelpositionThreshold, ], family = family) + 
        geom_text(aes_string(x = "labelx", y = "labely", 
                             label = "label"), size = donutLabelSize, data = df3[df3$ratio1 >= 
                                                                                   labelpositionThreshold, ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                          size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                          label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid.newpage()
    print(p1, vp = viewport(height = 1, width = 1))
    print(p3, vp = viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}


library(moonBook)
pdf("../Paper_Definitivo/Output/PieDonut_4G_IHC.pdf")
PieDonutC(df[,c(1:3)],mapping = aes(Status,Gene,count=Count),title = "",
          donutLabelSize = 4,r0 = 0.4,r2=1.2,r1=0.9,showPieName = F)
dev.off()


