library(maftools)
maf<-read.delim("./Data/maf_lof_or_pathogenic_2211.tsv")
tmb<-read_xlsx("./Data/df_b2m_tmb_metachr.xlsx")
clinical<-read.delim("./Data/clinical_history_cancer.tsv")
clinical$name<-gsub("_tumor","",clinical$tumor_sample)
mafpaired<-maf[maf$Tumor_Sample_Barcode%in%clinical$name,]


mafuro<-read.delim("./Data/Urothalial_maf_lof_or_pathogenic.tsv")
clinicaluro<-read.delim("./Data/Nouscom_LS-UC17Dec.txt")
clinicaluro<-clinicaluro[clinicaluro$Location%in%c("Bladder","UTUC"),]
uromapping<-read.delim("./Data/Patient_ID_sampleCegat_pair_tumor_normal.txt")
uromapping$sample_pair<-paste0("S000021_",uromapping$sample_pair)
selectedgenesuro<-merge(mafuro,uromapping,by.x="Tumor_Sample_Barcode",
                        by.y = "sample_pair")
selectedgenesuro$Patient_ID.x
names(selectedgenesuro)[1]<-"Tumorpair"
names(selectedgenesuro)[115]<-"Tumor_Sample_Barcode"
required_columns <- c(
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Variant_Classification",
  "Variant_Type",
  "Tumor_Sample_Barcode"
)

infosorted$Tumor_Sample_Barcode<-infosorted$name
infosorted$TMB<-infosorted$Tumor_mutation_burden
maftotale<-rbind(mafpaired[,required_columns],
                 selectedgenesuro[,required_columns])

maftotale$Tumor_Sample_Barcode%>%unique()
totmaf<-read.maf(rbind(mafpaired[,required_columns],
                       selectedgenesuro[,required_columns]),clinicalData = infosorted)

head(infosorted)


head(rbind(mafpaired[,required_columns],
           selectedgenesuro[,required_columns]))


genes <- c(
  "MSH3", "MSH6", "POLD1", "MLH3",  "POLD3", 
  "POLE", "ERCC5", "EXO1",  "POLR2A", "PARG", "PARP1",
  "POLG", "CDK7", "ERCC6", "DDB1", "LIG3", "MBD4", "PARP4",
  "POLR2C", "TDG"
)

names(infosorted)



svg("Oncoplot_multigenes_TMB.svg",width = 10,height = 6)
oncoplot(totmaf,genes = genes,
         topBarData = "TMB",
         clinicalFeatures = c("N_tumor","Tumor"),
         sortByAnnotation = T
           )

dev.off()

svg("Oncoplot_multigenes_NOUS-209.svg",width = 10,height = 6)
oncoplot(totmaf,genes = genes,
         topBarData = "NOUS209",
         clinicalFeatures = c("N_tumor","Tumor"),
         sortByAnnotation = T
)
dev.off()














library(maftools)

# Calculate mutation burden per sample
mut_burden <- getSampleSummary(totmaf)

# Get clinical annotations
anno_df <- totmaf@clinical.data

# Merge mutation data with annotations
combined_df <- merge(mut_burden, anno_df, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode")
dim(combined_df)
# Sort by annotation (e.g., Subtype) and mutation count
combined_df <- combined_df[order(combined_df$Tumor_Sample_Barcode, combined_df$total, decreasing = FALSE), ]

# Extract ordered sample list
ordered_samples <- combined_df$Tumor_Sample_Barcode

# Plot oncoplot with custom sample order and annotation
oncoplot(maf = totmaf,
         clinicalFeatures = "N_tumor",
         sampleOrder = ordered_samples,genes = genes)











