library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)

################# Fig 5 A 
follow<- read_excel("./Data/Followup_2_tumortypes.xlsx",
                    sheet = "mix")

View(follow)

merged_diff<-follow%>%
  group_by(patient_id)%>%
  arrange(age_at_colonoscopy)%>%
  mutate(diff=age_at_colonoscopy-min(age_at_colonoscopy))

merged_diff<-merged_diff%>%
  group_by(patient_id)%>%
  arrange(age_at_colonoscopy)%>%
  mutate(Range=max(age_at_colonoscopy)-min(age_at_colonoscopy))

ranking<-merged_diff[,c("Range","patient_id")]
ranking<-ranking[order(ranking$Range,decreasing = T),]
ranking<-unique(ranking)
ranking<-ranking[order(ranking$Range,decreasing = F),]
dim(ranking)
ranking$Position<-seq(1, 12)
merged_diff2<-merge(merged_diff,ranking,by="patient_id")


merged_diff2<-merged_diff2%>%
  group_by(patient_id)%>%
  mutate(diff_final=max(diff))



merged_diff2$Linea_color="gray"

for (pt in unique(c(197,798,651,237,392,507,767))){
  print(pt)
  ageprimary<-merged_diff2[merged_diff2$patient_id==pt&
                             merged_diff2$sample_id=="Primary"&
                             !is.na(merged_diff2$sample_id),]$age_at_colonoscopy
  agemeta<-merged_diff2[merged_diff2$patient_id==pt&
                          merged_diff2$sample_id=="Metachronous"&
                          !is.na(merged_diff2$sample_id),]$age_at_colonoscopy
  merged_diff2[merged_diff2$patient_id==pt&
                 merged_diff2$age_at_colonoscopy<agemeta&
                 merged_diff2$age_at_colonoscopy>=ageprimary,]$Linea_color="lightblue"
}



for (pt in unique(c(1013,1086))){
  print(pt)
  min_age=min(merged_diff2[merged_diff2$patient_id==pt&
                             !is.na(merged_diff2$sample_id)&
                             merged_diff2$sample_id=="Metachronous",]$age_at_colonoscopy)
  ageprimary<-merged_diff2[merged_diff2$patient_id==pt&
                             merged_diff2$sample_id=="Metachronous"&
                             !is.na(merged_diff2$sample_id)&
                             merged_diff2$age_at_colonoscopy==min_age,]$age_at_colonoscopy
  agemeta<-merged_diff2[merged_diff2$patient_id==pt&
                          merged_diff2$sample_id=="Metachronous"&
                          !is.na(merged_diff2$sample_id)&
                          merged_diff2$age_at_colonoscopy!=min_age,]$age_at_colonoscopy
  merged_diff2[merged_diff2$patient_id==pt&
                 merged_diff2$age_at_colonoscopy<agemeta&
                 merged_diff2$age_at_colonoscopy>=ageprimary,]$Linea_color="lightblue"
}




for (pt in unique(c(1122))){
  print(pt)
  
  agemeta<-merged_diff2[merged_diff2$patient_id==pt&
                          merged_diff2$Other=="UTUC/Bladder"&
                          !is.na(merged_diff2$Other),]$age_at_colonoscopy
  ageprimary<-merged_diff2[merged_diff2$patient_id==pt&
                             merged_diff2$sample_id=="Primary"&
                             !is.na(merged_diff2$sample_id),]$age_at_colonoscopy
  merged_diff2[merged_diff2$patient_id==pt&
                 merged_diff2$age_at_colonoscopy<agemeta&
                 merged_diff2$age_at_colonoscopy>=ageprimary,]$Linea_color="lightblue"
}



merged_diff2[merged_diff2$patient_id=="1122",]$Linea_color


# Version1
merged_diff2<-merged_diff2[merged_diff2$patient_id%in%
                             c(1122,1086,1013,798,767,651,392,237,197),]


## Version2
#merged_diff2<-merged_diff2[merged_diff2$patient_id%in%
#                            c(507,459),]


#merged_diff2<-merged_diff2[!(merged_diff2$patient_id==767&
#                 merged_diff2$age_at_colonoscopy>=46),]


merged_diff2[merged_diff2$No_available=="EC no available"&
               !is.na(merged_diff2$No_available),]$No_available<-"Other histology"
merged_diff2$patient_id<-factor(merged_diff2$patient_id,
                                levels = rev(c("1086",1013,798,197,1122,767,651,392,237)))

svg("./Output/Fig_5A_timeline.svg",width = 10,height = 6)
ggplot(merged_diff2,aes(x=diff, y=
                          factor(patient_id), group=patient_id))+
  geom_line(aes(color=Linea_color),size=2.5)+
  ## Primary-Metachr
  geom_point(merged_diff2[!is.na(merged_diff2$sample_id),],
             mapping=aes(x=diff, y=factor(patient_id), 
                         fill=sample_id), shape=21, stroke=1,size=5)+
  # Urinary
  geom_point(merged_diff2[!is.na(merged_diff2$Other),],
             mapping=aes(x = diff,
                         y = factor(patient_id),
                         fill=Other),shape=21,size=3,stroke=1,color="black")+
  # No available
  geom_point(merged_diff2[!is.na(merged_diff2$No_available),],
             mapping=aes(x=diff, y=factor(patient_id), 
                         col=No_available), shape=21, stroke=1,size=6)+
  geom_point(merged_diff2[!is.na(merged_diff2$Doublesample),],
             mapping=aes(x=diff, y=factor(patient_id), 
                         col=Doublesample), shape=21, stroke=2,size=6)+
  xlab("Years of follow up")+ylab("")+
  #scale_y_continuous(breaks =seq(1, 10, by = 2))+
  scale_shape_manual(values = c(0,3,24))+
  #scale_y_continuous(breaks =seq(1, 4, by = 1),
  #                  labels = c("Pt9","Pt4","Pt3","Pt7"))+
  ggprism::theme_prism(base_size = 15)+
  scale_color_manual(values = c("Synchronous Urinary"="tomato",
                                
                                "gray"="grey80",
                                "lightblue"="lightblue",
                                "CRC no available"="lightcoral",
                                "Other histology"="purple",
                                "Urinary no available"="khaki3"))+
  scale_fill_manual(values = c("Primary"="#1f77b4",
                               "Metachronous"="orange","UTUC/Bladder"="red"))+
  ylab("Patient ID")
dev.off()
ggsave("./Output/Fig_5A_timeline.pdf",width =8)


ggsave("./Output/Fig_5A_timeline.jpeg",width =8,dpi=300)





