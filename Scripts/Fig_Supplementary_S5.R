library(ggplot2)
library(ggprism)

 


plot1 = reac.csv('All_Prediction_data_MHC_II.txt', sep='\t')

output_dir = 'Plot/'
dir.create(output_dir)

barplot = ggplot(plot1[!plot1$Pt%in%c( '392_UM'),],aes(factor(Pt),Predictions,fill=Cond))+
  geom_bar(stat="identity",color="black",position = position_dodge())+
  ggprism::theme_prism()+
  xlab("")+
  ylab("# Good Binders")+
  #scale_y_continuous(breaks = seq(0,1600, by=200))+
  scale_fill_manual(values = c(Lost='#1f77b4',Kept='#2ca02c',Gained="#ff7f0e"))+
  theme(legend.position = "None")+ylim(c(0,1500))


png(paste0(output_dir, "Barplot_Class2.png"), res=300,width = 2400,height =900)
print(barplot)
dev.off()


pdf(paste0(output_dir, "Barplot_Class2.pdf"),width= 14, height = 6)
print(barplot)
dev.off()

svg(paste0(output_dir, "Barplot_Class2.svg"),width= 8, height = 4)
print(barplot)
dev.off()



library(ggsignif)

boxplot = ggplot(plot1[!plot1$Pt%in%c("392_UM"),],aes(Cond,Predictions,))+
  geom_boxplot(outliers = F)+
  geom_line(aes(group=factor(Pt)),linetype = "dashed",color="gray",linewidth=1)+
  geom_point(shape=21,size=3,aes(fill = factor(Pt)))+
  ggprism::theme_prism()+xlab("")+ylab("# Good Binders")+
  #theme(legend.position = "None")+
  scale_y_continuous(breaks = seq(0,2000, by=500), limits=c(0,2000))+
  geom_signif(comparisons = list(c("Lost","Kept"),c("Kept","Gained"),
                                 c("Lost","Gained") ),
              test.args = c(method = wilcox.test, paired=T),
              step_increase = 0.1)


png(paste0(output_dir, "Boxplot_Class2.png"),res=300,width = 1800,height =1500)
print(boxplot)
dev.off()



pdf(paste0(output_dir, "Boxplot_Class2.pdf"),width= 10, height = 8)
print(boxplot)
dev.off()


svg(paste0(output_dir, "Boxplot_Class2.svg"),width= 6, height = 4)
print(boxplot)
dev.off()