#Data transformation and analysis
library(xlsx)
library(plyr)
library(reshape2)
library(multcomp)
library(contrast)
library(car)

#Plotting and visualisation
library(gplots)
library(heatmap3)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(RColorBrewer)
library(plot3D)


library(gtools)

library(ellipse)

library(Unicode)

library(ggthemes)



devtools::install_github("PNorvaisas/PFun")
library(PFun)



library(grid)
library(gridExtra)





theme_Publication <- function(base_size=14) {
  
  (theme_foundation(base_size=base_size)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           #panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           #legend.position = "bottom",
           #legend.direction = "horizontal",
           #legend.key.size= unit(0.2, "cm"),
           #legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           #plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

theme_set(theme_Publication())


#theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/Metformin_project/Figures/Figure_2/")





odir<-'Joint_PCA'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")




getellipse<-function(x,y,sc=1) {
  as.data.frame(ellipse( cor(x, y),
                         scale=c(sd(x)*sc,sd(y)*sc),
                         centre=c( mean(x),mean(y)) ))
}



selcols<-c('PC1','PC2','Sample','Strain','Group','Metformin_mM','Dataset')


cel.FA.data<-read.csv('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_FA/PCA_data_OP50_OP50-MR.csv',row.names = 1)
cel.FA.data$Dataset<-'Lipidomics'
cel.FA.var<-read.csv('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_FA/PCA_variance_OP50_OP50-MR.csv',row.names = 1)
cel.FA.var$Dataset<-'Lipidomics'


cel.AA.data<-read.csv('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_AA_filtered/PCA_data_OP50_OP50-MR.csv',row.names = 1)
cel.AA.data$Dataset<-'Amino acids'
cel.AA.data<-rename(cel.AA.data,c('Drug_mM'='Metformin_mM'))

cel.AA.var<-read.csv('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_AA_filtered/PCA_variance_OP50_OP50-MR.csv',row.names = 1)
cel.AA.var$Dataset<-'Amino acids'




cel.RNAseq.PCA<-read.csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Results_1thrs_newannot/RNAseq_MDS_data.csv')


cel.RNAseq.info<-read.csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Pheno_data_new.csv')


cel.RNAseq.data<-cbind(cel.RNAseq.PCA,cel.RNAseq.info)
cel.RNAseq.data$Dataset<-'RNAseq'

datasets<-c('RNAseq','Amino acids','Lipidomics')


pcadata<-rbind(cel.RNAseq.data[,selcols],cel.FA.data[,selcols],cel.AA.data[,selcols])

pcadata$Metformin_mM<-as.factor(pcadata$Metformin_mM)

pcadata$Dataset<-factor(pcadata$Dataset,levels=datasets,labels=datasets)


AA.PC1<-round(cel.AA.var['Proportion of Variance',][[1]]*100,0)
AA.PC2<-round(cel.AA.var['Proportion of Variance',][[2]]*100,0)

FA.PC1<-round(cel.FA.var['Proportion of Variance',][[1]]*100,0)
FA.PC2<-round(cel.FA.var['Proportion of Variance',][[2]]*100,0)



ellipses<-ddply(pcadata,.(Dataset,Strain,Group,Metformin_mM), summarise, x=getellipse(PC1,PC2,0.75)$x,y=getellipse(PC1,PC2,0.75)$y ) 


ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  # xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  # ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  xlab('PC1')+
  ylab('PC2')+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=3,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(.~Dataset,scales = "free")

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Joint_PCA.pdf",sep=''),
             width=12,height=4, useDingbats=FALSE)




RNAseq<-ggplot(subset(pcadata,Dataset=='RNAseq'),aes(x=PC1,y=PC2,colour=Strain))+
  # xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  # ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  xlab('Leading logFC dimension 1')+
  ylab('Leading logFC dimension 2')+
  geom_path(data=subset(ellipses,Dataset=='RNAseq'), aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=3,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(.~Dataset,scales = "free")+
  theme(legend.position="none")


AA<-ggplot(subset(pcadata,Dataset=='Amino acids'),aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',AA.PC1,'% of variance',sep=''))+
  ylab(paste('PC2 - ',AA.PC2,'% of variance',sep=''))+
  geom_path(data=subset(ellipses,Dataset=='Amino acids'), aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=3,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  scale_y_continuous(breaks=seq(-20,20,by=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(.~Dataset,scales = "free")+
  theme(legend.position="none")

FA<-ggplot(subset(pcadata,Dataset=='Lipidomics'),aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',FA.PC1,'% of variance',sep=''))+
  ylab(paste('PC2 - ',FA.PC2,'% of variance',sep=''))+
  geom_path(data=subset(ellipses,Dataset=='Lipidomics'), aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=3,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  scale_y_continuous(breaks=seq(-20,20,by=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_grid(.~Dataset,scales = "free")+
  theme(legend.position="none")





grid.arrange(RNAseq,AA,FA,ncol=3)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Joint_PCA_withaxis.pdf",sep=''),
             width=12,height=4, useDingbats=FALSE)
