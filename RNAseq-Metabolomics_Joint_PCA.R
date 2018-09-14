#Data transformation and analysis
library(tidyverse)



#devtools::install_github("PNorvaisas/PFun")
library(PFun)



library(grid)
library(gridExtra)


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


setwd("~/Dropbox/Projects/2015-Metformin/RNAseq")

#save.image('Metformin_Fig2_Joint_PCA.RData')
#load('Metformin_Fig2_Joint_PCA.RData')


odir<-'Summary_Celegans_Metf_DR'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



#STcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
STcols <- c("red3","steelblue")#colorRampPalette(c("red", "blue4"))(6)
names(STcols) <- c("OP50-C","OP50-MR")
STlab<-'Strain'




dfiles<-c('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_FA/PCA_data_OP50_OP50-MR.csv',
          '~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_AA_filtered/PCA_data_OP50_OP50-MR.csv')
varfiles<-c('~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_FA/PCA_variance_OP50_OP50-MR.csv',
            '~/Dropbox/Projects/2015-Metformin/Metabolomics/Summary_Celegans_AA_filtered/PCA_variance_OP50_OP50-MR.csv')
tp<-c('Fatty acids','Amino acids')


cel.RNAseq.PCA<-read.csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Results_1thrs_newannot/RNAseq_MDS_data.csv')
cel.RNAseq.info<-read.csv('~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Pheno_data_new.csv')



datacols<-c('Dataset','Sample','Group','Replicate','Strain','Metformin_mM','PC1','PC2')

cel.RNAseq<-cbind(cel.RNAseq.PCA,cel.RNAseq.info) %>%
  mutate(Dataset='RNAseq') %>%
  select(datacols)



cel.var<-data.frame(File=varfiles,Dataset=tp) %>%
  group_by(File,Dataset) %>%
  do(read_csv(as.character(.$File))) %>%
  filter(X1=='Proportion of Variance') %>%
  mutate(X1=NULL) %>%
  ungroup %>%
  select(Dataset:PC2) %>%
  add_row(Dataset = "RNAseq") %>%
  mutate_at(c('PC1','PC2'),function(x) as.integer(x*100) )



datasets<-c('RNAseq','Amino acids','Fatty acids')

pcadata<-data.frame(File=dfiles,Dataset=tp) %>%
  group_by(File,Dataset) %>%
  do(read_csv(as.character(.$File)) ) %>%
  mutate(Metformin_mM=ifelse(is.na(Metformin_mM),Drug_mM,Metformin_mM),
         Drug_mM=NULL,
         X1=NULL) %>%
  ungroup %>%
  select(datacols) %>%
  rbind(cel.RNAseq) %>%
  mutate(Dataset=factor(Dataset,levels=datasets),
         PC1=ifelse(Dataset=='Amino acids',-PC1,PC1),
         Strain=recode(Strain,"OP50"="OP50-C")) %>%
  mutate_at(c('Strain','Group','Metformin_mM'),as.factor)



pcadata$Strain

glimpse(pcadata)



ellipses<-pcadata %>%
  group_by(Dataset,Group,Strain,Metformin_mM) %>%
  do(getellipse(.$PC1,.$PC2,sc = 1))





PCAplot<-function(data,ellipses) {
  data %>%
  ggplot(aes(x=PC1,y=PC2))+
    geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),colour=Strain,linetype=Metformin_mM),size=1)+ 
    geom_point(aes(fill= ifelse(Metformin_mM==0,as.character(Strain), NA ), colour=Strain ),size=3,stroke=1,shape=21)+
    scale_fill_manual(name = STlab,values =STcols,na.value=NA,guide=FALSE)+
    scale_color_manual(name = STlab,values =STcols,guide=guide_legend(order=1))+
    scale_linetype_manual("Metformin,\nmM",values=c("0"=1,"50"=2),
                          guide=guide_legend(override.aes = list(shape=c(21,21),
                                                                 size=1,
                                                                 linetype=c(1,3),
                                                                 colour="black",
                                                                 fill=c(1,NA)) ))+
    scale_x_continuous(breaks=seq(-20,20,by=1))+
    scale_y_continuous(breaks=seq(-20,20,by=1))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    facet_grid(.~Dataset,scales = "free")+
    theme(legend.position="none")
}

# +
#   theme(legend.position="none")

aa.var<-cel.var %>%
  filter(Dataset=='Amino acids')

fa.var<-cel.var %>%
  filter(Dataset=='Fatty acids')


pcadata %>%
  PCAplot(ellipses)+
  xlab('PC1')+
  ylab('PC2')


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Joint_PCA.pdf",sep=''),
             width=12,height=4, useDingbats=FALSE)


RNAseq<-pcadata %>%
  filter(Dataset=='RNAseq') %>%
  PCAplot(ellipses %>% filter(Dataset=='RNAseq'))+
  xlab('Leading logFC dimension 1')+
  ylab('Leading logFC dimension 2')


AA<-pcadata %>%
  filter(Dataset=='Amino acids') %>%
  PCAplot(ellipses %>% filter(Dataset=='Amino acids'))+
  xlab(paste('PC1 - ',aa.var$PC1,'% of variance',sep=''))+
  ylab(paste('PC2 - ',aa.var$PC2,'% of variance',sep=''))+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  scale_y_continuous(breaks=seq(-20,20,by=2))

FA<-pcadata %>%
  filter(Dataset=='Fatty acids') %>%
  PCAplot(ellipses %>% filter(Dataset=='Fatty acids'))+
  xlab(paste('PC1 - ',fa.var$PC1,'% of variance',sep=''))+
  ylab(paste('PC2 - ',fa.var$PC2,'% of variance',sep=''))+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  scale_y_continuous(breaks=seq(-20,20,by=2))



pcas<-grid.arrange(RNAseq,AA,FA,ncol=3)
pcas

#

ggsave(pcas,file=paste0(odir,"/Joint_PCA_withaxis.pdf"),
             width=115,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")






PCAplot_clean<-function(data,ellipses) {
  data %>%
    ggplot(aes(x=PC1,y=PC2))+
    geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),colour=Strain,linetype=Metformin_mM),size=1)+ 
    geom_point(aes(fill= ifelse(Metformin_mM==0,as.character(Strain), NA ), colour=Strain ),size=3,stroke=1,shape=21)+
    scale_fill_manual(name = STlab,values =STcols,na.value=NA,guide=FALSE)+
    scale_color_manual(name = STlab,values =STcols,guide=guide_legend(order=1))+
    scale_linetype_manual("Metformin,\nmM",values=c("0"=1,"50"=2),
                          guide=guide_legend(override.aes = list(shape=c(21,21),
                                                                 size=1,
                                                                 linetype=c(1,3),
                                                                 colour="black",
                                                                 fill=c(1,NA)) ))+
    scale_x_continuous(breaks=seq(-20,20,by=1))+
    scale_y_continuous(breaks=seq(-20,20,by=1))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


pcadata %>%
  filter(Dataset=='RNAseq') %>%
  PCAplot_clean(ellipses %>% filter(Dataset=='RNAseq'))+
  xlab('Leading logFC dimension 1')+
  ylab('Leading logFC dimension 2')

ggsave(RNAseql,file=paste0(odir,"/RNAseq_PCA_legend_larger.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


