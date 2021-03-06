#Figure numbering might have been changed
#Data transformation and analysis
library(tidyverse)
library(readxl)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)

#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")

#save.image('Metabolomics_Celegans_FA.RData')
#load('Metabolomics_Celegans_FA.RData')

odir<-'Summary_Celegans_FA'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


strains<-c('OP50','OP50-MR','OP50-crp','nhr-49')

met.names<-read_xlsx('Metformin FA worm data.xlsx',sheet='Fat_names') %>%
  rename(Metabolite=`Lipid name`)

removals<-read_xlsx('Metformin FA worm data.xlsx',sheet='Removals') %>%
  mutate(SampleMet=paste(Sample,Metabolite))


met<-read_xlsx('Metformin FA worm data.xlsx',sheet='Complete_table') %>%
  gather(Metabolite,Conc,everything(),-c(Name:Sample)) %>%
  group_by(Sample) %>%
  mutate(SampleMean=mean(Conc)) %>%
  ungroup %>%
  mutate(ConcNorm=Conc/(SampleMean/mean(SampleMean)),
         logConc=log2(ConcNorm),
         Metabolite=ifelse(Metabolite=="C18:2w6:cis","C18:2w6 cis",Metabolite),
         Strain=factor(Strain,levels=strains)) %>%
  group_by(Group,Metabolite) %>%
  mutate(logConc_group=mean(logConc)) %>%
  ungroup %>%
  mutate(logConc_fill=ifelse(paste(Sample,Metabolite) %in% removals$SampleMet, logConc_group,logConc)) %>%
  left_join(met.names) %>%
  mutate_at(c('Sample','Name','Metformin_mM','Group','Replicate','Metabolite','Common name','Clean name'),as.factor)


met %>%
  write_csv(paste0(odir,"/Raw_data.csv"))



metinfo<-met %>%
  group_by(Sample,Group,Strain,Metformin_mM) %>%
  summarise %>%
  data.frame

rownames(metinfo)<-metinfo$Sample

  
PCAres<-met %>%
  filter(Strain %in%  c('OP50','OP50-crp')) %>%
  PCAprep("Sample","Metabolite","logConc_fill",metinfo)


ellipses<-PCAres$Ellipses
pcadata<-PCAres$pcadata
PC1prc<-PCAres$PC1prc
PC2prc<-PCAres$PC2prc

# write.csv(pcadata,paste(odir,'/PCA_data_OP50_OP50-MR.csv',sep=''))
# write.csv(pcaresult,paste(odir,'/PCA_variance_OP50_OP50-MR.csv',sep=''))

ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste0(odir,"/PCA_OP50_OP50-crp.pdf"),
             width=12,height=9)


#Heatmap
heatshape<-met %>%
  filter(Strain %in%  c('OP50','OP50-MR')) %>%
  select(`Clean name`,Sample,logConc) %>%
  spread(Sample,logConc) %>%
  data.frame



rownames(heatshape)<-heatshape$Clean.name
heatshape$Clean.name<-NULL


#Order anotation by heatmap colnames
hanot<-metinfo[colnames(heatshape),] %>%
  select(Strain,Metformin_mM)

ha<-HeatmapAnnotation(df=hanot, col = list(Metformin_mM=c('50'='black','0'='white')))

heatshape.sc<-t(scale(t(heatshape)))

#make anotated heatmap
Heatmap(heatshape.sc,name = 'Z-score',
              #col=bgg4,
              column_names_side = 'top',
              clustering_method_rows='ward.D2',
              #clustering_method_columns ='ward.D2',
              top_annotation = ha,
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              row_names_max_width = unit(10, "cm"))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_OP50-OP50-MR.pdf",sep=''),
             width=8,height=9, useDingbats=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=12,height=9, useDingbats=FALSE)


ggplot(met,aes(x=`Clean name`,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
  geom_point()+
  geom_text(aes(label=Replicate),color='black',size=2)+
  scale_y_continuous(breaks=seq(-3,8,by=1) )+
  ylab('log 2 Concentration, nmol')+

  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)

ggplot(met,aes(x=`Clean name`,y=Conc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
  geom_point()+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Strain,labeller = labeller(Strain = label_both),ncol = 5)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_Conc_by_treatment.pdf",sep=''),
             width=25,height=10, useDingbats=FALSE)


ggplot(met,aes(x=Strain,y=logConc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
  geom_point()+
  scale_y_continuous(breaks=seq(-5,15,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log 2 Concentration, nmol')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Clean.name,ncol = 5)

#,labeller = labeller(Metablite = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=20,height=20, useDingbats=FALSE)


#Concentrations by groups

sum.c<-met %>%
  group_by(Metabolite,`Common name`,`Clean name`,Strain,Metformin_mM,Group) %>%
  summarise(Mean=mean(logConc),
            SD=sd(logConc)) %>%
  mutate(Index=paste(Group,Metabolite),
         VarPrc=ifelse(is.na(SD) ,Inf, (2^(SD)-1)*100 ) ) %>%
  #Reorder, then preserve order with factor levels
  arrange(VarPrc)%>%
  data.frame %>%
  mutate(Index=factor(Index, levels=Index,labels=Index)) 


sum.c %>%
  filter(VarPrc>300)


ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,1000,by=25),limits=c(0,100))+
  ylab('In-group variation, %')+
  coord_flip()


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=30, useDingbats=FALSE)



#Linear modelling
sel.groups<-unique(met$Group)
contrasts<-read.contrasts2('!Contrasts_Celegans_FA_metabolomics.xlsx')
contrasts$Contrasts.table

#Can be extended to include extra info that will be passed to results
contrasts.desc<-contrasts$Contrasts.table %>%
  select(-c(Target:Reference))

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix


#Need to define grouping and variables
results.all<-met %>%
  group_by(Metabolite,`Common name`,`Clean name`) %>%
  do(hypothesise2(.,"logConc~0+Group",contr.matrix)) %>%
  getresults(contrasts.desc)



results<-results.all$results
results.cast<-results.all$cast
results.castfill<-results.all$castfull

results.int<-results.all$multi


write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)

cnt.vals<-results %>%
  select(-one_of('Metabolite','Clean name','Common name') ) %>%
  colnames(.)


#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 2
cbrks<-seq(-amp,amp,by=1)
#gradcols<-c('black','purple','purple')
maincomp<-'Interaction strength'


gradcols<-c('blue4','blue','gray80','red','red4')


results.int %>%
  filter(x_Contrast=='C_T-C_C' & y_Contrast_type=='Treatment' & z_Contrast_type=='Interaction' & y_Strain==z_Strain) %>%
  mutate(y_Strain=factor(y_Strain,levels=c('OP50-MR','OP50-crp','nhr-49')) ) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=z_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor)+
  geom_point()+
  geom_text_repel(aes(label=`Clean name`),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  xlab('Metformin effect on OP50 (N2)')+
  ylab('Metformin effect on other strain')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  theme(panel.grid.minor = element_blank())+
  facet_wrap(~y_Strain)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_treatment.pdf',sep = ''),
             width=24,height=6,useDingbats=FALSE)


VolcanoPlot<-function(data) {
  data %>%
    ggplot(aes(x=logFC,y=logFDR,color=Strain))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
    geom_point()+
    ylim(0,15)+
    geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(`Clean name`),'')),size=2)+
    facet_wrap(~Description)
}

results %>%
  filter(Contrast_type=="Treatment") %>%
  VolcanoPlot+
  ggtitle('Metformin treatment effect in different strains')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)

results %>%
  filter(Contrast_type %in% c("Bacterial mutant","Worm mutant") ) %>%
  VolcanoPlot +
  ggtitle('Differences between bacterial/worm strains')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_mutant.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)


results %>%
  filter(Contrast_type %in% c("Interaction") ) %>%
  VolcanoPlot +
  ggtitle('Interaction between treatment and strain in comparison to OP50 (N2)')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_interaction_treatment-strain.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)




#Heatmap for summary


unique(as.character(results$Description))

comparisons<-c("Treatment effect on OP50","Treatment effect on OP50-MR","Treatment effect on OP50-crp","Treatment effect for nhr-49",
               "Mutant difference for OP50-MR","Mutant difference for OP50-crp","Mutant difference for nhr-49")




fig2conts<-c('C_T-C_C','MR_T-MR_C','nhr49_T-nhr49_C','MR_C-C_C','nhr49_C-C_C')

fig6conts<-c('C_T-C_C','MR_T-MR_C','crp_T-crp_C','MR_C-C_C','crp_C-C_C')

contsel<-fig6conts

heatsum<-results %>%
  filter(Contrast %in% contsel) %>%
  select(`Clean name`,Description,logFC) %>%
  spread(Description,logFC) %>%
  data.frame


rownames(heatsum)<-heatsum$Clean.name
heatsum$Clean.name<-NULL


d<-dist(as.matrix(heatsum),method = "euclidean")
h<-hclust(d,method="ward.D2")
ordmet<-rownames(heatsum[h$order,])



amp<-4

minv<- -amp
maxv<- amp

nstep<-maxv-minv

clrbrks<-seq(-amp,amp,by=2)

brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)

clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


results %>%
  filter(Contrast %in% contsel) %>%
  mutate(`Clean name`=factor(`Clean name`,levels=ordmet)) %>%
  ggplot(aes(x=Description,y=`Clean name`))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  #theme_minimal()+
  geom_text(aes(label=as.character(FDRStars)))+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
  xlab("Comparison")+
  ylab('Fatty acid')+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))



dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_heatmap.pdf',sep = ''),
             width=5,height=8,useDingbats=FALSE)





