#Data transformation and analysis

library(tidyverse)
library(readxl)


library(ggrepel)
library(RColorBrewer)

# library(ComplexHeatmap)
# library(circlize)

library(ggthemes)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_Publication())

#theme_set(theme_light())
# theme_update(panel.background = element_rect(colour = "black"),
#              axis.text = element_text(colour = "black"))



setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")

#load('Metabolomics_Ecoli_HMT.RData')
#save.image('Metabolomics_Ecoli_HMT.RData')

odir<-'Summary_Ecoli_HMT'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


# 
# .simpleCap <- function(x) {
#   s <- strsplit(x, " ")[[1]]
#   paste(toupper(substring(s, 1, 1)), substring(s, 2),
#         sep = "", collapse = " ")
# }




met.raw<-read_xlsx('HMT_Ecoli_complete/Raw_data.xlsx',sheet='Raw')


met.info<-read_xlsx('HMT_Ecoli_complete/Raw_data.xlsx',sheet='Samples')



met.TCA<-read_xlsx('HMT_Ecoli_complete/Raw_data.xlsx',sheet='TCA')




lmfill<-function(data,formula) {
  model<-lm(as.formula(formula),data=data)
  data$Prediction<-predict(model,data)
  return(data)
}




#Calculate fill concentrations from supplement data
met.TCAf<-met.TCA %>%
  rename(Metabolite=`Compound name`) %>%
  filter(Metabolite!="Isocitric acid") %>%
  gather(Sample_ID,Value,`OP50-1-C`:`CRP-OE-3-100`) %>%
  mutate(Value=as.numeric(Value),
         Metabolite=as.character(Metabolite),
         Metabolite=str_replace(Metabolite,"Malonyl CoA","Malonyl-CoA") ) %>%
  spread(ID,Value) %>%
  group_by(Metabolite) %>%
  do(lmfill(.,"Conc ~ RelArea"))

met.TCAf %>%
  filter(Prediction<0 ) 

met.TCAf %>%
  filter(Prediction>0 & is.na(Conc) )

#Plot supplement TCA cycle predictions
met.TCAf %>%
  gather(Measure,Value,Conc,Prediction) %>%
  ggplot(aes(x=RelArea,y=Value,color=Measure,shape=Measure))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  geom_point()+
  geom_line()+
  facet_wrap(~Metabolite,scale='free') 


met.TCAp<-met.TCAf %>%
  filter(Prediction>0 & is.na(Conc))%>%
  select(-Conc,-RelArea)

met.TCAp

grporder<-c("C_0_N","C_50_N",
            "C_0_G","C_50_G",
            "CRP_0_N","CRP_50_N",
            "C_0_I50","oeCRP_0_I50","oeCRP_0_I100")

strn.order<-c('OP50','CRP','oeCRP')
spl.order<-c('None','Glucose','IPTG50','IPTG100')

sgrp.order<-c('OP50 None','OP50 Glucose','CRP None','OP50 IPTG50','oeCRP IPTG50','oeCRP IPTG100')

met.all<-met.raw %>%
  rename(Metabolite=`Pathway Label`,
         Metabolite_full=`Compound name`,
         KEGG_ID=`KEGG ID`,
         HMDB_ID=`HMDB ID`) %>%
  gather(Sample_ID,Conc,`OP50-1-C`:`CRP-OE-3-100`) %>%
  left_join(met.info) %>%
  left_join(met.TCAp,by=c('Sample_ID','Metabolite')) %>%
  mutate(Conc=as.numeric(Conc),
         Conc=ifelse(is.na(Conc) & !is.na(Prediction),Prediction,Conc),
         Conc_log=log2(Conc),
         Metformin_mM=as.factor(Metformin_mM),
         SGroup=paste(Strain,Supplement),
         SGroup=factor(SGroup,levels=sgrp.order,labels=sgrp.order),
         Group=factor(Group,levels=grporder,labels=grporder),
         Strain=factor(Strain,levels=strn.order,labels=strn.order),
         Supplement=factor(Supplement,levels=spl.order,labels=spl.order) ) %>%
  select(Sample_ID,Replicate:Group,SGroup,Strain:Amount_ODmL,Metabolite_ID:HMDB_ID,everything(),-Prediction)
  
head(met.all)



#grepl('CoA',Metabolite)

met.all %>%
  filter(!is.na(Prediction))


write.csv(met.all,paste(odir,'/Raw_data_filled.csv',sep=''),row.names = FALSE)


met.clean<-met.all %>%
  tbl_df %>%
  filter(ID!='C_0_N_4') %>%
  filter(Metabolite_ID!='-') %>%
  #Remove missing metabolites
  group_by(Metabolite) %>%
  mutate(All_groups_missing=all(is.na(Conc))) %>%
  ungroup %>%
  filter(!All_groups_missing) %>%
  #Create filler
  group_by(Metabolite,Group) %>%
  mutate(All_obs_missing=all(is.na(Conc)),
         Fill_conc=mean(Conc,na.rm = TRUE),
         Fill_conc_log=mean(Conc_log,na.rm = TRUE)) %>%
  ungroup %>%
  mutate(Filled_conc=ifelse(is.na(Conc), Fill_conc,Conc),
         Filled_conc_log=ifelse(is.na(Conc_log), Fill_conc_log,Conc_log))


glimpse(met.clean)
head(met.clean)

View(met.clean)



Missing.mets<-met.clean %>%
  group_by(Metabolite) %>%
  summarise(Missing=sum(is.na(Conc)) ) %>%
  filter(Missing>0) %>%
  arrange(desc(Missing))


met.clean %>%
  group_by(Sample_ID) %>%
  summarise(Missing=sum(is.na(Conc)) ) %>%
  filter(Missing>0) %>%
  arrange(desc(Missing))


missing.obs<-met.clean %>%
  group_by(Group,Metabolite) %>%
  summarise(Missing=sum(is.na(Conc)) ) %>%
  group_by(Metabolite) %>%
  mutate(Missing_total=sum(Missing)) %>%
  ungroup %>%
  spread(Group,Missing)%>%
  arrange(desc(Missing_total))


View(Missing.mets)
head(missing.obs)
View(missing.obs)



write.csv(missing.obs,paste(odir,'/Missing_observations.csv',sep=''),row.names = FALSE)


#Remove missing metabolites

rm.mets<-missing.obs %>%
  filter(Missing_total>25)

#Remove metabolites that are missing
met.sel<-met.clean %>%
  filter(!Metabolite %in% rm.mets$Metabolite)


glimpse(met.sel)


bioinfo <- met.sel %>%
  group_by(ID,Group,Strain,Supplement,SGroup,Metformin_mM) %>%
  summarise %>%
  data.frame

rownames(bioinfo)<-bioinfo$ID


#PCA plots
# HC<-PCAres$HC
# pca<-PCAres$pca
# pcaloadings<-PCAres$Loadings
# pcashape=PCAres$pcashape

PCAplot<-function(PCAres) {
  ellipses<-PCAres$Ellipses
  pcadata<-PCAres$pcadata
  ggplot(pcadata,aes(x=PC1,y=PC2,colour=SGroup))+
    xlab(paste('PC1 - ',PCAres$PC1prc,'% of variance',sep=''))+
    ylab(paste('PC2 - ',PCAres$PC2prc,'% of variance',sep=''))+
    geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
    geom_point(aes(fill=factor( ifelse(Metformin_mM==0,SGroup, NA ) ) ),size=5,stroke=1,shape=21)+
    scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
    scale_fill_discrete(na.value=NA, guide="none")+
    guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
    geom_text_repel(aes(label=ID),size=2,color='black')+
    labs(colour='Strain & Supplement')+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


metcomplete<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc))) %>%
  ungroup %>%
  filter(!Filling_missing)

#Original - HMT linear scale filled
metcomplete %>%
  PCAprep('ID','Metabolite','Filled_conc',bioinfo) %>%
  PCAplot


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_All_HMT_linear_scale.pdf",sep=''),
             width=12,height=9)

#Vanilla filled - All
metcomplete %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_All.pdf",sep=''),
             width=12,height=9)



#No IPTG, oe filled
metcomplete %>%
  filter(!Group %in% c("C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_noIPTG_oe.pdf",sep=''),
             width=12,height=9)


#Control, IPTG, oe filled
metcomplete %>%
  filter(Group %in% c("C_0_N","C_50_N","C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_IPTG_oe.pdf",sep=''),
             width=12,height=9)

#No Glucose
metcomplete %>%
  filter(!Supplement=='Glucose') %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_NoGlucose.pdf",sep=''),
             width=12,height=9)


#What Filipe asked
met.sel %>%
  filter(!Supplement=='Glucose' & Strain %in% c('OP50','CRP')) %>%
  PCAdata %>%
  PCAplot

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_WhatFilipeAsked.pdf",sep=''),
             width=12,height=9)


heatcomplete<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Completeness=sum(!is.na(Conc_log))*100/n() )%>%
  ungroup %>%
  #At least 50% observations per metabolite
  filter(Completeness>50)


cols<-list(Metformin_mM = c('0' = "white", '50' = "Black"),
           Strain=c("OP50"="red3","CRP"="blue","oeCRP"="purple"),
           Supplement=c("None"="white","Glucose"="green2","IPTG50"="turquoise","IPTG100"="turquoise3") )
           
hinfo<-bioinfo %>%
  select(Strain, Metformin_mM,Supplement)

#Complete
heatcomplete %>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_All.pdf",sep=''),
             width=17,height=12, useDingbats=FALSE)

#Heatmap - No Glucose
heatcomplete %>%
  filter(!Supplement=='Glucose') %>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_noGlucose.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)

#Heatmap - no IPTG, oe
heatcomplete %>%
  filter(!Group %in% c("C_0_I50","oeCRP_0_I50","oeCRP_0_I100"))%>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_noIPTG_oe.pdf",sep=''),
             width=14,height=12, useDingbats=FALSE)


#Heatmap - IPTG, oe
heatcomplete %>%
  filter(Group %in% c("C_0_N","C_50_N","C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_IPTG_oe.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)
 

#Heatmap - WhatFilipeAsked
heatcomplete %>%
  filter(!Supplement=='Glucose' & Strain %in% c('OP50','CRP')) %>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_WhatFilipeAsked.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)










met.clean %>%
  #filter(Metabolite %in% c('Malonyl-CoA','Fumaric acid','Malic acid')) %>%
  ggplot(aes(x=SGroup,y=Conc_log,color=Metformin_mM))+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log2 Concentration, pmol/ODmL')+
  xlab('Group & Supplement')+
  labs(color='Metformin, mM')+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Metabolite,ncol = 4,scales='free_y')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite_TCA.pdf",sep=''),
             width=12,height=6, useDingbats=FALSE)



ggplot(met.clean,aes(x=SGroup,y=Conc,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  #scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('Concentration, pmol/ODmL')+
  xlab('Group & Supplement')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Metabolite,ncol = 4,scales='free_y')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_Conc_by_metabolite.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)



#Histograms
ggplot(met.clean,aes(x=Conc_log))+
  ggtitle('Distribution of log2 Concentration, pmol/ODmL')+
  geom_density(aes(y=..scaled..),fill='red',alpha=0.5)+
  geom_rug(aes(color=SGroup,linetype=Metformin_mM))+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  ylab('Scaled density')+
  xlab('log2 Concentration, pmol/ODmL')+
  labs(color='Strain & Supplement',
       linetype='Metformin, mM')+
  facet_wrap(~Metabolite,ncol = 4,scales='free_x')


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Density_logConc_by_metabolite.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)



ggplot(met.clean,aes(x=Conc))+
  ggtitle('Distribution of Concentration, pmol/ODmL')+
  geom_density(aes(y=..scaled..),fill='red',alpha=0.5)+
  geom_rug(aes(color=SGroup,linetype=Metformin_mM))+
  ylab('Scaled density')+
  xlab('Concentration, pmol/ODmL')+
  labs(color='Strain & Supplement',
       linetype='Metformin, mM')+
  facet_wrap(~Metabolite,ncol = 4,scales='free_x')

#aes(y=..scaled..),

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Density_Conc_by_metabolite.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)


#Summarise variance

sum.c<-met.clean %>%
  group_by(Metabolite,Group,Strain,Metformin_mM,Supplement) %>%
  summarise_at(vars(Conc,Conc_log),funs(SD=sd(.,na.rm = TRUE),Mean=mean(.,na.rm = TRUE))) %>%
  mutate(Index=paste(Group,Metabolite),
         VarPrc=ifelse(is.na(Conc_log_SD) ,Inf, (2^(Conc_log_SD)-1)*100 ) ) %>%
  arrange(VarPrc) %>%
  data.frame %>%
  mutate(Index=factor(Index, levels=Index,labels=Index))


sum.c %>%
  filter(VarPrc>300)


ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,1000,by=50),limits=c(0,400))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=10,height=150, useDingbats=FALSE)




#Linear modelling

contrasts<-read.contrasts('!Contrasts_Ecoli_HMT_metabolomics.xlsx')
strainlist<-c('OP50','CRP','oeCRP')

contrasts.desc<-contrasts$Contrasts.table %>%
  mutate(Strain=factor(Strain,levels=strainlist,labels=strainlist)) %>%
  select(-c(Target,Reference))
  

contr.matrix<-contrasts$Contrasts.matrix


contrasts.desc
contr.matrix

odir<-'Summary_Ecoli_HMT/NoGlucose_stats'
results.all<-met.clean %>%
  #filter(Metabolite=='Gly') %>%
  filter(Supplement!='Glucose') %>%
  group_by(Metabolite,Metabolite_ID,Metabolite_full,KEGG_ID,HMDB_ID) %>%
  do(hypothesise(.,"Conc_log~0+Group",contr.matrix)) %>%
  getresults(contrasts.desc)

#Separate results
results<-results.all$results

results.cast<-results.all$cast 
results.castfull<-results.all$castfull

results.multi<-results.all$multi


results.summary<-results %>%
  select(Metabolite,Contrast,logFC) %>%
  group_by(Metabolite,Contrast) %>%
  summarise(Estimate=!is.na(logFC)) %>%
  spread(Contrast,Estimate)

results.summary

results %>%
  group_by(Contrast) %>%
  summarise(Total=n(),
            Up=sum(logFC>0 & FDR<=0.05,na.rm=TRUE),
            Down=sum(logFC<0 & FDR<=0.05,na.rm=TRUE),
            All=sum(FDR<=0.05,na.rm=TRUE))


write.csv(results.summary,paste(odir,'/Missing_comparisons.csv',sep=''),row.names = FALSE)


#Consistency checks
met.check<-met.clean %>%
  group_by(Metabolite,Group) %>%
  summarise(Mean=mean(Conc_log,na.rm=TRUE),SD=sd(Conc_log,na.rm=TRUE))

met.check %>%
  filter(Metabolite=='Gly')


results %>%
  filter(Metabolite=='Gly') %>%
  select(Contrast,logFC,p.value,FDR)


write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)


# 
# 
# allcontrs<-unique(as.character(results$Contrast))
# allmets<-unique(as.character(results$Metabolite))



resultsmin<-results %>%
  select(Metabolite,Contrast,logFC:logFDR)


results.multi2<-multiplex(results,c("Metabolite","Metabolite_ID","Metabolite_full","KEGG_ID","HMDB_ID"),dim=2)




head(results.exp)
#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 8
cbrks<-seq(-amp,amp,by=2)
#gradcols<-c('black','purple','purple')
maincomp<-'Interaction strength'


gradcols<-c('blue4','blue','gray80','red','red4')


results.multi2 %>%
  filter(x_Contrast=='C_Metf' &
           y_Contrast %in% c('CRP_Metf','CGlu_Metf','dCRP', 'C_Glu','oeCRP50','oeCRP100')) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=y_logFC-x_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=ifelse(y_FDR<0.05 & x_FDR<0.05, as.character(Metabolite),"" ) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  ggtitle('Effects of various factors compared against metformin effect on OP50',
          subtitle = 'Metabolites with FDR<0.05 in both of the effects are shown')+
  xlab('Metformin effect on OP50')+
  ylab('Other effects')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~y_Description,ncol=2)+
  theme(panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_Various_vs_OP50_Meft.pdf',sep = ''),
             width=18,height=18,useDingbats=FALSE)



#Volcano plots
VolcanoPlot<-function(data,selmets=c()){
  plot<-ggplot(data,aes(x=logFC,y=logFDR,color=Strain))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
    geom_point()+
    scale_y_continuous(breaks=seq(0,20,by=2))+
    scale_x_continuous(breaks=seq(-10,10,by=2))+
    facet_wrap(~Description,ncol = 2)
  if (length(selmets)==0){
    plot<-plot+geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)
  } else {
    plot<-plot+geom_text_repel(aes(label=ifelse(Metabolite %in% selmets & FDR<0.05,as.character(Metabolite),'')),size=2)
  }
  
  return(plot)
}



results %>%
  filter(Contrast_type=='Treatment' ) %>%
  VolcanoPlot+
  ggtitle('Metformin treatment effect in different strains')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_treatment.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)


results %>%
  filter(Contrast_type=='Bacterial mutant' ) %>%
  VolcanoPlot+
  ggtitle('Mutant effects')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_mutant.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)



results %>%
  filter(!Contrast_type %in% c('Treatment') ) %>%
  VolcanoPlot+
  ggtitle('Other effects')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_other.pdf',sep = ''),
             width=16,height=16,useDingbats=FALSE)


results %>%
  filter(Contrast %in% c('C_Metf','CGlu_Metf','CRP_Metf','C_Glu','dCRP','oeCRP50','oeCRP100') ) %>%
  VolcanoPlot +
  ggtitle('Comparison of treatment, supplementation and mutant effects')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_MainEffects.pdf',sep = ''),
             width=16,height=16,useDingbats=FALSE)


results %>%
  filter(Contrast_type %in% c('Supplementation') ) %>%
  VolcanoPlot +
  ggtitle('Supplementation')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_supplementation.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)


results %>%
  filter(Contrast_type %in% c('Interaction') ) %>%
  VolcanoPlot+
  ggtitle('Interaction')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Interaction.pdf',sep = ''),
             width=12,height=6,useDingbats=FALSE)



sel.mets<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,FDR) %>%
  spread(Contrast,FDR) %>%
  filter(C_Metf<0.05  & CRP_Metf>0.05 & (oeCRP50<0.05 & oeCRP100<0.05) )%>%
  pull(Metabolite)

sel.mets2<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,logFC) %>%
  spread(Contrast,logFC) %>%
  filter(sign(C_Metf)==sign(oeCRP50) & sign(C_Metf)==sign(oeCRP100) ) %>%
  pull(Metabolite)

mets<-intersect(sel.mets2,sel.mets)



results %>%
  filter(Contrast %in% c('C_Metf','CRP_Metf','oeCRP50','oeCRP100') ) %>%
  VolcanoPlot(mets)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Publication_4.pdf',sep = ''),
             width=7,height=5,useDingbats=FALSE)





results %>%
  filter(Contrast %in% c('C_Metf','CRP_Metf','oeCRP50','oeCRP100')) %>%
  ggplot(aes(x=logFC,y=logFDR,color=Strain))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  scale_y_continuous(breaks=seq(0,20,by=2))+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=ifelse(Metabolite %in% mets  & FDR<0.05 & abs(logFC)>1,as.character(Metabolite),'')),size=2)+
  facet_wrap(~Description,ncol = 2)+
  ggtitle('Metformin & CRP')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_Metformin_and_CRP.pdf',sep = ''),
             width=12,height=8,useDingbats=FALSE)







#Venn diagram
results.ecocel %>%
  filter(MetaboliteU=="D-Glucosaminic Acid" & Contrast=="T-C")


venn.pass<-results %>%
  mutate(Up=logFC>0 & FDR<0.05,
         Down=logFC<0 & FDR<0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Up,Down,All) %>%
  filter(Pass) %>%
  group_by(Contrast,Type) %>%
  do(List=c(as.character(.$Metabolite)))


metmatch100<-venn.pass %>%
  filter(Contrast %in% c('C_Metf','oeCRP100') & Type!="All") %>%
  unite(OT,Contrast,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)


metmatch50<-venn.pass %>%
  filter(Contrast %in% c('C_Metf','oeCRP50') & Type!="All") %>%
  unite(OT,Contrast,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)

metmatch50<-metmatch50[c(4,3,2,1)]
metmatch100<-metmatch100[c(4,3,2,1)]

vcols<-c("red","red4","blue","blue4")


grid::grid.draw(VennDiagram::venn.diagram(metmatch50, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_HMT_oeCRP50.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)
dev.off()#



#Heatmap for summary
plotCHeatmap<-function(results,comparisons,mets){
  
  heatsum<-results %>%
    filter(Description %in% comparisons & Metabolite %in% mets) %>%
    select(Description,Metabolite,logFC) %>%
    spread(Description,logFC) %>%
    data.frame(check.names = FALSE,check.rows = FALSE)
  
  
  rownames(heatsum)<-heatsum$Metabolite
  heatsum$Metabolite<-NULL
  
  
  max(heatsum)
  min(heatsum)
  
  amp<-8
  
  minv<- -amp
  maxv<- amp
  
  nstep<-maxv-minv
  nstep<-8
  
  clrbrks<-seq(-amp,amp,by=2)
  clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)
  
  d<-dist(as.matrix(heatsum),method = "euclidean")
  h<-hclust(d,method="ward.D2")
  ordmet<-rownames(heatsum[h$order,])
  
  
  
  if (length(ordmet)!=length(unique(ordmet))){
    print("Non unique metabolites!")
  }
  
  
  results.sum<-results %>%
    filter(Metabolite %in% ordmet & Description %in% comparisons) %>%
    mutate(Metabolite=factor(Metabolite,levels=ordmet,labels=ordmet),
           Description=factor(Description,levels=comparisons,labels=comparisons),
           FDRstars=gtools::stars.pval(FDR))
  
  
  ggplot(results.sum,aes(x=Description,y=Metabolite))+
    geom_tile(aes(fill=logFC))+
    #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
    theme_minimal()+
    geom_text(aes(label=as.character(FDRstars)))+
    #scale_size_discrete(range = c(2,4))+#,breaks=brks
    scale_fill_gradientn(colours = clrscale,
                         breaks=clrbrks,limits=c(-amp,amp))+
    #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
    xlab("Comparison")+
    theme(axis.ticks=element_blank(),
          panel.border=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
}

unique(as.character(results$Description))





comparisons<-c("Treatment effect on OP50","Treatment effect on CRP",
               "Treatment effect on OP50+Glucose","Glucose effect on OP50","Mutant difference for CRP",
               "Mutant difference for oeCRP+IPTG50", "Mutant difference for oeCRP+IPTG100",
               "Interaction between metformin and Glucose","Interaction between metformin and CRP")



comparisons<-c("Treatment effect on OP50","Treatment effect on CRP",
               "Treatment effect on OP50+Glucose","Mutant difference for CRP","Glucose effect on OP50",
               "Mutant difference for oeCRP+IPTG50", "Mutant difference for oeCRP+IPTG100")


comparisons<-c("Treatment effect on OP50","Treatment effect on CRP",
               "Treatment effect on OP50+Glucose",
               "Mutant difference for oeCRP+IPTG50", "Mutant difference for oeCRP+IPTG100")



comparisons<-c("Treatment effect on OP50","Treatment effect on CRP",
               "Mutant difference for oeCRP+IPTG50", "Mutant difference for oeCRP+IPTG100")

#, "Mutant difference for oeCRP+IPTG100"


#"Treatment effect on OP50+Glucose",

#Filter metabolites for clean heatmap
sel.mets<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,FDR) %>%
  group_by(Metabolite) %>%
  filter(any(FDR<0.05) & n()>length(comparisons)*0.5 ) %>%
  ungroup %>%
  spread(Contrast,FDR)



#Chosen


comparisons<-c("Treatment effect on OP50","Treatment effect on CRP","Mutant difference for oeCRP+IPTG50","Mutant difference for oeCRP+IPTG100")

plotCHeatmap(results,comparisons,mets)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP50-100_stat_filter_conservative.pdf',sep = ''),
             width=3,height=5,useDingbats=FALSE)


#Generate table



dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Complete_tidy.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP_All.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP50-100_stat_filter.pdf',sep = ''),
             width=5,height=10,useDingbats=FALSE)




dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP50_stat_filter.pdf',sep = ''),
             width=4,height=8,useDingbats=FALSE)


