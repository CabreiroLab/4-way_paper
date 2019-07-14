#Figure numbering might have been changed
#Data transformation and analysis
library(tidyverse)
library(readxl)

library(ggrepel)
library(RColorBrewer)

# library(ComplexHeatmap)
# library(circlize)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


cwd<-"~/Dropbox/Projects/2015-Metformin/Metabolomics/"
setwd(cwd)


#load('Metabolomics_Ecoli_HMT.RData')
#save.image('Metabolomics_Ecoli_HMT.RData')




#odir<-'Summary_Ecoli_HMT'
odir<-'Summary_Ecoli_HMT/NoGlucose_stats'

dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



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


grporder<-c("C_0_N","C_50_N",
            "C_0_G","C_50_G",
            "CRP_0_N","CRP_50_N",
            "C_0_I50","oeCRP_0_I50","oeCRP_0_I100")

strn.order<-c('OP50-C','crp','oecrp')
spl.order<-c('None','Glucose','IPTG50','IPTG100')

sgrp.order<-c('OP50-C','OP50-C + Glucose','crp','OP50-C + IPTG50','oecrp + IPTG50','oecrp + IPTG100')

met.all<-met.raw %>%
  rename(Metabolite=`Pathway Label`,
         Metabolite_full=`Compound name`,
         KEGG_ID=`KEGG ID`,
         HMDB_ID=`HMDB ID`) %>%
  gather(Sample_ID,Conc,`OP50-1-C`:`CRP-OE-3-100`) %>%
  left_join(met.info) %>%
  left_join(met.TCAp,by=c('Sample_ID','Metabolite')) %>%
  mutate(Strain=recode(Strain,"OP50"="OP50-C","CRP"="crp","oeCRP"="oecrp"),
         Conc=as.numeric(Conc),
         Conc=ifelse(is.na(Conc) & !is.na(Prediction),Prediction,Conc),
         Conc_log=log2(Conc),
         Metformin_mM=as.factor(Metformin_mM),
         SGroup=ifelse(Supplement!="None",paste(Strain,Supplement,sep=" + "),Strain),
         SGroup=factor(SGroup,levels=sgrp.order,labels=sgrp.order),
         Group=factor(Group,levels=grporder,labels=grporder),
         Strain=factor(Strain,levels=strn.order,labels=strn.order),
         Supplement=factor(Supplement,levels=spl.order,labels=spl.order) ) %>%
  select(Sample_ID,Replicate:Group,SGroup,Strain:Amount_ODmL,Metabolite_ID:HMDB_ID,everything(),-Prediction)
  

write.csv(met.all,paste(odir,'/Raw_data_filled.csv',sep=''),row.names = FALSE)


met.clean<-met.all %>%
  tbl_df %>%
  filter(ID!='C_0_N_4') %>% # One sample removed
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

write.csv(missing.obs,paste(odir,'/Missing_observations.csv',sep=''),row.names = FALSE)


#Remove missing metabolites
rm.mets<-missing.obs %>%
  filter(Missing_total>25)

met.sel<-met.clean %>%
  filter(!Metabolite %in% rm.mets$Metabolite)


bioinfo <- met.sel %>%
  group_by(ID,Group,Strain,Supplement,SGroup,Metformin_mM) %>%
  summarise %>%
  data.frame

rownames(bioinfo)<-bioinfo$ID


#Clean data
metcomplete<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc))) %>%
  ungroup %>%
  filter(!Filling_missing)


#PEP
Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- c("0","50")
Metlab<-'Metformin,\nmM'


Pepyr<-metcomplete %>%
  filter(Metabolite %in% c('PEP','Pyruvic acid') & !is.na(Conc_log) & SGroup %in% c("OP50-C","OP50-C + Glucose","crp")) %>%
  mutate(Energy=ifelse(Metabolite=='PEP','PEP','Pyr'),
         Energy=factor(Energy,levels=c('Pyr','PEP'),
                       labels=c('Pyr','PEP')),
         SGroup=factor(SGroup,levels=c("OP50-C","OP50-C + Glucose","crp"),labels=c("OP50","OP50 + Glu","crp") ) )


pepyr.stats<-Pepyr %>%
  group_by(SGroup) %>%
  lmtest('Conc_log~Energy*Metformin_mM') %>%
  rename(Metabolite=Energy,Sample=SGroup)


pepyr.stats %>%
  filter(Contrast_type!='Metformin_mM') %>%
  select(-Metabolite) %>%
  write_csv(paste0(odir,"/PEP|Pyr_stats.csv"))

hj<-0.8
vj<-2
nx<--0.5

pepyr.stats %>%
  filter(Contrast_type=='Energy') %>%
  ggplot(aes(x=Sample,y=logFC,fill=Metformin_mM))+
  geom_hline(yintercept = 0,color="gray40")+
  geom_errorbar(aes(ymin=logFC-SE,ymax=logFC+SE),position = position_dodge(0.75),width = 0.2)+
  geom_col(position = "dodge",width=0.75)+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  scale_fill_manual(name = Metlab,values =Metcols)+
  scale_color_manual(name = Metlab,values =Metcols)+
  geom_text(data=filter(pepyr.stats,Contrast_type=="Energy"),aes(label=pStars,y=3.5,x=Sample,color=Metformin_mM),size=5,show.legend = FALSE, position = position_dodge(0.75))+
  geom_text(data=filter(pepyr.stats,Contrast_type=="Interaction"),aes(label=pStars,y=4,x=Sample),show.legend = FALSE,size=5,color="green4",position = position_dodge(0.75))+
  xlab("Conditions")+
  ylab('PEP/Pyr, logFC')


ggsave(file=paste0(odir,"/PEP|Pyr_ratio.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


Pepyroe<-metcomplete %>%
  filter(Metabolite %in% c('PEP','Pyruvic acid') & !is.na(Conc_log) & SGroup %in% c("OP50-C + IPTG50","oecrp + IPTG50","oecrp + IPTG100")) %>%
  mutate(Energy=ifelse(Metabolite=='PEP','PEP','Pyr'),
         Energy=factor(Energy,levels=c('Pyr','PEP')),
         Overexpression=case_when(Group=='C_0_I50' ~ "None",
                                  Group=="oeCRP_0_I50" ~ "50",
                                  Group=="oeCRP_0_I100" ~ "100"),
         Overexpression=factor(Overexpression, levels=c("None",'50',"100")))


pepyroe.stats<-Pepyroe %>%
  ungroup %>%
  lmtests('Conc_log~Energy*Overexpression') %>%
  rename(Metabolite=Energy)


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
    geom_point(aes(fill=factor( ifelse(Metformin_mM==0,SGroup, NA ) ) ),size=2,stroke=1,shape=21)+
    scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=3))+
    scale_fill_discrete(na.value=NA, guide="none")+
    guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
    #geom_text_repel(aes(label=ID),size=2,color='black')+
    labs(colour='Strain & Supplement')+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


#Original - HMT linear scale filled
metcomplete %>%
  PCAprep('ID','Metabolite','Filled_conc',bioinfo) %>%
  PCAplot


ggsave(file=paste0(odir,"/PCA_All_HMT_linear_scale.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


#Vanilla filled - All
metcomplete %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot

ggsave(file=paste0(odir,"/PCA_All.pdf"),
             width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")



#No IPTG, oe filled
metcomplete %>%
  filter(!Group %in% c("C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot


ggsave(file=paste(odir,"/PCA_noIPTG_oe.pdf",sep=''),
             width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")


#Control, IPTG, oe filled
metcomplete %>%
  filter(Group %in% c("C_0_N","C_50_N","C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot

ggsave(file=paste(odir,"/PCA_IPTG_oe.pdf",sep=''),
             width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")



#No Glucose
metcomplete %>%
  filter(!Supplement=='Glucose') %>%
  PCAprep('ID','Metabolite','Filled_conc_log',bioinfo) %>%
  PCAplot+
  labs(color="Condition")


ggsave(file=paste(odir,"/PCA_NoGlucose.pdf",sep=''),
             width=70,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


heatcomplete<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Completeness=sum(!is.na(Conc_log))*100/n() )%>%
  ungroup %>%
  #At least 50% observations per metabolite
  filter(Completeness>50)


cols<-list(Metformin_mM = c('0' = "white", '50' = "Black"),
           Strain=c("OP50-C"="red3","crp"="blue","oecrp"="purple"),
           Supplement=c("None"="white","Glucose"="green2","IPTG50"="turquoise","IPTG100"="turquoise3") )
           
hinfo<-bioinfo %>%
  select(Strain, Metformin_mM,Supplement)

#Complete
heatcomplete %>%
  HMap('ID','Metabolite','Conc_log',hinfo,cols=cols)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Heatmap_All.pdf",sep=''),
             width=8,height=6)


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
 


met.clean %>%
  #filter(Metabolite %in% c('Malonyl-CoA','Fumaric acid','Malic acid')) %>%
  ggplot(aes(x=SGroup,y=Conc_log,color=Metformin_mM))+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
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
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
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
strainlist<-c('OP50-C','crp','oecrp')

contrasts.desc<-contrasts$Contrasts.table %>%
  mutate(Strain=factor(Strain,levels=strainlist,labels=strainlist)) %>%
  select(-c(Target,Reference))
  

contr.matrix<-contrasts$Contrasts.matrix

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

results.summary %>%
  write.csv(paste(odir,'/Missing_comparisons.csv',sep=''),row.names = FALSE)


dif_summary<-results %>%
  group_by(Contrast,Description) %>%
  summarise(Total=n(),
            Up=sum(logFC>0 & FDR<=0.05,na.rm=TRUE),
            Down=sum(logFC<0 & FDR<=0.05,na.rm=TRUE),
            All=sum(FDR<=0.05,na.rm=TRUE))



dif_summary %>%
  write.csv(paste(odir,'/Stats_summary.csv',sep=''),row.names = FALSE)


#Consistency checks
met.check<-met.clean %>%
  group_by(Metabolite,Group) %>%
  summarise(Mean=mean(Conc_log,na.rm=TRUE),SD=sd(Conc_log,na.rm=TRUE))


write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)



results.multi2<-multiplex(results,c("Metabolite","Metabolite_ID","Metabolite_full","KEGG_ID","HMDB_ID"),dim=2)



#Plotting starts

erralpha<-0.6
errcolor<-'grey80'
lblsize<-2
amp<- 8
cbrks<-seq(-amp,amp,by=2)
#gradcols<-c('black','purple','purple')
maincomp<-'Difference,\nlogFC'


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


results.multi %>%
  filter(x_Contrast=='C_Metf' &
           y_Contrast %in% c('CRP_Metf','CGlu_Metf','oeCRP50','oeCRP100') &
           z_Contrast_type=='Interaction',
           z_Strain==y_Strain &
           z_Supplement==y_Supplement) %>%
  ggplot(aes(x=x_logFC,y=y_logFC,color=z_logFC))+
  geom_smooth(method="lm",fullrange=TRUE,se=FALSE, size=0.5,color='red')+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=x_NE,xmax=x_PE),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=y_NE,ymax=y_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=ifelse(z_FDR<0.01 & abs(z_logFC)>1 | Metabolite %in% mets, as.character(Metabolite),"" ) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+

  xlab('Metformin effect on OP50')+
  ylab('Other effects')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~y_Description,ncol=2)+
  theme(panel.grid.minor = element_blank())

ggsave(file=paste(odir,'/Scatter_Selection_vs_OP50_Meft.pdf',sep = ''),
       width=110,height=110,units='mm',scale=2,device=cairo_pdf,family="Arial")


scatter.fit<-results.multi %>%
  filter(x_Contrast=='C_Metf' &
           y_Contrast %in% c('CRP_Metf','CGlu_Metf','oeCRP50','oeCRP100') &
           z_Contrast_type=='Interaction',
         z_Strain==y_Strain &
           z_Supplement==y_Supplement) %>%
  group_by(y_Contrast) %>%
  do(broom::tidy(lm(y_logFC~x_logFC,data=.)))


scatter.r<-results.multi %>%
  filter(x_Contrast=='C_Metf' &
           y_Contrast %in% c('CRP_Metf','CGlu_Metf','oeCRP50','oeCRP100') &
           z_Contrast_type=='Interaction',
         z_Strain==y_Strain &
           z_Supplement==y_Supplement) %>%
  group_by(y_Contrast) %>%
  do(broom::glance(lm(y_logFC~x_logFC,data=.)))



#Volcano plots

amp<-6

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-10

clrbrks<-seq(-amp,amp,by=2)
#patclrscale <- colorRampPalette(c("purple", "gray50","green"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


VolcanoPlot<-function(data,selmets=c()){
  plot<-ggplot(data,aes(x=logFC,y=logFDR,colour=logFC))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
    geom_point()+

    scale_colour_gradientn(colours = clrscale,
                         breaks=clrbrks,limits=c(-amp,amp))+
    
    scale_y_continuous(breaks=seq(0,20,by=2))+
    scale_x_continuous(breaks=seq(-10,10,by=2))+
    facet_wrap(~Description,ncol = 2)
  if (length(selmets)==0){
    plot<-plot+geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)
  } else {
    plot<-plot+geom_text_repel(aes(label=ifelse(Metabolite %in% selmets & FDR<0.05,as.character(Metabolite),'')),size=5)
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



#Find metabolites which change in metformin treatment and oeCRP
sel.comp<-c("Treatment effect on OP50-C","Treatment effect on crp","Treatment effect on OP50-C + Glucose",
            "Mutant difference for oecrp+IPTG50", "Mutant difference for oecrp+IPTG100")


qthrs<-0.05
fcthrs<-1

sel.mets<-results %>%
  filter(Description %in% sel.comp) %>%
  select(Contrast,Metabolite,logFC,FDR) %>%
  gather(Stat,Value,logFC,FDR) %>%
  unite(CS,Contrast,Stat) %>%
  spread(CS,Value) %>%
  filter( (C_Metf_FDR<qthrs & abs(C_Metf_logFC)>fcthrs )  &
            (oeCRP50_FDR<qthrs & abs(oeCRP50_logFC)>fcthrs )  &
            (oeCRP100_FDR<qthrs & abs(oeCRP100_logFC)>fcthrs ) &
            !(CRP_Metf_FDR<qthrs & abs(CRP_Metf_logFC)>fcthrs ) &
            !(CGlu_Metf_FDR<qthrs & abs(CGlu_Metf_logFC)>fcthrs ) &
            sign(C_Metf_logFC)==sign(oeCRP50_logFC) & sign(C_Metf_logFC)==sign(oeCRP100_logFC)) %>% 
  pull(Metabolite)




mets<-sel.mets


results %>%
  group_by(Contrast) %>%
  arrange(FDR) %>%
  mutate(Show=row_number()==1) %>%
  filter(Contrast %in% c('C_Metf','CRP_Metf','C_Glu','CGlu_Metf','oeCRP50','oeCRP100') ) %>% # 
  VolcanoPlot(mets)+
  geom_text_repel(aes(label=ifelse(Show,as.character(Metabolite),'')),size=5)


#Glucose
ggsave(file=paste(odir,'/Volcano_Publication_4_logFCcolor_Glucose.pdf',sep = ''),
       width=110,height=110,units='mm',scale=2,device=cairo_pdf,family="Arial")


#mets,
results %>%
  filter(Contrast %in% c('C_Metf','CRP_Metf','C_Glu','CGlu_Metf','oeCRP50','oeCRP100') & Metabolite %in% c('Ornithine')) %>%
  View


#Venn diagram

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


#CRP match
metmatchcrp<-venn.pass %>%
  filter(Contrast %in% c('C_Metf','CRP_Metf') & Type!="All") %>%
  unite(OT,Contrast,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)

metmatchcrp<-metmatchcrp[c(4,3,2,1)]

vcols<-c("red","red4","blue","blue4")
dev.off()#
grid::grid.draw(VennDiagram::venn.diagram(metmatchcrp, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_HMT_crp.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)
dev.off()#



metmatchAll<-results %>%
  mutate(Up=logFC>0 & FDR<0.05,
         Down=logFC<0 & FDR<0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Up,Down,All) %>%
  filter(Pass) %>%
  filter(Contrast %in% c('C_Metf','oeCRP50','oeCRP100') & Type!='All') %>%
  group_by(Strain,Type) %>%
  do(List=c(as.character(.$Metabolite))) %>%
  unite(OT,Strain,Type) %>%
  select(OT,List) %>%
  spread(OT,List) %>%
  as.list %>%
  unlist(recursive=FALSE)


vcols<-c("red","red4","blue","blue4")
grid::grid.draw(VennDiagram::venn.diagram(metmatch50, col=vcols,cat.col=vcols ,NULL))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Venn_HMT_oeCRPAll.pdf",sep=''),
             width=3,height=3, useDingbats=FALSE)

dev.off()#



sel.comp<-c("Treatment effect on OP50-C","Treatment effect on crp",'Treatment effect on OP50-C + Glucose',"Mutant difference for oecrp+IPTG50","Mutant difference for oecrp+IPTG100")


#Filter metabolites for clean heatmap
sel.mets<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,FDR) %>%
  group_by(Metabolite) %>%
  filter(any(FDR<0.05) & n()>length(comparisons)*0.5 ) %>%
  ungroup %>%
  spread(Contrast,FDR)

#Chosen


amp<-6

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-6

clrbrks<-seq(-amp,amp,by=2)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


results %>%
  filter(Description %in% sel.comp & Metabolite %in% mets) %>%
  mutate(Description=factor(Description,levels=sel.comp,labels=sel.comp)) %>%
  clustorder("Metabolite","Description","logFC",descending = TRUE) %>%
  ggplot(aes(x=Description,y=Metabolite))+
  geom_tile(aes(fill=logFC))+
  geom_text(aes(label=as.character(FDRStars)))+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  xlab("Comparison")+
  theme_Heatmap()+
  theme(axis.text.x = element_text(angle=45))


ggsave(file=paste0(odir,'/Comparison_Heatmap_Treatments_and_Glucose_oeCRP50-100_stat_filter_conservative2.pdf'),
       width=55,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")




#Generate table

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Complete_tidy.pdf',sep = ''),
             width=100,height=80,units='mm',scale=2,device=cairo_pdf,family="Arial")


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP_All.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP50-100_stat_filter.pdf',sep = ''),
             width=5,height=10,useDingbats=FALSE)




dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP50_stat_filter.pdf',sep = ''),
             width=4,height=8,useDingbats=FALSE)




#For pathways

unique(metpath$PID)


patmets<-c("G1P","G6P","F6P","F1,6P","DHAP","Glyceraldehyde 3-phosphate","3-PG","2-PG","PEP","Pyruvic acid",
           "AcCoA","Malonyl-CoA","Citric acid","cis-Aconitic acid","2-OG","Succinic acid","Fumaric acid","Malic acid",
           "Asp","Asn","Arg","ArgSuccinate","Glu","Gln",
           "NAD+","NADP+","NADH","NADPH")



amp<-2

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-10

clrbrks<-seq(-amp,amp,by=1)
#patclrscale <- colorRampPalette(c("purple", "gray50","green"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)


results %>%
  filter(Contrast %in% c("C_Metf","CRP_Metf","oeCRP100") ) %>% #& Metabolite %in% patmets Metabolite=factor(Metabolite,levels=rev(patmets) ),
  mutate(logFC=ifelse(logFC> 2,2,logFC),
         logFC=ifelse(logFC< -2,-2,logFC)) %>%
  ggplot(aes(x=Description,y=Metabolite))+
  geom_tile(aes(fill=logFC))+
  #geom_text(aes(label=as.character(FDRStars)))+
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  xlab("Comparison")+
  theme_Heatmap()+
  theme(axis.text.x = element_text(angle=45))


ggsave(file=paste0(odir,'/Comparison_Heatmap_For_pathway_bluered.pdf'),
       width=55,height=110,units='mm',scale=2,device=cairo_pdf,family="Arial")


ggsave(file=paste0(odir,'/Comparison_Heatmap_For_pathway_bluered_all.pdf'),
       width=55,height=330,units='mm',scale=2,device=cairo_pdf,family="Arial")

#Enrichment

#Pathways to ECID
ecp<-read_csv("~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/EcoCyc_pathways_metabolites_2018-03-30.csv") %>%
  select(-X1) %>%
  rename(PID=`Pathway-ID`,Pathway=`Pathway-Name`,EC_ID=Compounds) %>%
  separate_rows(EC_ID,sep=";")


#ECID to KEGG & HMDB IDs
ecids<-read_csv("~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/EcoCyc_metabolites_KEGG_HMDB_2018-03-30.csv") %>%
  select(-X1) %>%
  rename(HMDB=`|HMDB|`,KEGG=`|LIGAND-CPD|`,Object_ID=Compound_ID,EC_ID=Instances) %>%
  gather(DB,ID,HMDB,KEGG) %>%
  #filter(!is.na(ID)) %>%
  separate_rows(EC_ID,sep=";")

ecids %>%
  filter(EC_ID=="|FRUCTOSE-6P|")


#Metabolites to KEGG & HMDB IDs
metids<-met.all %>%
  group_by(Metabolite_ID,Metabolite,Metabolite_full,KEGG_ID,HMDB_ID) %>%
  summarise %>%
  rename(KEGG=KEGG_ID,HMDB=HMDB_ID) %>%
  gather(DB,ID,KEGG,HMDB) %>%
  separate_rows(ID,sep=",") %>%
  filter(Metabolite_ID!="-")


metfix<-read_csv(paste0(odir,"/Missing_EcoCyc_pathway_links_manual.csv")) %>%
  rename(EC_IDfix=EC_ID)


#Metabolite - Pathway links
mecp<-metids %>%
  left_join(ecids) %>%
  left_join(metfix) %>%
  mutate(EC_ID=ifelse(is.na(EC_ID) & !is.na(EC_IDfix),EC_IDfix,EC_ID),
         EC_IDfix=NULL) %>%
  left_join(ecp) %>%
  filter(!is.na(PID) & !is.na(EC_ID)) %>%
  group_by(Metabolite_ID,Metabolite,PID) %>%
  filter(row_number() == 1) %>%
  select(Metabolite_ID,Metabolite,Metabolite_full,EC_ID,DB,ID,PID,Pathway)


#20 metabolites could not be assigned to pathwyas
length(unique(mecp$Metabolite))
length(unique(metids$Metabolite))

nopath<-setdiff(unique(metids$Metabolite),unique(mecp$Metabolite))
nopath

metids %>%
  filter(Metabolite %in% nopath) %>%
  group_by(Metabolite_ID,Metabolite,Metabolite_full) %>%
  summarise %>% 
  write_csv(paste0(odir,"/Missing_EcoCyc_pathway_links.csv"))


groups<-c("PID","Pathway")

metpath<-results %>%
  right_join(mecp %>% filter(!is.na(PID) & !is.na(EC_ID))) %>%
  filter(!is.na(FDR) & !is.na(PID) ) %>%
  group_by(Contrast,PID) %>%
  mutate(Pathway_size=n()) %>%
  ungroup
  
metpath %>%
  group_by(Contrast,PID,Pathway,Pathway_size) %>%
  summarise %>%
  group_by(Contrast,Pathway_size) %>%
  summarise(Count=n()) %>%
  View




enrtest<-results %>%
  group_by(Contrast,Description,Contrast_type,Strain,Supplement)


grouping<-group_vars(enrtest)


uniques %>%
  mutate(Up=logFC>0 & FDR <0.05,
         Down=logFC<0 & FDR <0.05,
         All=FDR<0.05) %>%
  #gather different thresholds
  gather(Type,Pass,Up,Down,All) %>%
  group_by_(.dots=c(grouping,"KEGG_ID","Pass","Type")) %>%
  summarise %>%
  group_by_(.dots=c(grouping,"Type")) %>%
  summarise(Unique_size=n(),
         Unique_pass=sum(Pass,na.rm = TRUE) ) %>%
  View


EcoCycenrich<-metpath %>%
  filter(Pathway_size>2) %>%
  group_by(Contrast,Description,Contrast_type,Strain,Supplement) %>%
  enrichment(groups,"Metabolite_ID",enrtype="regular")


EcoCycenrich %>%
  write_csv(paste0(odir,"/All_Enrichment_EcoCyc.csv"))



EcoCycenrich %>%
  filter(p.value<0.05) %>%
  View

EcoCycenrich %>%
  filter(Contrast %in% c("C_Metf","CRP_Metf","oeCRP50","oeCRP100") & Type %in% c("Up","Down")) %>%
  group_by(PID) %>%
  filter( any(p.value<0.05)) %>%
  mutate(logp=-log10(p.value) ) %>%
  ungroup %>%
  clustorder('Pathway',c("Contrast","Type"),'logp',descending=TRUE) %>%
  PlotEnrichment("Type","Pathway",fillval = "logpbin")+
  facet_grid(.~Description)+
  ylab("EcoCyc pathway")+
  labs(fill="p")+
  theme(strip.text.x = element_text(angle=90,size=7))

ggsave(file=paste(odir,'/Enrichment_heatmap_EcoCyc_p.pdf',sep = ''),
       width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")



#KEGG info
path2pathde<-limma::getKEGGPathwayNames('eco',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)
allpaths<-gsub('eco','',path2pathde$PathwayID)

write.csv(path2pathde,paste(odir,'/KEGG/KEGG_pathways.csv',sep=''))


#KEGG mapping

all.kegg.mappings<-read_csv(paste0(odir,"/All_KEGG_mappings_Complete.csv")) %>%
  
KEGGmets<-all.kegg.mappings %>%
  filter(all.mapped!="") %>%
  group_by(PathwayID,Description,all.mapped) %>%
  summarise %>%
  rename(KEGG_ID=all.mapped,Pathway=Description)


KEGGgroups<-c("PathwayID","Pathway")

KEGGenrich<-results.ecocel %>%
  separate_rows(KEGG_ID,sep=",") %>%
  left_join(KEGGmets) %>%
  
  filter(!is.na(KEGG_ID) & !is.na(PathwayID) & !is.na(FDR) & !PathwayID %in% c("eco01502")) %>%
  
  group_by(Contrast,Description,Contrast_type,Strain,Supplement) %>%
  enrichment(KEGGgroups,"KEGG_ID")


KEGGenrich %>%
  write_csv(paste0(odir,"/All_Enrichment_KEGG.csv"))


KEGGenrich %>%
  filter(Contrast %in% c("C_Metf","CRP_Metf","oeCRP50","oeCRP100") & Type %in% c("Up","Down")) %>%
  group_by(PathwayID) %>%
  filter( any(p.value<0.05)) %>%
  mutate(logp=-log10(p.value) ) %>%
  ungroup %>%
  clustorder('Pathway',c("Contrast","Type"),'logp',descending=TRUE) %>%
  PlotEnrichment("Type","Pathway",fillval = "logpbin")+
  facet_grid(.~Description)+
  ylab("KEGG pathway")+
  labs(fill="p")+
  theme(strip.text.x = element_text(angle=90,size=7))

ggsave(file=paste(odir,'/Enrichment_heatmap_KEGG_p.pdf',sep = ''),
       width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial")

