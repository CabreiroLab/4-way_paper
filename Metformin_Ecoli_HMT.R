#Data transformation and analysis
library(xlsx)
library(tidyverse)


# library(plyr)
# library(reshape2)
#library(multcomp)
#library(contrast)
#library(car)

#Plotting and visualisation
#library(gplots)
library(heatmap3)
library(ggplot2)
#library(ggbiplot)
library(ggrepel)
library(RColorBrewer)
#library(plot3D)


library(gtools)

library(ellipse)

#library(Unicode)



devtools::install_github("PNorvaisas/PFun")
library(PFun)


library(ggthemes)
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

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")
load('.RData')




odir<-'Summary_Ecoli_HMT'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}

getellipse<-function(x,y,sc=1) {
  as.data.frame(ellipse::ellipse( cor(x, y),
                         scale=c(sd(x)*sc,sd(y)*sc),
                         centre=c( mean(x),mean(y)) ))
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


met.raw<-read.xlsx2('HMT_Ecoli_complete/Raw_data.xlsx',sheetName='Raw',
                    header=TRUE,check.names=FALSE,
                    colClasses=c(rep("character",5),rep("numeric",34)) )

met.info<-read.xlsx2('HMT_Ecoli_complete/Raw_data.xlsx',sheetName='Samples',
                     header=TRUE,check.names=FALSE)


grporder<-c("C_0_N","C_50_N",
            "C_0_G","C_50_G",
            "CRP_0_N","CRP_50_N",
            "C_0_I50","oeCRP_0_I50","oeCRP_0_I100")

strn.order<-c('OP50','CRP','oeCRP')
spl.order<-c('None','Glucose','IPTG50','IPTG100')

sgrp.order<-c('OP50 None','OP50 Glucose','CRP None','OP50 IPTG50','oeCRP IPTG50','oeCRP IPTG100')

met.all<-met.raw %>%
  gather(Sample_ID,Conc,`OP50-1-C`:`CRP-OE-3-100`) %>%
  left_join(met.info) %>%
  mutate(Conc_log=log2(Conc),
         SGroup=paste(Strain,Supplement),
         Group=factor(Group,levels=grporder,labels=grporder),
         Strain=factor(Strain,levels=strn.order,labels=strn.order),
         Supplement=factor(Supplement,levels=spl.order,labels=spl.order),
         SGroup=factor(SGroup,levels=sgrp.order,labels=sgrp.order)) %>%
  rename(Metabolite=`Pathway Label`,
         Metabolite_full=`Compound name`,
         KEGG_ID=`KEGG ID`,
         HMDB_ID=`HMDB ID`) %>%
  select(Sample_ID,Replicate:Amount_ODmL,Metabolite_ID:HMDB_ID,everything())
  
head(met.all)


write.csv(met.all,paste(odir,'/Raw_data.csv',sep=''),row.names = FALSE)


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

met.sel<-met.clean %>%
  filter(!Metabolite %in% rm.mets$Metabolite)


glimpse(met.sel)

#Original - Full of holes... HMT linear scale
metslm<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc) %>%
  spread(Metabolite,Filled_conc)


#Vanilla filled - All
metslm<-met.sel %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc_log))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc_log) %>%
  spread(Metabolite,Filled_conc_log)



#No IPTG, oe filled
metslm<-met.sel %>%
  filter(!Group %in% c("C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc_log))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc_log) %>%
  spread(Metabolite,Filled_conc_log)


#Control, IPTG, oe filled
metslm<-met.sel %>%
  filter(Group %in% c("C_0_N","C_50_N","C_0_I50","oeCRP_0_I50","oeCRP_0_I100")) %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc_log))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc_log) %>%
  spread(Metabolite,Filled_conc_log)



#No Glucose
metslm<-met.sel %>%
  filter(!Supplement=='Glucose') %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc_log))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc_log) %>%
  spread(Metabolite,Filled_conc_log)


#What Filipe asked
metslm<-met.sel %>%
  filter(!Supplement=='Glucose' & Strain %in% c('OP50','CRP')) %>%
  group_by(Metabolite) %>%
  mutate(Filling_missing=any(is.na(Filled_conc_log))) %>%
  ungroup %>%
  filter(!Filling_missing) %>%
  select(Group,SGroup,ID,Strain,Metformin_mM,Supplement,Metabolite,Filled_conc_log) %>%
  spread(Metabolite,Filled_conc_log)


dim(metslm)
View(metslm)




#Find compounds with missing values


#All missing valus for compound
allmiss.cols<-apply(metslm, 2, function(x) all(is.na(x)))
metslm[,allmiss.cols]


anymiss.cols<-apply(metslm, 2, function(x) any(is.na(x)))
anymiss.rows<-apply(metslm, 1, function(x) any(is.na(x)))


metslm[,anymiss.cols]
metslm[anymiss.rows,]

missing.cols<-names(anymiss.cols[anymiss.cols==TRUE])
missing.rows<-rownames(metslm)[anymiss.rows==TRUE]

missing.cols
missing.rows




pca.dat<-metslm %>%
  select(- (Group:Supplement))



ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
plot(ir.pca,type='l')

pcadata<-data.frame(ir.pca$x)
pcadata[,c('ID','Group','SGroup','Strain','Metformin_mM','Supplement')]<-metslm[,c('ID','Group','SGroup','Strain','Metformin_mM','Supplement')]


pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)


ellipses<-pcadata %>%
  group_by(Strain,Group,SGroup,Metformin_mM,Supplement) %>%
  do(getellipse(.$PC1,.$PC2,1))


ggplot(pcadata,aes(x=PC1,y=PC2,colour=SGroup))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Group),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,SGroup, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  geom_text_repel(aes(label=ID),size=2,color='black')+
  labs(colour='Strain & Supplement')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_NoGlucose.pdf",sep=''),
             width=12,height=9)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_All_HMT_linear_scale.pdf",sep=''),
             width=12,height=9)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_WhatFilipeAsked.pdf",sep=''),
             width=12,height=9)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_All.pdf",sep=''),
             width=12,height=9)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_noIPTG_oe.pdf",sep=''),
             width=12,height=9)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_IPTG_oe.pdf",sep=''),
             width=12,height=9)



#Heatmap
met.heat<-met.sel

#Heatmap - No Glucose
met.heat<-met.sel %>%
  filter(!Supplement=='Glucose')



#Heatmap - no IPTG, oe
met.heat<-met.sel %>%
  filter(!Group %in% c("C_0_I50","oeCRP_0_I50","oeCRP_0_I100"))


#Heatmap - IPTG, oe
met.heat<-met.sel %>%
  filter(Group %in% c("C_0_N","C_50_N","C_0_I50","oeCRP_0_I50","oeCRP_0_I100"))
 


#Heatmap - WhatFilipeAsked
met.heat<-met.sel %>%
 filter(!Supplement=='Glucose' & Strain %in% c('OP50','CRP'))


heatshape<-met.heat %>%
  group_by(Metabolite) %>%
  mutate(Completeness=sum(!is.na(Conc_log))*100/n() )%>%
  ungroup %>%
  #At least 50% observations per metabolite
  filter(Completeness>50) %>%
  select(ID,Metabolite,Conc_log) %>%
  spread(ID,Conc_log) %>%
  data.frame

rownames(heatshape)<-heatshape$Metabolite
heatshape$Metabolite<-NULL

met.groups<-met.heat %>%
  group_by(ID,Strain,Metformin_mM,Supplement) %>%
  summarise %>%
  ungroup %>%
  mutate_at(vars(Strain:Supplement),funs(num=as.numeric(.) )) %>%
  mutate(Strain_col=rainbow(max(Strain_num))[Strain_num],
         Metformin_mM_col=ifelse(Metformin_mM=='0','white','black'),
         Supplement_col=brewer.pal(max(Supplement_num),'Set1')[Supplement_num] ) %>%
  data.frame %>%
  column_to_rownames('ID')


met.groups<-met.groups[colnames(heatshape),]


met.strain<-met.groups %>%
  group_by(Strain,Strain_col) %>%
  summarise

met.suppl<-met.groups %>%
  group_by(Supplement,Supplement_col) %>%
  summarise

met.metf<-met.groups %>%
  group_by(Metformin_mM,Metformin_mM_col) %>%
  summarise


clab<-cbind(met.groups$Metformin_mM_col,met.groups$Supplement_col,met.groups$Strain_col)
colnames(clab)<-c('Metformin, mM','Supplement','Strain')


hmap<-heatmap3(as.matrix(heatshape),
               scale = 'row',
               balanceColor = T,
               ColSideColors = clab,
               #hclustfun = function(x) hclust(dist(x,method="manhattan"),method="complete"),
               #distfun = function(x) dist(x,method="manhattan"),
               #labRow=FALSE,
               #key=TRUE,
               #ColSideLabs = 'Group',
               margins=c(10,10)
)

legend('topright',legend=met.strain$Strain,fill=met.strain$Strain_col, border=FALSE, bty="n",title='Strain')
legend('right',legend=met.suppl$Supplement,fill=met.suppl$Supplement_col, border=FALSE, bty="n",title='Supplement')
legend('bottomright',legend=met.metf$Metformin_mM,fill=met.metf$Metformin_mM_col, border=TRUE, bty="n",title='Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_noIPTG_oe.pdf",sep=''),
             width=14,height=12, useDingbats=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_IPTG_oe.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_All.pdf",sep=''),
             width=17,height=12, useDingbats=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_noGlucose.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Heatmap_WhatFilipeAsked.pdf",sep=''),
             width=15,height=12, useDingbats=FALSE)




ggplot(met.clean,aes(x=SGroup,y=Conc_log,color=Metformin_mM))+
  ggtitle('Comparison between Control and Treatment. Boxplot: +/-SD, Min/Max')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point()+
  scale_y_continuous(breaks=seq(-20,20,by=1))+
  geom_text(aes(label=Replicate),color='black',size=2)+
  ylab('log2 Concentration, pmol/ODmL')+
  xlab('Group & Supplement')+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_wrap(~Metabolite,ncol = 4,scales='free_y')

#,labeller = labeller(Metablite = label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Raw_logConc_by_metabolite.pdf",sep=''),
             width=15,height=80, useDingbats=FALSE)





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

#,labeller = labeller(Metablite = label_both)

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






#Summarise manually

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





#Metabolite summary
met.c<-met.clean %>%
  group_by(Metabolite) %>%
  summarise_at(vars(Conc,Conc_log),funs(SD=sd(.,na.rm = TRUE),Mean=mean(.,na.rm = TRUE)))
  

#Linear modelling


head(met.clean)




complete<-met.clean %>%
  filter(!Supplement %in% c('IPTG50','IPTG100')) %>%
  group_by(Metabolite_ID,Metabolite,Group) %>%
  summarise(Complete=sum(!is.na(Conc_log))>0 ) %>%
  spread(Group,Complete) %>%
  mutate(Complete=all(C_0_N:CRP_50_N))

complete %>% 
  filter(Metabolite %in% c('2-KIV','') )


all.groups<-as.character(unique(met.clean$Group))



contrasts<-read.contrasts('!Contrasts_Ecoli_HMT_metabolomics.xlsx','Contrasts_values',all.groups)


contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix


strainlist<-c('OP50','CRP','oeCRP')

contrasts.table$Strain<-factor(contrasts.table$Strain,levels=strainlist,labels=strainlist) #
contr.matrix


LA<-met.clean%>%
  filter(Metabolite=='Lactic acid') %>%
  group_by(Group) %>%
  summarise(Observations=sum(!is.na(Conc_log))) %>%
  data.frame


#New testing routine
contr.matrix %>%
  data.frame %>%
  rownames_to_column("Contrast") %>%
  gather(Group,Value,C_0_N:oeCRP_0_I100) %>%
  left_join(LA) %>%
  group_by(Contrast) %>%
  






#New linear modelling pipeline

lmdata<-met.clean %>%
  group_by(Metabolite,Metabolite_ID,Metabolite_full,KEGG_ID,HMDB_ID) %>%
  do(hypothesise2(.,'Conc_log~0+Group',contr.matrix)) %>%
  ungroup




results<-lmdata %>%
  group_by(Contrast) %>%
  mutate(FDR=p.adjust(p.value,method = 'fdr'),
         PE=logFC+SE,
         NE=logFC-SE,
         logFDR=-log10(FDR)) %>%
  ungroup %>%
  left_join(contrasts.table[,c('Contrast','Description','Contrast_type','Strain','Supplement')],by='Contrast') %>%
  left_join(met.mets) %>%
  mutate(Contrast=factor(Contrast,levels=contrasts.table$Contrast,labels=contrasts.table$Contrast),
         Description=factor(Description,levels=contrasts.table$Description,labels=contrasts.table$Description)) %>%
  select(Contrast,Description:Supplement,Metabolite_ID,Metabolite:HMDB_ID,everything())


results.summary<-results %>%
  select(Metabolite,Contrast,logFC) %>%
  group_by(Metabolite,Contrast) %>%
  summarise(Estimate=!is.na(logFC)) %>%
  spread(Contrast,Estimate)


results %>%
  filter(Metabolite=='Lactic acid')


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


results %>%
  filter(Metabolite=='2-Hydroxyglutaric acid') %>%
  select(Contrast,logFC,p.value,FDR)

results.old %>%
  filter(Metabolite=='2-Hydroxyglutaric acid') %>%
  select(Contrast,logFC,p.value,FDR)



# allresults<-hypothesise(lmshape,metid,contr.matrix,formula="0+Group")
# results.old<-allresults$All %>%
#   rename(Metabolite_ID=Variable) %>%
#   left_join(contrasts.table[,c('Contrast','Description','Contrast_type','Strain','Supplement')],by='Contrast') %>%
#   left_join(met.mets) %>%
#   select(Contrast,Description:Supplement,Metabolite_ID,Metabolite:HMDB_ID,everything())




results.m<-results %>%
  gather(Stat,Value,logFC:logFDR)
  

results.castfull<-results.m %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Metabolite_ID:HMDB_ID,CS,Value) %>%
  spread(CS,Value)
  


results.cast<-results.m %>%
  filter(Stat %in% c('logFC','FDR')) %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Metabolite_ID:HMDB_ID,CS,Value) %>%
  spread(CS,Value)
  


head(results.castfull)
head(results.cast)





write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)




allcontrs<-unique(as.character(results$Contrast))
allmets<-unique(as.character(results$Metabolite))




stats<-c('logFC','FDR','PE','NE')
contcombs<-apply(expand.grid(contrs, stats), 1, paste, collapse="_")



resultsmin<-results %>%
  select(Metabolite,Contrast,logFC:logFDR)



results.exp<-expand.grid(allmets,contrs, contrs) %>%
  rename(Metabolite=Var1,x=Var2,y=Var3) %>%
  left_join(results,by=c('Metabolite','y'='Contrast')) %>%
  left_join(resultsmin,by=c('Metabolite','x'='Contrast'),suffix = c("", "_x")) %>%
  filter(!is.na(Description) & !is.na(logFC_x))


results.exp %>%
  filter(is.na(Description))




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

results.exp %>%
  filter(x=='C_Metf' & y %in% c('CRP_Metf','CGlu_Metf',
                               'dCRP', 'C_Glu',
                               'oeCRP50','oeCRP100') ) %>%
  ggplot(aes(x=logFC_x,y=logFC,color=logFC-logFC_x))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbarh(aes(xmin=NE_x,xmax=PE_x),alpha=erralpha,color=errcolor,height=0)+
  geom_errorbar(aes(ymin=NE,ymax=PE),alpha=erralpha,color=errcolor,width=0)+
  geom_point()+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  geom_text_repel(aes(label=ifelse(FDR<0.05 | FDR_x<0.05, as.character(Metabolite),"" ) ),
                  size=lblsize,
                  force=2,
                  segment.colour=errcolor,
                  segment.alpha =erralpha)+
  ggtitle('Effects of various factors compared against metformin effect on OP50',
          subtitle = 'Metabolites with FDR<0.05 in either of effects are shown')+
  xlab('Metformin effect on OP50')+
  ylab('Other effects')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  facet_wrap(~Description,ncol=2)+
  theme(panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Scatter_Various_vs_OP50_Meft.pdf',sep = ''),
             width=18,height=18,useDingbats=FALSE)



#Volcano plots
VolcanoPlot<-function(data){
  plot<-ggplot(data,aes(x=logFC,y=logFDR,color=Strain))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
    geom_point()+
    scale_y_continuous(breaks=seq(0,20,by=1))+
    scale_x_continuous(breaks=seq(-10,10,by=1))+
    geom_text_repel(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),size=2)+
    facet_wrap(~Description,ncol = 3)
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








#Heatmap for summary

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



#Filter metabolites for clean heatmap
sel.mets<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,FDR) %>%
  group_by(Metabolite) %>%
  filter(any(FDR<0.05) & n()>length(comparisons)*0.5 ) %>%
  ungroup %>%
  spread(Contrast,FDR)



#Chosen
comparisons<-c("Treatment effect on OP50","Treatment effect on CRP","Mutant difference for oeCRP+IPTG50")

sel.mets<-results %>%
  filter(Description %in% comparisons) %>%
  select(Contrast,Metabolite,FDR) %>%
  spread(Contrast,FDR) %>%
  filter(C_Metf<0.05 & oeCRP50<0.05 & CRP_Metf>0.05 )








#Generate table
heatsum<-results %>%
  filter(Description %in% comparisons & Metabolite %in% sel.mets$Metabolite) %>%
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
h<-hclust(d)



#rownames(heatsum[h$order,])
# brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))
# bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)
#reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
#rownames(heatsum)
# hm<-heatmap3(as.matrix(heatsum),key=TRUE,Colv=FALSE,trace='none',col=bgg,
#              xlab='Comparison',Rowv=TRUE,breaks = brks,dendrogram="row",scale="none")
#ordmet<-rownames(heatsum[hm$rowInd,])

ordmet<-rownames(heatsum[h$order,])


if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}


results.sum<-results %>%
  filter(Metabolite %in% ordmet & Description %in% comparisons) %>%
  mutate(Metabolite=factor(Metabolite,levels=ordmet,labels=ordmet),
         Description=factor(Description,levels=comparisons,labels=comparisons),
         FDRstars=stars.pval(FDR))


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


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Complete_tidy.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)




dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Treatments_and_oeCRP_tidy.pdf',sep = ''),
             width=6,height=16,useDingbats=FALSE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_Complete.pdf',sep = ''),
             width=8,height=16,useDingbats=FALSE)



dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Treatment_OP50_oeCRP_vs_CRP.pdf',sep = ''),
             width=6,height=8,useDingbats=FALSE)
