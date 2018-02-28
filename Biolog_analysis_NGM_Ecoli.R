library(tidyverse)
library(ggrepel)

library(ComplexHeatmap)
library(circlize)

#library(RColorBrewer)


#library(limma)

#devtools::install_github("PNorvaisas/PFun")
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



pole<-function(x,y) {
  r<-sqrt(x^2+y^2)
  o<- -2*atan(y/(x))/(pi/2)
  return(o) 
}

pole2<-function(x,y) {
  o<- (atan(y/(x)))/(2*pi)
  return(o) 
}







cwd<-"~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/"
setwd(cwd)

#load('Ecoli.RData')
#save.image('Ecoli.RData')


odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



info<-read_csv('../Biolog/Biolog_metabolites_EcoCyc_Unique_PM1-PM5.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))

head(info)


# info<-info %>%
#   mutate(Metabolite=trimws(as.character(Metabolite)),
#          Name=trimws(as.character(Name)),
#          MetaboliteU=ifelse(Metabolite=='Negative Control','NGM',Metabolite),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM3B'),
#                                 paste(Metabolite,'N',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='P-Source',
#                                 paste(Metabolite,'P',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='S-Source',
#                                 paste(Metabolite,'S',sep='_'),MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM1-A1','NGM_C1',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM2A-A1','NGM_C2',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM3B-A1','NGM_N',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM4A-A1','NGM_P',MetaboliteU),
#          MetaboliteU=ifelse(Index=='PM4A-F1','NGM_S',MetaboliteU)) %>%
#   rename(Metabolite_class=Description)
# 
# 
# head(info)



#Get data

data<-read_csv('Data/Summary.csv') %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A')) %>%
  mutate(Sample=paste(Strain,ifelse(Type=='Control','C','T'),sep='_'),
         SampleID=paste(Sample,Replicate,sep='_'),
         Sample=factor(Sample,
                        levels=c("OP50Sens_C","OP50Sens_T"),
                        labels=c("OP50Sens_C","OP50Sens_T"))) %>%
  left_join(info[,c('Index','MetaboliteU')],by='Index') %>%
  rename(Metabolite=Name,G=Int_750nm_log,GR=a_log) %>%
  select(File:Metformin_mM,Replicate,Well,Index,Sample:MetaboliteU,Metabolite:Group,G,GR) %>%
  gather(key=Measure,value=Value,G,GR) %>%
  group_by(Plate,Type,Replicate,Group,Measure) %>%
  mutate(Value_ref=Value[Metabolite=='Negative Control'],
         Value_norm=Value-Value_ref) %>%
  ungroup %>%
  mutate_at(c('SampleID','Sample','Strain','Metformin_mM'),as.factor)


data.nc<-data %>%
  filter(Metabolite!='Negative Control' & Plate !='PM5')

PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))




bioinfo<-data.nc %>%
  group_by(SampleID,Sample,Strain,Metformin_mM) %>%
  summarise %>%
  ungroup %>%
  data.frame

rownames(bioinfo)<-bioinfo$SampleID



pcashape<- data.nc %>%
  filter(Measure=='G') %>%
  select(SampleID,Index,Value_norm) %>%
  spread(Index,Value_norm) %>%
  data.frame


rownames(pcashape)<-pcashape$SampleID
pcashape$SampleID<-NULL


hca_sample<-hclust(dist(pcashape,method="euclidean"),method="ward.D2")
#hca_variable<-hclust(dist(t(log_data),method="manhattan"),method="complete")

plot(hca_sample, labels=samples,
     hang=-1, main="Cluster Dendrogram", xlab="", sub="", cex=1)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hierarchical_Clustering.pdf",sep=''),
             width=9,height=6)

ir.pca <- prcomp(clean_data,
                 center = TRUE,
                 scale. = TRUE)

write.csv(ir.pca[2],paste(odir,"/PCA_loadings.csv",sep=''))

plot(ir.pca,type='l')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_components.pdf",sep=''),
             width=9,height=6)

pcadata<-data.frame(ir.pca$x) %>%
  rownames_to_column('SampleID') %>%
  left_join(bioinfo)





pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)



ellipses<-pcadata %>%
  group_by(Sample,Strain,Metformin_mM) %>%
  do(getellipse(.$PC1,.$PC2,1) ) %>%
  data.frame
  
# 
# pcadata$Metformin_mM<-as.factor(pcadata$Metformin_mM)
# ellipses$Metformin_mM<-as.factor(ellipses$Metformin_mM)

ggplot(pcadata,aes(x=PC1,y=PC2,colour=Strain))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_path(data=ellipses, aes(x=x, y=y,group=interaction(Sample),linetype=Metformin_mM),size=1)+ 
  geom_point(aes(fill=factor( ifelse(Metformin_mM==0,Strain, NA ) ) ),size=5,stroke=1,shape=21)+
  scale_linetype_manual("Metformin, mM",values=c("0"=1,"50"=2))+
  scale_fill_discrete(na.value=NA, guide="none")+
  guides(linetype = guide_legend(override.aes = list(shape=c(21,21),size=1,linetype=c(1,3),colour='black',fill=c(1,NA))))+
  #geom_text(aes(label=Sample))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)



head(data.n)

#Heatmap
heatshape<- data.nc %>%
  filter(Measure=='G') %>%
  select(SampleID,MetaboliteU,Value_norm) %>%
  spread(SampleID,Value_norm) %>%
  data.frame


rownames(heatshape)<-heatshape$MetaboliteU

heatshape$MetaboliteU<-NULL
head(heatshape)



#Order anotation by heatmap colnames
hanot<-bioinfo[colnames(heatshape),] %>%
  select(Metformin_mM)

ha<-HeatmapAnnotation(df=hanot, col = list(Metformin_mM=c('50'='black','0'='white')))

heatshape.sc<-t(scale(t(heatshape)))

rownames(heatshape.sc)<-NULL

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
             file=paste(odir,"/Heatmap.pdf",sep=''),
             width=6,height=9, useDingbats=FALSE)





#Get LM ready
#Based on selected data
#Get negative control wells
head(data.nc)




contrasts<-read.contrasts2('!Contrasts.xlsx')

contrasts.desc<-contrasts$Contrasts.table %>%
  select(Description:Strain)

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix



results.all<-data.nc %>%
  group_by(Measure,Group,Plate,Well,Index,Metabolite,MetaboliteU,EcoCycID,KEGG_ID) %>%
  do(hypothesise2(.,'Value_norm~0+Sample',contr.matrix)) %>%
  getresults(contrasts.desc)




logFDRbreaks<-c(-1,1.3,2,3,14)
logFDRbins<-c('N.S.','p<0.05','p<0.01','p<0.001')



results<-results.all$results
results.cast<-results.all$cast %>% filter(Measure=='G')
results.castfull<-results.all$castfull %>%
  filter(Measure=='G') %>%
  mutate(`T-C_logFDR_bin`=cut(`T-C_logFDR`, breaks=logFDRbreaks,labels=logFDRbins),
         Pole=ifelse(`C_logFC`>0,
                     pole2(`C_logFC`,`T_logFC`),
                     pole2(`C_logFC`,`T_logFC`)+0.5),
         Pole=ifelse(Pole>0.625,Pole-1,Pole),
         Pole360=Pole*360)

results.castfull.c<-results.all$castfull

results.multi<-results.all$multi %>% filter(Measure=='G')


write.csv(results,paste(odir,'/Ecoli_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/Ecoli_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/Ecoli_results_sidebyside_full.csv',sep=''),row.names = FALSE)



results.celr<-read_csv('Celegans/Summary/Celegans_results_withNGM.csv') %>%
  filter(! Index %in% c('Controls-A1','Controls-B1')) %>%
  select(Index,logFC:logFDR_bin)

celvars<-base::setdiff(colnames(results.celr),'Index')

results.cels<-results.celr %>%
  rename_(.dots=setNames(celvars,paste0('Cel_',celvars)) )


#Main data
results.castcomb<-results.castfull %>%
  left_join(results.cels) %>%
  arrange(desc(`T-C_logFC`)) %>%
  mutate(MetaboliteU=factor(MetaboliteU,labels=MetaboliteU))


#With growth rate
results.castcomb.c<-results.castfull.c %>%
  left_join(results.cels)



#View(results.castcomb)
write.csv(results.castcomb,paste(odir,'/Combined_results_sidebyside_full.csv',sep=''),row.names = FALSE)





#Pick data to show

#Only filtered
#Carboxylic acids to remove
coxy<-c('Itaconic Acid','Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')
selmets<-info %>%
  filter(!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy) ) %>%
  select(Group,Plate,Well,Metabolite,MetaboliteU)


selectcast<-subset(results.castcomb,MetaboliteU %in% selmets$MetaboliteU )
selectcast.c<-subset(results.castcomb.c,MetaboliteU %in% selmets$MetaboliteU )


dim(selectcast)

metorder<-as.character(results.castcomb$MetaboliteU)
data.nc<-data.nc %>%
  mutate(MetaboliteU=factor(MetaboliteU,levels=metorder))

data.nc %>%
  filter(Measure=='G') %>%
  ggplot(aes(x=MetaboliteU,y=Value_norm,color=Type))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity") +
  geom_point()+
  ylab('Growth logFC vs NGM')+
  xlab('Metabolite')+
  coord_flip()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ecoli_logFC_tall_raw.pdf",sep=''),
             width=9,height=50)



#Volcano plots

ggplot(selectcast,aes(x=`T-C_logFC`,y=`T-C_logFDR`,color=-Cel_logFC))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=`T-C_NE`,xmax=`T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(`Cel_logFC`) ))+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  ylab('-logFDR')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  labs(color='Significance (FDR)')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_Ecoli_Growth.pdf",sep=''),
             width=12,height=9)



# ggplot(selectcast.c,aes(x=`GR_T-C_logFC`,y=`GR_T-C_logFDR`,color=-Cel_logFC))+
#   geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
#   geom_errorbarh(aes(xmin=`GR_T-C_NE`,xmax=`GR_T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
#   geom_point(aes(size=abs(`Cel_logFC`) ))+
#   scale_colour_gradientn(colours = gradcols,
#                          breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
#   labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
#   scale_size(range = c(0.25, 7),name=maincomp)+
#   ylab('-logFDR')+
#   xlab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
#   labs(color='Significance (FDR)')
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Volcano_Ecoli_GrowthRate.pdf",sep=''),
#              width=12,height=9)




#Growth rate analysis

fit<-lm(GR_T_logFC~GR_C_logFC,selectcast.c)
lmeq<-lm_eqn(fit)
a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]

thres<-0.05

ggplot(selectcast.c,aes(y=GR_T_logFC,x=GR_C_logFC,color=-Cel_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=GR_T_NE,ymax=GR_T_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=GR_C_NE,xmax=GR_C_PE),alpha=erralpha,color=errcolor,height=0)+
  ggtitle(paste('Scatterplot of metformin and metabolite supplementation effects ',grp,sep=''),
          subtitle = paste('Metabolites with FDR<',thres,' are marked',sep='') )+
  
  geom_point(size=3)+#aes(size=abs(Cel_logFC))
  #coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  scale_x_continuous(breaks=seq(-10,10,by=1))+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  
  xlab('Doubling rate - Control, log2(OD)/h')+
  ylab('Doubling rate - +50mM Metformin, log2(OD)/h')+
  #eval(parse(text = intfdr)) < thres & abs( eval(parse(text = intvar)) ) > 0.75 
  # geom_text_repel(aes(label=ifelse(Cel_FDR<thres, as.character(Metabolite),'')),
  #                 size=lblsize,nudge_y = 0.3,
  #                 force=1,
  #                 segment.colour=errcolor,
  #                 segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  #scale_size(range = c(0.25, 7),name=maincomp)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_Control-Treatment_Complete_GrowthRate_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)



#Growth rate and growth comparison

head(selectcast.c)

ggplot(selectcast.c,aes(x=`G_T-C_logFC`,y=`GR_T-C_logFC`,color=-Cel_logFC) )+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  #geom_abline(aes(intercept=b,slope=a),alpha=1,color='red')+
  #geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=`GR_T-C_NE`,ymax=`GR_T-C_PE`),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=`G_T-C_NE`,xmax=`G_T-C_PE`),alpha=erralpha,color=errcolor,height=0)+
  geom_point(aes(size=abs(Cel_logFC)))+
  geom_text_repel(aes(label=ifelse(Cel_FDR<thres, as.character(Metabolite),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ylab('Doubling rate difference (Treatment-Control), log2(OD)/h')+
  xlab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  scale_size(range = c(0.25, 7),name=maincomp)+
  coord_cartesian(xlim=c(-3,3),ylim = c(-2,2))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_GrowthRate_vs_Growth_vs_Celegans.pdf",sep=''),
             width=13,height=9, useDingbats=FALSE)


#Driving effects for interaction
ggplot(selectcast,aes(y=`T-C_logFC`,x=Pole360,color=-Cel_logFC))+
  geom_point()+
  geom_point(aes(size=abs(Cel_logFC)))+
  geom_text_repel(aes(label=ifelse(Cel_FDR<0.05, as.character(Metabolite),'')),
                  size=lblsize,nudge_y = 0.3,
                  force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)+
  ggtitle('Breakdown of driving components for metabolite-metformin interaction')+
  #scale_x_discrete(labels=c("-90"="-T","-45"="-TC","0" = "C", "45" = "TC","90" = "T","135"="T-C","180"="-C","225"="-T-C","270"="-T"))+
  scale_x_continuous(breaks=seq(-135,225,by=45),labels=c("-T-C","-T","-TC","C","TC","T","T-C","-C","-T-C")  )+
  scale_y_continuous(breaks=seq(-4,4,by=1),limits = c(-3,3) )+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp))+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  xlab('Driving component')+
  ylab('Metabolite - metformin interaction as growth logFC vs NGM')+
  scale_size(range = c(0.25, 7),name=maincomp)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Celegans_Ecoli_Driving_effects_for_interaction_CelSignificant.pdf",sep=''),
             width=16,height=12, useDingbats=FALSE)





#Enrichment
#Generate tables for KEGG enrichment


dir.create(paste(odir,"/Biolog_enrichment/",sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

enrt<-subset(data.nc,Group %in% c('C-Source','N-Source','P-Source','S-Source') & 
               !Name %in% c('Negative Control',"D,L-Malic Acid","beta-Methyl-D-Galactoside","O-Phospho-D-Tyrosine") &
               !KEGG_ID %in% c("") &
               !( Plate %in% c('PM3B','PM4A') & Name %in% PM1PM2$Name ) &
               !( Plate %in% c('PM3B','PM4A') & KEGG_ID %in% PM1PM2$KEGG_ID ) & 
               !( Plate=='PM4A' & KEGG_ID=='C00097' ) )[,c('SampleID','Sample','Plate','Index','Group','Name','KEGG_ID','Measure','Value_norm')]
head(enrt)

enrc<-dcast(enrt,SampleID+Sample+Measure~KEGG_ID,value.var = 'Value_norm')

head(enrc)

KEGGS<-as.character(unique(enrc$KEGG_ID))

#Find duplicates
#Works when aggregate function is length
dKEGGS<-KEGGS[apply(enrc[,KEGGS],2, function(x) any(x>1))]
dKEGGS
#subset(enrs,KEGG_ID %in% dKEGGS[2])
length(unique(enrt$KEGG_ID))


KEGGbac<-data.frame(unique(enrt$KEGG_ID))

write.csv(KEGGbac,paste(odir,"/Biolog_enrichment/Background.csv",sep=''),row.names = FALSE,col.names = NA,quote=FALSE)


for (gr in c('C-Source','N-Source','P-Source','S-Source','Complete')){
  print(gr)
  if (gr=='Complete') {
    enrs<-subset(enrt,Measure=='G')
  } else {
    enrs<-subset(enrt,Group==gr & Measure=='G')
  }
  enrc<-dcast(enrs,SampleID+Sample~KEGG_ID,value.var = 'Value_norm')
  write.csv(enrc,paste(odir,"/Biolog_enrichment/Biolog_Enrichment_",gr,"_T-C_",strn,".csv",sep=''),row.names = FALSE)
}





#KEGG info
library(limma)

library(org.EcK12.eg.db)
library(pathview)

#allpaths<-c('01110')




path2pathde<-getKEGGPathwayNames('eco',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)
allpaths<-gsub('eco','',path2pathde$PathwayID)

write.csv(path2pathde,paste(odir,'/KEGG/KEGG_pathways.csv',sep=''))




#Read KEGG enrichment results
allqea<-read.csv(paste(odir,'/Biolog_enrichment/pathway_results.csv',sep=''),sep=',')

allqea<-rename(allqea,c('X'='KEGG pathway'))
allqea$logFDR<- -log10(as.numeric(allqea$FDR))
allqea$Impact<-as.numeric(allqea$Impact)
allqea$Hits<-as.numeric(allqea$Hits)
allqea$`Total.Cmpd`<-as.numeric(allqea$`Total.Cmpd`)
allqea$Ratio<-allqea$Hits/allqea$`Total.Cmpd`

allqea.a<-merge(path2pathde,allqea,by.x='Description',by.y='KEGG pathway',all.y=TRUE)
subset(allqea.a,is.na(PathwayID))

head(allqea.a)

thres<-0.05
ggplot(allqea.a,aes(x=Impact,y=logFDR,size=Ratio,color=Hits/`Total.Cmpd`))+
  geom_point()+
  ylab('-log(FDR)')+
  xlab('Pathway impact')+
  geom_text_repel(aes(label=ifelse(FDR < thres & Impact > 0.25 , as.character(`Description`),'') ),
                  size=2,nudge_y = 0.3,force=1,
                  segment.colour=errcolor,
                  segment.alpha =segalpha)



#KEGG mapping

#Only unique values with preference for PM1 PM2A
all.results.rcp<-subset(results.castfull,
                        !Metabolite %in% c('Negative Control',"D,L-Malic Acid","beta-Methyl-D-Galactoside","O-Phospho-D-Tyrosine") &
                          !KEGG_ID %in% c("") &
                          !is.na(KEGG_ID) &
                          !( Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Name ) &
                          !( Plate %in% c('PM3B','PM4A') & KEGG_ID %in% PM1PM2$KEGG_ID ) & 
                          !( Plate=='PM4A' & KEGG_ID=='C00097'))

keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/KEGG_pathways'

pathcomp<-c('OP50Sens_T-OP50Sens_C_logFC')


mapkey<-'KEGG_ID'
nrow(all.results.rcp)
gdata<-all.results.rcp
#gdata<-all.results.rcp
nrow(gdata)

dupl<-duplicated(gdata[,mapkey])
print('KEGG ID duplicates')
print(table(dupl))
gdataf<-gdata[!dupl,]
rownames(gdataf)<-gdataf[,mapkey]

print(paste('Total KEGG pathways to plot:',length(allpaths)))

keggdir<-paste(cwd,odir,'/KEGG',sep='')
dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd(keggdir)
#comp.gr<-'Main'

limitcpd<-3

# pathcomp<-comp.list[[comp.gr]]
# outsuffx<-comp.gr
pv.out <- pathview(cpd.data = gdataf[,pathcomp,drop=FALSE],
                   pathway.id = allpaths,
                   cpd.idtype = "kegg",
                   species = "eco",
                   out.suffix = 'Complete',
                   kegg.dir = keggxml,
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   min.nnodes = 0,
                   limit=list(gene=2,cpd=limitcpd),node.sum = "mean",
                   low = list(gene = "blue", cpd = "blue"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="red"))
setwd(cwd)

all.kegg.mappingst<-data.frame()
paths<-names(pv.out)
for (pth in paths) {
  result<-pv.out[[pth]][['plot.data.cpd']]
  result$PathwayID<-pth
  all.kegg.mappingst<-rbind(all.kegg.mappingst,result)
}
all.kegg.mappingst$mol.col<-NULL

all.kegg.mappings<-merge(path2pathde,all.kegg.mappingst,by='PathwayID',all.y=TRUE)


write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings_Complete.csv',sep = ''),row.names = FALSE)



