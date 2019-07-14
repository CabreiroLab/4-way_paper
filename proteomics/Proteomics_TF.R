library(gplots)
library(ggplot2)


library(xlsx)
library(plyr)
library(reshape2)
library(ggbiplot)
library(ggrepel)
library(car)
library(heatmap3)
library(gtools)

library(RColorBrewer)
library(grid)
library(gridExtra)

# library(rgl)
# library(pca3d)

library(multcomp)
library(contrast)
library(scales)


library(pathview)
library(org.EcK12.eg.db)
library(limma)

theme_set(theme_light())
theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Proteomics/")
getreplicates<-function(df,uniqs) {
  df$Replicate<-1
  repeat {
    if(!anyDuplicated(df[,c(uniqs,'Replicate')]))
    { break }
    df$Replicate<-ifelse(duplicated(df[,c(uniqs,'Replicate')]),df$Replicate+1,df$Replicate)
  }
  return(df$Replicate)
}


map<-read.csv('UniProt_gene_name.txt',sep='\t')
colnames(map)<-c('UniProtID','Gene')

prot<-read.xlsx2('Data Filipe - Metformin Full.xlsx',
                 sheetName = 'Combined with previous analysis',endRow = 165,
                 colClasses = c('numeric','character','character',rep('numeric',9),rep('character',3)))
prot$FDR_TvC<-p.adjust(prot$ANOVA_TvC,method='fdr')
prot$FDR_RvS<-p.adjust(prot$ANOVA_RvS,method='fdr')


prota<-merge(prot,map,by.x='Uniprot',by.y='UniProtID',all.x=TRUE)




annot<-read.csv('../Annotations/Ecoli/UniPort_GeneID_20170719.csv',header = TRUE)

# Loop through the dataframe looking for duplicate pairs of
# Runs and Indices and increment the index if it's a duplicate



prota$Replicate<-getreplicates(prota,c('Gene'))

colnames(prota)

colnames(prota)[c(5:9,11:13,16:17)]

meas<-colnames(prota)[c(5:9,11:13,16:17)]


prota.m<-melt(prota,id.vars = c('Uniprot','Gene','Replicate'),measure.vars = meas,variable.name = 'Comparison',value.name = 'Value')
prota.m$Uniprot<-as.character(prota.m$Uniprot)
prota.m$Uniprot<-ifelse(prota.m$Uniprot=='B1XGK9','P77581',prota.m$Uniprot)



prota.m$Stat<-ifelse(grepl('ANOVA',prota.m$Comparison),'ANOVA','logFC')
prota.m$Comparison<-gsub('ANOVA_','',prota.m$Comparison)
prota.m$Stat<-ifelse(grepl('FDR',prota.m$Comparison),'FDR',prota.m$Stat)
prota.m$Comparison<-gsub('FDR_','',prota.m$Comparison)


prota.ma<-merge(prota.m,annot[,c('b_ID','Gene_name')],by.x='Gene',by.y='Gene_name',all.x=TRUE)
#prota.ms<-ddply(prota.m,.(Gene,Uniprot,Stat,Comparison),summarise,Value=mean(Value))

prota.c<-dcast(prota.ma,Gene+b_ID+Uniprot+Replicate+Comparison~Stat,value.var = 'Value')

prota.c<-subset(prota.c,!is.na(ANOVA) & !is.na(logFC))


#Order by significance
prota.c<-prota.c[order(prota.c$Gene,prota.c$Comparison,prota.c$FDR,decreasing = FALSE),]
subset(prota.c,Comparison=='TvC')

#Remove less significant duplicated values
prota.u<-prota.c[!duplicated(prota.c[,c('Gene','Comparison')]),]

prota.um<-melt(prota.u,id.vars = c('Gene','b_ID','Uniprot','Comparison'),measure.vars = c('ANOVA','FDR','logFC'),
               variable.name = 'Stat',value.name = 'Value')

prota.um



prota.uc<-dcast(prota.um,Gene+b_ID+Uniprot~Comparison+Stat,value.var = 'Value')
prota.uc<-subset(prota.uc,Gene!='<NA>')


prota.uc$TvC_Rank<-ifelse(prota.uc$TvC_logFC<0,'-1','3')
prota.uc$RvS_Rank<-ifelse(prota.uc$RvS_logFC<0,'-1','3')

head(prota.uc)

length(prota.uc$Gene)

length(intersect(prota.uc$Gene,annot$Gene_name))
length(intersect(prota.uc$Gene,annot$Synonym_1))
length(intersect(prota.uc$Gene,annot$Synonym_2))
length(intersect(prota.uc$Gene,annot$Synonym_3))
length(intersect(prota.uc$Gene,annot$Synonym_4))


write.csv(prota.uc,'Gene_proteomics.csv',quote = FALSE,row.names = FALSE)



TF<-read.csv('../Annotations/Ecoli/Transcription_Factors/network_tf_gene.txt',comment.char = '#',sep = '\t',header=FALSE)
colnames(TF)<-c('TF','Gene','Effect','Evidence','Support','NAs')
TF$NAs<-NULL




TFa<-merge(TF,prota.uc,by='Gene',all=TRUE)

TFa<-subset(TFa,!is.na(TF))



Gene_counts<-data.frame(table(TFa$Gene))
Gene_counts.u<-subset(Gene_counts,Freq==1)



TFa$Unique<-TFa$Gene %in% Gene_counts.u$Var1
TFa<-TFa[,c('TF','Gene','Effect','Evidence','Support','Unique',"TvC_ANOVA","TvC_FDR","TvC_logFC","RvS_ANOVA","RvS_FDR","RvS_logFC" )]

write.csv(TFa,'TF_gene.csv',quote = FALSE)

TFa.sel<-TFa


#TFa.sel<-subset(TFa,Unique)

length(TFa$Gene)
length(TFa.sel$Gene)


TFau<-TFa.sel[!duplicated(TFa.sel[,c('Gene','TF')]),]


allgenes<-length(unique(TFau$Gene))

TFsum<-ddply(TFau,.(TF),summarise,TF_total=length(Gene), TvC_sig=sum(TvC_ANOVA<0.05,na.rm = TRUE),RvS_sig=sum(RvS_ANOVA<0.05,na.rm = TRUE))

TvC_sig<-length(unique(subset(TFau,TvC_ANOVA<0.05)$Gene))
RvS_sig<-length(unique(subset(TFau,RvS_ANOVA<0.05)$Gene))

TFsum$TvC_total<-TvC_sig
TFsum$RvS_total<-RvS_sig

allgenes
TvC_sig
RvS_sig

TFsum$TvC_p<-phyper(TFsum$TvC_sig-1,TFsum$TvC_total,allgenes-TvCsig,TFsum$TF_total,lower.tail =FALSE)
TFsum$RvS_p<-phyper(TFsum$RvS_sig-1,TFsum$RvS_total,allgenes-RvSsig,TFsum$TF_total,lower.tail = FALSE)


TFsum$TvC_FDR<-p.adjust(TFsum$TvC_p,method='fdr')
TFsum$RvS_FDR<-p.adjust(TFsum$RvS_p,method='fdr')
TFsum<-TFsum[order(TFsum$TvC_p),]

subset(TFsum,TF=='FNR')
TvC_sig/allgenes
30/513

TFsum<-TFsum[,c('TF','TF_total','TvC_total','RvS_total','TvC_sig','RvS_sig','TvC_p','TvC_FDR','RvS_p','RvS_FDR')]
write.csv(TFsum,'TF_enrichment.csv')

all<-merge(TFa,TFsum,by='TF',all.x=TRUE)
all$Evidence<-NULL

write.csv(all,'TF_gene_all.csv',quote = FALSE)