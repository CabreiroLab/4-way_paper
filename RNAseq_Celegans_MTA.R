options(java.parameters = "-Xmx4096m") 

#That's for annotations
# library(org.Ce.eg.db)
# library(celegans.db)
# library(GO.db)
library(biomaRt)

# library(ballgown)
# library(edgeR)

#library(genefilter)

#library(dplyr)

library(devtools)



library(tidyverse)

library(reshape2)



library(heatmap3)
library(gplots)



library(grid)
library(gridExtra)




library(xlsx)


library(ggrepel)

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




cwd<-"~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/"
setwd(cwd)


odir<-'Results_1thrs_newannot/MTA'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





#Ensembl annotation
celegans87 = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl",version=87)
listDatasets(celegans87)
listFilters(celegans87)

listAttributes(celegans87)

celegans.annotation <- getBM(attributes=c('wormbase_gene_seq_name','ensembl_gene_id','external_gene_name','description'),
                             mart = celegans87)
#
head(celegans.annotation)
dim(celegans.annotation)



celegans.annotation<-celegans.annotation[!duplicated(celegans.annotation[,c('wormbase_gene_seq_name','ensembl_gene_id')]),]


head(celegans.annotation)
dim(celegans.annotation)







MTAfiles<-c('bootstrap_MTA_SM_to_S.xlsx','bootstrap_MTA_SM_to_RM.xlsx','bootstrap_MTA_RM_to_R_fixed.xlsx')
conts<-c('SM_S','SM_RM','RM_R')

trans<-c('SM_to_S','SM_to_RM','RM_to_R')

MTAinfo<-data.frame('File'=MTAfiles,'Contrast'=conts,'Transition'=trans) %>%
  mutate(Contrast=as.character(Contrast),
         Transition=as.character(Transition))






MTAall<-data.frame()

for (fli in 1:nrow(MTAinfo)){
  
  fl<-as.character(MTAinfo[fli,'File'])
  print(fl)
  dat<-read.xlsx(paste0('~/Dropbox/Projects/OP50_Celegans_holobiont/Data/RNAseq_Bootstrapping_separate_4_removed_1thrs/results/',fl),
                 sheetName = paste0('bootstrap_MTA_',MTAinfo[fli,'Transition']),check.names=FALSE) %>%
    rename(Gene=Var1) %>%
    gather(Bootstrap,Value,contains('_')) %>%
    mutate(Contrast=as.character(MTAinfo[fli,'Contrast']))
  
  
  # dat<-dat.r[,c('Gene number','Gene ID','NA','Mean')]
  # dat<-plyr::rename(dat,c('NA'='Description') )
  # dat$Contrast<-MTAinfo[fli,'Contrast']
  
  MTAall<-rbind(MTAall,dat)
  
}


cuts<-MTAall%>%
  filter(Gene=='cutoff') %>%
  rename(Cutoff=Value) %>%
  mutate(Gene=NULL)


MTAvals<-MTAall%>%
  filter(Gene!='cutoff') %>%
  left_join(cuts) %>%
  mutate(Pass=Value>Cutoff*1.01) %>%
  group_by(Contrast,Bootstrap,Pass) %>%
  mutate(Rank=rank(Value),
         Percentile=rank(Value)*100/length(Value),
         N=length(Value)) %>%
  ungroup %>%
  mutate(Percentile=ifelse(Pass,Percentile,0),
         Rank=ifelse(Pass,Rank,0),
         Complete=ifelse(Bootstrap=='None_None','Complete','Trimmed')) %>%
  gather(PR,Meas,Percentile,Rank)

  
  

  
  
MTAsum<-MTAvals%>% 
  #Summarise over all bootstraps
  group_by(Gene,Contrast,Complete,PR) %>%
  summarise(Mean=mean(Meas),
            SD=sd(Meas),
            Median=median(Meas),
            Q25=quantile(Meas,0.25),
            Q75=quantile(Meas,0.75)) %>%
  ungroup %>%
  gather(Stat,Value,Median,Mean,SD,Q25,Q75) %>%
  spread(Complete,Value) %>%
  mutate(Score=ifelse(grepl('_Mean|_Median',Stat),0.2*Complete+0.8*Trimmed,Trimmed),
         PR=ifelse(PR=='Percentile','P','R')) %>%
  mutate(ensembl_gene_id=celegans.annotation[match(as.character(Gene),celegans.annotation$external_gene_name),'ensembl_gene_id'],
         ensembl_gene_id=ifelse(is.na(ensembl_gene_id),
                                celegans.annotation[match(as.character(Gene),celegans.annotation$wormbase_gene_seq_name),'ensembl_gene_id'],
                                ensembl_gene_id)) %>%
  left_join(celegans.annotation)


MTA<-MTAsum %>%
  unite(Contrast_PR_Stat,Contrast,PR,Stat,remove = TRUE) %>%
  select(-Complete,-Trimmed) %>%
  spread(Contrast_PR_Stat,Score)
  





MTAvals.sel<-subset(MTAvals,PR=='Percentile' & Gene %in% as.character(subset(MTA,SM_S_P_Median>80)$Gene) )
MTAsum.sel<-subset(MTAsum,PR=='P' & Gene %in% as.character(subset(MTA,SM_S_P_Median>80)$Gene) & Stat=='Median')


ggplot(MTAvals.sel,aes(x=Gene,y=Meas))+
  geom_violin()+
  geom_point()+
  geom_point(data=MTAsum.sel,aes(x=Gene,y=Score),color='red')+
  ylab('Percentile')+
  facet_grid(.~Contrast)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




head(MTAsum)
View(MTA)

subset(MTA,Contrast=='SM_S')

head(MTA)
  



head(MTAval)



length(colnames(MTAall))

length(intersect(colnames(MTAall),colnames(dat)))




subset(MTA,`Gene ID`=='C01B4.6')

subset(MTA,`Gene ID`=='T20B3.1')


head(MTA)




head(celegans.annotation)



write.csv(MTA,paste(odir,'/MTA_results.csv',sep=''),row.names = FALSE)



erralpha<-1
errcolor<-'grey80'


q<-0.95
qx<-quantile(MTA$SM_S_P_Mean,q,na.rm=TRUE) 
qy<-quantile(MTA$RM_R_P_Mean,q,na.rm=TRUE)

ggplot(MTA,aes(x=SM_S_P_Mean,y=RM_R_P_Mean))+
  geom_errorbar(aes(ymin=RM_R_P_Mean-RM_R_P_SD,ymax=RM_R_P_Mean+RM_R_P_SD),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=SM_S_P_Mean-SM_S_P_SD,xmax=SM_S_P_Mean+SM_S_P_SD),alpha=erralpha,color=errcolor,height=0)+
  geom_vline(xintercept = qx,color='red',alpha=0.5)+
  geom_hline(yintercept = qy,color='red',alpha=0.5)+
  geom_point()+
  geom_text_repel(aes(label=ifelse( SM_S_P_Mean>qx & RM_R_P_Mean>qy,as.character(Gene),'') ) )+
  scale_x_continuous(breaks=seq(0,100,by=10))+
  scale_y_continuous(breaks=seq(0,100,by=10))+
  coord_cartesian(xlim=c(0,100),ylim=c(0,100))+
  xlab('SM-S average percentile')+
  ylab('RM-R average percentile')
#  geom_text(aes(label=ifelse( SM_S_Percentile>95& RM_R_Percentile>95,as.character(`Gene ID`),'') ) )

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Scatter_SM_S_RM_R_effects.pdf",sep=''),width=6,height=6)






q<-0.95
qx<-quantile(MTA$SM_S_P_Median,q,na.rm=TRUE) 
qy<-quantile(MTA$RM_R_P_Median,q,na.rm=TRUE)

ggplot(MTA,aes(x=SM_S_P_Median,y=RM_R_P_Median))+
  geom_errorbar(aes(ymin=RM_R_P_Q25,ymax=RM_R_P_Q75),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=SM_S_P_Q25,xmax=SM_S_P_Q75),alpha=erralpha,color=errcolor,height=0)+
  geom_vline(xintercept = qx,color='red',alpha=0.5)+
  geom_hline(yintercept = qy,color='red',alpha=0.5)+
  geom_point()+
  geom_text_repel(aes(label=ifelse( SM_S_P_Median>qx & RM_R_P_Median>qy,as.character(Gene),'') ) )+
  scale_x_continuous(breaks=seq(0,100,by=10))+
  scale_y_continuous(breaks=seq(0,100,by=10))+
  coord_cartesian(xlim=c(0,100),ylim=c(0,100))+
  xlab('SM-S median percentile')+
  ylab('RM-R median percentile')




dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Scatter_SM_S_RM_R_effects_Median.pdf",sep=''),width=6,height=6)




q<-0.95
qx<-quantile(MTA.a$SM_S_Mean,q,na.rm=TRUE) 
qy<-quantile(MTA.a$SM_RM_Mean,q,na.rm=TRUE)

ggplot(MTA.a,aes(x=SM_S_Mean,y=SM_RM_Mean))+
  geom_vline(xintercept = qx,color='red',alpha=0.5)+
  geom_hline(yintercept = qy,color='red',alpha=0.5)+
  scale_x_continuous(breaks=seq(0,1000,by=100),limits=c(0,1000))+
  scale_y_continuous(breaks=seq(0,1000,by=100),limits=c(0,1000))+
  geom_point()+
  geom_text_repel(aes(label=ifelse( SM_S_Percentile> 95& SM_RM_Percentile>95 ,as.character(`Gene ID`),'') ) )

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Scatter_SM_S_SM_RM_effects.pdf",sep=''),width=6,height=6)





head(MTA.a)

#Venn
pr<-95
SM_S_5<-as.character(subset(MTA.a,SM_S_Percentile>95 )$ensembl_gene_id)
RM_R_5<-as.character(subset(MTA.a,RM_R_Percentile>95)$ensembl_gene_id)
SM_RM_5<-as.character(subset(MTA.a,SM_RM_Percentile>95)$ensembl_gene_id)


metf<-list('SM-S top 5%'=SM_S_5,
                  'RM-R top 5%'=RM_R_5)

library('Vennerable')
plot(Venn(metf),show = list(Faces = FALSE),
     doWeights = FALSE)#,type='ellipses'

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_SM_S_RM_R_effects.pdf",sep=''),width=6,height=6)


inter<-list('SM-S top 5%'=SM_S_5,
            'SM-RM top 5%'=SM_RM_5)

plot(Venn(inter),show = list(Faces = FALSE),
     doWeights = FALSE)#,type='ellipses'

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_SM_S_SM_RM_effects.pdf",sep=''),width=6,height=6)



# source("https://bioconductor.org/biocLite.R")
# biocLite("RDAVIDWebService")

library("RDAVIDWebService")



david<-DAVIDWebService$new(email="povilas.norvaisas.14@ucl.ac.uk", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

BG <- addList(david, unique(MTA.a$ensembl_gene_id), idType="ENSEMBL_GENE_ID", listName="All", listType="Background")

getCurrentBackgroundListPosition(david)

setCurrentBackgroundPosition(david, 1)


setTimeOut(david, 500000)


comps<-c('SM_S','SM_RM','RM_R')


DAVID<-data.frame()
for (ctr in comps) {
  for (thr in c(95,90)) {
    print(paste(ctr,thr))
    genes<-MTA.a[MTA.a[,paste0(ctr,'_Percentile')]>thr,'ensembl_gene_id']
    print(length(genes))
    
    FG <- addList(david, genes, idType="ENSEMBL_GENE_ID", listName=paste(ctr,thr,sep='_'), listType="Gene")
    
    setCurrentBackgroundPosition(david, 1)
    
    print("Sleeping 5s")
    Sys.sleep(5)
    
    getFunctionalAnnotationChartFile(david, paste0('Results_1thrs_newannot/MTA/',paste0('DAVID_',ctr,'_',thr,'.tsv')))
    
    FuncAnnotChart <- getFunctionalAnnotationChart(david)
    
    FuncAnnotChart$Contrast<-ctr
    FuncAnnotChart$Threshold<-thr
    
    DAVID<-rbind(DAVID,FuncAnnotChart)

  }
}


DAVID<-subset(DAVID, Term!='cel01130:Biosynthesis of antibiotics')


write.csv(DAVID,paste(odir,'/DAVID_results.csv',sep=''),row.names = FALSE)


head(DAVID)


# DAVIDfiles<-list.files('Results_1thrs_newannot/MTA/')
# DAVIDfiles<-DAVIDfiles[grepl('DAVID_',DAVIDfiles)&grepl('.txt',DAVIDfiles)]
# 
# 
# DAVID<-data.frame()
# for (fl in DAVIDfiles) {
#   print(fl)
#   contr<-gsub('DAVID_','',fl)
#   contr<-gsub('.txt','',contr)
#   print(contr)
#   dv<-read.table(paste0('Results_1thrs_newannot/MTA/',fl),sep='\t',header = TRUE,check.names = FALSE)
#   dv$Contrast<-contr
#   DAVID<-rbind(DAVID,dv)
# }
# 
# 
# head(DAVID)
# 
# DAVID<-subset(DAVID, Term!='cel01130:Biosynthesis of antibiotics')
# 
# 
# write.csv(DAVID,paste(odir,'/DAVID_results.csv',sep=''),row.names = FALSE)



DAVID.c<-dcast(subset(DAVID,Threshold==90),Category+Term~Contrast,value.var = 'FDR')



head(DAVID.c)
dim(DAVID.c)


DAVID.f<-DAVID.c[apply(DAVID.c[,c('SM_S','RM_R','SM_RM')],1,function(x) any(x<0.05,na.rm = TRUE)),]
dim(DAVID.f)
rownames(DAVID.f)<-paste(DAVID.f$Category,DAVID.f$Term)





reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 5)
enbrks<-c(0,-log(0.05,10),2,3,4,5)





hdata<-DAVID.f[,comps]



hdata<- -log10(hdata)

hdatafill<-hdata
hdatafill[is.na(hdatafill)]<-0
hdatafill[hdatafill==Inf]<-30

#Use the dummy dataset to get ordering of rows (enrichment terms)
#Graphical settings here are irrelevent for the final plot
#Data needs to be supplied as a matrix - table with row/column names defined implicitly and only numeric values
hmap<-heatmap.2(data.matrix(hdatafill),
                key=TRUE,
                Colv=FALSE,
                trace='none',
                col=gyrs,
                xlab='Comparison',
                Rowv=TRUE,
                dendrogram="row",
                scale="none",
                na.color="white", 
                symkey=FALSE,
                breaks=enbrks,
                reorderfun=reorderfun_mean,
                cexRow=0.8,
                cexCol=0.5)

#Use ordering with the real data
heatmap.2(data.matrix(hdata),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=hmap$rowDendrogram, # Here comes in the ordering
          trace='none',
          col=gyrs, # Colour scale
          breaks=enbrks, # Colour breaks
          xlab='Comparison',
          dendrogram='none', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))

dev.copy2pdf(device=cairo_pdf,file=paste0(odir,'/Enrichment_heatmap_DAVID_top10.pdf'),
             width=10,height=7,useDingbats=FALSE)




View(DAVID.c)
