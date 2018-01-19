options(java.parameters = "-Xmx4096m") 

#That's for annotations
library(org.Ce.eg.db)
library(org.EcK12.eg.db)
library(celegans.db)
library(GO.db)
library(biomaRt)

#library(ballgown)
library(edgeR)

library(RSkittleBrewer)
library(genefilter)

#library(dplyr)

library(devtools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(ggbiplot)


library(heatmap3)

#library(ggrepel)
library(AnnotationDbi)
library(grid)
library(gridExtra)

library(devtools)
library(data.table)

library(dtplyr)
library(pathview)

library(xlsx)


library(ggrepel)

devtools::install_github("PNorvaisas/PFun")
library(PFun)




#install.packages('dtplyr')

#install.packages('data.table')
# install_github('alyssafrazee/RSkittleBrewer')
# install_github("vqv/ggbiplot")

#library(plotly)
# 
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("org.EcK12.eg.db"))

# 
# install.packages("grid","gridExtra")


theme_set(theme_light())



write.xl<-function(data,explanations,ofile,mode='Full') {
  print(paste('Writing:',ofile))
  explst<-subset(explanations,Column %in% colnames(data))
  explst<-explst[match(colnames(data),explst$Column),]
  if (mode=='Readme'){
    write.xlsx2(explst, file=ofile, sheetName="Readme",row.names = FALSE,showNA=FALSE)
  } else {
    write.xlsx2(explst, file=ofile, sheetName="Readme",row.names = FALSE,showNA=FALSE)
    write.xlsx2(data, file=ofile, sheetName="Data", append=TRUE,row.names = FALSE)#showNA=FALSE
  }
}


mymerge<-function(all.results,results) {
  if (dim(all.results)[[1]]==0) {
    all.results<-results
  } else {
    all.results<-merge(all.results,results,all.x=TRUE,all.y=TRUE)
  }
  return(all.results)
}



#Order by expression in transcripts select data for each dataset



cwd<-"~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_Rotenone/"
keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Celegans/KEGG_pathways/'
setwd(cwd)
odir<-'Results'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


# Generate nice Excel tables for publication
# explanations<-read.xlsx2('Readme.xlsx',sheetName = 'Columns',
#                          stringsAsFactors = TRUE,
#                          header=TRUE)



pheno_data = read.csv("Pheno_data.csv")





#Ensembl annotation
celegans87 = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl",version=87)
listDatasets(celegans87)
listFilters(celegans87)
listAttributes(celegans87)

celegans.annotation <- getBM(attributes=c('wormbase_gene_seq_name','ensembl_gene_id','external_gene_name','entrezgene','description'),
                             mart = celegans87)
celegans.annotation<-rename(celegans.annotation,c('transcript_length'='length'))
head(celegans.annotation)
dim(celegans.annotation)

table(duplicated(celegans.annotation$wormbase_gene_seq_name))
#59286

#dupl.entrez<-duplicated(celegans.annotation$entrezgene)
#celegans.annotation<-celegans.annotation[!dupl.entrez,]

counts<-read.xlsx2('GSE46051_ce_22samples_DMSO_Rotenone_counts_and_RPKM.xls',sheetName = 'gene_counts',header = TRUE,colClasses = c(rep('character',3),rep('numeric',22)))
head(counts)
dim(counts)
counts$genename<-as.character(counts$genename)

length(intersect(counts$gene,celegans.annotation$wormbase_gene_seq_name))

nomatch<-setdiff(counts$gene,celegans.annotation$wormbase_gene_seq_name)
nomatch


expression<-merge(celegans.annotation,counts,by.x='wormbase_gene_seq_name',by.y='gene',all.y=TRUE)
subset(expression,is.na(ensembl_gene_id))
expression$genename

expression$external_gene_name<-ifelse(is.na(expression$external_gene_name),expression$genename,expression$external_gene_name)
expression$genename<-NULL
expression$chromosome<-NULL

dim(counts)
dim(expression)

table(duplicated(expression$wormbase_gene_seq_name))

expression<-expression[!duplicated(expression$wormbase_gene_seq_name),]
rownames(expression)<-expression$wormbase_gene_seq_name
dim(expression)

#Split table for edgeR

ancols.2<-setdiff(colnames(expression),c(as.character(pheno_data$ids)) )
ancols.2

gene.info<-expression[,ancols.2]
gene.data<-expression[,as.character(pheno_data$ids)]


head(gene.info)
head(gene.data)

pheno_data$Group<-factor(pheno_data$Group,levels=c('1_C','1_R','5_C','5_R','10_C','10_R','20_C','20_R'))
pheno_data$Group

pheno_data

dt<-DGEList(counts=gene.data,
            genes=gene.info,
            group=factor(pheno_data$Group))


#How to choose right threshold

dim(dt)
apply(dt$counts, 2, sum)
#At least one count per million (cpm) reads in at least 4 consistent samples (one group)




cpms<-cpm(dt)

ugroups<-unique(as.character(pheno_data$Group))
ugroups


thrs<- 0.25
keep<-FALSE
for (ug in ugroups){
  consis<-length(subset(pheno_data,Group==ug)$ids)
  keep<-rowSums(cpms[,grepl(ug,as.character(pheno_data$Group))]>thrs )>= consis | keep
}

dim(dt)
d <- dt[keep,]
dim(d)


cpms[!keep,]


#Normalisation
dim(dt)
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d$samples



#Visualise filtered counts
dtc<-as.data.frame(cpm(dt$counts,log=TRUE,prior.count = 1))
dtc$Feature<-rownames(dtc)

dtcm<-melt(dtc, id.vars = 'Feature',variable.name='Sample',value.name = 'logCPM')
dtcm$Sample<-as.factor(dtcm$Sample)

raw<-ggplot(dtcm,aes(x=logCPM,color=Sample))+
  geom_density()+
  xlim(-10,15)+
  ggtitle('Raw')
#raw
filtered<-ggplot(dtcm[keep,],aes(x=logCPM,color=Sample))+
  geom_density()+
  xlim(-10,15)+
  ggtitle('Filtered')
#filtered
grid.arrange(raw,filtered,ncol=2,
             top='logCPM before and after filtering')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_logCPM_filtering.pdf",sep=''),
             width=12,height=6, useDingbats=FALSE)



#method="bcv",
plotMDS(d,col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:8, pch=10)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)


#method="bcv",
plotMDS(cpms,col=as.numeric(d$samples$group))
legend("bottomright", as.character(unique(d$samples$group)), col=1:8, pch=10)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV_nonorm_nofilt.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)




#Heatmap
logCPM <- cpm(d, log=TRUE, prior.count=1)

heatshape<-as.matrix(logCPM)
colnames(heatshape)==pheno_data$ids

day.fac<-pheno_data$Day
treat.fac<-pheno_data$Treatment

day.col<-revalue(as.character(day.fac),c("1"="plum","5"="red","10"="red4","20"="black"))
treat.col<-revalue(as.character(treat.fac),c("Control"="white","Rotenone"="black"))

#Create legend
clab<-cbind(treat.col,day.col)
colnames(clab)<-c('Treatment','Day')


hmap<-heatmap3(as.matrix(heatshape),
               scale = 'row',
               #hclustfun = function(x) hclust(dist(x,method="manhattan"),method="complete"),
               #distfun = function(x) dist(x,method="manhattan"),
               labRow=FALSE,
               #key=TRUE,
               #ColSideLabs = 'Group',
               ColSideColors = clab)

legend('topright',legend=c('1','5','10','20'),fill=c("plum","red","red4","black"), border=TRUE, bty="n",title='Day')
legend('right',legend=c('Control','Rotenone'),fill=c('white','black'), border=TRUE, bty="n",title='Treatment')


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Heatmap_logCPM.pdf",sep=''),
             width=12,height=9, useDingbats=FALSE)




#GLM
design.mat <- model.matrix(~ 0 + d$samples$group) #
colnames(design.mat) <- c(levels(d$samples$group))
design.mat

# priorn<-5
# priordf<-55
# d2<-estimateDisp(d,design.mat,prior.df = priordf)
# print(paste('Prior.n=',d2$prior.n,' Residual.df=',d2$prior.df/d2$prior.n,' Prior.df=',d2$prior.df,sep=''))


d2<-estimateDisp(d,design.mat)

plotBCV(d2)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Dispersion.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#All design matrices
design.mat
fit <- glmFit(d2, design.mat)



contrasts.all<-read.contrasts('!Contrasts.xlsx','Contrasts_values',colnames(design.mat))



contrasts<-contrasts.all$Contrasts.matrix[,colnames(design.mat)]
contrasts


allcomparisons<-contrasts.all$Contrasts.table
allcomparisons

allcomp<-rownames(contrasts)



gp.KEGG<-subset(getGeneKEGGLinks('cel',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) & PathwayID!='path:cel01100')
#gp.GO<-subset(getGeneGOLinks('cel',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) )





#GLM testing
all.results.None<-data.table()
all.results.1FC<-data.table()
all.results<-data.table()

all.KEGGenrichment<-data.table()
all.GOenrichment<-data.table()

summary.all<-data.frame()

for (cid in 1:nrow(allcomparisons)) {
  comp<-allcomparisons[cid,'Contrast']
  desc<-allcomparisons[cid,'Description']
  #contr<-contrasts[[as.character(comp)]]
  contr<-contrasts[comp,]
  
  print(as.character(comp))
  print(as.character(desc))
  print(contr)
  
  thre<-'None'
  enrcoef<-'logFC'
  print(paste('Threshold:',thre,sep=' '))
  lrt <- glmLRT(fit, contrast=contr)
  
  print('...Differential gene expression')
  de.glm <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
  print(summary(de.glm))
  
  sumr<-t(summary(de.glm))
  rownames(sumr)<-comp
  summary.all<-rbind(summary.all,sumr)
  
  de2tags <- rownames(d2)[as.logical(de.glm)]
  cairo_pdf(paste(odir,'/MVA_',comp,'_Threshold-',thre,'.pdf',sep = ''),width=9,height=6)
  plotSmear(lrt, de.tags=de2tags)
  abline(h = c(-1, 1), col = "blue")
  title(paste(desc,comp,sep='\n'))
  dev.off()
  
  result<-as.data.frame(topTags(lrt,n=dim(d2$genes)[1]))
  result$Comparison<-comp
  result$Threshold<-thre
  
  all.results<-rbind(all.results,data.table(result))
  
  if (summary(de.glm)[3]>10 | summary(de.glm)[1] > 10) {
    print('...KEGG enrichment')
    enrichment.KEGG<-kegga(lrt,coef=enrcoef,geneid = 'entrezgene',species.KEGG='cel',gene.pathway = gp.KEGG,convert=TRUE,FDR=0.05,trend=TRUE)#,trend='length'
    
    enrichment.KEGG$ID<-rownames(enrichment.KEGG)
    enrichment.KEGG<-rename(enrichment.KEGG,c("Pathway"="Description"))
    enrichment.KEGG$Comparison<-comp
    enrichment.KEGG$Threshold<-thre
    all.KEGGenrichment<-rbind(all.KEGGenrichment,data.table(enrichment.KEGG))
    
    print('...GO enrichment')
    enrichment.GO<-goana(lrt,species='Ce',coef=enrcoef,geneid='entrezgene',convert=TRUE,FDR=0.05,trend=TRUE)#,trend='length'
    
    enrichment.GO$ID<-rownames(enrichment.GO)
    enrichment.GO<-rename(enrichment.GO,c("Term"="Description"))
    enrichment.GO$Comparison<-comp
    enrichment.GO$Threshold<-thre
    all.GOenrichment<-rbind(all.GOenrichment,data.table(enrichment.GO))
    
  } else {
    print("No DE genes!")
  }
  print('---------')

}


colnames(summary.all)<-c('Contrast','Change','Count')

summary.counts<-dcast(summary.all,Contrast~Change,value.var='Count')
summary.counts

write.csv(summary.counts,paste(odir,'/Contrast_summary.csv',sep = ''),row.names = FALSE)


all.GOenrichment<-data.frame(all.GOenrichment)
dim(all.GOenrichment)

all.KEGGenrichment<-data.frame(all.KEGGenrichment)
dim(all.KEGGenrichment)

all.results$Comparison<-factor(all.results$Comparison,
                               levels=allcomp,
                               labels=allcomp)
all.results<-data.frame(all.results)
dim(all.results)


#Get gene pathway mappings
path2eg<-getGeneKEGGLinks('cel',convert = TRUE)
path2eg<-subset(path2eg,!is.na(GeneID) & !is.na(PathwayID) & PathwayID!='path:cel01100')
path2eg$PathID<-gsub('path:','',path2eg$PathwayID)
path2eg$EntrezGene<-path2eg$GeneID
#head(path2eg)


path2pathde<-getKEGGPathwayNames('cel',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)

go2eg<-toTable(org.Ce.egGO2EG)
go2eg<-subset(go2eg,!is.na(gene_id) & !is.na(go_id) )
go2eg$GO_id<-go2eg$go_id
go2eg$EntrezGene<-go2eg$gene_id
#head(go2eg)


allKEGGinfo.genes<-ddply(path2eg,.(PathID),summarise,Genes=paste(EntrezGene,collapse=';'))
allKEGGinfo.paths<-ddply(path2eg,.(EntrezGene),summarise,Pathways=paste(PathID,collapse=';'))


allGOinfo.genes<-ddply(go2eg,.(GO_id),summarise,Genes=paste(EntrezGene,collapse=';'))
allGOinfo.terms<-ddply(go2eg,.(EntrezGene),summarise,GO_terms=paste(GO_id,collapse=';'))


all.results.k<-merge(all.results,allKEGGinfo.paths,by.x='entrezgene',by.y='EntrezGene',all.x=TRUE)
all.results.kg<-merge(all.results.k,allGOinfo.terms,by.x='entrezgene',by.y='EntrezGene',all.x=TRUE)

options(scipen=5)




write.csv(all.results.kg,paste(odir,'/All_results.csv',sep = ''),row.names = FALSE)
#write.xl(all.results.kg,explanations,paste(odir,'/All_results.xlsx',sep = ''))

write.csv(all.KEGGenrichment,paste(odir,'/All_results_KEGG.csv',sep = ''),row.names = FALSE)
#write.xl(all.KEGGenrichment,explanations,paste(odir,'/All_results_KEGG.xlsx',sep = ''))

write.csv(all.GOenrichment,paste(odir,'/All_results_GO.csv',sep = ''),row.names = FALSE)
#write.xl(all.GOenrichment,explanations,paste(odir,'/All_results_GO.xlsx',sep = ''))




head(all.results.kg)

change_sum<-ddply(all.results.kg,.(Comparison),summarise,FDR_UP=sum(FDR<0.05 & logFC>0),FDR_None=sum(FDR>0.05),FDR_DOWN=sum(FDR<0.05 & logFC<0),
      P_UP=sum(PValue<0.05 & logFC>0),P_None=sum(PValue>0.05),P_DOWN=sum(PValue<0.05 & logFC<0))

write.csv(change_sum,paste(odir,'/DGE_counts.csv',sep = ''))

dim(all.results.kg)





#Heatmaps ans volcano plots
annotations<-c('KEGG','GO')


reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 5)
pearson<-function(x) as.dist((1-cor(t(x)))/2)

enbrks<-c(0,-log(0.05,10),2,3,4,5)



#Volcano plots

txtsize<-2
nudgex<-0.1
nudgey<-0.1
erralpha<-1
errcolor<-'grey80'
baralpha<-0.2
barcolor<-'red'

mapkey<-'entrezgene'

sp<-1
tp<-25


#For pathways
allcomp



#For volcano
volc.grps<-allcomp


#For heatmaps
heat.all<-allcomp
heat.write<-allcomp

heat.genh<-c('GenR','Age_general')
heat.age<-c('Age5_C','Age10_C','Age20_C','GenAge5','GenAge10','GenAge10v5')#,'Age_general'
heat.rot<-c('Rot_1','Rot_5','Rot_10','Rot_20','GenR')
heat.int<-c('Int5','Int10','Int20')

heat.grps<-list('All'=heat.all,'Write'=heat.write,'General'=heat.genh,'Aging'=heat.age,'Rotenone'=heat.rot,'Interaction'=heat.int)

names(heat.grps)




#For KEGG mappings
KEGG.main<-paste(setdiff(allcomp,c('Age27_eat2','Age_eat2','DR_general','Age_general','DR_27')),'logFC',sep='_')
KEGG.gen<-paste(c('GenR','GenAge5','Int5') ,'logFC',sep='_')
KEGG.list<-list('Main'=KEGG.main,'General'=KEGG.gen)




add.path<-c('cel04213','cel04212')
allpaths<-c()

setwd(cwd)


head(all.results.kg)




ancols.3<-c('Comparison','wormbase_gene_seq_name','ensembl_gene_id','external_gene_name','entrezgene','description','Pathways','GO_terms')

measure.vars<-c('logFC','FDR')

# 
# t_name+gene_name+entrezgene+ensembl_transcript_id+ensembl_gene_id+external_transcript_name+description+Transcripts+ExpressedTranscripts+Pathways+GO_terms+
#   chr+strand+start+end+num_exons+length



thres<-'None'
enrcoef<-'logFC'
all.results.r<-all.results.kg
#all.results.r$unshrunk.logFC<-NULL
all.results.r$Threshold<-NULL

#Melt
all.results.rm<-melt(all.results.r,id.vars = ancols.3,
                     measure.vars=measure.vars,
                     variable.name = 'Stat',value.name = 'Value')

all.results.rmc<-subset(all.results.rm,Stat %in% c(enrcoef,'FDR'))
all.results.rcp<-dcast(all.results.rmc,wormbase_gene_seq_name+ensembl_gene_id+external_gene_name+entrezgene+description+Pathways+GO_terms
                         ~Comparison+Stat,value.var = 'Value')

head(all.results.rcp)

#Reorder by significance
signsum<-apply(all.results.rcp[,grep('FDR',colnames(all.results.rcp))],1,function(x) sum(-log10(x)))
all.results.rcp<-all.results.rcp[order(signsum,decreasing = TRUE),]

nrow(all.results.rcp)

gdata<-subset(all.results.rcp,! is.na(entrezgene))
nrow(gdata)


#Do selection by gene id, because gene_id is associated with function
dupl<-duplicated(gdata[,mapkey])
print('Entrez Gene duplicates')
print(table(dupl))


#nrow(gdataf)

write.csv(all.results.rcp,paste(odir,'/All_results_sidebyside_Threshold-',thres,'.csv',sep = ''),row.names = FALSE)
#write.xl(all.results.rcp,explanations,paste(odir,'/All_results_sidebyside_Threshold-',thres,'.xlsx',sep = ''),'Readme')


head(all.results.r)

#Volcano plots
all.plot<-subset(all.results.r,Comparison %in% volc.grps)

head(all.plot)
vol<-ggplot(all.plot,aes(x=logFC,y=-log10(FDR)))+
  geom_hline(yintercept = -log(0.05,10),color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = -1,color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = 1,color=barcolor,alpha=baralpha)+
  ggtitle(paste('Volcano plot for differentially expressed genes. Threshold: ',thres,sep=''),
          subtitle='Labels are shown for some of the significantly affected genes')+
  geom_point(alpha=0.9,size=1)+
  xlim(-15,15)+
  #   geom_label_repel(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp,as.character(gene_name),'')),
  #             size=txtsize,colour = "red")+
  geom_text(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp & FDR<0.05,as.character(external_gene_name),'')),
            hjust=1, vjust=-0.5,size=txtsize,colour = "red")+
  facet_wrap(~Comparison,ncol = 3)
fname<-paste(odir,'/Volcano_Threshold-',thres,'.pdf',sep = '')
print(fname)
cairo_pdf(fname,width=30,height=60)
print(vol)
dev.off()




#Enrichment heatmaps
for (ent in annotations) {
  print(paste('Enrichment:',ent),sep=' ')
  if (ent=='KEGG') {
    all.enrichment<-subset(all.KEGGenrichment,Threshold==thres)
    idvars<-c('ID','Description','N','Comparison')
  } else {
    all.enrichment<-subset(all.GOenrichment,Threshold==thres)
    idvars<-c('ID','Ont','Description','N','Comparison')
  }
  all.enrichment$ID<-gsub('path:','',all.enrichment$ID)
  all.enrichment$logP_Up<- -log10(all.enrichment$P.Up)
  all.enrichment$logP_Down<- -log10(all.enrichment$P.Down)
  all.enrichment$Comparison<-factor(all.enrichment$Comparison,
                                    levels=allcomp,
                                    labels=allcomp)
  
  
  enrichmentm<-melt(all.enrichment,id.vars = idvars,variable.name = 'Stat',value.name = 'Value')
  #enrichmentm$Value<-as.numeric(enrichmentm$Value)
  enrichmentfe<-subset(enrichmentm,Stat %in% c('P.Up','P.Down') & Comparison %in% heat.write)
  enrichmentf<-subset(enrichmentm,Stat %in% c('logP_Up','logP_Down'))
  
  
  if (ent=='KEGG') {
    #For writing
    enrichmentce<-dcast(enrichmentfe,ID+Description+N~Comparison+Stat,value.var = 'Value',fill = NA)
    enrichmentcea<-merge(enrichmentce,allKEGGinfo.genes,by.x='ID',by.y='PathID',all.x=TRUE)
    #For plotting
    #enrichmentf<-subset(enrichmentm,Stat %in% c('logP_Up','logP_Down') & Comparison %in% comp.heat)
    enrichmentc<-dcast(enrichmentf,ID+Description+N~Comparison+Stat,value.var = 'Value',fill = NA)
    rownames(enrichmentc)<-paste(enrichmentc$ID,enrichmentc$Description,sep=' ')
    ontologies<-c(1)
    #nrow(gdataf)
  } else {
    #For writing
    enrichmentce<-dcast(enrichmentfe,ID+Ont+Description+N~Comparison+Stat,value.var = 'Value',fill = NA)
    enrichmentcea<-merge(enrichmentce,allGOinfo.genes,by.x='ID',by.y='GO_id',all.x=TRUE)
    #For plotting
    enrichmentc<-dcast(enrichmentf,ID+Ont+Description+N~Comparison+Stat,value.var = 'Value',fill = NA)
    rownames(enrichmentc)<-paste(enrichmentc$Ont,enrichmentc$ID,enrichmentc$Description,sep=' ')
    ontologies<-c('BP','MF','CC')
  }
  
  write.csv(enrichmentcea,paste(odir,'/All_results_sidebyside_',ent,'_Threshold-',thres,'.csv',sep = ''),row.names = FALSE)
  #write.xl(enrichmentcea,explanations,paste(odir,'/All_results_sidebyside_',ent,'_Threshold-',thres,'.xlsx',sep = ''),'Readme')
  
  for (ont in ontologies) {
    for (comp.grh in names(heat.grps)) {
      comp.heat<-heat.grps[[comp.grh]]
      if (ent=='GO') {
        enrichmentco<-subset(enrichmentc,Ont==ont)
        fname<-paste(odir,'/Enrichment_heatmap_',ent,'-',ont,'_Threshold-',thres,'_',comp.grh,'.pdf',sep = '')
      } else {
        enrichmentco<-enrichmentc
        fname<-paste(odir,'/Enrichment_heatmap_',ent,'_Threshold-',thres,'_',comp.grh,'.pdf',sep = '')
      }
      
      #hdata<-enrichmentco[,!colnames(enrichmentco) %in% c('ID','Description','N','Ont')]
      #Filter just selected columns
      hdata<-enrichmentco[,grep(paste(comp.heat,collapse="|"),colnames(enrichmentco)) ]
      dim(hdata)
      hdata<-hdata[!apply(hdata, 1, function(x) {all(as.numeric(x) < -log(0.05,10),na.rm=TRUE)}),]
      dim(hdata)
      
      #allpaths<-as.character(rownames(hdata))
      
      hdatafill<-hdata
      hdatafill[is.na(hdatafill)]<-0
      hdatafill[hdatafill==Inf]<-30
      
      print(paste('Enrichment terms to plot:',nrow(hdata)))
      hmap<-heatmap.2(data.matrix(hdatafill),key=TRUE,Colv=FALSE,trace='none',col=gyrs,
                      xlab='Comparison',Rowv=TRUE,
                      dendrogram="row",scale="none",na.color="white",
                      cexRow=0.8,cexCol=0.5,symkey=FALSE,
                      breaks=enbrks,
                      reorderfun=reorderfun_mean)
      print(fname)
      hgh<-dim(hdata)[1]/6
      wdh<-7
      #17
      cairo_pdf(fname,width=wdh,height=hgh+wdh/2)
      heatmap.2(data.matrix(hdata),key=FALSE,Colv=FALSE,
                trace='none',col=gyrs,
                xlab='Comparison',
                Rowv=hmap$rowDendrogram,key.xlab='-log10(FDR)',
                dendrogram='row',scale="none",na.color="white",
                cexRow=0.7,cexCol=0.7,margin=c(10,20),
                lwid=c(0.2,0.8),symkey=FALSE,lhei=c(0.05,0.95),
                breaks=enbrks)
      #margin=c(16,16),
      dev.off()
    }
  }
}



#KEGG plotting
nrow(all.results.rcp)
gdata<-subset(all.results.rcp,! is.na(entrezgene))
nrow(gdata)
dupl<-duplicated(gdata[,mapkey])
print('Entrez Gene duplicates')
print(table(dupl))
gdataf<-gdata[!dupl,]
rownames(gdataf)<-gdataf[,mapkey]

allpaths<-gsub('cel','',allKEGGinfo.genes$PathID)

print(paste('Total KEGG pathways to plot:',length(allpaths)))


keggdir<-paste(odir,'/KEGG',sep='')
dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd(keggdir)
comp.gr<-'General'
pathcomp<-KEGG.list[[comp.gr]]
pathcomp
outsuffx<-comp.gr
pv.out <- pathview(gene.data = gdataf[,pathcomp,drop=FALSE],
                   pathway.id = allpaths,
                   species = "cel",
                   out.suffix = outsuffx,
                   kegg.dir = keggxml,
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   min.nnodes = 0,
                   limit=list(gene=1,cpd=1),node.sum = "mean",
                   low = list(gene = "blue", cpd = "green"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="yellow"))
setwd(cwd)

all.kegg.mappingst<-data.frame()
paths<-names(pv.out)
for (pth in paths) {
  print(pth)
  result<-pv.out[[pth]][['plot.data.gene']]
  if(!is.null(result) ){
    result$PathwayID<-pth
    all.kegg.mappingst<-mymerge(all.kegg.mappingst,result)
  }
}


all.kegg.mappingst$mol.col<-NULL

all.kegg.mappings<-merge(path2pathde,all.kegg.mappingst,by='PathwayID',all.y=TRUE)




write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings.csv',sep = ''),row.names = FALSE)
#write.xl(all.kegg.mappings,explanations,paste(odir,'/All_KEGG_mappings.xlsx',sep=''))





#Venn
library('Vennerable')
#Might need to deactive reshape and Vennerable to use melt in other code


DR_3.change<-subset(all.results.rcp,`DR_3_FDR`<0.05 )[,'gene_id']
DR_15.change<-subset(all.results.rcp,`DR_15_FDR`<0.05 )[,'gene_id']
DR_27.change<-subset(all.results.rcp,`DR_27_FDR`<0.05 )[,'gene_id']

venn.DR<-list('DR 3 days'=DR_3.change,
              'DR 15 days'=DR_15.change,
              'DR 27 days'=DR_27.change)


plot(Venn(venn.DR),show = list(Faces = FALSE),
     doWeights = FALSE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_DR.pdf",sep=''),width=6,height=6)



Age_N2.change<-subset(all.results.rcp,`Age_N2_FDR`<0.05 )[,'gene_id']
Age_eat2.change<-subset(all.results.rcp,`Age_eat2_FDR`<0.05 )[,'gene_id']
Age27_eat2.change<-subset(all.results.rcp,`Age27_eat2_FDR`<0.05 )[,'gene_id']

venn.age<-list('Age N2'=Age_N2.change,
               'Age eat-2'=Age_eat2.change,
               'Age eat-2 (27)'=Age27_eat2.change)


plot(Venn(venn.age),show = list(Faces = FALSE),
     doWeights = FALSE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_Age.pdf",sep=''),width=6,height=6)



Age_N2.up<-subset(all.results.rcp,`Age_N2_FDR`<0.05 & `Age_N2_logFC`>0)[,'gene_id']
Age_N2.down<-subset(all.results.rcp,`Age_N2_FDR`<0.05 & `Age_N2_logFC`<0)[,'gene_id']

Age_eat2.up<-subset(all.results.rcp,`Age_eat2_FDR`<0.05 & `Age_eat2_logFC`>0)[,'gene_id']
Age_eat2.down<-subset(all.results.rcp,`Age_eat2_FDR`<0.05 & `Age_eat2_logFC`<0)[,'gene_id']

venn.age.dir<-list('Age eat-2 down'=Age_eat2.down,'Age N2 down'=Age_N2.down,'Age eat-2 up'=Age_eat2.up,'Age N2 up'=Age_N2.up)


plot(Venn(venn.age.dir),show = list(Faces = FALSE),
     doWeights = FALSE,type='ellipses')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_Age_direction.pdf",sep=''),width=6,height=6)



