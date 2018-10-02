options(java.parameters = "-Xmx4096m") 

#That's for annotations
library(org.Ce.eg.db)
library(org.EcK12.eg.db)
library(celegans.db)
library(GO.db)
library(biomaRt)

library(ballgown)
library(edgeR)

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




# mymerge<-function(all.results,results) {
#   if (dim(all.results)[[1]]==0) {
#     all.results<-results
#   } else {
#     all.results<-merge(all.results,results,all.x=TRUE,all.y=TRUE)
#   }
#   return(all.results)
# }



#Order by expression in transcripts select data for each dataset

cwd<-"~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/"
keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Celegans/KEGG_pathways/'
setwd(cwd)
#load(".RData")

odir<-'Results_1thrs_newannot'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





pheno_data = read.csv("Pheno_data.csv")




cel = ballgown(dataDir = "~/Dropbox/Projects/Metformin_RNAseq/Ballgown_Metformin", samplePattern = "", pData=pheno_data)


#Get counts
full_table <- texpr(cel, 'all')
rownames(full_table)<-full_table$t_name


head(full_table)

cnts.g<-read.table('Counts/BN_EN_WB235_gene_count_matrix.csv',sep=',',
                   header = TRUE,row.names = 1)
cnts.t<-read.table('Counts/BN_EN_WB235_transcript_count_matrix.csv',sep=',',
                   header = TRUE,row.names = 1)

coverage<-colnames(full_table)[grep("cov.", colnames(full_table))]
FPKM<-colnames(full_table)[grep("FPKM.", colnames(full_table))]

dim(cnts.g)
dim(cnts.t)
dim(full_table)




#Ensembl annotation
celegans87 = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl",version=87)
listDatasets(celegans87)
listFilters(celegans87)

listAttributes(celegans87)

celegans.annotation <- getBM(attributes=c('wormbase_gene_seq_name','ensembl_transcript_id','ensembl_gene_id','external_gene_name','entrezgene','description'),
                             mart = celegans87)
#
head(celegans.annotation)
dim(celegans.annotation)
#59286

#dupl.entrez<-duplicated(celegans.annotation$entrezgene)
#celegans.annotation<-celegans.annotation[!dupl.entrez,]



annotation.j<-merge(full_table,cnts.t,by.x='t_name',by.y='row.names',all.x=TRUE)



head(annotation.j)
dim(annotation.j)

annotation.j$ensembl_transcript_id<-ifelse(grepl('transcript',annotation.j$t_name), gsub('transcript:','',annotation.j$t_name),NA  ) 
subset(annotation.j,ensembl_transcript_id=='R10E11.2')


annotation<-merge(annotation.j,celegans.annotation,by='ensembl_transcript_id',all.x=TRUE,all.y=FALSE)


dim(annotation)

tdupl<-annotation[duplicated(annotation$t_name),'t_name']

subset(annotation,t_name %in% tdupl)



subset(celegans.annotation,ensembl_transcript_id=='R10E11.2')
subset(annotation,ensembl_transcript_id=='R10E11.2')



#trans.counts<-ddply(annotation,.(gene_id),summarise,Transcripts=length(t_id))
#annotation$Transcripts<-trans.counts$Transcripts[match(annotation$gene_id,trans.counts$gene_id)]



ancols<-setdiff(colnames(annotation),c(as.character(pheno_data$ids),coverage,FPKM) )


annot.ordered<-annotation[order(annotation$Transcripts,annotation$gene_id,decreasing = TRUE),]
head(subset(annot.ordered,gene_name!='.'))


get.expression<-function(data,annotations,samples) {
  print('Counting transcripts...')
  trans.counts<-ddply(data,.(gene_id),summarise,Transcripts=length(t_id))
  data$Transcripts<-trans.counts$Transcripts[match(data$gene_id,trans.counts$gene_id)]
  annotations<-c(annotations,'Transcripts')
  print('Ordering data by expression level...')
  #Orde by normalised expression level
  
  #data<-data[order(apply(data[,samples],1,sum)/data[,'length'],decreasing = TRUE),c(annotations,samples)]
  
  coverage<-paste('cov',samples,sep='.')
  data<-data[order(apply(data[,coverage],1,sum),decreasing = TRUE),c(annotations,samples)]
  
  
  #Remove non-expressed transcripts
  expressed<-data[apply(data[,samples],1,sum)>0,]
  #Remove possible duplicate annotation
  expressed<-expressed[!duplicated(expressed$t_id),]
  
  print('Counting expressed transcripts...')
  #Get number of expressed transcripts
  extrans.counts<-ddply(expressed,.(gene_id),summarise,ExpressedTranscripts=length(t_id),UniqueGene=length( setdiff(unique(as.character(ensembl_gene_id) ),c('NA') ))==1 )
  
  #Remove annotations for low expression to take annotation for greatest expression at gene level
  expressed.gene<-expressed[!duplicated(expressed$gene_id),]
  
  
  #Put gene level data
  print('Estimating read counts per gene...')
  cnts<-ddply(expressed,"gene_id",function(x) apply(x[,samples],2,sum))
  expressed.gene[,samples]<-cnts[ match(expressed.gene$gene_id,cnts$gene_id),samples]
  
  expressed.gene[,c('ExpressedTranscripts','UniqueGene')]<-extrans.counts[match(expressed.gene$gene_id,extrans.counts$gene_id),c('ExpressedTranscripts','UniqueGene')]
  expressed[,c('ExpressedTranscripts','UniqueGene')]<-extrans.counts[match(expressed$gene_id,extrans.counts$gene_id),c('ExpressedTranscripts','UniqueGene')]
  data[,c('ExpressedTranscripts','UniqueGene')]<-extrans.counts[match(data$gene_id,extrans.counts$gene_id),c('ExpressedTranscripts','UniqueGene')]
  
  
  #Add info about number of expressed transcripts
  print('Finalizing...')
  
  annotations<-c(annotations,'ExpressedTranscripts','UniqueGene')
  
  expressed.gene<-expressed.gene[,c(annotations,samples)]
  expressed<-expressed[,c(annotations,samples)]
  data<-data[,c(annotations,samples)]
  return(list('Transcript'=expressed,'Gene'=expressed.gene,'Annotation'=data) )
}



#Get expression at gene level
expression.all<-get.expression(annotation,ancols,as.character(pheno_data$ids))


expression<-expression.all$Gene
rownames(expression)<-expression$gene_id

head(expression)

subset(expression,gene_id=='MSTRG.16015')
cnts.g['MSTRG.16015',]
#View(expression)



subset(annotation,gene_id=='MSTRG.7505')
subset(expression,gene_id=='MSTRG.7505')
cnts.g['MSTRG.7505',]




#Sanity checks
cnts.g.f<-cnts.g[apply(cnts.g,1,sum)>0,]

dim(cnts.g.f)

cnt.matches<-data.frame(cbind(apply(expression[,as.character(pheno_data$ids)],1,sum),
      apply( cnts.g.f[ match(expression$gene_id,rownames(cnts.g.f)),as.character(pheno_data$ids) ],1,sum)))

colnames(cnt.matches)<-c('Expression','Raw')
cnt.matches$Equal<-cnt.matches$Expression==cnt.matches$Raw

subset(cnt.matches,Equal==FALSE)




write.csv(expression.all$Transcript,paste(odir,"/Raw_data_for_transcripts.csv",sep=''),row.names=FALSE)
write.csv(expression,paste(odir,"/Raw_data_for_genes.csv",sep=''),row.names=FALSE)









expression<-read.csv('Results_1thrs_newannot/Raw_data_for_genes.csv',sep=',',header = TRUE)


#Split table for edgeR

ancols.2<-setdiff(colnames(expression),c(as.character(pheno_data$ids)) )

gene.info<-expression[,ancols.2]
gene.data<-expression[,as.character(pheno_data$ids)]


head(gene.info)
head(gene.data)


dim(gene.data)




groups<-c('C','C','C','C','M','M','M','M','R','R','R','R','RM','RM','RM','RM')


dt<-DGEList(counts=gene.data,
            genes=gene.info,
            group=factor(pheno_data$Group))


#How to choose right threshold

dim(dt)
apply(dt$counts, 2, sum)
#At least one count per million (cpm) reads in at least 4 consistent samples (one group)

thrs<- 1
consis<- 4


keep<-group.filter(dt,1)$Keep


table(keep)



d <- dt[keep,]
dim(d)

#Consis 4 >0
#22501->14193

#Consis 4 >1
#22501->11184


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
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=10)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#Batch removal
logCPM <- cpm(d, log=TRUE, prior.count=1)
logCPMc <- removeBatchEffect(logCPM,as.character(pheno_data$Batch) )



MDSdata<-plotMDS(logCPMc, col=as.numeric(d$samples$group),plot=TRUE)
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=1)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV_batch-adjusted.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



MDSxy<-data.frame(cbind(MDSdata$x,MDSdata$y))

colnames(MDSxy)<-c('PC1','PC2')


write.csv(MDSxy,paste(odir,"/RNAseq_MDS_data.csv",sep=''),row.names=FALSE)





#Remove numbers from sample IDs
# groups<-gsub('[[:digit:]]+', '', colnames(logCPM))
# annData<-factor(groups)
# ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]
# #terrain.colors
# 
# #Heatmap raw
# heatmap3(as.matrix(logCPM),
#          scale = 'row',
#          ColSideLabs = 'Group',
#          balanceColor=TRUE,
#          ColSideColors=ColSideColors,
#          labRow=FALSE)
# 
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/RNAseq_Heatmap_logCPM.pdf",sep=''),
#              width=9,height=9, useDingbats=FALSE)




strains<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$id)),'Strain'])
metf<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$id)),'Treatment'])


strains.fac<-factor(strains)
metf.fac<-factor(metf)

strains.col<-ifelse(strains.fac=='C','white','black')
metf.col<-ifelse(metf.fac=='Control','green','red')

clab<-cbind(metf.col,strains.col)
colnames(clab)<-c('Metformin, mM','Strain')



#Heatmap adjusted
heatmap3(as.matrix(logCPMc),
         scale = 'row',
         #ColSideLabs = 'Group',
         balanceColor=TRUE,
         ColSideColors=clab,
         labRow=FALSE)


legend('topright',legend=c('OP50','OP50-MR'),fill=c('white','black'), border=TRUE, bty="n",title='Strain')
legend('right',legend=c('0','50'),fill=c('green','red'), border=FALSE, bty="n",title='Metformin, mM')



dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Heatmap_logCPMc_batchadj.pdf",sep=''),
             width=9,height=9, useDingbats=FALSE)







#GLM tests
batch
d$samples

design.mat <- model.matrix(~ 0 + d$samples$group+as.character(pheno_data$Batch)) #
colnames(design.mat) <- c(levels(d$samples$group),'Batch')
design.mat

priorn<-5
priordf<-55
d2<-estimateDisp(d,design.mat,prior.df = priordf)


print(paste('Prior.n=',d2$prior.n,' Residual.df=',d2$prior.df/d2$prior.n,' Prior.df=',d2$prior.df,sep=''))

plotBCV(d2)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Dispersion.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)











#All design matrices


design.mat
fit <- glmFit(d2, design.mat)


# gofsum<-gof(fit,plot=TRUE,adjust='fdr',pcutoff = 0.05)
# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/RNAseq_GOF_qq-plot.pdf",sep=''),
#              width=9,height=6, useDingbats=FALSE)

##Outlier genes
# table(gofsum$outlier)
# 
# gofsum$outlier[gofsum$outlier==TRUE]
# 
# #Significant outliers
# table(gofsum$gof.pvalues<0.05)





contrasts<-list('SM-S'=c(-1,1,0,0,0),
                'RM-R'=c(0,0,-1,1,0),
                'R-S'=c(-1,0,1,0,0),
                'SM-S-(RM-R)'=c(-1,1,1,-1,0),
                'T-C_general'=c(-1,1,-1,1,0),
                'R-S_general'=c(-1,-1,1,1,0),
                'S-R_general'=c(1,1,-1,-1,0))


allcomp<-c('SM-S','RM-R','R-S','SM-S-(RM-R)','T-C_general','R-S_general','S-R_general')



comparisons<-c('SM-S','RM-R','R-S','SM-S-(RM-R)',
               'T-C_general','R-S_general','S-R_general')

oldids<-c('M-C','RM-R','R-C','M-C-(RM-R)',
          'T-C_general','R-C_general','C-R_general')



description<-c('Treatment effect while on sensitive strain',
               'Treatment effect while on resistant strain',
               'Difference between strains in control (R-S)',
               'Difference between strains in treatment\n(R adjusted)',
               'General treatment over control',
               'General resistant over sensitive',
               'General sensitive over resistant')



allcomparisons<-data.frame(comparison=comparisons,
                           ID=oldids,
                           description=description)

allcomparisons


gp.KEGG<-subset(getGeneKEGGLinks('cel',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) & PathwayID!='path:cel01100')
#gp.GO<-subset(getGeneGOLinks('cel',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) )


all.results.None<-data.table()
all.results.1FC<-data.table()
all.results<-data.table()

all.KEGGenrichment<-data.table()
all.GOenrichment<-data.table()


for (cid in 1:nrow(allcomparisons)) {
  comp<-allcomparisons[cid,'comparison']
  oldid<-allcomparisons[cid,'ID']
  desc<-allcomparisons[cid,'description']
  contr<-contrasts[[as.character(comp)]]
  
  print(as.character(comp))
  print(as.character(desc))
  print(contr)
  for (thre in c('None')) { #,'1FC'
    enrcoef<-'logFC'
    print(paste('Threshold:',thre,sep=' '))
    if(thre=='None') {
      lrt <- glmLRT(fit, contrast=contr)
    } else {
      #enrcoef<-'unshrunk.logFC'
      lrt <- glmTreat(fit, contrast=contr,lfc=1) 
    }
    
    print('...Differential gene expression')
    de.glm <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
    print(summary(de.glm))
    de2tags <- rownames(d2)[as.logical(de.glm)]
    cairo_pdf(paste(odir,'/MVA_',comp,'_Threshold-',thre,'.pdf',sep = ''),width=9,height=6)
    plotSmear(lrt, de.tags=de2tags)
    abline(h = c(-1, 1), col = "blue")
    title(paste(desc,comp,sep='\n'))
    dev.off()
    
    result<-as.data.frame(topTags(lrt,n=dim(d2$genes)[1]))
    result$Comparison<-comp
    result$Threshold<-thre
    
    if('LR' %in% colnames(result)){
      result$unshrunk.logFC<-as.double(NA)
    } else {
      result$LR<-as.double(NA)
    }
    
    all.results<-rbind(all.results,data.table(result))

    
    if (summary(de.glm)[3]>10 | summary(de.glm)[1] > 10) {
      print('...KEGG enrichment')
      enrichment.KEGG<-kegga(lrt,coef=enrcoef,geneid = 'entrezgene',species.KEGG='cel',gene.pathway = gp.KEGG,convert=TRUE,FDR=0.05,trend='length')
      
      enrichment.KEGG$ID<-rownames(enrichment.KEGG)
      enrichment.KEGG<-rename(enrichment.KEGG,c("Pathway"="Description"))
      enrichment.KEGG$Comparison<-comp
      enrichment.KEGG$Threshold<-thre
      all.KEGGenrichment<-rbind(all.KEGGenrichment,data.table(enrichment.KEGG))
      
      print('...GO enrichment')
      enrichment.GO<-goana(lrt,species='Ce',coef=enrcoef,geneid='entrezgene',convert=TRUE,FDR=0.05,trend='length')
      
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
}



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


write.csv(all.KEGGenrichment,paste(odir,'/All_results_KEGG.csv',sep = ''),row.names = FALSE)


write.csv(all.GOenrichment,paste(odir,'/All_results_GO.csv',sep = ''),row.names = FALSE)










# 
# cwd<-"~/Projects/2015-Metformin/RNAseq/Celegans_metformin/"
# keggxml<-'~/Projects/2015-Metformin/Annotations/Celegans/KEGG_pathways/'
# setwd(cwd)
# odir<-'Results'
# dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
# 
# 
# 
# all.results.kg<-read.table(paste(odir,'/All_results.csv',sep = ''),header = TRUE,sep=',')
# all.KEGGenrichment<-read.table(paste(odir,'/All_results_KEGG.csv',sep = ''),header = TRUE,sep=',')
# all.GOenrichment<-read.table(paste(odir,'/All_results_GO.csv',sep = ''),header = TRUE,sep=',')


dim(all.results.kg)

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
comp.main<-c('SM-S_logFC','RM-R_logFC','SM-S-(RM-R)_logFC')
comp.gen<-c('T-C_general_logFC','S-R_general_logFC')

#For heatmaps
comp.all<-c('SM-S','RM-R','SM-S-(RM-R)','R-S','T-C_general','S-R_general')
comp.write<-c('SM-S','RM-R','SM-S-(RM-R)','R-S','T-C_general','R-S_general','S-R_general')
comp.mainh<-c('SM-S','RM-R','SM-S-(RM-R)')

comp.grps<-list('All'=comp.all,'Write'=comp.write,'Main'=comp.mainh)

comp.list<-list('Main'=comp.main,'General'=comp.gen)

#grep('SM-S_logFC',colnames(gdataf))

add.path<-c('cel04213','cel04212')
allpaths<-c()

setwd(cwd)


head(all.results.kg)

ancols.3<-c('Comparison','gene_id','ensembl_transcript_id','wormbase_gene_seq_name','ensembl_gene_id','external_gene_name','entrezgene','description','Transcripts','ExpressedTranscripts','UniqueGene',
            'chr','strand','start','end','num_exons','length','Pathways','GO_terms')

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


head(all.results.rm)

all.results.rmc<-subset(all.results.rm,Stat %in% c(enrcoef,'FDR'))

all.results.rcp<-dcast(all.results.rmc,gene_id+ensembl_transcript_id+wormbase_gene_seq_name+ensembl_gene_id+external_gene_name+entrezgene+description+Transcripts+ExpressedTranscripts+UniqueGene+
                         Pathways+GO_terms+chr+strand+start+end+num_exons+length~Comparison+Stat,value.var = 'Value')

head(all.results.rcp)

#Reorder by significance
signsum<-apply(all.results.rcp[,grep('FDR',colnames(all.results.rcp))],1,function(x) sum(-log10(x)))
all.results.rcp<-all.results.rcp[order(signsum,decreasing = TRUE),]

nrow(all.results.rcp)

gdata<-subset(all.results.rcp,! is.na(entrezgene))
nrow(gdata)
#   nrow(gdata)
#   o<-order(rowSums(gdata[,grep('_FDR',colnames(gdata))]),decreasing=FALSE)
#   gdata<-gdata[o,]

#Do selection by gene id, because gene_id is associated with function
dupl<-duplicated(gdata[,mapkey])
print('Entrez Gene duplicates')
print(table(dupl))



#nrow(gdataf)

write.csv(all.results.rcp,paste(odir,'/All_results_sidebyside_Threshold-',thres,'.csv',sep = ''),row.names = FALSE)



#Volcano plots
all.plot<-subset(all.results.r,Comparison %in% c('SM-S','RM-R','SM-S-(RM-R)'))
vol<-ggplot(all.plot,aes(x=logFC,y=-log10(FDR)))+
  geom_hline(yintercept = -log(0.05,10),color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = -1,color=barcolor,alpha=baralpha)+
  geom_vline(xintercept = 1,color=barcolor,alpha=baralpha)+
  ggtitle(paste('Volcano plot for differentially expressed genes. Threshold: ',thres,sep=''),
          subtitle='Labels are shown for some of the significantly affected genes')+
  geom_point(alpha=0.9,size=1)+
  xlim(-10,10)+
  #   geom_label_repel(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp,as.character(gene_name),'')),
  #             size=txtsize,colour = "red")+
  geom_text(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp & FDR<0.05,as.character(gene_name),'')),
            hjust=1, vjust=-0.5,size=txtsize,colour = "red")+
  facet_grid(.~Comparison)
fname<-paste(odir,'/Volcano_Threshold-',thres,'.pdf',sep = '')
print(fname)
cairo_pdf(fname,width=18,height=9)
print(vol)
dev.off()
#Enrichment

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
  enrichmentfe<-subset(enrichmentm,Stat %in% c('P.Up','P.Down') & Comparison %in% comp.write)
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

  
  for (ont in ontologies) {
    for (comp.grh in c('All','Main')) {
      comp.heat<-comp.grps[[comp.grh]]
      if (ent=='GO') {
        enrichmentco<-subset(enrichmentc,Ont==ont)
        fname<-paste(odir,'/Enrichment_heatmap_',ent,'-',ont,'_Threshold-',thres,'_',comp.grh,'.pdf',sep = '')
      } else {
        enrichmentco<-enrichmentc
        fname<-paste(odir,'/Enrichment_heatmap_',ent,'_Threshold-',thres,'_',comp.grh,'.pdf',sep = '')
      }
      
      #hdata<-enrichmentco[,!colnames(enrichmentco) %in% c('ID','Description','N','Ont')]
      #Filter just selected columns
      hdata<-enrichmentco[,grep(paste(comp.heat,collapse="|"),colnames(enrichmentco))]
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
comp.gr<-'Main'
pathcomp<-comp.list[[comp.gr]]
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
    #all.kegg.mappingst<-mymerge(all.kegg.mappingst,result)
    all.kegg.mappingst<-bind_rows(all.kegg.mappingst,result)
  }
}


all.kegg.mappingst$mol.col<-NULL

all.kegg.mappings<-merge(path2pathde,all.kegg.mappingst,by='PathwayID',all.y=TRUE)

write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings.csv',sep = ''),row.names = FALSE)







all.results.rcp<-read.csv(paste(odir,'/All_results_sidebyside_Threshold-None.csv',sep = ''),check.names = FALSE)





head(all.results.rcp)



#Venn
library('Vennerable')
#Might need to deactive reshape and Vennerable to use melt in other code






SM_S.up<-subset(all.results.rcp,`SM-S_FDR`<0.05 & `SM-S_logFC`>0)[,'gene_id']
SM_S.down<-subset(all.results.rcp,`SM-S_FDR`<0.05 & `SM-S_logFC`<0)[,'gene_id']
RM_R.up<-subset(all.results.rcp,`RM-R_FDR`<0.05 & `RM-R_logFC`>0)[,'gene_id']
RM_R.down<-subset(all.results.rcp,`RM-R_FDR`<0.05 & `RM-R_logFC`<0)[,'gene_id']

venn.groups<-list('RM-R.down'=RM_R.down,
                  'SM-S.down'=SM_S.down,
                  'RM-R.up'=RM_R.up,
                  'SM-S.up'=SM_S.up)

ups<-intersect(SM_S.down,RM_R.up)
downs<-intersect(SM_S.up,RM_R.down)

write.csv(subset(all.results.rcp,gene_id %in% downs),paste(odir,'/Flipping_downregulated.csv',sep = ''),row.names = FALSE)
write.csv(subset(all.results.rcp,gene_id %in% ups),paste(odir,'/Flipping_upregulated.csv',sep = ''),row.names = FALSE)



S.change<-subset(all.results.rcp,`SM-S_FDR`<0.05 )[,'gene_id']
R.change<-subset(all.results.rcp,`RM-R_FDR`<0.05 )[,'gene_id']
L.change<-subset(all.results.rcp,`SM-S-(RM-R)_FDR`<0.05 )[,'gene_id']

venn.change<-list('SM-S'=S.change,
            'RM-R'=R.change)




plot(Venn(venn.groups),show = list(Faces = FALSE),
     doWeights = FALSE,type='ellipses')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram.pdf",sep=''),width=6,height=6)

plot(Venn(venn.change),show = list(Faces = FALSE),
     doWeights = FALSE)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_general_change.pdf",sep=''),width=6,height=6)


venn.longevity<-list('SM-S'=S.change,
                  'RM-R'=R.change,
                  'Longevity'=L.change)

plot(Venn(venn.longevity),show = list(Faces = FALSE),
     doWeights = TRUE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_longevity.pdf",sep=''),width=6,height=6)



detach("package:Vennerable", unload=TRUE)
detach("package:reshape", unload=TRUE)
detach("package:reshape2", unload=TRUE)

library(reshape2)





#Comparison scatter





fit<-lm(`RM-R_logFC`~`SM-S_logFC`,data=all.results.rcp)
res<-summary(fit)
res
res$r.squared


dim(subset(all.results.rcp,sqrt( (`SM-S_logFC`*1.5)^2+(`RM-R_logFC`)^2 ) >5 &
             (`SM-S_FDR`<0.05 |
             `RM-R_FDR`<0.05) ))



library(ellipse)


ellipsoid<-as.data.frame(ellipse( 0.45,
                                  scale=c(1,2),
                                  centre=c( 0,0) ))

ggplot()+geom_path(data=ellipsoid, aes(x=x, y=y),size=1)



all.results.rcp$Change<-'None'
all.results.rcp$Change<-ifelse(all.results.rcp$`SM-S_FDR`<0.05,'Sensitive',all.results.rcp$Change)
all.results.rcp$Change<-ifelse(all.results.rcp$`RM-R_FDR`<0.05,'Resistant',all.results.rcp$Change)
all.results.rcp$Change<-ifelse(all.results.rcp$`RM-R_FDR`<0.05 & all.results.rcp$`SM-S_FDR`<0.05,'Both',all.results.rcp$Change)

all.results.rcp$Change<-factor(all.results.rcp$Change,levels=c('None','Sensitive','Resistant','Both'))




allS<-subset(all.results.rcp,Change=='Sensitive')
allR<-subset(all.results.rcp,Change=='Resistant')
allB<-subset(all.results.rcp,Change=='Both')

topS<-allS[order(allS$`SM-S_FDR`) ,'t_id'][1:10]
topR<-allR[order(allR$`RM-R_FDR`) ,'t_id'][1:10]
topB<-allB[order(allB$`SM-S_FDR`,allB$`RM-R_FDR`) ,'t_id'][1:10]


topAll<-c(51908)#topS,topR,topB,topR,



all.results.rcp[order(abs(all.results.rcp$`SM-S-(RM-R)_logFC`),decreasing = TRUE),]


amp<-5
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue2','blue2','gray80','red2','red2')



subset(all.results.rcp,gene_name=='acs-2')

ggplot(all.results.rcp,aes(x=`SM-S_logFC`,y=`RM-R_logFC`))+
  geom_vline(xintercept = 0,alpha=0.2)+
  geom_hline(yintercept = 0,alpha=0.2)+
  geom_abline(slope=1,intercept = 0,alpha=0.2)+
  geom_abline(slope = fit$coefficients[[2]],intercept = fit$coefficients[[1]],color='red',alpha=0.5)+
  geom_point(data=all.results.rcp,aes(color=`SM-S-(RM-R)_logFC` ),stroke=0,size=1)+
  labs(color='Significant change')+
  xlab('Metformin effect when on OP50, logFC')+
  ylab('Metformin effect when on OP50MR, logFC')+
  scale_colour_gradientn(colors = gradcols,limits=c(-amp,amp),name='Longevity effect')+
  #scale_size(range=c(0,8),breaks=seq(0,12,by=2),name='Amplitude of change')+
  scale_x_continuous(breaks=seq(-12,12,by=2),limits=c(-12,10))+
  scale_y_continuous(breaks=seq(-10,10,by=2),limits=c(-10,10))+
  #geom_text_repel(aes(label=ifelse( t_id %in% topAll,as.character(gene_name),'') ))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_SensRes_longevity_new.pdf",sep=''),
             width=7,height=4,useDingbats=FALSE)






all<-1529+754+953

(1529+953)/all

13/all

head(all.results.rcp)




#Read REVIGO results
revigo.all<-data.frame()

for (fl in list.files(paste(odir,'REVIGO',sep='/')) ) {
  print(fl)
  if ( grepl('_Up',fl) ) {
    tp<-'Up'
  } else {
    tp<-'Down'
  }
  contrast.t<-gsub(paste0('_',tp,'.csv'),'',fl )
  contrast<-gsub('REVIGO-','',contrast.t )
  
  revigo<-read.table(paste(odir,'REVIGO',fl,sep='/'),sep=',',header=TRUE)
  revigo$Type<-tp
  revigo$Contrast<-contrast
  revigo.all<-rbind(revigo.all,revigo)
}

revigo.all$FDR<-10^revigo.all$log10.p.value
revigo.all$logFDR<- -revigo.all$log10.p.value

contrasts<-unique(revigo.all$Contrast)

contrasts.ord<-c('SM-S','RM-R','SM-S-(RM-R)','T-C_general')

revigo.all$Contrast<-factor(revigo.all$Contrast,levels=contrasts.ord,labels=contrasts.ord)
revigo.all$Type<-factor(revigo.all$Type,levels=c('Up','Down'),labels=c('Up','Down'))

write.csv(revigo.all,paste(odir,'REVIGO_all.csv',sep = '/'))




head(revigo.all)


revigo.cast<-dcast(subset(revigo.all, plot_X!='null'),term_ID+description~Contrast+Type,value.var='FDR')
write.csv(revigo.cast,paste(odir,'REVIGO_side-by-side.csv',sep = '/'))



reorderfun_mean = function(d,w) { reorder(d, w, agglo.FUN = mean) }
gyrs<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = 7)

enbrks<-c(0,-log(0.05,10),2,3,4,5,10,100)


hdata<-dcast(subset(revigo.all,plot_X!='null' & Contrast!='T-C_general'),term_ID+description~Contrast+Type,value.var='logFDR') #


encols<-setdiff(colnames(hdata),c('term_ID','description','Index'))
encols<-encols[grepl('SM-S',encols,fixed = TRUE)]

head(hdata)
dim(hdata)

#hdata<-hdata[!apply(hdata, 1, function(x) {any(as.numeric(x) < -log(0.0001,10),na.rm=TRUE)}),]

hdata<-hdata[apply(hdata[,encols], 1, function(x) {any(as.numeric(x) > 4,na.rm=TRUE)}),]


#50 terms
dim(hdata)


sel_terms<-as.character(hdata$term_ID)

  

hdata2.r<-dcast(subset(revigo.all,Contrast!='T-C_general' & term_ID %in% sel_terms),term_ID+description~Contrast+Type,value.var='logFDR') #

hdata2.r$Index<-paste(hdata2.r$term_ID,hdata2.r$description,sep=' - ')
rownames(hdata2.r)<-hdata2.r$Index


write.csv(hdata2.r,paste0(odir,'/Enrichment_heatmap_REVIGO_GO.csv'))


hdata2<-hdata2.r[,encols]

head(hdata2)


#Create a dummy dataset where missing values are filled with 0 and Infinite values (p=0, log(0)=-Inf) with large numbers
#This is needed for the clustering/rearangement algorithm to work
hdatafill<-hdata2
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
heatmap.2(data.matrix(hdata2),
          key=FALSE, #Colour key needs to be provided separately
          Colv=FALSE, #Columns should not be reordered
          Rowv=hmap$rowDendrogram, # Here comes in the ordering
          trace='none',
          col=gyrs, # Colour scale
          breaks=enbrks, # Colour breaks
          xlab='Comparison',
          dendrogram='row', #Row dendogram, but should be changed to none, as dendrogram represents our data with filled-in values
          scale="none", #Should values be normalised in rows or columns - No
          na.color="white", # What colour to use with not missing values
          symkey=FALSE, #Provided colour scale is not symetrical
          cexRow=0.7, #Some figure scaling parameters. Works only after a lot of experimentation
          cexCol=0.7,
          margin=c(10,20),
          lwid=c(0.2,0.8),
          lhei=c(0.05,0.95))



dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Enrichment_heatmap_REVIGO_high_significance_summary_logp_4.pdf",sep=''),
             width=7,height=15,useDingbats=FALSE)




