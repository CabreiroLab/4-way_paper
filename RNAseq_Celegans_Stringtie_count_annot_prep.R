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




#devtools::install_github("PNorvaisas/PFun")
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



cwd<-"~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/"
keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Celegans/KEGG_pathways/'
setwd(cwd)
#load(".RData")

odir<-'Results_1thrs_newannot'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


pheno_data = read.csv("Pheno_data.csv")


#get Ballgown output files
cel = ballgown(dataDir = "~/Dropbox/Projects/Metformin_RNAseq/Ballgown_Metformin", samplePattern = "", pData=pheno_data)
#This is necessary for the annotation of transcripts

#Get counts
full_table <- texpr(cel, 'all')
rownames(full_table)<-full_table$t_name


head(full_table)

#Get read counts
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

#This is for cleanup of Stringtie output
#Stringtie predicts expression of unannotated transcripts, which is interesting, but complicates the analysis at gene level
#In the future I would reccomend the use of HTSeq
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



