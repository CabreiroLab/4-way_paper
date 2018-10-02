options(java.parameters = "-Xmx4096m") 

#That's for annotations
library(org.Ce.eg.db)
library(org.EcK12.eg.db)
library(celegans.db)
library(GO.db)


library(ballgown)
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


library(biomaRt)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



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

length(rownames(cnts.g))
length(rownames(cnts.t))
length(rownames(full_table))

cnts.ga<-merge(cnts.g,full_table,by.x='row.names',by.y='gene_id',all.x=TRUE)
cnts.ga<-rename(cnts.ga,c('Row.names'='gene_id'))


head(cnts.g)

table(grepl('transcript:',cnts.ga$gene_id))
table(grepl('transcript:',rownames(cnts.g)))



cnts.ga$ENSEMBL_transcript<-ifelse(grepl('transcript',cnts.ga$t_name), gsub('transcript:','',cnts.ga$t_name),NA  )
cnts.ga$ENSEMBL<-ifelse(grepl('WBGene',cnts.ga$gene_id),gsub('gene:','',cnts.ga$gene_id),NA)

table(!is.na(cnts.ga$ENSEMBL))
table(!is.na(cnts.ga$ENSEMBL_transcript))

head(cnts.ga)

table(grepl('transcript:',cnts.ga$ENSEMBL))
table(grepl('transcript:',cnts.ga$ENSEMBL))


#Fix ENSEMBL_transcript ID for gene level
#This part is not really releveant anymore, because annotations are later acquired using BioMart
#There should be a better way without a loop
prcl<-0
for (ri in 1:nrow(cnts.ga)){
  ep<-unlist(strsplit(cnts.ga[ri,'ENSEMBL_transcript'],paste("([\\.])",sep=''), perl=TRUE))
  epc<-ep
  if (length(epc)>1 & !is.na(epc)) {
    if (grepl('t',epc[2])!=1 ) {
      x <- gregexpr("[0-9]+", epc[2])
      epc[2] <- unlist(regmatches(epc[2], x))
    }
    cnts.ga[ri,'ENSEMBL_gene']<-paste(epc[1],epc[2],sep='.')
  } else {
    cnts.ga[ri,'ENSEMBL_gene']<-NA
  }
  prc<-ri*100/nrow(cnts.ga)
  if (prc-prcl>=5){
    prcl<-prc
    print(paste(round(prc),'%',sep=''))
  }
}


cnts.ga$entrezgene<-NA

#By ENSEMBL transcript/protein IDs
egENSEMBLP<-toTable(org.Ce.egENSEMBLPROT2EG)
head(egENSEMBLP)
dim(egENSEMBLP)

mENP<-match(cnts.ga$ENSEMBL_transcript,egENSEMBLP$prot_id)
table(!is.na(mENP))
cnts.ga$entrezgene<-ifelse(is.na(cnts.ga$entrezgene),egENSEMBLP$gene_id[mENP],cnts.ga$entrezgene)
table(!is.na(cnts.ga$entrezgene))
#31305

egENSEMBL<-toTable(org.Ce.egENSEMBL2EG)
head(egENSEMBL)
dim(egENSEMBL)
#ENSEMBL to entrezgene
dim(subset(egENSEMBL,is.na(gene_id) | is.na(ensembl_id)))
me<-match(cnts.ga$ENSEMBL,egENSEMBL$ensembl_id)
table(!is.na(me))
cnts.ga$entrezgene<-ifelse(is.na(cnts.ga$entrezgene),egENSEMBL$gene_id[me],cnts.ga$entrezgene)
table(!is.na(cnts.ga$entrezgene))
#31434

egWB<-toTable(org.Ce.egWORMBASE)
head(egWB)
dim(egWB)
#Wormbase to entrezgene
mWB<-match(cnts.ga$ENSEMBL,egWB$wormbase_id)
table(!is.na(mWB))
cnts.ga$entrezgene<-ifelse(is.na(cnts.ga$entrezgene),egWB$gene_id[mWB],cnts.ga$entrezgene)
table(!is.na(cnts.ga$entrezgene))
#31469


egALIAS<-toTable(org.Ce.egALIAS2EG)
head(egALIAS)
dim(egALIAS)

#Missing ALIAS to entrezgene
dim(subset(egALIAS,is.na(gene_id) | is.na(alias_symbol)))
mAL<-match(cnts.ga$gene_name,egALIAS$alias_symbol)
table(!is.na(mAL))
cnts.ga$entrezgene<-ifelse(is.na(cnts.ga$entrezgene),egALIAS$gene_id[mAL],cnts.ga$entrezgene)
table(!is.na(cnts.ga$entrezgene))
#58013




subset(cnts.ga,is.na(entrezgene))
table(!is.na(cnts.ga$entrezgene))

#Based on entrezgene
table(!is.na(cnts.ga$ENSEMBL))

#Fill ENSEMBL ids by Entrez Gene
# mENS<-match(cnts.ga$entrezgene,egENSEMBL$gene_id)
# table(!is.na(mENS))
# cnts.ga$ENSEMBL<-ifelse(is.na(cnts.ga$ENSEMBL),egENSEMBL$ensembl_id[mENS],cnts.ga$ENSEMBL)
# table(!is.na(cnts.ga$ENSEMBL))
#
#
# #Fill ENSEMBL ids by Entrez Gene using WB
# mENS_WB<-match(cnts.ga$entrezgene,egWB$gene_id)
# table(!is.na(mENS_WB))
# cnts.ga$ENSEMBL<-ifelse(is.na(cnts.ga$ENSEMBL),egWB$wormbase_id[mENS_WB],cnts.ga$ENSEMBL)
# table(!is.na(cnts.ga$ENSEMBL))



# egSYM2EG<-toTable(org.Ce.egSYMBOL2EG)
# m<-match(cnts.gaf$gene_name,egSYM2EG$symbol)
# cnts.gaf$entrezgene<-egSYM2EG$gene_id[m]

# egREFSEQ<-toTable(org.Ce.egREFSEQ2EG)
# m<-match(cnts.gaf$entrezgene,egREFSEQ$gene_id)
# cnts.gaf$RefSeqID<-egREFSEQ$accession[m]

# egGENEDESC<-toTable(org.Ce.egGENENAME)
# mgd<-match(cnts.gaf$entrezgene,egGENEDESC$gene_id)
# cnts.gaf$Description<-egGENEDESC$gene_name[mgd]
#
# msym<-match(cnts.gaf$entrezgene,egSYMBOL$gene_id)
# cnts.gaf$Symbol<-egSYMBOL$symbol[msym]
#




length(unique(cnts.ga$gene_id))
length(unique(cnts.ga$t_id))

length(unique(cnts.ga$ENSEMBL))
length(unique(cnts.ga$ENSEMBL))


codes<-c('C1','C2','C3','C4',
         'M1','M2','M3','M4',
         'R1','R2','R3','R4',
         'RM1','RM2','RM3','RM4')


#Remove non-expressed transcripts
dim(cnts.ga)
expressed<-rowSums(cnts.ga[, coverage])>0
cnts.gae<-cnts.ga[expressed,]
dim(cnts.gae)
#62270->36128


#get all gene names
# allgenes<-data.frame(table(cnts.gae$gene_id))
#
# #Find duplicated genes - multiple transcripts per gene
# Duplicated.g<-subset(allgenes,Freq>1)$Var1
# cnts.ag<-subset(cnts.gae,gene_id %in% Duplicated.g)
# dim(cnts.ag)
#19709 has duplicates

#Mark genes which have duplicates
# unq<-ddply(cnts.gae,.(gene_id),summarise,Unique=length(unique(ENSEMBL))==1)
# table(unq$Unique)
#22501 TRUE - All unique by ENSEMBL

# unq<-ddply(cnts.gae,.(gene_id),summarise,Unique=length(unique(entrezgene))==1)
# table(unq$Unique)
#2422 20079

unq<-ddply(cnts.gae,.(gene_id),summarise,Unique=length(unique(gene_name))==1)
table(unq$Unique)
#2445 20056


head(cnts.gae)

ma<-match(cnts.gae$gene_id,unq$gene_id)
cnts.gae$Unique<-unq$Unique[ma]
table(cnts.gae$Unique)
#26716 unique gene names
#9412 has dupicates or more gene names
#ifelse(cnts.gae$gene_id %in% Duplicated.g,TRUE,FALSE)

#Unique checks for gene regions with multiple possible functions before filtering


#Expressed but in different levels, select most highly expressed
#Select gene annotation for one of duplicates which is most expressed
nrow(cnts.gae)
o<-order(rowSums(cnts.gae[,c(coverage)]),decreasing=TRUE)
cnts.gae<-cnts.gae[o,]




#Do selection by gene id, because gene_id is unique in Ballgown gene summary
dupl<-duplicated(cnts.gae$gene_id)
table(dupl)

#Filter out duplicates
cnts.gaf<-cnts.gae[!dupl,]
nrow(cnts.gaf)

#Top expressed
head(cnts.gaf,n=100)
table(is.na(cnts.gaf$ENSEMBL))
table(is.na(cnts.gaf$entrezgene))

#36128->22501


head(cnts.gaf)

cnts.gaf


#Missing annotations
missing<-subset(cnts.gaf,is.na(entrezgene) & gene_name!='.')[,c('t_id','t_name','gene_id','gene_name','entrezgene','ENSEMBL')]
dim(missing)
#231
#missing

dots<-subset(cnts.gaf,is.na(entrezgene) & gene_name=='.')[,c('t_id','t_name','gene_id','gene_name','entrezgene','ENSEMBL')]
dim(dots)
#1023

dots


cols.cov<-colnames(cnts.gaf)[grep('cov',colnames(cnts.gaf))]
cols.FPKM<-colnames(cnts.gaf)[grep('FPKM',colnames(cnts.gaf))]
colorder<-c('t_name','t_id','entrezgene','ENSEMBL','ENSEMBL_gene','ENSEMBL_transcript','gene_id','gene_name',
            'chr','strand','start','end','length','num_exons',
            codes,cols.cov,cols.FPKM)


cnts.gaf.w<-cnts.gaf[,colorder]
write.csv(cnts.gaf.w,paste(odir,"/Raw_data_for_genes.csv",sep=''),row.names=FALSE)



cnts.gae.w<-cnts.gae[,setdiff(colorder,c('entrezgene'))]
write.csv(cnts.gae.w,paste(odir,"/Raw_data_for_genes_Unfiltered.csv",sep=''),row.names=FALSE)
