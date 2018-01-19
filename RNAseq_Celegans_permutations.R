options(java.parameters = "-Xmx4096m") 

#That's for annotations
library(org.Ce.eg.db)
library(org.EcK12.eg.db)
#library(celegans.db)
library(GO.db)


#library(ballgown)
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

library(gtools)

theme_set(theme_light())

Statprep<-function(df) {
  dfn= arrange(df,pval)
  dfn$logFC<-log2(dfn$fc)
  dfn$logqval<--log10(dfn$qval)
  dfn$logpval<--log10(dfn$pval)
  return(dfn)
}



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
odir<-'Results_1thrs_Bootstraping'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



expression<-read.csv("Results_1thrs_newannot/Raw_data_for_genes.csv")
rownames(expression)<-expression$gene_id

head(expression)



pheno_data = read.csv("Pheno_data.csv")
rownames(pheno_data)<-pheno_data$ids
codes<-as.character(pheno_data$ids)
ugroups<-as.character(unique(pheno_data$Group))

ancols.2<-setdiff(colnames(expression),c(as.character(pheno_data$ids)) )



gene.info<-expression[,ancols.2]
gene.data<-expression[,as.character(pheno_data$ids)]





# codes<-c('C1','C2','C3','C4',
#          'M1','M2','M3','M4',
#          'R1','R2','R3','R4',
#          'RM1','RM2','RM3','RM4')



# groups<-c('C','C','C','C','M','M','M','M','R','R','R','R','RM','RM','RM','RM')

# 
# batch.lvl<-c('1','1','2','2',
#              '1','2','2','2',
#              '1','2','2','2',
#              '1','1','2','2')
# sample.info<-data.frame(Codes=codes,Groups=groups,Batches=batch.lvl)
# rownames(sample.info)<-sample.info$Codes
# 
# sample.info


#pheno_data = read.csv("Pheno_data.csv")


#Set comparisons
contrasts<-list('SM-S'=c(-1,1,0,0,0),
                'RM-R'=c(0,0,-1,1,0),
                'R-S'=c(-1,0,1,0,0),
                'RM-SM'=c(0,-1,0,1,0),
                'SM-R'=c(0,1,-1,0,0),
                'RM-S'=c(-1,0,0,1,0))



comparisons<-c('SM-S','RM-R','R-S','RM-SM','SM-R','RM-S')

oldids<-c('M-C','RM-R','R-C','RM-M','M-R','RM-C')

description<-c('Treatment effect while on sensitive strain',
               'Treatment effect while on resistant strain',
               'Difference between strains in control (R-S)',
               'Difference between strains in treatment (RM-SM)',
               'Diagonal difference (SM-R)',
               'Diagonal difference (RM-S)')



allcomparisons<-data.frame(comparison=comparisons,
                           ID=oldids,
                           description=description)





#t(combn(1:16,2))

#t(combn(4,2,simplify = TRUE))



#All unique sample pairs
pairs.t<-t(combn(ugroups,2))


pairs<-rbind(c('None','None'),pairs.t)

pairs
samp.rem<-expand.grid(1:4, 1:4)


#How many removals?
rm.no<-5

#Expression threshold
thrs<-1

#All design matrices
#design.mat

all.results<-data.table()

for (p in 1:nrow(pairs) ){
  pair<-pairs[p,]
  print(pair)
  
  if (all(pair=='None')) {
    rm.noi<-1
  } else {
    rm.noi<-rm.no
  }
  
  rnd.smp1<-sample(1:16,rm.noi,replace = FALSE)
  rnd.smp2<-sample(1:9,rm.noi,replace = FALSE)
  
  rnd.smp<-cbind(rnd.smp1,rnd.smp2)
  
  for (sri in 1:nrow(rnd.smp) ) {
    sri1<-rnd.smp[sri,1]
    sri2<-rnd.smp[sri,2]
    
    sp1<-samp.rem[sri1,]
    
    samp.rem2<-expand.grid(setdiff(1:4,sp1[,1]),setdiff(1:4,sp1[,2]))
    
    sp2<-samp.rem2[sri2,]
    
    #rmsp<-paste(pair,sp1,sep='')
    if (all(pair=='None')) {
      rmsp<-pair
    } else {
      rmsp<-c(paste(pair,sp1,sep=''),paste(pair,sp2,sep=''))
    }
    
    rmsp.string<-paste(rmsp,collapse="_")
    print(rmsp)
    
    sample.info.i<-pheno_data[setdiff(codes,rmsp),]
    
    #print(sample.info.i)
    
    group.freq<-data.frame(table(sample.info.i$Group))
    rownames(group.freq)<-group.freq$Var1
    codes.i<-as.character(sample.info.i$ids)
    
    gene.data<-expression[,codes.i]
    
    dt<-DGEList(counts=gene.data,
                genes=gene.info,
                group=as.character(sample.info.i$Group) )
    
    
    keep<- rowSums(cpm(dt)[,grep('^C[[:digit:]]',colnames(dt))]>thrs)>= group.freq['C','Freq'] |
      rowSums(cpm(dt)[,grep('^M[[:digit:]]',colnames(dt))]>thrs)>= group.freq['M','Freq'] |
      rowSums(cpm(dt)[,grep('^R[[:digit:]]',colnames(dt))]>thrs)>= group.freq['R','Freq'] | 
      rowSums(cpm(dt)[,grep('^RM[[:digit:]]',colnames(dt))]>thrs)>= group.freq['RM','Freq']
    
    d <- dt[keep,]
    d$samples$lib.size <- colSums(d$counts)
    d <- calcNormFactors(d)
    

    #GLM tests
    design.mat <- model.matrix(~ 0 + d$samples$group+as.character(sample.info.i$Batch)) #
    colnames(design.mat) <- c(levels(d$samples$group),'Batch')
    
    priorn<-5
    priordf<-55
    d2<-estimateDisp(d,design.mat,prior.df = priordf)
    

    fname<-paste(odir,"/RNAseq_Dispersion_",rmsp.string,".pdf",sep='')
    print(fname)
    cairo_pdf(fname,width=9,height=6)
    plotBCV(d2)
    dev.off()
    
    maxgenes<-dim(d2$genes)[1]
    
    fit <- glmFit(d2, design.mat)
    
    
    for (cid in 1:nrow(allcomparisons)) {
      comp<-allcomparisons[cid,'comparison']
      oldid<-allcomparisons[cid,'ID']
      desc<-allcomparisons[cid,'description']
      contr<-contrasts[[as.character(comp)]]
      
      print(as.character(comp))
      print(as.character(desc))
      print(contr)
      
      thre<-'None'
      #for (thre in c('None')) {
      #Do just "None" Threshold
      #,'1FC'
      enrcoef<-'logFC'
      print(paste('Threshold:',thre,sep=' '))
      lrt <- glmLRT(fit, contrast=contr)

      
      print('...Differential gene expression')
      de.glm <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
      print(summary(de.glm))
      # de2tags <- rownames(d2)[as.logical(de.glm)]
      # cairo_pdf(paste(odir,'/MVA_',comp,'_Threshold-',thre,'_',rmsp.string,'.pdf',sep = ''),width=9,height=6)
      # plotSmear(lrt, de.tags=de2tags)
      # abline(h = c(-1, 1), col = "blue")
      # title(paste(desc,comp,sep='\n'))
      # dev.off()
      
      result<-as.data.frame(topTags(lrt,n=maxgenes))
      result$Comparison<-comp
      result$Threshold<-thre
      
      if('LR' %in% colnames(result)){
        result$unshrunk.logFC<-as.double(NA)
      } else {
        result$LR<-as.double(NA)
      }
      
      result.dt<-data.table(result)
      result.dt$PairRM<-paste(pair,collapse='_')
      result.dt$SamplesRM<-rmsp.string
      
      all.results<-rbind(all.results,result.dt)
      
      print('---------')
    
    }
    
  }
}



all.results$Comparison<-factor(all.results$Comparison,
                               levels=comparisons,
                               labels=comparisons)

#

all.results<-data.frame(all.results)
dim(all.results)


head(all.results)



all.results.r<-subset(all.results,Threshold=='None')
all.results.r$Threshold<-NULL
all.results.r$unshrunk.logFC<-NULL

head(all.results.r)




ancols.ex<-setdiff(colnames(all.results.r), ancols.2)

all.results.rmc<-melt(all.results.r,id.vars = c('Comparison',ancols.2,'PairRM','SamplesRM'),
                     measure.vars=c('logFC','FDR'),
                     variable.name = 'Stat',value.name = 'Value')


head(all.results.rmc)


#For metabolic modeling
#Choose only the necessary comparisons
all.results.model<-subset(all.results.rmc,Comparison %in% c('SM-S','RM-R','R-S','RM-SM','RM-S','SM-R') & !is.na(external_gene_name))

all.results.model$Comparison<-as.character(all.results.model$Comparison)
all.results.model$Stat<-as.character(all.results.model$Stat)


all.results.model$Comparison<-gsub('-','x',all.results.model$Comparison)
all.results.model$Comparison<-revalue(all.results.model$Comparison,c('SMxS'='SxSM',
                                                                     'RMxR'='RxRM',
                                                                     'RxS'='SxR',
                                                                     'RMxSM'='SMxRM',
                                                                     'RMxS'='SxRM',
                                                                     'SMxR'='RxSM'))

unique(all.results.model$Comparison)


all.results.model$Stat<-gsub('logFC','fc',all.results.model$Stat)
all.results.model$Stat<-gsub('FDR','padj',all.results.model$Stat)

all.results.model$Stat<-factor(all.results.model$Stat,levels=c('padj','fc'),labels=c('padj','fc'))



head(all.results.model)
all.results.modelc<-dcast(all.results.model,ensembl_gene_id+wormbase_gene_seq_name+ensembl_transcript_id+external_gene_name+PairRM+SamplesRM~Comparison+Stat,value.var = 'Value')




dim(all.results.modelc)
head(all.results.modelc)

colnames(all.results.modelc)

colnames(all.results.modelc)[7:18]

colnames(all.results.modelc)[7:18]<-gsub('_','',colnames(all.results.modelc)[7:18])


head(all.results.modelc)

dim(all.results.modelc)

sample<-read.csv('~/Dropbox/Projects/OP50_Celegans_holobiont/Celegans/MTA_by_RNAseq/sample-data.csv',sep='\t')
sample$external_gene_id<-as.character(sample$external_gene_id)
sample$ensembl_gene_id<-as.character(sample$ensembl_gene_id)
dim(sample)


sample$external_gene_id
sample$ensembl_gene_id


dim(all.results.modelc)



length(intersect(sample$external_gene_id,all.results.modelc$external_gene_name))

length(intersect(sample$ensembl_gene_id,all.results.modelc$external_gene_name))


length(intersect(sample$ensembl_gene_id,all.results.modelc$wormbase_gene_seq_name))


all.results.modelc$CK_external_gene_id<-NA


#Match by Ex Ex
match.CK_Ex.Ex<-match(all.results.modelc$external_gene_name,sample$external_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CK_Ex.Ex],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))

#Match by En Ex
match.CK_En.Ex<-match(all.results.modelc$ensembl_transcript_id,sample$external_gene_id )
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CK_En.Ex],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))



#Match by Ex En
match.CK_Ex.En<-match(all.results.modelc$external_gene_name,sample$ensembl_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CK_Ex.En],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))


#Match by En t En
match.CK_Ent.En<-match(all.results.modelc$ensembl_transcript_id,sample$ensembl_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CK_Ent.En],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))



#Match by En g En
match.CK_Eng.En<-match(all.results.modelc$wormbase_gene_seq_name,sample$ensembl_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CK_Eng.En],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))



head(all.results.modelc)


head(all.results.modelc)
colnames(all.results.modelc)
all.results.modelc<-all.results.modelc[,c('CK_external_gene_id',colnames(all.results.modelc)[1:18])]
head(all.results.modelc)


# 
# 
# 
# gnid<-read.table('../../Metabolic_models/Celegans/genelist.txt',header=FALSE)
# colnames(gnid)<-c('Gene_id')
# 
# colnames(all.results.modelc)
# 
# 
# length(intersect(gnid$Gene_id,all.results.modelc$CK_external_gene_id))
# length(intersect(gnid$Gene_id,all.results.modelc$gene_name))
# length(intersect(gnid$Gene_id,all.results.modelc$EntrezGene))
# 
# intersect(gnid$Gene_id,all.results.modelc$EntrezProtClean)
# 
# 
# View(subset(all.results.modelc,EntrezProt %in% intersect(gnid$Gene_id,all.results.modelc$EntrezProt)))

# 
# write.table(matches[,c('CK_external_gene_id','gene_name','EntrezGene','EntrezProt')],'../../Metabolic_models/Celegans/Fixed_matches.csv',sep='\t',row.names = FALSE,quote = FALSE)


write.table(all.results.modelc,'~/Dropbox/Projects/OP50_Celegans_holobiont/Data/Celegans_metformin_RNAseq_Bootstrapping_complete_4_removed_1thrs.csv',sep='\t',row.names = FALSE,quote = FALSE)

all.bs<-unique(as.character(all.results.modelc$SamplesRM))

for (bs in all.bs) {
  all.results.bs<-subset(all.results.modelc,SamplesRM==bs)
  all.results.bs$PairRM<-NULL
  all.results.bs$SamplesRM<-NULL
  fname<-paste('~/Dropbox/Projects/OP50_Celegans_holobiont/Data/RNAseq_Bootstrapping_separate_4_removed_1thrs/Celegans_metformin_RNAseq_Bootstrapping_',bs,'.csv',sep='')
  write.table(all.results.bs,fname,sep='\t',row.names = FALSE,quote = FALSE)
}





