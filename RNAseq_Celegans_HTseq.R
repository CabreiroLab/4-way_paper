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



cwd<-"~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/"
keggxml<-'~/Dropbox/Projects/2015-Metformin/Annotations/Celegans/KEGG_pathways/'
setwd(cwd)
#load(".RData")

odir<-'Results_1thrs_HTseq'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")




pheno_data = read.csv("Pheno_data_new.csv")

pheno_data


#HTseq genes

gfiles<-list.files('HTseq')
gfiles<-gfiles[grepl('HTseq_counts_gene_',gfiles)]

HTgr<-data.frame()
for (fl in gfiles){
  print(fl)
  sample<-gsub('HTseq_counts_gene_','',fl)
  sample<-gsub('.txt','',sample)
  rd<-read.csv(paste0('HTseq/',fl),sep='\t',header = FALSE,comment.char = '-')
  rownames(rd)<-rd$V1
  rd$V1<-NULL
  colnames(rd)<-sample
  if (dim(HTgr)[[1]]==0) {
    HTgr<-rd
  } else {
    HTgr<-cbind(HTgr,rd)
  }
}




rownames(HTgr)<-gsub('gene:','',rownames(HTgr))


HTg.stats<-HTgr[grepl('__',rownames(HTgr)),]

HTg<-HTgr[!(grepl('__',rownames(HTgr)) | grepl('transcript',rownames(HTgr))),]




dim(HTgr)

#Some lost trasncripts?
HTgt<-HTgr[grepl('transcript',rownames(HTgr)),]

apply(HTgt,2,sum)

HTg[order(apply(HTg,1,sum),decreasing = TRUE),]

dim(HTg)









#Ensembl annotation
celegans87 = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl",version=87)
listDatasets(celegans87)
listFilters(celegans87)

listAttributes(celegans87)

celegans.annotation <- getBM(attributes=c('ensembl_gene_id','external_gene_name','wormbase_gene_seq_name','entrezgene','description'),
                             mart = celegans87)


#
head(celegans.annotation)
dim(celegans.annotation)
#47065



#dupl.entrez<-duplicated(celegans.annotation$entrezgene)
#celegans.annotation<-celegans.annotation[!dupl.entrez,]



#Remove duplicated entries with multiple entrez gene ids
dupl.ensembl_g_id<-duplicated(celegans.annotation$ensembl_gene_id)
celegans.annotation.c<-celegans.annotation[!dupl.ensembl_g_id,]

dim(celegans.annotation.c)


#47065
dim(HTg)
length(intersect(rownames(HTg),unique(celegans.annotation$ensembl_gene_id)))


annotation<-merge(celegans.annotation.c,HTg,by.x='ensembl_gene_id',by.y='row.names',all.y=TRUE)
head(annotation)
dim(annotation)

expression<-annotation[apply(annotation[,as.character(pheno_data$ids)],1,sum)>0,]

dim(expression)
#22653


rownames(expression)<-expression$ensembl_gene_id
head(expression)


#write.csv(expression.all$Transcript,paste(odir,"/Raw_data_for_transcripts.csv",sep=''),row.names=FALSE)

write.csv(expression,"Raw_data_for_genes_HTseq.csv")
#write.xl(cnts.gaf.w,explanations,paste(odir,"/Raw_data_for_genes.xlsx",sep=''))









expression<-read.csv('Raw_data_for_genes_HTseq.csv',check.names = FALSE,row.names=1)




ancols<-setdiff(colnames(expression),c(as.character(pheno_data$ID)) )
ancols

pheno_data

#Split table for edgeR


gene.info<-expression[,ancols]


gene.data<-expression[,as.character(pheno_data$ID)]


head(gene.info)


colnames(gene.data)<-as.character(pheno_data$Sample)


head(gene.data)


dt<-DGEList(counts=gene.data,
            genes=gene.info,
            group=factor(pheno_data$Group))


#How to choose right threshold

dim(dt)
apply(dt$counts, 2, sum)
#At least one count per million (cpm) reads in at least 4 consistent samples (one group)

thrs<- 1
consis<- 4



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

#HTseq
#22501->11513



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




strains<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$Sample)),'Strain'])
metf<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$Sample)),'Metformin_mM'])


strains.fac<-factor(strains)
metf.fac<-factor(metf)

strains.col<-ifelse(strains.fac=='OP50','white','black')
metf.col<-ifelse(metf.fac=='0','green','red')

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
#d2<-estimateDisp(d,design.mat,prior.df = priordf)



d2<-estimateDisp(d,design.mat)
print(paste('Prior.n=',d2$prior.n,' Residual.df=',d2$prior.df/d2$prior.n,' Prior.df=',d2$prior.df,sep=''))


plotBCV(d2)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Dispersion.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)




#All design matrices
design.mat
fit <- glmFit(d2, design.mat)


gofsum<-gof(fit,plot=TRUE,adjust='fdr',pcutoff = 0.05)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_GOF_qq-plot.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#Outlier genes
table(gofsum$outlier)

gofsum$outlier[gofsum$outlier==TRUE]

#Significant outliers
table(gofsum$gof.pvalues<0.05)





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
  enrcoef<-'logFC'
  #print(paste('Threshold:',thre,sep=' '))
  lrt <- glmLRT(fit, contrast=contr)
  
  
  print('...Differential gene expression')
  de.glm <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
  print(summary(de.glm))
  
  sumr<-t(summary(de.glm))
  rownames(sumr)<-comp
  summary.all<-rbind(summary.all,sumr)
  
  de2tags <- rownames(d2)[as.logical(de.glm)]
  cairo_pdf(paste(odir,'/MVA_',comp,'.pdf',sep = ''),width=9,height=6)
  plotSmear(lrt, de.tags=de2tags)
  abline(h = c(-1, 1), col = "blue")
  title(paste(desc,comp,sep='\n'))
  dev.off()
  
  result<-as.data.frame(topTags(lrt,n=dim(d2$genes)[1]))
  result$Comparison<-comp
  #result$Threshold<-thre
  
  
  all.results<-rbind(all.results,data.table(result))
  
  
  if (summary(de.glm)[3]>10 | summary(de.glm)[1] > 10) {
    print('...KEGG enrichment')
    enrichment.KEGG<-kegga(lrt,coef=enrcoef,geneid = 'entrezgene',species.KEGG='cel',gene.pathway = gp.KEGG,convert=TRUE,FDR=0.05)#,trend='length'
    
    enrichment.KEGG$ID<-rownames(enrichment.KEGG)
    enrichment.KEGG<-rename(enrichment.KEGG,c("Pathway"="Description"))
    enrichment.KEGG$Comparison<-comp
    #enrichment.KEGG$Threshold<-thre
    all.KEGGenrichment<-rbind(all.KEGGenrichment,data.table(enrichment.KEGG))
    
    print('...GO enrichment')
    enrichment.GO<-goana(lrt,species='Ce',coef=enrcoef,geneid='entrezgene',convert=TRUE,FDR=0.05)#,trend='length'
    
    enrichment.GO$ID<-rownames(enrichment.GO)
    enrichment.GO<-rename(enrichment.GO,c("Term"="Description"))
    enrichment.GO$Comparison<-comp
    #enrichment.GO$Threshold<-thre
    all.GOenrichment<-rbind(all.GOenrichment,data.table(enrichment.GO))
    
  } else {
    print("No DE genes!")
  }
  print('---------')
}



colnames(summary.all)<-c('Contrast','Change','Count')

summary.all$Change<-revalue(summary.all$Change,c('0'='None','-1'='Down','1'='Up'))

summary.counts<-dcast(summary.all,Contrast~Change,value.var='Count')
summary.counts

summary.counts$Prc_change<-apply(summary.counts[,c('Down','Up')],1,sum)*100/apply(summary.counts[,c('Down','None','Up')],1,sum)

write.csv(summary.counts,paste(odir,'/Contrast_summary.csv',sep = ''),row.names = FALSE)




all.GOenrichment<-data.frame(all.GOenrichment)
dim(all.GOenrichment)

all.KEGGenrichment<-data.frame(all.KEGGenrichment)
dim(all.KEGGenrichment)



all.results<-data.frame(all.results)
dim(all.results)
all.results$unshrunk.logFC<-NULL
all.results$Threshold<-NULL



all.results$Comparison<-factor(all.results$Comparison,
                               levels=allcomp,
                               labels=allcomp)



all.results.a<-merge(all.results,allcomparisons[,c('Contrast','Description')],by.x='Comparison',by.y='Contrast',all.x=TRUE)


head(all.results.a)





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


all.results.k<-merge(all.results.a,allKEGGinfo.paths,by.x='entrezgene',by.y='EntrezGene',all.x=TRUE)
all.results.kg<-merge(all.results.k,allGOinfo.terms,by.x='entrezgene',by.y='EntrezGene',all.x=TRUE)

options(scipen=5)



head(all.results.kg)
dim(all.results.kg)


colnames(all.results.kg)
oldorder<-colnames(all.results.kg)
all.results.kg<-all.results.kg[,c('Comparison','Description','ensembl_gene_id','external_gene_name','entrezgene','wormbase_gene_seq_name','description','logFC','logCPM','LR','PValue','FDR','Pathways','GO_terms')]
print('All columns taken?')
length(colnames(all.results.kg))==length(oldorder)


head(all.results.kg)

write.csv(all.results.kg,paste(odir,'/All_results.csv',sep = ''),row.names = FALSE)
#write.xl(all.results.kg,explanations,paste(odir,'/All_results.xlsx',sep = ''))

write.csv(all.KEGGenrichment,paste(odir,'/All_results_KEGG.csv',sep = ''),row.names = FALSE)
#write.xl(all.KEGGenrichment,explanations,paste(odir,'/All_results_KEGG.xlsx',sep = ''))

write.csv(all.GOenrichment,paste(odir,'/All_results_GO.csv',sep = ''),row.names = FALSE)
#write.xl(all.GOenrichment,explanations,paste(odir,'/All_results_GO.xlsx',sep = ''))







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
allcomp


add.path<-c('cel04213','cel04212')
allpaths<-c()

setwd(cwd)


head(all.results.kg)

mcols<-c('ensembl_gene_id','external_gene_name','entrezgene','description','Pathways','GO_terms','Comparison')

measure.vars<-c('logFC','FDR')


#Melt
all.results.rm<-melt(all.results.kg,id.vars = mcols,
                     measure.vars=measure.vars,
                     variable.name = 'Stat',value.name = 'Value')

all.results.rmc<-subset(all.results.rm,Stat %in% c(enrcoef,'FDR'))
all.results.rcp<-dcast(all.results.rmc,ensembl_gene_id+external_gene_name+entrezgene+description+Pathways+GO_terms~Comparison+Stat,value.var = 'Value')

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



write.csv(all.results.rcp,paste(odir,'/All_results_sidebyside.csv',sep = ''),row.names = FALSE)
#write.xl(all.results.rcp,explanations,paste(odir,'/All_results_sidebyside_Threshold-',thres,'.xlsx',sep = ''),'Readme')



#Volcano plots


comp.vgrps<-list('Treatment'=c('OP_D','MR_D','Longevity'),
                 'General'=c('D','MROP_G'))



for (cvg in names(comp.vgrps)) {
  #cvg<-'Ef_adaptation_3h'
  
  all.plot<-subset(all.results,Comparison %in% comp.vgrps[[cvg]])
  
  
  vol<-ggplot(all.plot,aes(x=logFC,y=-log10(FDR)))+
    geom_hline(yintercept = -log(0.05,10),color=barcolor,alpha=baralpha)+
    geom_vline(xintercept = -1,color=barcolor,alpha=baralpha)+
    geom_vline(xintercept = 1,color=barcolor,alpha=baralpha)+
    ggtitle(paste0('Volcano plot for differentially expressed genes: ',cvg),
            subtitle='Labels are shown for some of the significant changes')+
    geom_point(alpha=0.9,size=1)+
    scale_x_continuous(breaks=seq(-10,10,by=1),limits=c(-5,5))+
    scale_y_continuous(breaks=seq(0,20,by=2),limits=c(0,18))+
    #scale_y_continuous(breaks=seq(0,20,by=2),limits=c(0,10))+
    #   geom_label_repel(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp,as.character(gene_name),'')),
    #             size=txtsize,colour = "red")+
    geom_text_repel(aes(label=ifelse( -log10(FDR)>4 & abs(logFC)>1 ,as.character(external_gene_name),'')),
                    size=txtsize,colour = "red")+
    # geom_text_repel(aes(label=ifelse( FDR<0.05 & abs(logFC)>1 ,as.character(external_gene_name),'')),
    #                 size=txtsize,colour = "red")+
    facet_wrap(~Comparison,ncol=6)
  vol
  
  fname<-paste(odir,'/Volcano_plot_',cvg,'.pdf',sep = '')
  print(fname)
  cairo_pdf(fname,width=24,height=8)
  print(vol)
  dev.off()
  
}



#Enrichment

#For heatmaps
comp.all<-allcomp

comp.write<-allcomp





comp.grps<-list('All'=allcomp,
                'Treatment'=c('OP_D','MR_D','Longevity'),
                'General'=c('D','MROP_G'))


for (ent in annotations) {
  print(paste('Enrichment:',ent),sep=' ')
  if (ent=='KEGG') {
    all.enrichment<-all.KEGGenrichment #,Threshold==thres)
    idvars<-c('ID','Description','N','Comparison')
  } else {
    all.enrichment<-all.GOenrichment #,Threshold==thres)
    idvars<-c('ID','Ont','Description','N','Comparison')
  }
  all.enrichment$ID<-gsub('path:','',all.enrichment$ID)
  all.enrichment$logP_Up<- -log10(all.enrichment$P.Up)
  all.enrichment$logP_Down<- -log10(all.enrichment$P.Down)
  all.enrichment$Comparison<-factor(all.enrichment$Comparison,
                                    levels=allcomp,
                                    labels=allcomp)
  
  enrichmentm<-reshape2::melt(data.frame(all.enrichment),id.vars = idvars,variable.name = 'Stat',value.name = 'Value')
  
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
  
  write.csv(enrichmentcea,paste(odir,'/All_results_sidebyside_',ent,'.csv',sep = ''),row.names = FALSE)
  #write.xl(enrichmentcea,explanations,paste(odir,'/All_results_sidebyside_',ent,'_Threshold-',thres,'.xlsx',sep = ''),'Readme')
  
  for (ont in ontologies) {
    for (comp.grh in names(comp.grps)) {
      comp.heat<-comp.grps[[comp.grh]]
      if (ent=='GO') {
        enrichmentco<-subset(enrichmentc,Ont==ont)
        fname<-paste0(odir,'/Enrichment_heatmap_',ent,'-',ont,'_',comp.grh,'.pdf')
      } else {
        enrichmentco<-enrichmentc
        fname<-paste0(odir,'/Enrichment_heatmap_',ent,'_',comp.grh,'.pdf')
      }
      
      #hdata<-enrichmentco[,!colnames(enrichmentco) %in% c('ID','Description','N','Ont')]
      #Filter just selected columns
      hdata<-enrichmentco[,grep(paste(paste0(comp.heat,'_logP'),collapse="|"),colnames(enrichmentco)) ]
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

pathcomp<-c('OP_D_logFC','MR_D_logFC','Longevity_logFC')
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

SM_S.up<-as.character(subset(all.results.rcp,`OP_D_FDR`<0.05 & `OP_D_logFC`>0)[,'ensembl_gene_id'])
SM_S.down<-as.character(subset(all.results.rcp,`OP_D_FDR`<0.05 & `OP_D_logFC`<0)[,'ensembl_gene_id'])
RM_R.up<-as.character(subset(all.results.rcp,`MR_D_FDR`<0.05 & `MR_D_logFC`>0)[,'ensembl_gene_id'])
RM_R.down<-as.character(subset(all.results.rcp,`MR_D_FDR`<0.05 & `MR_D_logFC`<0)[,'ensembl_gene_id'])

S.change<-as.character(subset(all.results.rcp,`OP_D_FDR`<0.05 )[,'ensembl_gene_id'])
R.change<-as.character(subset(all.results.rcp,`MR_D_FDR`<0.05 )[,'ensembl_gene_id'])
L.change<-as.character(subset(all.results.rcp,`Longevity_FDR`<0.05 )[,'ensembl_gene_id'])



venn.groups<-list('OP50-MR Down'=RM_R.down,
                  'OP50-MR Up'=RM_R.up,
                  'OP50 Down'=SM_S.down,
                  'OP50 Up'=SM_S.up)

ups<-intersect(SM_S.down,RM_R.up)
downs<-intersect(SM_S.up,RM_R.down)

write.csv(subset(all.results.rcp,gene_id %in% downs),paste(odir,'/Flipping_downregulated.csv',sep = ''),row.names = FALSE)
write.csv(subset(all.results.rcp,gene_id %in% ups),paste(odir,'/Flipping_upregulated.csv',sep = ''),row.names = FALSE)




plot(Venn(venn.groups),show = list(Faces = FALSE),
     doWeights = FALSE,type='ellipses')

dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram.pdf",sep=''),width=6,height=6)





venn.longevity<-list('OP50 Metformin'=S.change,
                  'OP50-MR Metformin'=R.change,
                  'Longevity'=L.change)

plot(Venn(venn.longevity),show = list(Faces = FALSE),
     doWeights = TRUE)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,"/Venn_diagram_longevity.pdf",sep=''),width=6,height=6)



detach("package:Vennerable", unload=TRUE)
detach("package:reshape", unload=TRUE)
detach("package:reshape2", unload=TRUE)

library(reshape2)







#Comparison scatter


fit<-lm(`MR_D_logFC`~`OP_D_logFC`,data=all.results.rcp)
res<-summary(fit)
res
res$r.squared


dim(subset(all.results.rcp,sqrt( (`OP_D_logFC`*1.5)^2+(`MR_D_logFC`)^2 ) >5 &
             (`OP_D_FDR`<0.05 |
             `MR_D_FDR`<0.05) ))



library(ellipse)


ellipsoid<-as.data.frame(ellipse( 0.45,
                                  scale=c(1,2),
                                  centre=c( 0,0) ))

ggplot()+geom_path(data=ellipsoid, aes(x=x, y=y),size=1)



all.results.rcp$Change<-'None'
all.results.rcp$Change<-ifelse(all.results.rcp$`OP_D_FDR`<0.05,'Sensitive',all.results.rcp$Change)
all.results.rcp$Change<-ifelse(all.results.rcp$`MR_D_FDR`<0.05,'Resistant',all.results.rcp$Change)
all.results.rcp$Change<-ifelse(all.results.rcp$`MR_D_FDR`<0.05 & all.results.rcp$`OP_D_FDR`<0.05,'Both',all.results.rcp$Change)

all.results.rcp$Change<-factor(all.results.rcp$Change,levels=c('None','Sensitive','Resistant','Both'))




allS<-subset(all.results.rcp,Change=='Sensitive')
allR<-subset(all.results.rcp,Change=='Resistant')
allB<-subset(all.results.rcp,Change=='Both')

topS<-allS[order(allS$`OP_D_FDR`) ,'t_id'][1:10]
topR<-allR[order(allR$`MR_D_FDR`) ,'t_id'][1:10]
topB<-allB[order(allB$`OP_D_FDR`,allB$`MR_D_FDR`) ,'t_id'][1:10]


topAll<-c(51908)#topS,topR,topB,topR,



all.results.rcp[order(abs(all.results.rcp$`Longevity_logFC`),decreasing = TRUE),]


amp<-5
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue2','blue2','gray80','red2','red2')



subset(all.results.rcp,gene_name=='acs-2')


max(all.results.rcp$OP_D_logFC)
min(all.results.rcp$OP_D_logFC)

max(all.results.rcp$MR_D_logFC)
min(all.results.rcp$MR_D_logFC)


ggplot(all.results.rcp,aes(x=`OP_D_logFC`,y=`MR_D_logFC`))+
  geom_vline(xintercept = 0,alpha=0.2)+
  geom_hline(yintercept = 0,alpha=0.2)+
  geom_abline(slope=1,intercept = 0,alpha=0.2)+
  geom_abline(slope = fit$coefficients[[2]],intercept = fit$coefficients[[1]],color='red',alpha=0.5)+
  geom_point(data=all.results.rcp,aes(color=`Longevity_logFC` ),stroke=0,size=1)+
  labs(color='Significant change')+
  xlab('Metformin effect when on OP50, logFC')+
  ylab('Metformin effect when on OP50-MR, logFC')+
  scale_colour_gradientn(colors = gradcols,limits=c(-amp,amp),name='Longevity effect')+
  scale_size(range=c(0,8),breaks=seq(0,12,by=2),name='Amplitude of change')+
  scale_x_continuous(breaks=seq(-12,12,by=2),limits=c(-10,10))+
  scale_y_continuous(breaks=seq(-10,10,by=2),limits=c(-4,8))+
  #geom_text_repel(aes(label=ifelse( t_id %in% topAll,as.character(gene_name),'') ))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_SensRes_longevity.pdf",sep=''),
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





#Rubbish
# compare (group 2 - group 1) to 0:
# this is equivalent to comparing group 1 to group 2
lrt12 <- glmLRT(fit, contrast=contrasts[['SM-S']])#Specific treatment effect


#lrt12 <- glmLRT(fit, contrast=c(-1,1,1,-1,0))#Lifespan

de12.glm <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
summary(de12.glm)
de2tags12 <- rownames(d2)[as.logical(de12.glm)]
plotSmear(lrt12, de.tags=de2tags12)
abline(h = c(-1, 1), col = "blue")

View(topTags(lrt12, n=30))


#subset(full_table,gene_id %in% c('MSTRG.18297'))

# getGeneKEGGLinks('cel',convert = TRUE)
# getKEGGPathwayNames('cel',remove=TRUE)

gp<-subset(getGeneKEGGLinks('cel',convert = TRUE),!is.na(GeneID) & ! is.na(PathwayID) & PathwayID!='path:cel01100')

#lrt12$genes$EntrezGene

keg12<-kegga(lrt12,geneid = 'EntrezGene',coef='logFC',convert=TRUE,
             species.KEGG='cel',FDR=0.05,gene.pathway=gp.KEGG,trend='length')

topKEGG(keg12)

keg12$Comparison<-'SM/S'


#writre script to compare gene numbers from enrichment and in table


go12<-goana(lrt12,species='Ce',geneid='EntrezGene',convert=TRUE,FDR=0.05)


View(topGO(go12,n=50))

#df residual needs to be 11!

# 
biocLite(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db"))
biocLite("pathview")
#biocLite("gage")



# source("http://bioconductor.org/biocLite.R")
# biocLite(c("pathview"))






#Pathview



gdata<-subset(all.results.rcp,Comparison=='SM-S' & Threshold=='None' & ! is.na(EntrezGene))

nrow(all.results.rcp)
gdata<-subset(all.results.rcp,! is.na(EntrezGene))


nrow(gdata)
o<-order(gdata[,'SM-S_FDR'],decreasing=FALSE)
gdata<-gdata[o,]


head(gdata)

mapkey<-'EntrezGene'

#Do selection by gene id, because gene_id is associated with function
dupl<-duplicated(gdata[,mapkey])
table(dupl)

gdataf<-gdata[!dupl,]
rownames(gdataf)<-gdataf[,mapkey]
nrow(gdataf)

#gdatam<-as.matrix(gdataf[,'logFC',drop=FALSE])

head(gdataf)


i<-1

setwd(cwd)

odir<-paste('/Pathview_test',sep='')

ecoxml<-'~/Projects/2015-Metformin/Annotations/Ecoli/KEGG_pathways/'


download.kegg(pathway.id = allpaths,species = 'eco',kegg.dir = keggxml)





allpaths<-gsub('cel','',allKEGGinfo.genes$PathID)




pv.out <- pathview(gene.data = gdataf[,c('SM-S_logFC','RM-R_logFC','SM-S-(RM-R)_logFC'),drop=FALSE],
                   pathway.id = "00010",species = "cel",
                   kegg.dir = "KEGG_pathways_Celegans",
                   out.suffix = "comparison",
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   low = list(gene = "blue", cpd = "green"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="yellow"),
                   limit=list(gene=1,cpd=1),node.sum = "mean",
                   trans.fun = list(gene=NULL,cpd=NULL))




# keggdir<-paste("~/Projects/2015-Metformin/Sequencing/CElegans/",odir,sep='')
# dir.create(keggdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# setwd(keggdir)




pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = T)


pv.out <- pathview(gene.data = gdataf[,'logFC',drop=FALSE], pathway.id = "cel00020",
                   species = "cel", out.suffix = "test",kegg.native = T,map.symbol=TRUE,
                   new.signature=FALSE,plot.col.key=TRUE)#, ,gene.idtype='ENSEMBLPROT',gene.idtype='entrez'



pv.out <- pathview(gene.data = gdataf[,c('SM-S_logFC','RM-R_logFC','SM-S-(RM-R)_logFC'),drop=FALSE],
                   pathway.id = "04212",species = "cel", out.suffix = "comparison",
                   same.layer=FALSE,kegg.native = T,map.symbol=FALSE,
                   map.null=FALSE,new.signature=FALSE, plot.col.key=TRUE,
                   res=300,
                   low = list(gene = "blue", cpd = "green"),
                   mid=list(gene = "gray", cpd = "gray"),
                   high = list(gene = "red", cpd ="yellow"),
                   limit=list(gene=2,cpd=2),node.sum = "mean",
                   trans.fun = list(gene=NULL,cpd=NULL))

pv.out




pv.out

#pv.out


min(gdataf[,c('SM-S_logFC','RM-R_logFC','SM-S-(RM-R)_logFC'),drop=FALSE])
max(gdataf[,c('SM-S_logFC','RM-R_logFC','SM-S-(RM-R)_logFC'),drop=FALSE])

#,,gene.idtype="entrez"
#gene.annotpkg='org.Ce.eg.db'

subset(gdataf,ForKEGG=='CELE_F20H11.3')
subset(all.results.kg,EntrezProt=='F54H12.1')

pv.out

data(gene.idtype.list)
gene.idtype.list



path2eg<-getGeneKEGGLinks('cel',convert = TRUE)
path2eg<-subset(path2eg,!is.na(GeneID))

subset(path2eg,is.na(GeneID))

path2eg$PathID<-gsub('path:','',path2eg$PathwayID)
path2eg$EntrezGene<-path2eg$GeneID
head(path2eg)


dim(path2eg)

length(unique(all.results.kg$EntrezGene))

length(intersect(unique(all.results.kg$EntrezGene),path2eg$EntrezGene))

head(path2eg)


allTCA<-subset(path2eg,PathID=='cel01200')
dim(allTCA)



