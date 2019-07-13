options(java.parameters = "-Xmx4096m") 

#That's for annotations
# library(org.Ce.eg.db)
# library(org.EcK12.eg.db)
# library(celegans.db)
# library(GO.db)



library(edgeR)

# library(RSkittleBrewer)
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


#devtools::install_github("PNorvaisas/PFun")
library(PFun)

#load('Celegans_RNAseq.RData')
#save.image('Celegans_RNAseq.RData')



theme_set(theme_light())


#In code should be replaced with rowbind or similar alternative
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
odir<-'Results_1thrs_newannot'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


pheno_data = read.csv("Pheno_data.csv")




#Ensembl annotation
celegans87 = useEnsembl(biomart="ensembl", dataset="celegans_gene_ensembl",version=87)
listDatasets(celegans87)
listFilters(celegans87)

listAttributes(celegans87)

celegans.annotation <- getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id','external_gene_name','description'),
                             mart = celegans87)




head(celegans.annotation)
dim(celegans.annotation)



#Get clean read counts
cnts.gaf.o<-read.table('Results/Raw_data_for_genes.csv',sep = ',',header = TRUE)
cnts.gaf.o$EntrezProt<-NULL
cnts.gaf.o$Symbol<-NULL
cnts.gaf.o$gene_name<-NULL
cnts.gaf.o<-rename(cnts.gaf.o,c('EntrezGene'='entrezgene'))


#Use ENSEMBL transcript IDs (WormBase IDs) for annotations
cnts.gaf.o$ensembl_transcript_id<-ifelse(grepl('transcript',cnts.gaf.o$t_name), gsub('transcript:','',cnts.gaf.o$t_name),NA  ) 


head(cnts.gaf.o)
dim(cnts.gaf.o)

cnts.gaf<-merge(cnts.gaf.o,celegans.annotation,by='ensembl_transcript_id',all.x=TRUE,all.y=FALSE)


dim(cnts.gaf)
head(cnts.gaf)

ancols<-c('ensembl_transcript_id','ensembl_gene_id','t_id','gene_id','external_gene_name','entrezgene','description','chr','strand','start','end','num_exons','length')

#Split table for edgeR
gene.info<-cnts.gaf[,ancols]
gene.data<-cnts.gaf[,as.character(pheno_data$ids)]

raw.data<-cnts.gaf[,c("gene_id","ensembl_gene_id",as.character(pheno_data$ids))]


colnames(cnts.gaf)

colnames(raw.data)



ids<-as.character(pheno_data$ids)



#Split raw data by sample
for (id in ids) {
  tbl<-raw.data[,c("gene_id","ensembl_gene_id",id)]
  write.table(tbl,file = paste0('Raw_counts/',id,"_gene_read_counts.txt"),sep = '\t',row.names = FALSE)
}


write.csv(raw.data,paste0(odir,"/Raw_clean_gene_data.csv"))







table(!is.na(gene.info$entrezgene))


head(gene.info)
head(gene.data)

groups<-c('C','C','C','C','M','M','M','M','R','R','R','R','RM','RM','RM','RM')



dt<-DGEList(counts=gene.data,
            genes=gene.info,
            group=as.character(pheno_data$Group) )





#How to choose right threshold

dim(dt)
apply(dt$counts, 2, sum)
#At least one count per million (cpm) reads in at least 4 consistent samples (one group)

thrs<- 0
consis<- 4

keep<- rowSums(cpm(dt)[,grep('^C[[:digit:]]',colnames(dt))]>thrs)>= consis |
  rowSums(cpm(dt)[,grep('^M[[:digit:]]',colnames(dt))]>thrs)>= consis |
  rowSums(cpm(dt)[,grep('^R[[:digit:]]',colnames(dt))]>thrs)>= consis | 
  rowSums(cpm(dt)[,grep('^RM[[:digit:]]',colnames(dt))]>thrs)>= consis

#keep <- rowSums(cpm(dt)>1) >= samplethreshold


d <- dt[keep,]
dim(d)


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

#Normalisation
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d$samples


#method="bcv",
plotMDS(d,col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=10)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#Batch removal
logCPM <- cpm(d, log=TRUE, prior.count=1)
batch<-c('1','1','2','2',
         '1','2','2','2',
         '1','2','2','2',
         '1','1','2','2')
logCPMc <- removeBatchEffect(logCPM, batch)

plotMDS(logCPMc, col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=1)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_BCV_batch-adjusted.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)





#Remove numbers from sample IDs
groups<-gsub('[[:digit:]]+', '', colnames(logCPM))
annData<-factor(groups)
ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]
#terrain.colors

#Heatmap raw
heatmap3(as.matrix(logCPM),
         scale = 'row',
         ColSideLabs = 'Group',
         balanceColor=TRUE,
         ColSideColors=ColSideColors,
         labRow=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Heatmap_logCPM.pdf",sep=''),
             width=9,height=9, useDingbats=FALSE)


groups<-gsub('[[:digit:]]+', '', colnames(logCPMc))
annData<-factor(groups)
ColSideColors<-rainbow(length(unique(annData)))[as.numeric(annData)]
#terrain.colors


#Heatmap adjusted
heatmap3(as.matrix(logCPMc),
         scale = 'row',
         ColSideLabs = 'Group',
         ColSideColors=ColSideColors,
         balanceColor=TRUE,
         labRow=FALSE)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Heatmap_logCPMc_batchadj.pdf",sep=''),
             width=9,height=9, useDingbats=FALSE)




#GLM tests
# batch.lvl<-c('1','1','2','2',
#                     '1','2','2','2',
#                     '1','2','2','2',
#                     '1','1','2','2')
design.mat <- model.matrix(~ 0 + d$samples$group+pheno_data$Batch) #
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


contrasts<-list('SM-S'=c(-1,1,0,0,0),
                'RM-R'=c(0,0,-1,1,0),
                'R-S'=c(-1,0,1,0,0),
                'RM-SM'=c(0,-1,0,1,0),
                'SM-R'=c(0,1,-1,0,0),
                'RM-S'=c(-1,0,0,1,0),
                'T-C_general'=c(-1,1,-1,1,0),
                'R-S_general'=c(-1,-1,1,1,0),
                'S-R_general'=c(1,1,-1,-1,0),
                'SM-S-(RM-R)'=c(-1,1,1,-1,0))


allcomp<-c('SM-S','RM-R','SM-S-(RM-R)','R-S','RM-SM','SM-R','RM-S','T-C_general','R-S_general','S-R_general')



comparisons<-c('SM-S','RM-R','R-S','RM-SM','SM-R','RM-S',
               'T-C_general','R-S_general','S-R_general','SM-S-(RM-R)')

oldids<-c('M-C','RM-R','R-C','RM-M','M-R','RM-C',
          'T-C_general','R-C_general','C-R_general','M-C-(RM-R)')

description<-c('Treatment effect while on sensitive strain',
               'Treatment effect while on resistant strain',
               'Difference between strains in control (R-S)',
               'Difference between strains in treatment (RM-SM)',
               'Diagonal difference (SM-R)',
               'Diagonal difference (RM-S)',
               'General treatment over control',
               'General resistant over sensitive',
               'General sensitive over resistant',
               'Difference between strains in treatment\n(R adjusted)')



allcomparisons<-data.frame(comparison=comparisons,
                           ID=oldids,
                           description=description)




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
  for (thre in c('None')) {
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
    
    all.results<-bind_rows(all.results,data.table(result))

    
    if (summary(de.glm)[3]>10 | summary(de.glm)[1] > 10) {
      print('...KEGG enrichment')
      enrichment.KEGG<-kegga(lrt,coef=enrcoef,geneid = 'entrezgene',species.KEGG='cel',gene.pathway = gp.KEGG,convert=TRUE,FDR=0.05,trend='length')
      enrichment.KEGG$ID<-rownames(enrichment.KEGG)
      enrichment.KEGG<-rename(enrichment.KEGG,c("Pathway"="Description"))
      enrichment.KEGG$Comparison<-comp
      enrichment.KEGG$Threshold<-thre
      all.KEGGenrichment<-bind_rows(all.KEGGenrichment,data.table(enrichment.KEGG))
      
      print('...GO enrichment')
      enrichment.GO<-goana(lrt,species='Ce',coef=enrcoef,geneid='entrezgene',convert=TRUE,FDR=0.05,trend='length')
      enrichment.GO$ID<-rownames(enrichment.GO)
      enrichment.GO<-rename(enrichment.GO,c("Term"="Description"))
      enrichment.GO$Comparison<-comp
      enrichment.GO$Threshold<-thre
      all.GOenrichment<-bind_rows(all.GOenrichment,data.table(enrichment.GO))
      
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

# 
# all.results$Comparison<-factor(all.results$Comparison,
#                                levels=c('SM-S','RM-R','SM-S-(RM-R)','R-S','T-C_general','R-S_general','S-R_general'),
#                                labels=c('SM-S','RM-R','SM-S-(RM-R)','R-S','T-C_general','R-S_general','S-R_general'))



all.results<-data.frame(all.results)
dim(all.results)


#Get gene pathway mappings
path2eg<-getGeneKEGGLinks('cel',convert = TRUE)
path2eg<-subset(path2eg,!is.na(GeneID) & !is.na(PathwayID) & PathwayID!='path:cel01100')
path2eg$PathID<-gsub('path:','',path2eg$PathwayID)
path2eg$entrezgene<-path2eg$GeneID
#head(path2eg)


path2pathde<-getKEGGPathwayNames('cel',remove.qualifier = TRUE)
path2pathde$PathwayID<-gsub('path:','',path2pathde$PathwayID)

go2eg<-toTable(org.Ce.egGO2EG)
go2eg<-subset(go2eg,!is.na(gene_id) & !is.na(go_id) )
go2eg$GO_id<-go2eg$go_id
go2eg$entrezgene<-go2eg$gene_id
#head(go2eg)


allKEGGinfo.genes<-ddply(path2eg,.(PathID),summarise,Genes=paste(entrezgene,collapse=';'))
allKEGGinfo.paths<-ddply(path2eg,.(entrezgene),summarise,Pathways=paste(PathID,collapse=';'))


allGOinfo.genes<-ddply(go2eg,.(GO_id),summarise,Genes=paste(entrezgene,collapse=';'))
allGOinfo.terms<-ddply(go2eg,.(entrezgene),summarise,GO_terms=paste(GO_id,collapse=';'))


all.results.k<-merge(all.results,allKEGGinfo.paths,by='entrezgene',all.x=TRUE)
all.results.kg<-merge(all.results.k,allGOinfo.terms,by='entrezgene',all.x=TRUE)

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



c('ensembl_transcript_id','ensembl_gene_id','t_id','gene_id','external_gene_name','entrezgene','description','chr','strand','start','end','num_exons','length')


idvariables<-c('Comparison',ancols,'Pathways','GO_terms')
measvariables<-c('logFC','FDR')


for (thres in c('None')) {
  enrcoef<-'logFC'
  all.results.r<-subset(all.results.kg,Threshold==thres)
  #all.results.r$unshrunk.logFC<-NULL
  all.results.r$Threshold<-NULL
  
  #Melt
  all.results.rm<-melt(all.results.r,
                       id.vars = idvariables,
                       measure.vars=measvariables,
                       variable.name = 'Stat',value.name = 'Value')

  all.results.rmc<-subset(all.results.rm,Stat %in% c(enrcoef,'FDR'))
  all.results.rcp<-dcast(all.results.rmc,t_id+ensembl_transcript_id+ensembl_gene_id+external_gene_name+entrezgene+Pathways+GO_terms+
                           chr+strand+start+end+num_exons+length~Comparison+Stat,value.var = 'Value')

  #head(all.results.rcp)
  
  #Reorder by significance
  signsum<-apply(all.results.rcp[,grep('FDR',colnames(all.results.rcp))],1,function(x) sum(-log10(x)))
  all.results.rcp<-all.results.rcp[order(signsum,decreasing = TRUE),]
  
  nrow(all.results.rcp)
  
  gdata<-subset(all.results.rcp,! is.na(entrezgene))
  nrow(gdata)
  nrow(gdata)
  o<-order(rowSums(gdata[,grep('_FDR',colnames(gdata))]),decreasing=FALSE)
  gdata<-gdata[o,]
  
  #Do selection by gene id, because gene_id is associated with function
  dupl<-duplicated(gdata[,mapkey])
  print('Entrez Gene duplicates')
  print(table(dupl))
  gdataf<-gdata[!dupl,]
  rownames(gdataf)<-gdataf[,mapkey]
  
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
    geom_text(aes(label=ifelse(-log10(FDR)>-sp*logFC^2+tp & FDR<0.05,as.character(external_gene_name),'')),
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
}



all.results.r<-subset(all.results.kg,Threshold=='None')
all.results.r$Threshold<-NULL


all.results.rm<-melt(all.results.r,id.vars = c('Comparison','gene_id','gene_name','entrezgene','ENSEMBL','ENSEMBL_gene','ENSEMBL_transcript',
                                               'Unique','chr','strand','start','end','num_exons','length','t_id','Pathways','GO_terms'),
                     variable.name = 'Stat',value.name = 'Value')

all.results.rmc<-subset(all.results.rm,Stat %in% c('logFC','FDR'))
all.results.rcp<-dcast(all.results.rmc,t_id+gene_name+entrezgene+ENSEMBL+ENSEMBL_gene+ENSEMBL_transcript+Unique+Pathways+GO_terms+
                         chr+strand+start+end+num_exons+length~Comparison+Stat,value.var = 'Value')





#For metabolic modeling
#Choose only the necessary comparisons
all.results.model<-subset(all.results.rmc,Comparison %in% c('SM-S','RM-R','R-S','RM-SM','RM-S','SM-R') & gene_name!='.')

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


all.results.modelc<-dcast(all.results.model,gene_name+entrezgene+EntrezProt~Comparison+Stat,value.var = 'Value')


dim(all.results.modelc)

colnames(all.results.modelc)[4:15]<-gsub('_','',colnames(all.results.modelc)[4:15])

head(all.results.modelc)


sample<-read.csv('../../Metabolic_models/Celegans/sample-data.csv',sep='\t')
sample$external_gene_id<-as.character(sample$external_gene_id)
sample$ensembl_gene_id<-as.character(sample$ensembl_gene_id)
dim(sample)


all.results.modelc$EntrezProt



all.results.modelc$EntrezProtClean


length(intersect(sample$external_gene_id,all.results.modelc$gene_name))
length(intersect(sample$ensembl_gene_id,all.results.modelc$gene_name))
length(intersect(sample$ensembl_gene_id,all.results.modelc$EntrezProt))



all.results.modelc$CK_external_gene_id<-NA

#Match by gene name
match.gene<-match(all.results.modelc$gene_name,sample$external_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.gene],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))

#Match by CK ensemble
match.CKen<-match(all.results.modelc$gene_name,sample$ensembl_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CKen],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))

#Match by CK ensemble and Entrez
match.CKEP<-match(all.results.modelc$EntrezProt,sample$ensembl_gene_id)
all.results.modelc$CK_external_gene_id<-ifelse(is.na(all.results.modelc$CK_external_gene_id),sample$external_gene_id[match.CKEP],all.results.modelc$CK_external_gene_id)
table(is.na(all.results.modelc$CK_external_gene_id))

table(is.na(all.results.modelc$CK_external_gene_id))


subset(all.results.modelc,EntrezProt=='D1009.1a')



length(colnames(all.results.modelc))

all.results.modelc<-all.results.modelc[,c('CK_external_gene_id',colnames(all.results.modelc)[1:15])]

gnid<-read.table('../../Metabolic_models/Celegans/genelist.txt',header=FALSE)
colnames(gnid)<-c('Gene_id')

colnames(all.results.modelc)


length(intersect(gnid$Gene_id,all.results.modelc$CK_external_gene_id))
length(intersect(gnid$Gene_id,all.results.modelc$gene_name))
length(intersect(gnid$Gene_id,all.results.modelc$entrezgene))

intersect(gnid$Gene_id,all.results.modelc$EntrezProtClean)


View(subset(all.results.modelc,EntrezProt %in% intersect(gnid$Gene_id,all.results.modelc$EntrezProt)))


write.table(matches[,c('CK_external_gene_id','gene_name','entrezgene','EntrezProt')],'../../Metabolic_models/Celegans/Fixed_matches.csv',sep='\t',row.names = FALSE,quote = FALSE)


write.table(all.results.modelc,'../../Metabolic_models/Celegans/Celegans_metformin_RNAseq.csv',sep='\t',row.names = FALSE,quote = FALSE)



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
  result<-pv.out[[pth]][['plot.data.gene']]
  result$PathwayID<-pth
  all.kegg.mappingst<-bind_rows(all.kegg.mappingst,result)
}
all.kegg.mappingst$mol.col<-NULL

all.kegg.mappings<-merge(path2pathde,all.kegg.mappingst,by='PathwayID',all.y=TRUE)




write.csv(all.kegg.mappings,paste(odir,'/All_KEGG_mappings.csv',sep = ''),row.names = FALSE)









#Venn
library('Vennerable')
#Might need to deactive reshape and Vennerable to use melt in other code

SM_S.up<-subset(all.results.rcp,`SM-S_FDR`<0.05 & `SM-S_logFC`>0)[,'t_id']
SM_S.down<-subset(all.results.rcp,`SM-S_FDR`<0.05 & `SM-S_logFC`<0)[,'t_id']
RM_R.up<-subset(all.results.rcp,`RM-R_FDR`<0.05 & `RM-R_logFC`>0)[,'t_id']
RM_R.down<-subset(all.results.rcp,`RM-R_FDR`<0.05 & `RM-R_logFC`<0)[,'t_id']

venn.groups<-list('RM-R.up'=RM_R.up,
                  'RM-R.down'=RM_R.down,
                  'SM-S.up'=SM_S.up,
                  'SM-S.down'=SM_S.down)

ups<-intersect(SM_S.down,RM_R.up)
downs<-intersect(SM_S.up,RM_R.down)

write.csv(subset(all.results.rcp,t_id %in% downs),paste(odir,'/Flipping_downregulated.csv',sep = ''),row.names = FALSE)
write.csv(subset(all.results.rcp,t_id %in% ups),paste(odir,'/Flipping_upregulated.csv',sep = ''),row.names = FALSE)



S.change<-subset(all.results.rcp,`SM-S_FDR`<0.05 )[,'t_id']
R.change<-subset(all.results.rcp,`RM-R_FDR`<0.05 )[,'t_id']
L.change<-subset(all.results.rcp,`SM-S-(RM-R)_FDR`<0.05 )[,'t_id']

venn.change<-list('SM-S'=S.change,
            'RM-R'=R.change)


subset(all.results.rcp,t_id %in% downs)[,'gene_name']

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


amp<-14
cbrks<-seq(-amp,amp,by=3.5)
gradcols<-c('blue2','blue2','gray80','red2','red2')



subset(all.results.rcp,gene_name=='acs-2')

ggplot(all.results.rcp,aes(x=`SM-S_logFC`,y=`RM-R_logFC`))+
  geom_vline(xintercept = 0,alpha=0.2)+
  geom_hline(yintercept = 0,alpha=0.2)+
  geom_abline(slope=1,intercept = 0,alpha=0.2)+
  geom_abline(slope = fit$coefficients[[2]],intercept = fit$coefficients[[1]],color='red',alpha=0.5)+
  geom_point(data=all.results.rcp,aes(size=abs(`SM-S-(RM-R)_logFC`),
                                      color=`SM-S-(RM-R)_logFC` ),stroke=0)+
  labs(color='Significant change')+
  xlab('Metformin effect when on OP50, logFC')+
  ylab('Metformin effect when on OP50MR, logFC')+
  scale_colour_gradientn(colors = gradcols,limits=c(-amp,amp),name='Longevity effect')+
  scale_size(range=c(0,8),breaks=seq(0,12,by=2),name='Amplitude of change')+
  scale_x_continuous(breaks=seq(-12,12,by=2),limits=c(-12,10))+
  scale_y_continuous(breaks=seq(-10,10,by=2),limits=c(-10,10))+
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




#Garbage
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

#lrt12$genes$entrezgene

keg12<-kegga(lrt12,geneid = 'entrezgene',coef='logFC',convert=TRUE,
             species.KEGG='cel',FDR=0.05,gene.pathway=gp.KEGG,trend='length')

topKEGG(keg12)

keg12$Comparison<-'SM/S'


#writre script to compare gene numbers from enrichment and in table


go12<-goana(lrt12,species='Ce',geneid='entrezgene',convert=TRUE,FDR=0.05)


View(topGO(go12,n=50))

#df residual needs to be 11!

# 
biocLite(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db"))
biocLite("pathview")
#biocLite("gage")



# source("http://bioconductor.org/biocLite.R")
# biocLite(c("pathview"))






#Pathview



gdata<-subset(all.results.rcp,Comparison=='SM-S' & Threshold=='None' & ! is.na(entrezgene))

nrow(all.results.rcp)
gdata<-subset(all.results.rcp,! is.na(entrezgene))


nrow(gdata)
o<-order(gdata[,'SM-S_FDR'],decreasing=FALSE)
gdata<-gdata[o,]


head(gdata)

mapkey<-'entrezgene'

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
path2eg$entrezgene<-path2eg$GeneID
head(path2eg)


dim(path2eg)

length(unique(all.results.kg$entrezgene))

length(intersect(unique(all.results.kg$entrezgene),path2eg$entrezgene))

head(path2eg)


allTCA<-subset(path2eg,PathID=='cel01200')
dim(allTCA)



