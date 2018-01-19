
#Ballgown
#Remove transcripts with low expression:
cel_filt = subset(cel,"rowVars(texpr(cel)) >1",genomesubset=TRUE)

#Get raw data

filt_table <- texpr(cel_filt , 'all')

ftm<-melt(filt_table,
          id.vars = c('t_id','chr','strand','start','end','t_name','num_exons','length','gene_id','gene_name'),
          variable.name = 'Var',
          value.name = 'Val')

ftm = transform(ftm, ID=colsplit(ftm$Var, "\\.", c('Stat', 'Sample')) )
ftm$Stat<-ftm$ID.Stat
ftm$Sample<-ftm$ID.Sample
ftm$ID.Stat<-NULL
ftm$ID.Sample<-NULL

head(ftm)

ftms<-dcast(ftm,t_id+chr+strand+start+end+t_name+num_exons+length+gene_id+gene_name+Sample~Stat,
            value.var = 'Val')

ftms$logFPKM<-log2(ftms$FPKM+1)

ftms$Group<-unlist(lapply(ftms$Sample,function(x) as.character(ifelse(nchar(x)==2,substr(x,1,1),substr(x,1,2)))))


flm<-dcast(ftms,Sample+Group~t_id,value.var = 'logFPKM',
           fill = as.numeric(NA),drop=TRUE)

head(flm[,c(1:20)])

miss<-apply(flm, 2, function(x) any(is.na(x)))
missing<-names(miss[miss==TRUE])

pca.group<-flm$Group

pca.dat<-flm[,!colnames(flm) %in% c('Sample','Group',missing)]

ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE,
                 na.action=na.omit) 
summary(ir.pca)
plot(ir.pca,type='l')
generalpca <- ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       groups = pca.group,
                       ellipse = TRUE,
                       circle = TRUE,
                       var.axes = 0)+
  scale_color_discrete(name = '')+
  theme(legend.direction = 'vertical',legend.position = 'right')
generalpca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_PCA.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#Figure out whether the gene_ids perfectly match gene names

#Strain comparison
#Get transcripts
strain.results_transcripts.t = stattest(cel_filt,
                                        feature="transcript",covariate="Strain",adjustvars =
                                          c("Treatment",'Batch'), getFC=TRUE, meas="FPKM")
strain.results_transcripts.n =
  data.frame(gene_name=ballgown::geneNames(cel_filt),
             geneIDs=ballgown::geneIDs(cel_filt), strain.results_transcripts.t)

#Get genes
strain.results_genes.t = stattest(cel_filt, feature="gene",
                                  covariate="Strain", adjustvars = c("Treatment",'Batch'), getFC=TRUE,
                                  meas="FPKM")

strain.results_genes.n<-merge(strain.results_genes.t,full_table[,c('gene_id','gene_name')],by.x='id',by.y='gene_id')

strain.results_genes<-Statprep(strain.results_genes.n)
strain.results_transcripts<-Statprep(strain.results_transcripts.n)

head(strain.results_transcripts)
head(strain.results_genes)


# custom model matrices:
### create example data:
# set.seed(43)
# strain = sample(c('M','F'), size=nrow(pData(cel_filt)), replace=TRUE)
# age = sample(21:52, size=nrow(pData(cel_filt)), replace=TRUE)
# 
# ### create design matrices:
# mod = model.matrix(~ sex + age)
# mod0 = model.matrix(~ pData(bg)$group + pData(bg)$time)
# 
# ### build model:
# adjusted_results = stattest(bg, feature='transcript', meas='FPKM',
#                             mod0=mod0, mod=mod)


#Treatment comparison
#Get transcripts
treat.results_transcripts.t = stattest(cel_filt,
                                       feature="transcript",covariate="Treatment",adjustvars =
                                         c("Strain",'Batch'), getFC=TRUE, meas="FPKM")
treat.results_transcripts.n =
  data.frame(gene_name=ballgown::geneNames(cel_filt),
             geneIDs=ballgown::geneIDs(cel_filt), treat.results_transcripts.t)

#Get genes
treat.results_genes.t = stattest(cel_filt, feature="gene",
                                 covariate="Treatment", adjustvars = c("Strain",'Batch'), getFC=TRUE,
                                 meas="FPKM")

treat.results_genes.n<-merge(treat.results_genes.t,full_table[,c('gene_id','gene_name')],by.x='id',by.y='gene_id')


treat.results_genes<-Statprep(treat.results_genes.n)
treat.results_transcripts<-Statprep(treat.results_transcripts.n)

head(treat.results_transcripts)
head(treat.results_genes)

dim(subset(treat.results_genes,qval<0.05 & abs(logFC>1)))



hist(treat.results_genes$qval, main='Ballgown q-values: treat comparison', col="grey",
     xlab='Range of q-values for the more significant transcripts')




#Comparison made Treatment/Control
subset(full_table,gene_id==treat.results_genes[1,'id'])

subset(full_table,gene_id==strain.results_genes[1,'id'])


`CM = subset(cel_filt, 'Strain == 'C'', genomesubset=FALSE)`

CM = subset(cel_filt,"texpr(cel)Strain ==C",genomesubset=TRUE)


fpkm = texpr(cel,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$Treatment),las=2,ylab='log2(FPKM+1)')

ballgown::geneNames(cel_filt)[12]


ggplot(treat.results_genes,aes(x=logFC,y=logqval))+
  geom_vline(xintercept = c(-1,1),color='red')+
  geom_hline(yintercept =-log10(0.05),color='red')+
  geom_point(color='grey20')+
  xlim(-7,7)+
  xlab('logFC Treatment')+
  ylab('-log10(q)')+
  geom_text(aes(label=ifelse(logqval>-log10(0.05) & abs(logFC)>1,as.character(gene_name),'')),hjust=1, vjust=-0.5,size=3,colour = "red")

ggplot(treat.results_transcripts,aes(x=logFC,y=logqval))+
  geom_vline(xintercept = c(-1,1),color='red')+
  geom_hline(yintercept =-log10(0.05),color='red')+
  geom_point(color='grey20')+
  xlim(-7,7)+
  xlab('logFC Treatment')+
  ylab('-log10(q)')+
  geom_text(aes(label=ifelse(logqval>-log10(0.05) & abs(logFC)>1,as.character(gene_name),'')),hjust=1, vjust=-0.5,size=3,colour = "red")




treat.g.sign<-subset(treat.results_genes,abs(logFC)>1 & qval<0.05)
strain.g.sign<-subset(strain.results_genes,abs(logFC)>1 & qval<0.05)

treat.t.sign<-subset(treat.results_transcripts,abs(logFC)>1 & qval<0.05)
strain.t.sign<-subset(strain.results_transcripts,abs(logFC)>1 & qval<0.05)

head(treat.t.sign,n=100)





genenr<-1
plotTranscripts(ballgown::geneIDs(cel_filt)[genenr],
                cel_filt, main=c(paste('Gene',ballgown::geneIDs(cel_filt)[genenr],'in sample C1')), sample=c('C1'))


head(treat.results_transcripts,n=50)
head(treat.results_genes,n=50)


