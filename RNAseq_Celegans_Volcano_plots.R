library(tidyverse)
library(PFun)
library(broom)
library(viridis) #Nice colors
library(readxl)
library(ggrepel)


library(PFun)

library(edgeR)


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau



cnames<-c('SM-S'='OP50-C treatment','RM-R'='OP50-MR treatment','R-S'='Strain difference','SM-S-(RM-R)'='Longevity effect')



setwd("~/Dropbox/Projects/2015-Metformin/RNAseq/Celegans_metformin/Results_1thrs_newannot")

data <- read_csv('All_results.csv') %>%
  filter(Comparison %in% c('SM-S','RM-R','R-S','SM-S-(RM-R)')) %>%
  group_by(Comparison) %>%
  arrange(FDR) %>%
  mutate(logFDR=-log10(FDR),
         Description=cnames[Comparison],
         Description=factor(Description,levels=as.vector(cnames)),
         Name=ifelse(row_number()<=10 & FDR<0.05,gene_name,NA ) )


unique(data$Comparison)

data %>%
  ggplot(aes(x=logFC,y=logFDR))+
  geom_hline(yintercept = -log10(0.05),color='red4',alpha=0.5)+
  geom_point(size=1,color='gray50' )+
  geom_text_repel(aes(label=Name),color='red2')+
  ylab('-log10(FDR)')+
  facet_wrap(~Description,ncol = 2)

ggsave(file='Volcano_plot.pdf',
       width=100,height=100,units='mm',scale=2,device=cairo_pdf,family="Arial")





pheno_data<-read_csv("../Pheno_data.csv")


expression<-read_csv('Raw_data_for_genes.csv')


#Split table for edgeR

ancols.2<-setdiff(colnames(expression),c(as.character(pheno_data$ids)) )

gene.info<-expression[,ancols.2]
gene.data<-expression[,as.character(pheno_data$ids)]


head(gene.info)
head(gene.data)


dim(gene.data)






dt<-DGEList(counts=gene.data,
            genes=gene.info,
            group=factor(pheno_data$Group))



groups<-pheno_data$Group



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



#Batch removal
logCPM <- cpm(d, log=TRUE, prior.count=1)
# batch<-c('1','1','2','2',
#          '1','2','2','2',
#          '1','2','2','2',
#          '1','1','2','2')
logCPMc <- removeBatchEffect(logCPM,as.character(pheno_data$Batch) )

#MDS

MDSdata<-plotMDS(logCPMc, col=as.numeric(d$samples$group),plot=TRUE)
legend("bottomleft", as.character(unique(d$samples$group)), col=1:4, pch=1)




#Heatmap


strains<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$id)),'Strain'])
metf<-as.character(pheno_data[match(colnames(logCPMc),as.character(pheno_data$id)),'Treatment'])


strains.fac<-factor(strains)
metf.fac<-factor(metf)

strains.col<-ifelse(strains.fac=='C','white','black')
metf.col<-ifelse(metf.fac=='Control','green','red')

clab<-cbind(metf.col,strains.col)
colnames(clab)<-c('Metformin, mM','Strain')


library(heatmap3)
#Heatmap adjusted
heatmap3(as.matrix(logCPMc),
         scale = 'row',
         #ColSideLabs = 'Group',
         balanceColor=TRUE,
         ColSideColors=clab,
         labRow=FALSE)


legend('topright',legend=c('OP50','OP50-MR'),fill=c('red3','blue3'), border=TRUE, bty="n",title='Strain')
legend('right',legend=c('0','50'),fill=c('white','black'), border=FALSE, bty="n",title='Metformin, mM')



#

bioinfo<-pheno_data %>%
  mutate(Metformin_mM=ifelse(Treatment=='Control',0,50) %>% as.factor,
         Strain=ifelse(Strain=='C','OP50-C','OP50-MR')%>% as.factor) %>%
  data.frame

rownames(bioinfo)<-bioinfo$ids
bioinfo$ids<-NULL

hinfo<-bioinfo %>%
  select(Strain, Metformin_mM)


hinfo

CPMdata<-logCPMc %>%
  data.frame %>%
  mutate(ID=row_number()) %>%
  gather(ids,logCPM,everything(),-ID) %>%
  left_join(bioinfo %>%
              mutate(ids = rownames(bioinfo)) ) 


bioinfo %>%
  mutate(ids = rownames(bioinfo))

glimpse(CPMdata)



cols<-list(Metformin_mM = c('0' = "white", '50' = "black"),
           Strain=c("OP50-C"="red","OP50-MR"="blue") )

#Complete
CPMdata %>%
  HMap('ids','ID','logCPM',hinfo,cols=cols,scalesel = 'row',cmc = 'ward.D2',cmr='ward.D2')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/RNAseq_Heatmap_logCPMc_batchadj.pdf",sep=''),
             width=9,height=9, useDingbats=FALSE)





