library(gplots)
library(ggplot2)
library(xlsx)
library(plyr)
library(reshape2)
library(ggbiplot)
library(ggrepel)
library(car)
library(heatmap3)

library(RColorBrewer)
library(grid)
library(gridExtra)

library(plot3D)

# library(rgl)
# library(pca3d)

library(multcomp)
library(contrast)

theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Metabolomics/")



odir<-'Summary_nucleotides'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

read.and.clean<-function(fname,sheet='') {
  print(fname)
  #print(grepl('.xlsx',fname))
  #print(grepl('.csv',fname))
  if (grepl('xlsx',fname)) {
    dft<-read.xlsx2(fname,sheetName = sheet,
                    stringsAsFactors = FALSE,
                    header=FALSE)
  } else if (grepl('csv',fname)) {
    dft<-read.csv(fname,sep=',',quote = "\"",stringsAsFactors = FALSE,
                  header=FALSE)
    dft[1,which(dft[1,] == "")] <- "Metabolite Set"
  }
  #print(length(colnames(dft)))
  
  colnames(dft)<-dft[1,]
  df<-dft[-1,]
  
  return(df)
}



nuc.raw<-read.and.clean('Ecoli_nucleotides.xlsx','Raw')
annot<-read.and.clean('Nucleotide_annotation.xlsx','Annotation')

nuc.m<-melt(nuc.raw,id.vars = 'Metabolite',variable.name = 'Sample',value.name = 'Conc')


nuc<-merge(nuc.m,annot,by.x = 'Metabolite',by.y='Compound',all.x = TRUE)




nuc$Metabolite<-trimws(nuc$Metabolite)

nuc$Conc<-as.numeric(nuc$Conc) 
nuc$logConc<-log(nuc$Conc,2)
nuc$Replicate<-as.numeric(gsub("\\D", "", nuc$Sample))
nuc$Type<-as.character(gsub("[[:digit:]]", "", nuc$Sample))
nuc$Class<-as.factor(nuc$Class)
nuc$Type<-factor(nuc$Type,levels=c('C','CM','R','RM'),labels=c('C','CM','R','RM'))
nuc$Metabolite<-as.factor(nuc$Metabolite)


metslm<-dcast(nuc,Sample+Type+Replicate~Metabolite,value.var = c('logConc'),fill = as.numeric(NA),drop=TRUE)
cols<-setdiff(colnames(metslm),c('Sample','Type','Replicate'))



#Find compounds with missing values
miss<-apply(metslm, 2, function(x) any(is.na(x)))
missing<-names(miss[miss==TRUE])

missing

pca.group<-metslm$Type

pca.dat<-metslm[,cols]

ir.pca <- prcomp(pca.dat,
                 center = TRUE,
                 scale. = TRUE) 
summary(ir.pca)
plot(ir.pca,type='l')

generalpca <- ggbiplot(ir.pca, obs.scale = 1,
                       var.scale = 1,
                       groups = pca.group,
                       ellipse = TRUE,
                       circle = TRUE,
                       var.axes = 0)+
  theme(legend.direction = 'vertical',legend.position = 'right')
generalpca

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA.pdf",sep=''),
             width=12,height=9)

allmets<-unique(as.character(nuc$Metabolite))

cleanmets<-setdiff(allmets,missing)






#Use data with deviant growth removed
lmshape<-dcast(nuc,Sample+Type+Replicate+Class+Group~Metabolite,value.var = 'logConc',drop=TRUE)

contrasts<-list('CM-C'=c(-1,1,0,0),
                  'RM-R'=c(0,0,-1,1),
                  'R-C'=c(-1,0,1,0),
                  'T-C_general'=c(-1,1,-1,1),
                  'R-C_general'=c(-1,-1,1,1),
                  'CM-C-(RM-R)'=c(-1,1,1,-1),
                  'C'=c(1,0,0,0),
                  'CM'=c(0,1,0,0),
                  'R'=c(0,0,1,0),
                  'RM'=c(0,0,0,1))

contrasts
#comparisons<-names(contrasts)

descriptions<-list('CM-C'='Treatment effect on OP50',
                   'RM-R'='Treatment effect on OP50MR',
                   'R-C'='OP50MR over OP50 in control',
                   'T-C_general'='General Treatment over Control',
                   'R-C_general'='General OP50MR over OP50',
                   'CM-C-(RM-R)'='Longevity effect',
                   'C'='OP50 control',
                   'CM'='OP50 treatment',
                   'R'='OP50MR control',
                   'RM'='OP50MR treatment')

cont<-ldply(contrasts)
rownames(cont)<-cont$.id
cont$.id<-NULL
colnames(cont)<-c('C','CM','R','RM')
cont
#contr<-as.matrix(cont)

contr.matrix<-as.matrix(cont[,c('C','CM','R','RM')])

contr.matrix


desc<-ldply(descriptions[rownames(cont)])
rownames(desc)<-desc$.id
desc$.id<-NULL
colnames(desc)<-c('Description')

desc

design<-cbind(desc,cont)
design

write.csv(design,paste(odir,'/!Design.csv',sep=''),row.names = TRUE)




allresults.t<-data.frame()
prec<-0
met.len<-length(cleanmets)

for (metid in 1:length(cleanmets)) {
  met<-as.character(cleanmets[metid])
  precn<-metid*100/length(cleanmets)
  if (precn-prec>5) {
    print(paste(round(precn,digits=0),'%',sep=''))
    prec<-precn
  }
  model<-lm(paste("`",met,"`~0+Type",sep=""),lmshape)
  lmod_glht <- glht(model, linfct = contr.matrix)
  result<-summary(lmod_glht,test=adjusted("none"))
  res<-ldply(result$test[c('coefficients','sigma','tstat','pvalues')])
  res$Metabolite<-met
  allresults.t<-rbind(allresults.t,res)
}

result




head(allresults.t)
allresults.m<-melt(allresults.t,id.vars = c('.id','Metabolite'),variable.name = 'Comparison',value.name = 'Value')
allresults<-dcast(allresults.m,Metabolite+Comparison~`.id`,value.var = 'Value')




results<-rename(allresults,c('coefficients'='logFC','sigma'='SE','pvalues'='p.value','tstat'='t.value'))
results$FDR<-p.adjust(results$p.value,method = 'fdr')
results$SD<-results$SE*sqrt(4)
results$PE<-results$logFC+results$SE
results$NE<-results$logFC-results$SE
results$PD<-results$logFC+results$SD
results$ND<-results$logFC-results$SD
results$logFDR<--log10(results$FDR)

results.exp<-results[,c('Comparison','Metabolite','logFC','SE','t.value','p.value','FDR')]

head(results)
write.csv(results.exp,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)



results.m<-melt(results,id.vars = c('Comparison','Metabolite'),variable.name = 'Stats',value.name = 'Value')
results.cast<-dcast(results.m,Metabolite~Comparison+Stats,value.var = 'Value')

results.exp<-subset(results.m,Stats %in% c('logFC','FDR'))
results.castexp<-dcast(results.exp,Metabolite~Comparison+Stats,value.var = 'Value')

head(results.cast)
head(results.castexp)

write.csv(results.castexp,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)






erralpha<-1
errcolor<-'grey80'

brks<-seq(0,8,by=1)
gradcols<-c('black','purple','purple')


ggplot(subset(results,!Comparison %in% c('C','CM','R','RM')),aes(x=logFC,y=logFDR))+
  geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=erralpha,color=errcolor,height=0)+
  geom_point()+
  ggtitle('Metabolic changes in OP50 (C) and OP50MR (R) treated with metformin')+
  geom_text(aes(label=ifelse(FDR <0.05,as.character(Metabolite),'')),nudge_y = 0.1,size=2)+
  facet_wrap(~Comparison)

dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Volcano_all.pdf',sep = ''),
             width=12,height=9,useDingbats=FALSE)




xs<-c('CM-C','R-C','RM-R','CM-C-(RM-R)',
      'T-C_general','R-C_general')


maincomp<-'logFC'

for (x in xs){
  xvar<-paste('`',x,'_logFC`',sep='')
  yvar<-paste('`',x,'_logFDR`',sep='')
  xPE<-paste('`',x,'_PE`',sep='')
  xNE<-paste('`',x,'_NE`',sep='')
  
  volcano<-ggplot(results.cast,aes_string(x=xvar,y=yvar))+
    geom_hline(yintercept = -log10(0.05),color='red',alpha=0.5,linetype='longdash')+
    # geom_vline(xintercept = 1,color='red',alpha=0.5,linetype='longdash')+
    # geom_vline(xintercept = -1,color='red',alpha=0.5,linetype='longdash')+
    geom_errorbarh(aes_string(xmin=xNE,xmax=xPE),
                   alpha=erralpha,color=errcolor,height=0)+
    geom_point(aes(size=abs(eval(parse(text = xvar))),
                   color=abs(eval(parse(text = xvar)))))+
    scale_x_continuous(breaks=seq(-14,14,by=1),limits=c(-8,8))+
    scale_y_continuous(breaks=seq(0,7,by=1),limits=c(0,8))+
    scale_colour_gradientn(colours = gradcols,
                           breaks=brks,limits=c(0,7),name=maincomp)+
    scale_size(range = c(0.25, 4),name=maincomp)+
    geom_text(aes(label=ifelse(eval(parse(text = yvar)) > -log10(0.05),as.character(Metabolite),'')),nudge_y = 0.2,size=2)+
    # geom_text_repel(aes(label=ifelse(eval(parse(text = yvar)) > -log10(0.05),
    #                                  as.character(Metabolite),'')),
    #                 size=4,
    #                 nudge_y = 0.3,
    #                 segment.colour = errcolor)+
    coord_cartesian(ylim=c(0, 7),xlim=c(-7,7))
  
  fname<-paste(odir,"/Volcano_",x,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=9,height=6)
  print(volcano)
  dev.off()
}



subset(results,Metabolite=='cAMP')





