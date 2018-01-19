library(xlsx)
library(plyr)
library(reshape2)
library(tidyr)
library(ggplot2)

library(gplots)
library(gtools)
library(rgl)
library(cluster)
library(heatmap3)

library(multcomp)
library(contrast)
library(RColorBrewer)

library(grid)
library(gridExtra)
library(ggrepel)


devtools::install_github("PNorvaisas/PFun")
library(PFun)

theme_set(theme_light())

getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}


mycolsplit<-function(column,cols,sep=':'){
  elems<-unlist(strsplit(column,paste("([\\",sep,"])",sep=''), perl=TRUE))
  return(as.data.frame(matrix(elems,ncol=cols,byrow=TRUE)))
}


MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}

mymerge<-function(all.results,results) {
  if (dim(all.results)[[1]]==0) {
    all.results<-results
  } else {
    all.results<-merge(all.results,results,all.x=TRUE,all.y=TRUE)
  }
  return(all.results)
}



library(ggthemes)

theme_Publication <- function(base_size=14) {
  
  (theme_foundation(base_size=base_size)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           #panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           #legend.position = "bottom",
           #legend.direction = "horizontal",
           #legend.key.size= unit(0.2, "cm"),
           #legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           #plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

theme_set(theme_Publication())


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir<-'Summary_Filipe_D1_D5_D8'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



rawdata<-read.xlsx2('Filipe_bac_growth/Bacterial_growth_curve_D1_D5_D8.xlsx',
                    sheetName = 'Clean',header = TRUE,check.names = FALSE,
                    colClasses = c(rep('character',3),rep('numeric',217)))
head(rawdata)
dim(rawdata)


colnames(rawdata)[4:24]

rawdata.adj<-rawdata

rawdata.adj[,4:220]<-rawdata[,4:220]-apply(rawdata[,4:24],1,mean)

rawdata[,1:20]
rawdata.adj[,1:24]


data.r<-melt(rawdata,id.vars = c('Strain','Day','Replicate'),variable.name = 'Time_min',value.name='OD_raw')
data.r$Time_min<-as.numeric(as.character(data.r$Time_min))
data.r$Time_h<-data.r$Time_min/60




baseline<-ddply(subset(data.r,Time_min<=100),.(Strain,Day,Replicate),summarise,Baseline=mean(OD_raw))
data<-merge(data.r,baseline,by=c('Day','Strain','Replicate'),all.x=TRUE)
data$OD<-data$OD_raw-data$Baseline


head(data)


ggplot(data,aes(x=Time_h,y=OD))+
  geom_line(aes(group=interaction(Replicate)))+
  xlab('Time, h')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_grid(Strain~Day,labeller = labeller(.rows = label_both, .cols = label_both))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=9,height=6)



data.sum<-ddply(data,.(Day,Strain,Time_min,Time_h),summarise,OD_Mean=mean(OD),OD_SD=sd(OD))

head(data.sum)

ggplot(data.sum,aes(x=Time_h,y=OD_Mean,color=Strain,fill=Strain))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5)+
  xlab('Time, h')+
  ylab('OD')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_grid(Day~.,labeller=labeller(.rows = label_both))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_by_day.pdf",sep=''),
             width=6,height=6)






data.int<-ddply(data,.(Day,Strain,Replicate),summarise,OD_Int=sum(OD)*5/60)
data.int$logOD_Int<-log2(data.int$OD_Int)

data.int$Sample<-paste(data.int$Strain,data.int$Day,sep='_')
data.int$Sample<-gsub('-','',data.int$Sample)

head(data.int)

sel.samples<-unique(data.int$Sample)




contrasts<-read.contrasts('!Contrasts_Filipe_Bac.xlsx','Contrasts_values',sel.samples,variables = c('Day'))

contrasts.table<-contrasts$Contrasts.table
contr.matrix<-contrasts$Contrasts.matrix
contr.matrix


contrasts.table



allresults<-hypothesise(data.int,c('OD_Int','logOD_Int'),contr.matrix,formula='0+Sample')$All


results<-merge(contrasts.table[,c('Contrast','Description','Contrast_type','Strain','Day')],allresults,by='Contrast',all.x=TRUE)

results$pStars<-stars.pval(results$p.value)
results$pStars<-gsub('\\.','',as.character(results$pStars))
results$pStars<-gsub(' ','',as.character(results$pStars))
results$Day<-as.factor(results$Day)

head(results)


write.csv(results,paste(odir,'/Growth_results.csv',sep=''),row.names = FALSE)


results.m<-melt(results,id.vars = c('Contrast','Contrast_type','Description','Strain','Day','Variable'),
                measure.vars = c('logFC','p.value','SE','t.value','PE','NE','FDR','logFDR'),
                variable.name = 'Stat',value.name = 'Value')


results.castfull<-dcast(subset(results.m,Stat %in% c('logFC','SE','FDR','NE','PE','logFDR') ),Variable~Contrast+Stat,value.var = 'Value')
results.cast<-dcast(subset(results.m,Stat %in% c('logFC','FDR')),Variable~Contrast+Stat,value.var = 'Value')

head(results.cast)

head(results)

dwidth<-0.5

ODstars_day<-subset(results,Variable=='OD_Int' & Contrast_type=='Strain')[,c('Day','pStars')]
ODstars_day$Strain<-'N2'

head(data.int)

ggplot(data.int, aes(x=Day,y=OD_Int,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('OD integral, OD*h')+
  # geom_point(shape = 21,position =
  #              position_jitterdodge(jitter.width = 0.2, jitter.height=0,
  #                                   dodge.width= dwidth))+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=ODstars_day,aes(x=Day,y=4.5,label=as.character(pStars)))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_OD_int_by_day.pdf",sep=''),
             width=6,height=4)


logODstars_day<-subset(results,Variable=='logOD_Int' & Contrast_type=='Strain')[,c('Day','pStars')]
logODstars_day$Strain<-'N2'


ggplot(data.int, aes(x=Day,y=logOD_Int,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('log2 OD integral, log2(OD*h)')+
  # geom_point(shape = 21,position =
  #              position_jitterdodge(jitter.width = 0.2, jitter.height=0,
  #                                   dodge.width= dwidth))+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=logODstars_day,aes(x=Day,y=2.25,label=as.character(pStars)))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_logOD_int_by_day.pdf",sep=''),
             width=6,height=4)



#By strain

ODstars_strain<-subset(results,Variable=='OD_Int' & Contrast_type=='Time' & !Contrast %in% c('N2_8-5','Phm2_8-5') )
ODstars_strain
head(data.int)

ggplot(data.int, aes(x=Day,y=OD_Int,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('OD integral, OD*h')+
  # geom_point(shape = 21,position =
  #              position_jitterdodge(jitter.width = 0.2, jitter.height=0,
  #                                   dodge.width= dwidth))+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=ODstars_strain,aes(x=Day,y=4.5,label=as.character(pStars)))+
  facet_grid(~Strain)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_OD_int_by_strain.pdf",sep=''),
             width=6,height=4)





logODstars_strain<-subset(results,Variable=='OD_Int' & Contrast_type=='Time' & !Contrast %in% c('N2_8-5','Phm2_8-5') )


ggplot(data.int, aes(x=Day,y=logOD_Int,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('log2 OD integral, log2(OD*h)')+
  # geom_point(shape = 21,position =
  #              position_jitterdodge(jitter.width = 0.2, jitter.height=0,
  #                                   dodge.width= dwidth))+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=logODstars_strain,aes(x=Day,y=2.25,label=as.character(pStars)))+
  facet_grid(~Strain)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_logOD_int_by_strain.pdf",sep=''),
             width=6,height=4)



