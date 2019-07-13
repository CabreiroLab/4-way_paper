library(gplots)
library(ggplot2)
library(xlsx)

library(plyr)
library(dplyr)


library(reshape2)
library(ggbiplot)
library(ggrepel)
library(car)
library(heatmap3)
library(gtools)

library(RColorBrewer)
library(grid)
library(gridExtra)

library(plot3D)

# library(rgl)
# library(pca3d)

library(multcomp)
library(contrast)
library(scales)
library(Rtsne)
library(ggrepel)




devtools::install_github("PNorvaisas/PFun")
library(PFun)



MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
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



#theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))





setwd("~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/Celegans/")


load('Celegans.RData')

#save.image("Celegans.RData")

odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



info<-read.table('../../Biolog/Biolog_metabolites_EcoCyc.csv',sep=',',quote = '"',header = TRUE)
info<-info %>% filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))


dupmet<-data.frame(table(trimws(as.character(info$Metabolite))))
dupls<-dupmet[dupmet$Freq>1,'Var1']
dupls

info.b<-info %>%
  mutate(Metabolite=trimws(as.character(Metabolite)),
         Name=trimws(as.character(Name)),
         Metabolite_Uniq=ifelse(Metabolite=='Negative Control','NGM',Metabolite),
         Metabolite_Uniq=ifelse(Metabolite %in% dupls & Plate %in% c('PM3B'),
                                 paste(Metabolite,'N',sep='_'),Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='P-Source',
                                      paste(Metabolite,'P',sep='_'),Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Metabolite %in% dupls & Plate %in% c('PM4A') & Group=='S-Source',
                                 paste(Metabolite,'S',sep='_'),Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Index=='PM1-A1','NGM_C1',Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Index=='PM2A-A1','NGM_C2',Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Index=='PM3B-A1','NGM_N',Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Index=='PM4A-A1','NGM_P',Metabolite_Uniq),
         Metabolite_Uniq=ifelse(Index=='PM4A-F1','NGM_S',Metabolite_Uniq))
  

subset(info.b,Metabolite %in% dupls)

head(info)
subset(info.b,Metabolite=='Negative Control')
head(info)

dupmet<-data.frame(table(as.character(info.b$Metabolite_Uniq)))

dupls<-dupmet[dupmet$Freq>1,'Var1']
dupls

subset(info.b,Metabolite_Uniq %in% dupls)


ref.cA<-data.frame(Plate='Controls',Well=paste0('A',seq(1,8)),Index='Controls-A1',Group='Controls',Metabolite='NGM',Metabolite_Uniq='NGM')
ref.cB<-data.frame(Plate='Controls',Well=paste0('B',seq(1,8)),Index='Controls-B1',Group='Controls',Metabolite='Negative Control',Metabolite_Uniq='NGM_C')

ref.c<-rbind(ref.cA,ref.cB)


info<-merge(info.b,ref.c,all=TRUE)

View(info)



artefacts<-read.xlsx2('Artefacts_new.xlsx',sheetName = 'Artefacts')


data.raw<-read.csv('Summary_all_new_auto_threshold.csv',sep='\t',header=TRUE,check.names = FALSE)



data.ctr<-read.csv('Summary_controls_All_auto_threshold.csv',sep='\t',header=TRUE,check.names = FALSE)


unique(as.character(data.ctr$Well))
data.ctr<-data.ctr %>%
  mutate(Well=ifelse(grepl('NoMetf',File),
                     paste0('A',as.numeric(gsub('[[:alpha:]]|_|\\.','',File))),
                     paste0('B',as.numeric(gsub('[[:alpha:]]|_|\\.','',File)))),
         Plate='Controls') %>%
  filter(!Replicate %in% c('Rep1','Rep2','Rep3'))



data<-rbind(data.raw,data.ctr)


for (a in 1:nrow(artefacts)){
  data<-subset(data,!(Replicate==paste('Rep',as.character(artefacts[a,]$Replicate),sep='') &
                        Plate==as.character(artefacts[a,]$Plate) &
                        Well==as.character(artefacts[a,]$Well) &
                        Worm==as.character(artefacts[a,]$Worm)) )
}






subset(data,Plate=='Controls')[,c('Replicate','Plate','Well','Worm')]


table(data$Replicate)


head(data)


#colnames(data)<-gsub('X','',colnames(data))
#Remove zero pixels

data.c<-data[,!grepl('0.000',colnames(data))]



#All brightness levels
brights<-setdiff(colnames(data.c),c('Replicate','Plate','Well','File','Worm'))



data.ca<-data.c %>%
  merge(info[,c('Plate','Well','Index','Metabolite','Metabolite_Uniq','Group')],by=c('Plate','Well'),all.x=TRUE) %>%
  mutate(Row=as.character(gsub("[[:digit:]]", "", Well)),
         Col=as.numeric(gsub("[[:alpha:]]", "", Well)))


head(data.ca)

data.ca$Index

subset(data.ca,is.na(Index))[,c('Replicate','Plate','Well','Index','File','Worm','Metabolite','Metabolite_Uniq')]

subset(data.ca,grepl('NGM',File))[,c('Replicate','Plate','Well','Index','File','Worm','Metabolite','Metabolite_Uniq')]

table(data.ca$Metabolite_Uniq)



#Melt to long table
data.l<-data.ca %>%
  melt(measure.vars = brights,variable.name='Brightness',value.name = 'Frequency') %>%
  mutate(Brightness=as.numeric(as.character(Brightness)),
         LogBrightness=log2(Brightness),
         Brightness_sum=cut(Brightness, breaks=seq(0,1.0,0.01),labels=seq(0.01,1.0,0.01) ),
         LogBrightness_sum=cut(LogBrightness, breaks=seq(-6,0,0.1),labels=seq(-5.9,0,0.1)) ) %>%
  filter(LogBrightness>= -6) 



subset(data.l,Plate=='Controls')

min(data.l$Brightness)



#Make explicit long table from Freq table
#Can take some time
data.m<-data.l %>%
  ddply(.(Replicate,Plate,Well,Row,Col,Index,File,Worm,Metabolite,Metabolite_Uniq,Group),plyr::summarise,Brightness=rep(Brightness,Frequency) ) %>%
  mutate(BrightnessLog=log2(Brightness))

dim(data.m)
head(data.m)


table(data.m$Replicate)



ggplot(data.m,aes(x=BrightnessLog))+
  geom_histogram()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Overall_logBrightness.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



#Make summary with quantiles

rows<-c('A','B','C','D','E','F','G','H')
summary.a<-ddply(data.m,.(Replicate,Plate,Well,Index,Row,Col,File,Worm,Metabolite,Metabolite_Uniq,Group),summarise,
               Mean=mean(BrightnessLog),
               SD=sd(BrightnessLog),
               Median=median(BrightnessLog),
               Sum=sum(BrightnessLog),
               N=length(BrightnessLog),
               Q80=quantile(BrightnessLog,0.8),
               Q90=quantile(BrightnessLog,0.9),
               Q95=quantile(BrightnessLog,0.95)) %>%
  mutate(Row=factor(Row,levels=rows,labels=rows),
         Col=as.factor(Col),
         Metabolite=as.factor(Metabolite),
         Metabolite=relevel(Metabolite,ref='Negative Control'),
         Metabolite_Uniq=as.factor(Metabolite_Uniq),
         Worm=as.factor(Worm))



#Summary over all measures
summary.m<-melt(summary.a,id.vars = c('Replicate','Plate','Well','Index','Row','Col','File','Worm','Metabolite','Metabolite_Uniq','Group','N'),
                measure.vars = c('Mean','SD','Sum','Median','N','Q80','Q90','Q95'),
                variable.name = 'Stat',value.name = 'Raw')


ref<-subset(summary.m,Metabolite=='Negative Control') #Bad controls
ref<-plyr::rename(ref,c('Raw'='Ref'))




#Summarise over all worms
ref.rp<-ddply(ref,.(Replicate,Plate,Group,Stat),plyr::summarise,Ref=mean(Ref))

ref.rp


head(ref.rp)
subset(ref.rp,Replicate=='Rep5')


refstat<-subset(summary.m, Stat=='Q90' & Metabolite %in% c('Negative Control','NGM'))

ggplot(refstat,aes(x=Plate,y=Raw))+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
  geom_point(aes(color=Replicate))+
  ylab('log2 Brightness')+
  ylim(-3,0)+
  ggtitle('Negative control log2 Brightness')+
  facet_grid(.~Metabolite)



dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Negative_Controls_Allplates_Q90.pdf",sep=''),
             width=18,height=6)



#Manual scoring
scores<-read.xlsx2(file='Scores/Scores.xlsx',sheetName = 'Scores')
scores[scores=='']<- NA
scores.m<-melt(scores,measure.vars=c('Rep1','Rep2','Rep3'),variable.name='Replicate',value.name = 'Score')



#Use reference for each plate and replicate

wormsize<-5314

summary.norm<-summary.m %>%
  merge(ref.rp[,c('Plate','Replicate','Group','Stat','Ref')],
        by=c('Plate','Replicate','Group','Stat'),all.x=TRUE) %>%
  mutate(Norm=Raw-Ref,
         Worms=round(N/wormsize),
         Worms=ifelse(Worms==0,1,Worms)) %>%
  melt(id.vars = c('Replicate','Plate','Well','Index','Row','Col','File','Worm','Metabolite','Metabolite_Uniq','Group','N','Worms','Stat'),
        variable.name = 'Type',value.name = 'Value')
  
summary.cast<-summary.norm %>% 
  dcast(Replicate+Plate+Well+Index+Row+Col+File+Worm+Metabolite+Metabolite_Uniq+Group+N+Worms~Type+Stat,
        value.var = 'Value',drop = TRUE) %>%
  merge(scores.m[,c('Replicate','Plate','Well','Score')],by=c('Replicate','Plate','Well'),all.x=TRUE)
  



#Get approximate worm counts
wormsize<-median(subset(summary.cast,N<8500)$N )
  
  

# 
# summary.ref<-merge(summary.m,ref.rp[,c('Plate','Replicate','Group','Stat','Ref')],
#                    by=c('Plate','Replicate','Group','Stat'),all.x=TRUE)
# #Subtract reference measure
# summary.ref$Norm<-summary.ref$Raw-summary.ref$Ref
# 
# 
# summary.refm<-melt(summary.ref,
#                    id.vars = c('Replicate','Plate','Well','Index','Row','Col','File','Worm','Metabolite','Metabolite_Uniq','Group','N','Stat'),
#                    variable.name = 'Type',value.name = 'Value')
# 
# table(summary.refm$Replicate)
# head(summary.refm)
# 
# summary.cast.t<-dcast(summary.refm,Replicate+Plate+Well+Index+Row+Col+File+Worm+Metabolite+Metabolite_Uniq+Group+N~Type+Stat,
#                       value.var = 'Value',drop = TRUE)
# 
# 
# table(summary.cast.t$Replicate)
# 
# summary.cast<-summary.cast.t %>%
#   merge(scores.m[,c('Replicate','Plate','Well','Score')],by=c('Replicate','Plate','Well'),all.x=TRUE) %>%
#   mutate(Worms=round(N/wormsize),
#          Worms=ifelse(Worms==0,1,Worms))
         

summary.cast$UniqRep=makereplicates(summary.cast[,c('Replicate','Plate','Index')])






min(summary.cast$Norm_Q90)
max(summary.cast$Norm_Q90)

min(summary.cast$Norm_Q90)-max(summary.cast$Norm_Q90)

head(summary.m)



#Worm size distribution
ggplot(summary.cast,aes(x=N))+
  geom_histogram(binwidth=500)+
  scale_x_continuous(breaks=seq(0,30000,by=1000))+
  xlab('Pixels, N')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Worm_size_distribution.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)









dim(summary.cast)

head(summary.cast)


#View(summary.cast)


subset(summary.cast,Plate=='Controls')


write.csv(summary.cast,paste(odir,"/Celegans_raw_summary.csv",sep=''),row.names = FALSE)


#Linear modelling
#library(mblm)

#Generate template to overlay data on

head(summary.cast)
usestat<-'Norm_Q90'

cleanmets<-c(unique(as.character(info$Index)))




#Calculate SE
#se <- function(x) sd(x)/sqrt(length(x))

summary.sum<-summary.norm %>%
  filter(!(Replicate=='Rep1' & (Plate=='PM1'| Plate=='PM2A'))  ) %>%
  ddply(.(Plate,Index,Metabolite,Metabolite_Uniq,Group,Stat,Type),summarise,
        Count=length(Worms),
        Count_norm=sum(Worms),
        Size_Mean=mean(N/Worms),
        Size_SD=sd(N/Worms),
        Mean=mean(Value),
        SD=sd(Value),
        SE=SD/sqrt(length(Value))) %>%
  melt(measure.vars=c('Mean','SD','SE'),variable.name='SumStat',value.name='SumValue') %>%
  dcast(Plate+Index+Metabolite+Metabolite_Uniq+Group+Count+Count_norm+Size_Mean+Size_SD~Type+Stat+SumStat,value.var = 'SumValue')


head(summary.sum)


# summary.sum<-ddply(summary.sel,.(Plate,Index,Metabolite,Metabolite_Uniq,Group),summarise,
#                    Count=length(Worms),
#                    Count_norm=sum(Worms),
#                    Size_Mean=mean(N/Worms),
#                    Size_SD=sd(N/Worms))




summary.sel<-subset(summary.cast,!(Replicate=='Rep1' & (Plate=='PM1'| Plate=='PM2A'))  )

lmshape.meas<-dcast(summary.sel,Replicate+UniqRep~Index,value.var = usestat)
lmshape.meas$Type<-'Measure'
dim(lmshape.meas)
colnames(lmshape.meas)


#For references:
temp<-info[!duplicated(info$Index),c('Plate','Well','Index','Group')]

refs<-subset(summary.sel,Metabolite %in% c('Negative Control') )
ref.m<-merge(temp,refs[,c('Replicate','Plate','Worm','UniqRep','Group',usestat)],by=c('Plate','Group'),all.y=TRUE) 




#Use NGM as reference
# refs<-subset(summary.sel,Metabolite %in% c('NGM') )
# ref.m<-merge(temp,refs[,c('Replicate','Worm','UniqRep',usestat)],all.x=TRUE,all.y=TRUE)


subset(ref.m,Plate=='PM1')



lmshape.ref<-dcast(ref.m,Replicate+UniqRep~Index,value.var = usestat)
lmshape.ref$Type<-'Reference'
dim(lmshape.ref)

#View(lmshape.ref)

lmshape.ref[,c('Replicate','UniqRep','Controls-A1','Controls-B1')]



#Ref average = 0?
apply(lmshape.ref[,cleanmets],2,mean,na.rm=TRUE)


#Measurement average!=0 ?
apply(lmshape.meas[,cleanmets],2,mean,na.rm=TRUE)




lmshape<-rbind(lmshape.meas,lmshape.ref)
lmshape$Type<-factor(lmshape$Type,levels=c('Reference','Measure'),labels=c('Reference','Measure'))




allresults.t<-data.frame()
prec<-0
met.len<-length(cleanmets)


#met<-'PM1_A2'

for (metid in 1:length(cleanmets)) {
  met<-as.character(cleanmets[metid])
  precn<-metid*100/length(cleanmets)
  if (precn-prec>5) {
    print(paste(round(precn,digits=0),'%',sep=''))
    prec<-precn
  }
  model<-lm(paste("`",met,"`~Type",sep=""),lmshape) #0+Type when NGM is used
  result<-summary(model)
  res<-as.data.frame(result$coefficients)
  res$Index<-met
  res$Comparison<-rownames(res)
  rownames(res)<-NULL
  allresults.t<-rbind(allresults.t,res[res$Comparison=='TypeMeasure',])
}

allresults.t$Comparison<-NULL
head(allresults.t)



logFDRbreaks<-c(-1,1.3,2,3,50)

results.in<-allresults.t %>%
  merge(info[!duplicated(info$Index),],by='Index',all.x=TRUE) %>%
  rename(logFC=Estimate,SE=`Std. Error`,p.value=`Pr(>|t|)`,t.value=`t value`) %>%
  mutate(FDR=p.adjust(p.value,method = 'fdr'),
         Bonferroni=p.adjust(p.value,method = 'bonferroni'),
         PE=logFC+SE,
         NE=logFC-SE,
         logFDR=-log10(FDR)) %>%
  merge(summary.sum[,c('Index','Count','Count_norm','Size_Mean','Size_SD','Raw_Q90_Mean','Raw_Q90_SD','Raw_Q90_SE','Norm_Q90_Mean','Norm_Q90_SD','Norm_Q90_SE')],by='Index',all.x=TRUE) %>%
  mutate(logFDR_bin=cut(logFDR, breaks=logFDRbreaks,labels=c('N.S.','p<0.05','p<0.01','p<0.001')),
         Metabolite_Uniq=factor(Metabolite_Uniq,
                                levels=Metabolite_Uniq[order(logFC)],
                                labels=Metabolite_Uniq[order(logFC)]),
         Name=NA) %>%
  dplyr::select(Plate,Well,Index,Metabolite:Metabolite_Uniq,EcoCycID:Description,logFC:p.value,FDR:logFDR_bin)




head(summary.sum)

head(results.in)

dim(results.in)



write.csv(results.in,paste0(odir,'/Celegans_results_withNGM.csv'),row.names = FALSE)



#Continue

results.in<-read.csv(paste0(odir,'/Celegans_results_withNGM.csv'))

#results.in<-read.csv(paste0(odir,'/Celegans_results_against_NGM.csv'))





head(results.exp)


psig<-subset(results.in,FDR<0.05)
psig[order(psig$logFC),c('Index','Metabolite_Uniq','logFC','SE','p.value','FDR','Bonferroni','Count_norm')]



subset(summary.cast,Replicate=='Rep2' & Plate=='PM1' & Well=='C4')




data.m.sel<-subset(data.m,!(Replicate=='Rep1' & (Plate=='PM1' | Plate=='PM2A')) )

for (plate in c('PM1','PM2A','PM3B','PM4A','Controls')) {
  platedata<-subset(data.m.sel,Plate==plate)
  platestats<-subset(summary.sel,Plate==plate )
  
  cumdistlog<-ggplot(platedata,aes(x=BrightnessLog,color=Replicate))+
    #geom_histogram(position='identity',binwidth=0.05)+
    stat_ecdf(aes(group=interaction(Replicate,Worm)))+
    scale_x_continuous(breaks=seq(-5,0,by=1))+
    scale_y_continuous(labels = percent)+
    geom_vline(data=platestats,aes(xintercept=Raw_Q90,color=Replicate,group=interaction(Replicate,Worm)),linetype=1,alpha=0.5)+
    #geom_vline(data=platestats,aes(xintercept=Raw_Median,color=Replicate),linetype=3,alpha=0.5)+
    # geom_vline(data=platestats,aes(group=interaction(File,Replicate),
    #                                xintercept=BrightnessLog_Mean),color='black')+
    ggtitle(paste('Cumulative brightness distribution in log2 scale: ',plate,sep=''),subtitle='Q90 - solid line')+
    ylab('Cumulative frequency, %')+
    xlab('log2 Brightness')+
    #geom_density() + 
    #stat_bin(aes(y=cumsum(..count..)),geom="line")+
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = -90))+
    facet_grid(Row~Col)
  
  fname<-paste(odir,"/Cumulative_distribution_log2_",plate,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=18,height=9)
  print(cumdistlog)
  dev.off()
}



#Scoring plots

scored<-subset(summary.cast,Replicate %in% c('Rep2','Rep3'))
scored$Uniq_Index<-paste(scored$Replicate,scored$Index,sep='-')
scored$Score<-as.numeric(scored$Score)
scored$ScoreFac<-as.factor(as.character(scored$Score))


selmets<-c('PM1-G2',
           'PM1-B12',
           'PM1-G12',
           'PM1-H1',
           'PM1-D8',
           'PM1-F11',
           'PM2A-A11',
           'PM2A-B1',
           'PM2A-B3',
           'PM2A-B12',
           'PM2A-H9',
           'PM2A-E5',
           'PM2A-E12',
           'PM2A-H10',
           'PM4A-A10',
           'PM3B-D7',
           'PM3B-A5',
           'PM3B-E8',
           'PM3B-C5',
           'PM3B-F11',
           'PM3B-H11',
           'PM3B-A12',
           'PM4A-G6',
           'PM4A-H6',
           'PM4A-E6',
           'PM4A-A12',
           'PM4A-H9',
           'PM4A-H12')



ggplot(scored,aes(x=Score,y=Raw_Q90))+
  geom_boxplot(aes(group=interaction(Plate,Replicate,ScoreFac )) )+
  geom_point()+
  scale_x_continuous(limits=c(-1,5),breaks=seq(-5,5,by=1))+
  #geom_text(aes(label=Well),color='gray30')+
  # geom_text(aes(label=ifelse(Index %in% selmets,paste(as.character(Well),as.character(Metabolite),sep=' ' ),'')),
  #           color='red',size=2,nudge_y=0.1)+
  geom_text_repel(aes(label=ifelse(Index %in% selmets,paste(as.character(Well),as.character(Metabolite_Uniq),sep=' ' ),'')),
                  color='red',size=2)+
  facet_grid(Plate~Replicate)

# dev.copy2pdf(device=cairo_pdf,
#              file=paste(odir,"/Scoring_vs_Raw_Q90.pdf",sep=''),
#              width=12,height=15)


ggplot(scored,aes(x=Score,y=Norm_Q90))+
  geom_point()+
  geom_text(aes(label=Well),color='gray30',size=4,nudge_y=0.5)+
  # geom_text_repel(aes(label=ifelse(Uniq_Index %in% selmets,as.character(Metabolite),'')),
  #                 color='gray30')+
  facet_grid(Plate~Replicate)










#Significance plots

ggplot(results.in,aes(x=logFC,y=logFDR,color=Count))+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=0.5)+
  geom_point()+
  labs(color='Worms, N')+
  xlab('logFC vs Negative Control')+
  geom_text_repel(aes(label=ifelse(FDR<0.05,as.character(Metabolite_Uniq),"") ) )

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_plot.pdf",sep=''),
             width=20,height=12)



ggplot(results.in,aes(x=logFC,y=logFDR,color=Count))+
  geom_errorbarh(aes(xmin=NE,xmax=PE),alpha=0.5)+
  geom_point()+
  labs(color='Worms, N')+
  xlab('logFC vs Negative Control')+
  geom_text_repel(aes(label=ifelse((FDR<0.05 & logFC<0) | Metabolite=='Negative Control',as.character(Metabolite_Uniq),"") ) )

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Volcano_plot_against_NGM.pdf",sep=''),
             width=20,height=12)



ggplot(summary.sel,aes(x=N,y=Worms,color=Replicate))+
  geom_point()+
  xlab('Selection size, pixels')
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Worm_count_estimation.pdf",sep=''),
             width=9,height=6)



ggplot(results.in,aes(x=logFC,y=Size_Mean))+
  geom_errorbar(aes(ymin=Size_Mean-Size_SD,ymax=Size_Mean+Size_SD))+
  geom_errorbarh(aes(xmin=NE,xmax=PE))+
  xlab('Worm brightness logFC vs Negative control')+
  ylab('Worm size, pixels')+
  geom_point()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Worm_size_brightness_correlation.pdf",sep=''),
             width=9,height=6)




ggplot(results.in,aes(x=Metabolite_Uniq,y=logFC,color=logFDR_bin))+
  geom_hline(yintercept = 0,color='red',alpha=0.5)+
  geom_errorbar(aes(ymin=NE,ymax=PE))+
  geom_point()+
  xlab('Metabolite')+
  ylab('Brightness logFC vs Negative control')+
  labs(color='Significance (FDR)')+
  scale_colour_manual(values = c("gray40", "red4", "red3",'red'))+
  coord_flip()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Worm_logFC_tall_new_against_NGM.pdf",sep=''),
             width=9,height=50)


ggplot(results.in,aes(x=logFC,y=Count_norm,color=Plate))+
  #geom_errorbar(aes(ymin=Size_Mean-Size_SD,ymax=Size_Mean+Size_SD))+
  geom_errorbarh(aes(xmin=NE,xmax=PE))+
  geom_text_repel(aes(label=ifelse(logFC>0.75 & Count_norm<20,as.character(Metabolite_Uniq),'') ),
                  color='black',size=3)+
  xlab('Worm brightness logFC vs Negative control')+
  ylab('Estimated number of worms, N')+
  geom_point()
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Worm_number_vs_logFC.pdf",sep=''),
             width=12,height=8)


ggplot(results.in,aes(x=logFC))+
  geom_histogram(binwidth=0.2)







#Comparisons based on brightness distribution

head(data.l)


summary<-TRUE



#Normalisation for each metabolite in each worm
#This will take a lot of time!!!
data.la.sel<-subset(data.l,!(Replicate=='Rep1' & (Plate=='PM1'| Plate=='PM2A' )) )

data.norm<-data.la.sel %>%
  ddply(.(Replicate,Plate,Well,Row,Col,Index,Worm,Group,Metabolite,Metabolite_Uniq,LogBrightness_sum),summarise,Frequency=sum(Frequency)) %>%
  ddply(.(Replicate,Plate,Well,Row,Col,Index,Worm,Group,Metabolite,Metabolite_Uniq),mutate,N=sum(Frequency)) %>%
  mutate(Norm=Frequency/N)

head(data.norm)



#Get distributions for controls in each plate
ref<-data.la.sel %>%
  filter(Metabolite=='Negative Control') %>%
  ddply(.(Replicate,Plate,Group,LogBrightness_sum),summarise,Frequency=sum(Frequency)) %>%
  ddply(.(Replicate,Plate,Group),mutate,N=sum(Frequency)) %>%
  mutate(Reference=Frequency/N)

head(ref)

rows<-c('A','B','C','D','E','F','G','H')

data.n<-data.norm %>%
  rename(Absolute=Norm) %>%
  merge(ref[,c('Replicate','Plate','Group','LogBrightness_sum','Reference')],by=c('Replicate','Plate','Group','LogBrightness_sum'),all.x=TRUE) %>%
  mutate(Relative=Absolute-Reference) %>%
  melt(measure.vars = c('Absolute','Relative'),variable.name = 'Measure',value.name = 'Value')
  

data.Wsum<-data.n %>%
  ddply(.(Plate,Well,Row,Col,Index,Group,Metabolite,Metabolite_Uniq,LogBrightness_sum,Measure),summarise,
        Mean=mean(Value),SD=sd(Value)) %>%
  mutate(LogBrightness_num=as.numeric(as.character(LogBrightness_sum)),
         Row=factor(Row,levels=rows,labels=rows),
         Col=as.factor(Col))


data.Isum<-data.n %>%
  ddply(.(Index,Group,Metabolite,Metabolite_Uniq,LogBrightness_sum,Measure),summarise,
        Mean=mean(Value),SD=sd(Value)) %>%
  mutate(LogBrightness_num=as.numeric(as.character(LogBrightness_sum)),
         Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         Metabolite=as.character(Metabolite),
         Metabolite=ifelse(Metabolite=='NGM','Negative Control',Metabolite),
         Metabolite=as.factor(Metabolite),
         Metabolite=relevel(Metabolite,ref='Negative Control'))



dim(data.Isum)


#saveRDS(data.Isum,"Celegans_fluorescence_distributions.rds")



head(data.nsum)
#

results.in[grepl('Acetoacetic',results.in$Metabolite_Uniq),]

indxsel<-c('Controls-A1','PM1-A1','PM1-A2','PM4A-E6','PM1-G7')



res.sel<-results.in %>%
  filter(Index %in% indxsel) %>%
  mutate(Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         Metabolite=as.character(Metabolite),
         Metabolite=ifelse(Index=='Controls-A1','Negative Control',Metabolite),
         Metabolite=as.factor(Metabolite),
         Metabolite=relevel(Metabolite,ref='Negative Control'))


data.Isel<-data.Isum %>% filter(Index %in% indxsel & Measure=='Absolute')



ggplot(res.sel,aes(color=Type))+
  geom_vline(xintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean'],color='grey80',linetype='longdash')+
  geom_rect(aes(xmin=Raw_Q90_Mean-Raw_Q90_SE,xmax=Raw_Q90_Mean+Raw_Q90_SE,fill=Type),ymin=-Inf,ymax=Inf,alpha=0.2,color=NA)+
  geom_vline(aes(xintercept=Raw_Q90_Mean,color=Type))+
  geom_ribbon(data=data.Isel,aes(x=LogBrightness_num,ymin=(Mean-SD)*100,ymax=(Mean+SD)*100,fill=Type),alpha=0.5,color=NA)+
  geom_line(data=data.Isel,aes(x=LogBrightness_num,y=Mean*100,color=Type))+
  #geom_point()+
  ylab('Frequency, %')+
  xlab('log2 Brightness')+
  scale_x_continuous(breaks=seq(-20,20,by=1))+
  coord_cartesian(ylim=c(0,15))+
  facet_grid(.~Metabolite)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Brightness_distributions_example_new.pdf",sep=''),
             width=9,height=3)



dwidth<-0.5

ggplot(res.sel,aes(x=Metabolite,y=Raw_Q90_Mean,color=Type))+
  geom_hline(yintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean'],color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=Raw_Q90_Mean-Raw_Q90_SE,ymax=Raw_Q90_Mean+Raw_Q90_SE),position = position_dodge(width = dwidth),width=0.5)+
  geom_point(position = position_dodge(width = dwidth))+
  theme(axis.text.x = element_text(angle = 90,hjust=1))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Brightness_distributions_comparison_new.pdf",sep=''),
             width=4,height=4)



ggplot(res.sel,aes(x=Metabolite,y=logFC,color=Type))+
  geom_hline(yintercept=0,color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,4))+
  geom_errorbar(aes(ymin=logFC-SE,ymax=logFC+SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth))+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Brightness_distributions_comparison_logFC.pdf",sep=''),
             width=6,height=2)




for (plate in c('PM1','PM2A','PM3B','PM4A','Controls')) {
  platedata<-subset(data.Wsum,Plate==plate)
  
  distlog<-ggplot(platedata,
                  aes(x=LogBrightness_num,y=Mean,color=Measure,fill=Measure))+
    geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.5)+
    geom_line()+
    #geom_point()+
    scale_x_continuous(breaks=seq(-6,0,by=1))+
    scale_y_continuous(labels = percent)+
    ggtitle(paste('Brightness distribution in log2 scale: ',plate,sep=''))+
    ylab('Frequency, %')+
    xlab('log2 Brightness')+
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = -90))+
    facet_grid(Row~Col)
  
  fname<-paste(odir,"/Distribution_log2_",plate,".pdf",sep = '')
  print(fname)
  cairo_pdf(fname,width=18,height=9)
  print(distlog)
  dev.off()
}


#For PCA and tSNE
data.sc<-dcast(subset(data.nsum,Measure=='Rel'),Plate+Well+Index+Metabolite+Metabolite_Uniq~LogBrightness_sum,value.var = 'Mean')
data.sc<-subset(data.sc,Metabolite!='Negative Control')


dim(data.sc)

head(data.sc)

brightsc<-setdiff(colnames(data.sc),c('Replicate','Plate','Well','Worm','Index','Metabolite','Metabolite_Uniq'))

#data.sc$Sum<-apply(data.sc[, brightsc],1,sum)
#data.sc[,brightsc]<-data.sc[,brightsc]/data.sc$Sum
#data.sc$Index<-paste(data.sc$Plate,data.sc$Well,sep='_')
rownames(data.sc)<-data.sc$Index


head(data.sc)

#data.sca<-merge(data.sc,info[,c('Plate','Well','Metabolite')],by=c('Plate','Well'))



pcashape<-data.sc[,brightsc]
metabolites<-data.sc$Metabolite_Uniq
replicates<-data.sc$Replicate
plates<-data.sc$plate


hca_sample<-hclust(dist(pcashape,method="euclidean"),method="complete")

plot(hca_sample, labels=metabolites,
     hang=-1, main="Cluster Dendrogram", xlab="", sub="", cex=1)
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hierarchical_Clustering_all.pdf",sep=''),
             width=60,height=10)



ir.pca <- prcomp(pcashape,
                 center = FALSE,
                 scale. = FALSE)


pcadata<-merge(data.frame(ir.pca$x),data.sc[,c('Plate','Well','Index','Metabolite','Metabolite_Uniq')],by = 0)
pcadata$Metabolite<-as.factor(pcadata$Metabolite)
pcadata$Metabolite_Uniq<-as.factor(pcadata$Metabolite_Uniq)


#Get results
pcadata.m<-merge(pcadata,results.in[,c('Plate','Well','Index','logFC','FDR','SE')],by=c('Plate','Well','Index'))

pcaresult<-summary(ir.pca)$importance
PC1prc<-round(pcaresult['Proportion of Variance',][[1]]*100,0)
PC2prc<-round(pcaresult['Proportion of Variance',][[2]]*100,0)


amp<-2
nstep<- 8
minv<- -amp
maxv<- amp
brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue4","blue4", "gray90", "red4", "red4"))(n = nstep)


ggplot(pcadata.m,aes(x=PC1,y=PC2,color=logFC))+
  xlab(paste('PC1 - ',PC1prc,'% of variance',sep=''))+
  ylab(paste('PC2 - ',PC2prc,'% of variance',sep=''))+
  geom_text(aes(label=Metabolite_Uniq),nudge_y=0.005)+
  #geom_point(size=1,stroke=1,shape=21)+
  scale_colour_gradientn(colours = bgg,
                         breaks=brks,limits=c(minv,maxv))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/PCA_by_metabolite.pdf",sep=''),
             width=15,height=15, useDingbats=FALSE)




colors = rainbow(length(metabolites))
names(colors) = unique(metabolites)


## Executing the algorithm on curated data
tsne <- Rtsne(pcashape, dims = 2, perplexity=100, verbose=TRUE, max_iter = 1000)
tsnedata<-as.data.frame(tsne$Y)


tsnedata[,c('Plate','Well','Metabolite')]<-data.sc[,c('Plate','Well','Metabolite')]
tsnedata.m<-merge(tsnedata,results.in[,c('Plate','Well','logFC','FDR','SE','Metabolite_Uniq')],by=c('Plate','Well'))


#Plotting
erralpha<-1
errcolor<-'grey80'
eqsize<-3



amp<-2
nstep<- 8
minv<- -amp
maxv<- amp
brks<-seq(minv,maxv,by=(maxv-minv)/(nstep))


ggplot(tsnedata.m,aes(x=V1,y=V2,color=logFC))+
  geom_text(aes(label=Metabolite_Uniq),size=2,segment.alpha = 0)+
  xlab('X')+
  ylab('Y')+
  scale_colour_gradientn(colours = bgg,
                         breaks=brks,limits=c(minv,maxv))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/tSNE_by_metabolite.pdf",sep=''),
             width=15,height=15, useDingbats=FALSE)






#install.packages('emdist')
#library('emdist')

# A <- matrix(1:6 / sum(1:6), 2)
# B <- matrix(c(0, 0, 0, 0, 0, 1), 2)

#emd2d(A, B, dist="manhattan")




library(dplyr)
library(reshape2)



#Make distance matrix
freq.sel<-data.f.sum[,c('Index','Metabolite','LogBrightness_sum','Mean','SD')]

dim(data.f.sum)

pairwise<-full_join(freq.sel,freq.sel,by=c('LogBrightness_sum'),all=TRUE)

dim(pairwise)


emds<-ddply(pairwise,.(LogBrightness_sum,Index.x,Index.y,Metabolite.x,Metabolite.y),summarise,EMD=emd2d(as.matrix(Mean.x),as.matrix(Mean.y),dist="manhattan"))

emds[order(emds$EMD,decreasing = TRUE),]


distmat<-dcast(emds,Index.x~Index.y,value.var = 'EMD')




emds<-ddply(data.ns,.(Plate,Well,Metabolite,Group),summarise,EMD=emd2d(as.matrix(Norm),as.matrix(Ref),dist="manhattan"))

emds[order(emds$EMD,decreasing = TRUE),]


A2<-subset(data.ns,Plate=='PM1' & Well=='A2')

emd2d(as.matrix(A2$Norm),as.matrix(A2$Ref))

#data.sum<-ddply(data.ns,.(Plate,Well,Metabolite,LogBrightness_sum),summarise,P=sum(Rel))




# 
# 
# fit<-lm(Norm_Q95~Index,summary.cast)
# mod.results<-summary(fit)
# 
# allresults.comb<-as.data.frame(mod.results$coefficients)
# allresults.comb$Index<-gsub('Index','',rownames(allresults.comb))
# allresults.ct<-subset(allresults.comb,Index!='(Intercept)')
# 
# results.comb<-resultsprep(allresults.ct,info)
# 
# 
# head(results.comb)
# write.csv(results.comb,paste(odir,'/All_results_Combined_LM_Q95.csv',sep=''),row.names = FALSE)
# 
# subset(results.comb,Metabolite=='D-Arabinose')
# subset(results.comb,FDR<0.05)
# 
# 
# 
# 
# 
# 
# fitm<-mblm(NormQ90~Metabolite,platedata,repeated=TRUE)
# results<-summary(fitm)
# results
# 
# 
# 
# 
# 
# 
# fit<-lm(Mean_log~Type,subset(filesummary,Supplement=='Glucose'))
# results<-summary(fit)
# results
# 
# 
# fit<-lm(Mean_log~Type,subset(filesummary,Supplement=='L-Proline'))
# results<-summary(fit)
# results
# 

# cumdist<-ggplot(platedata,aes(x=Brightness,color=Replicate))+
#   #geom_histogram(position='identity',binwidth=0.05)+
#   stat_ecdf(aes(group=interaction(File,Replicate)))+
#   scale_x_continuous(breaks=seq(0,1,by=0.25))+
#   geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#                                  xintercept=Raw_Brightness_Q90),color='red')+
#   geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#                                  xintercept=Raw_Brightness_Median),color='blue')+
#   # geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#   #                                xintercept=BrightnessLog_Mean),color='black')+
#   ggtitle(paste('Cumulative brightness distribution in linear scale: ',plate,sep=''),subtitle='Estimates: Q90 - red, Median - blue')+
#   ylab('Cumulative frequency, %')+
#   xlab('Brightness')+
#   #geom_density() + 
#   #stat_bin(aes(y=cumsum(..count..)),geom="line")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -90))+
#   facet_grid(Row~Col)
# 
# fname<-paste(odir,"/Cumulative_distribution_",plate,".pdf",sep = '')
# print(fname)
# cairo_pdf(fname,width=12,height=6)
# print(cumdist)
# dev.off()
# 
# dist<-ggplot(platedata,aes(x=Brightness,color=Replicate))+
#   #geom_histogram(position='identity',binwidth=0.05)+
#   geom_density(aes(group=interaction(File,Replicate)),col=2) +
#   scale_x_continuous(breaks=seq(0,1,by=0.25))+
#   geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#                                  xintercept=Raw_Brightness_Q90),color='red')+
#   geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#                                  xintercept=Raw_Brightness_Median),color='blue')+
#   ggtitle(paste('Brightness distribution in linear scale: ',plate,sep=''),subtitle='Estimates: Q90 - red, Median - blue')+
#   xlab('Brightness')+
#   ylab('Frequency, %')+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -90))+
#   facet_grid(Row~Col)
# 
# fname<-paste(odir,"/Distribution_",plate,".pdf",sep = '')
# print(fname)
# cairo_pdf(fname,width=12,height=6)
# print(dist)
# dev.off()



# distlog<-ggplot(platedata,aes(x=Brightness,color=Replicate))+
#   #geom_histogram(position='identity',binwidth=0.05)+
#   geom_density(col=2) +
#   scale_x_continuous(breaks=seq(-5,0,by=1))+
#   geom_vline(data=platestats,aes(xintercept=Raw_Q90,color=Replicate),linetype=1,alpha=0.5)+
#   geom_vline(data=platestats,aes(xintercept=Raw_Median,color=Replicate),linetype=3,alpha=0.5)+
#   ggtitle(paste('Brightness distribution in log2 scale: ',plate,sep=''),subtitle='Estimates: Q90 - solid, Median - dashed')+
#   xlab('log2 Brightness')+
#   ylab('Frequency, %')+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = -90))+
#   facet_grid(Row~Col)
# 
# fname<-paste(odir,"/Distribution_log2_",plate,".pdf",sep = '')
# print(fname)
# cairo_pdf(fname,width=12,height=6)
# print(distlog)
# dev.off()
# 
# 
# normdist<-ggplot(platestats,aes(y=Norm_Q90,x=0))+
#   #geom_histogram(position='identity',binwidth=0.05)+
#   #stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "identity") +
#   geom_boxplot(color='gray50',width=0.5)+
#   geom_hline(yintercept = 0,color='gray50')+
#   geom_point(aes(color=Replicate)) +
#   scale_y_continuous(breaks=seq(-3,3,by=1))+
#   scale_x_continuous(breaks=c())+
#   ylab('Normalised log2 brightness')+
#   xlab('')+
#   # geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#   #                                xintercept=Raw_BrightnessLog_Q90),linetype=1)+
#   # geom_vline(data=platestats,aes(group=interaction(File,Replicate),
#   #                                xintercept=Raw_BrightnessLog_Median),linetype=4)+
#   ggtitle(paste('Relative fluorescence: ',plate,sep=''),
#           subtitle='')+
#   theme(panel.grid.major = element_blank())+
#   facet_grid(Row~Col)
# 
# fname<-paste(odir,"/Normalised_distribution_log2_",plate,".pdf",sep = '')
# print(fname)
# cairo_pdf(fname,width=12,height=6)
# print(normdist)
# dev.off()
