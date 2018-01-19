library(tidyverse)

library(gplots)
library(gtools)
library(ggrepel)

library(reshape2)
#library(tidyr)


library(heatmap3)




#This should go with my package
#library(multcomp)
#library(contrast)




library(RColorBrewer)


#library(grid)
#library(gridExtra)


library(xlsx)

library(ggthemes)

library(limma)

devtools::install_github("PNorvaisas/PFun")
library(PFun)

#remove(hypothesise)


MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


getinfo<-function(cof) {
  df<-data.frame(cof)
  df$Comparisons<-rownames(df)
  return(df)
}

pole<-function(x,y) {
  r<-sqrt(x^2+y^2)
  o<- -2*atan(y/(x))/(pi/2)
  return(o) 
}


lm_eqn = function(m) {
  fres<-summary(m)
  l <- list(a = format(abs(coef(m)[2]), digits = 2),
            b = format(coef(m)[1], digits = 2),
            r2 = format(summary(m)$r.squared, digits = 2),
            p2 = format(pf(fres$fstatistic[1], fres$fstatistic[2], fres$fstatistic[3],lower.tail = FALSE)[[1]], digits = 2,scientific=TRUE));
  
  if (coef(m)[2] >= 0)  {
    cof <- substitute(italic(y) == b + a %.% italic(x),l)
    full <- substitute(italic(y) == b + a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  } else {
    cof <- substitute(italic(y) == b - a %.% italic(x),l) 
    full <- substitute(italic(y) == b - a %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  }
  
  stat<-substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~p2,l)
  return(list('Coef'=as.character(as.expression(cof)),
              'Stat'=as.character(as.expression(stat)),
              'Full'=as.character(as.expression(full)),
              'Atop'=as.character(as.expression(paste(cof,'\n',stat,sep='') ))))                 
}


getellipse<-function(x,y,sc=1) {
  as.data.frame(ellipse::ellipse( cor(x, y),
                                  scale=c(sd(x)*sc,sd(y)*sc),
                                  centre=c( mean(x),mean(y)) ))
}



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


pole<-function(x,y) {
  r<-sqrt(x^2+y^2)
  o<- -2*atan(y/(x))/(pi/2)
  return(o) 
}


pole2<-function(x,y) {
  o<- (atan(y/(x)))/(2*pi)
  return(o) 
}




pearson<-function(x) as.dist((1-cor(t(x)))/2)

theme_set(theme_Publication())
#theme_set(theme_light())



cwd<-"~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/"
setwd(cwd)



odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")




info<-read.table('../Biolog/Biolog_metabolites_EcoCyc.csv',sep=',',quote = '"',header = TRUE) %>%
  filter(Plate %in% c('PM1','PM2A','PM3B','PM4A'))


dupmet<-data.frame(table(trimws(as.character(info$Metabolite))))
dupls<-as.character(dupmet[dupmet$Freq>1,'Var1'])
dupls

info<-info %>%
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
         Metabolite_Uniq=ifelse(Index=='PM4A-F1','NGM_S',Metabolite_Uniq)) %>%
  rename(Metabolite_class=Description)


head(info)

PM1PM2<-subset(info,Plate %in% c('PM1','PM2A'))


results.in<-read.csv('Celegans/Summary/Celegans_results_withNGM.csv')

results.castcomb<-read.csv('Summary/Combined_results_sidebyside_full.csv')


data.Isum<-readRDS('Celegans/Celegans_fluorescence_distributions.rds')





#Time series
ectimem<-read.table('Data/Data.csv',sep=',',quote = '"',header = TRUE,check.names = FALSE) %>%
  filter(Data=='750nm_f') %>%
  gather(Time_sec,OD,ends_with(".0")) %>%
  mutate(Time_sec=as.numeric(as.character(Time_sec)),
         Time_min=Time_sec/60,
         Time_h=Time_min/60)

tsum<-ectimem %>%
  group_by(Index,Type,Time_sec,Time_h) %>%
  summarise(Mean=mean(OD),SD=sd(OD),SE=SD/sqrt(length(OD))) %>%
  left_join(info)

  





#Pick data to show

#Carboxylic acids to remove
coxy<-c('Itaconic Acid','Caproic Acid','Capric Acid','4-Hydroxy Benzoic Acid','2-Hydroxy Benzoic Acid')


selectcast<-subset(results.castcomb,!(Plate %in% c('PM3B','PM4A') & Metabolite %in% PM1PM2$Metabolite ) & !Metabolite %in% c('Negative Control',coxy) )



ecolisum<-results.castcomb %>%
  select(Index,Metabolite,Metabolite_Uniq,contains('T_G'),contains('C_G')) %>%
  gather(Type,Value,contains('T_G'),contains('C_G')) %>%
  separate(Type,c('Type','Measure','Stat'),sep='_',remove=TRUE) %>%
  spread(Stat,Value)




head(selectcast)
#Scatter plots



maincomp<-'C. elegans\nphenotype rescue\n(acs-2 GFP)'


xmin<- -3
xmax<- 4
ymin<- -3
ymax<- 4


lblsize<-4

erralpha<-1
errcolor<-'grey80'
segalpha=0.5

amp<-2
cbrks<-seq(-amp,amp,by=1)
gradcols<-c('blue4','blue','gray80','red','red4')



ncontrol<-subset(results.castcomb,Metabolite=='Negative Control' )

subset(selectcast,C_logFC<0 & C_FDR<0.05)


grp<-'Complete'

fit<-lm(T_logFC~C_logFC,selectcast)
lmeq<-lm_eqn(fit)
a<-fit$coefficients[[2]]
b<-fit$coefficients[[1]]


xmin<- -3
xmax<- 4
ymin<- -3
ymax<- 4

thres<-0.05

gradcols<-c('blue4','blue','gray80','red','red4')
#gradcols<-c('green4','green','gray80','red','red4')
#subset( ,Cel_logFC<0)
#,C_logFC<0 & C_FDR<0.05 & T_logFC<0 & T_FDR<0.05)
ggplot(selectcast,aes(y=T_logFC,x=C_logFC,color=-Cel_logFC))+
  geom_vline(xintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 0,color='gray70',alpha=0.5,linetype='longdash')+
  geom_abline(aes(intercept=1,slope=1),color='blue',alpha=0.5,linetype='longdash')+
  geom_abline(aes(slope=1,intercept=0),color='grey',linetype='longdash',size=0.5)+
  geom_errorbar(aes(ymin=T_NE,ymax=T_PE),alpha=erralpha,color=errcolor,width=0)+
  geom_errorbarh(aes(xmin=C_NE,xmax=C_PE),alpha=erralpha,color=errcolor,height=0)+
  #ggtitle(paste('Scatterplot of metformin and metabolite supplementation effects ',grp,sep='') )+
  ggtitle('Scatterplot of metformin and metabolite supplementation effects')+
  #,subtitle = paste('Metabolites with FDR<',thres,' are marked',sep='')
  geom_point(size=3)+
  # geom_point(aes(size=abs(`Cel_logFC`) ))+
  coord_cartesian(xlim=c(xmin,xmax),ylim = c(ymin,ymax))+
  # scale_x_continuous(breaks=seq(-10,10,by=1))+
  # scale_y_continuous(breaks=seq(-10,10,by=1))+
  xlab('Growth logFC vs NGM - Control')+
  ylab('Growth logFC vs NGM - +50mM Metformin')+
  #eval(parse(text = intfdr)) < thres & abs( eval(parse(text = intvar)) ) > 0.75 
  geom_text(aes(label=ifelse(Cel_logFC>0 & Cel_FDR<0.05, as.character(Metabolite_Uniq),"" ))  )+
  # geom_text_repel(aes(label=ifelse(`G_T-C_FDR`<thres & abs(`G_T-C_logFC`)>0.5 & `G_T_logFC`>1, as.character(Metabolite_Uniq),'')),
  #                 size=lblsize,nudge_y = 0.3,
  #                 force=1,
  #                 segment.colour=errcolor,
  #                 segment.alpha =segalpha)+
  scale_colour_gradientn(colours = gradcols,
                         breaks=cbrks,limits=c(-amp,amp),name=maincomp)+
  labs(color='C. elegans\nphenotype rescue\n(acs-2 GFP)')+
  # scale_size(range = c(0.25, 7),name=maincomp)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Scatter_Control-Treatment_Complete_Growth_Celegans_clean_refline_Cel<0.pdf",sep=''),
             width=10,height=6, useDingbats=FALSE)








indxsel<-c('PM1-A1','PM1-A2','PM4A-E6','PM1-G7') #'Controls-A1',










indxsel<-c('Controls-A1','PM1-A1','PM1-A2','PM1-G7','PM4A-E6')

mlevels<-c('Negative Control','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')
mlabels<-c('NGM','L-Arabinose','Acetoacetic Acid','Phosphono Acetic Acid')






res.sel<-results.in %>%
  filter(Index %in% indxsel) %>%
  mutate(Type=as.factor(ifelse(Index=='Controls-A1','Control','Treatment')),
         Metabolite=as.character(Metabolite),
         Metabolite=ifelse(Index=='Controls-A1','Negative Control',Metabolite),
         Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))





data.Isel<-data.Isum %>%
  filter(Index %in% indxsel & Measure=='Absolute') %>%
  mutate(Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))



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
             file=paste(odir,"/Summary_Celegans_brightness_distribution.pdf",sep=''),
             width=9,height=3)


tsel<-tsum %>%
  filter(Index %in% indxsel) %>%
  mutate( Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))


ggplot(tsel,aes(x=Time_h,y=Mean,fill=Type,color=Type))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),color=NA,alpha=0.2)+
  geom_line()+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  ylab("OD")+
  xlab("Time, h")+
  labs(fill="Type")+
  facet_grid(~Metabolite)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Summary_Ecoli_growth_curves.pdf",sep=''),
             width=9,height=3)









dwidth<-0.5

ggplot(res.sel,aes(x=Metabolite,y=Raw_Q90_Mean,color=Type))+
  geom_hline(yintercept=res.sel[res.sel$Index=='PM1-A1','Raw_Q90_Mean'],color='grey80',linetype='longdash')+
  scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=Raw_Q90_Mean-Raw_Q90_SE,ymax=Raw_Q90_Mean+Raw_Q90_SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=3)+
  ylab('log2 Fluoresecence brightness')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Celegans_brightness_comparison.pdf",sep=''),
             width=5,height=5)






data.sel<-ecolisum %>%
  filter(Index %in% indxsel) %>%
  mutate(Type=ifelse(Type=='C','Control','Treatment'),
         Metabolite=factor(Metabolite,levels=mlevels,labels=mlabels))
  
  
  
ggplot(data.sel,aes(x=Metabolite,y=Mean,color=Type))+
  geom_hline(yintercept=data.sel[data.sel$Index=='PM1-A1' & data.sel$Type=='Control','Mean'],color='grey80',linetype='longdash')+
  #scale_y_continuous(breaks=seq(-20,20,by=1),limits=c(-4,0))+
  geom_errorbar(aes(ymin=Mean-SE,ymax=Mean+SE),position = position_dodge(width = dwidth),width=0.25)+
  geom_point(position = position_dodge(width = dwidth),size=3)+
  ylab('log2 Growth AUC')+
  theme(axis.text.x = element_text(angle = 90,hjust=1))


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Summary_Ecoli_growth_comparison.pdf",sep=''),
             width=5,height=5)




