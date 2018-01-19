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


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir<-'Summary_mutants'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


data<-read.table('Clean/2017-Rosie_Mutants/Summary.csv',sep=',',quote = '"',header = TRUE)


strainlist<-c('OP50Sens','OP50MR','crp','crr','cyaA','mlc','cpdA','ptsG')


data$Strains<-factor(data$Strains, levels=strainlist,labels=strainlist)
data$Metformin<-as.factor(data$Metformin_conc_mM)

contr<-subset(data,Strains=='OP50Sens')

contr<-rename(contr,c('Int_600nm_f'='Int_Sens','Int_600nm_log'='Int_Sens_log','a_log'='a_Sens'))
contr0<-subset(contr,Metformin==0)
contr0<-rename(contr0,c('Int_Sens'='Int_Sens0','Int_Sens_log'='Int_Sens0_log','a_Sens'='a_Sens0'))

datac<-merge(data,contr[,c('Plate','Replicate','Metformin','Int_Sens','Int_Sens_log','a_Sens')],by=c('Plate','Replicate','Metformin'),all.y=TRUE)

datac<-merge(datac,contr0[,c('Plate','Replicate','Int_Sens0','Int_Sens0_log','a_Sens0')],by=c('Plate','Replicate'),all.y=TRUE)



datac$Int_600nm_rel<-datac$Int_600nm_f/datac$Int_Sens
datac$Int_600nm_rel_log<-datac$Int_600nm_log-datac$Int_Sens_log
datac$Int_600nm_rel0_log<-datac$Int_600nm_log-datac$Int_Sens0_log


datac$a_rel<-datac$a_log-datac$a_Sens
datac$a_rel0<-datac$a_log-datac$a_Sens0




datacm<-melt(datac,id.vars = c('Strains','Metformin','Metformin_conc_mM'),measure.vars = c('Int_600nm_f','Int_600nm_log',
                                                                                           'Int_600nm_rel','Int_600nm_rel_log','Int_600nm_rel0_log',
                                                                                           'a_log','a_Sens','a_Sens0'),
             variable.name = 'Estimate',value.name = 'Value')


summaryt<-ddply(datacm,.(Strains,Metformin,Metformin_conc_mM,Estimate),
             summarise,Mean=mean(Value,na.rm=TRUE),SD=sd(Value,na.rm=TRUE))


summarym<-melt(summaryt,id.vars = c('Strains','Metformin','Metformin_conc_mM','Estimate'),measure.vars = c('Mean','SD'),
             variable.name = 'Stat',value.name = 'Value')

summary<-dcast(subset(summarym,!Stat %in% c('Int_600nm_rel0_log')),
               Strains+Metformin+Metformin_conc_mM~Estimate+Stat,mean,value.var = c('Value'),
               fill = as.numeric(NA),drop=TRUE)


write.csv(summary,
          paste(odir,"/Summary.csv",sep=''))


# fit2<-lm(Int_600nm_rel_log~0+Strains*Metformin,datac)
# summary(fit2)


lmshape<-dcast(datac,Strains+Replicate+Plate~Metformin,value.var = 'Int_600nm_rel_log',drop = TRUE)


#Comparison done within each metformin range, thus restricting variation
#Making general linear model leads to skewed SE estimation

metfc<-c('0','25','50','75','100','150')

relresults<-data.frame()
for (mc in metfc) {
  formula<-as.formula(paste("`",mc,"`","~0+Strains",sep=''))
  fit<-lm(formula,lmshape)
  results<-summary(fit)
  res<-getinfo(results$coefficients)
  res$Metformin<-mc
  relresults<-mymerge(relresults,res)
}



#fit<-lm(Int_600nm_rel_log~0+Strains:Metformin,datac)
#results<-summary(fit)
#results


#relresults<-getinfo(results$coefficients)
relresults$Comparisons<-gsub('Strains','',relresults$Comparisons)
relresults$Strains<-relresults$Comparisons
#relresults$Comparisons<-gsub('Metformin','',relresults$Comparisons)


#relresults[,c('Strains','Metformin')]<-mycolsplit(relresults$Comparisons,2,':')
rownames(relresults)<-NULL
relresults$Comparisons<-NULL
relresults$Strains<-factor(relresults$Strains, levels=strainlist,labels=strainlist)
relresults$Metformin<-factor(relresults$Metformin,levels=c(0,25,50,75,100,150),labels=c(0,25,50,75,100,150))
relresults$Metformin_conc_mM<-as.numeric(as.character(relresults$Metformin))
relresults<-rename(relresults,c('Estimate'='logFC','Std..Error'='SE','Pr...t..'='p.value'))
relresults$FDR<-p.adjust(relresults$p.value,method = 'fdr')

relresults$Prc<-(2^relresults$logFC)*100-100
relresults$Prc_min<-(2^(relresults$logFC-relresults$SE))*100-100
relresults$Prc_max<-(2^(relresults$logFC+relresults$SE))*100-100
relresults$Stars<-as.character(stars.pval(relresults$FDR))
relresults$Stars<-ifelse(relresults$Stars==".","",relresults$Stars)
head(relresults)


write.csv(relresults[,c('Strains','Metformin_conc_mM','logFC','SE','Prc','Prc_min','Prc_max','t.value','p.value','FDR')],
          paste(odir,"/Stats.csv",sep=''))




ggplot(datac,aes(x=Metformin,y=Int_600nm_log,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Absolute_growth_log.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)

ggplot(datac,aes(x=Metformin,y=Int_600nm_f,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Absolute_growth.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)


ggplot(datac,aes(x=Metformin,y=a_log,color=Strains))+
  #geom_boxplot()+
  ylab('Growth rate, Doubling/h')+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Absolute_growth_rate.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



ggplot(datac,aes(x=Metformin,y=(2^Int_600nm_rel0_log)*100-100,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Relative growth to OP50Sens, %')+
  ylim(-100,50)+
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)



ggplot(datac,aes(x=Metformin,y=a_rel0,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Growth rate (relative to OP50Sens), doubling/h')+
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_rate.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)





ggplot(datac,aes(x=Metformin,y=(2^Int_600nm_rel0_log)*100-100,color=Strains))+
  geom_hline(yintercept=0,colour='red',alpha=0.5)+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Relative growth change in comparison to OP50Sens, %')+
  ylim(-100,50)+
  xlab('Metformin, mM')+
  facet_grid(.~Strains)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_spread.pdf",sep=''),
             width=12,height=5, useDingbats=FALSE)



ggplot(datac,aes(x=Metformin,y=a_rel0,color=Strains))+
  geom_hline(yintercept=0,colour='red',alpha=0.5)+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Growth rate (relative to OP50Sens), doubling/h')+
  #ylim(-100,50)+
  xlab('Metformin, mM')+
  facet_grid(.~Strains)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_rate_spread.pdf",sep=''),
             width=12,height=5, useDingbats=FALSE)







ggplot(datac,aes(x=Metformin,y=(2^Int_600nm_rel_log)*100-100,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Relative growth change in comparison to OP50Sens (by metformin), %')+
  ylim(-100,200)+
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_by_metformin.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)




ggplot(datac,aes(x=Metformin,y=a_rel,color=Strains))+
  #geom_boxplot()+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  ylab('Growth rate (relative to OP50Sens), doubling/h')+
  #ylim(-100,200)+
  xlab('Metformin, mM')

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_rate_by_metformin.pdf",sep=''),
             width=9,height=6, useDingbats=FALSE)








ggplot(datac,aes(x=Metformin,y=(2^Int_600nm_rel_log)*100-100,color=Strains))+
  geom_hline(yintercept=0,colour='red',alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  #geom_point(data=relresults,aes(x=Metformin,y=(2^logFC)*100-100))+
  #geom_text(data=relresults,aes(label=stars.pval(p.value),x=Metformin,y=(2^logFC)*100-100 ),position='dodge',size=5,color='black')+
  #geom_boxplot()+
  ylab('Relative growth change in comparison to OP50Sens (by metformin), %')+
  ylim(-150,300)+
  xlab('Metformin, mM')+
  facet_grid(.~Strains)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_by_metformin_spread.pdf",sep=''),
             width=12,height=5, useDingbats=FALSE)





ggplot(datac,aes(x=Metformin,y=a_rel,color=Strains))+
  geom_hline(yintercept=0,colour='red',alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMMax, geom="boxplot",position = "dodge") +
  #geom_point(data=relresults,aes(x=Metformin,y=(2^logFC)*100-100))+
  #geom_text(data=relresults,aes(label=stars.pval(p.value),x=Metformin,y=(2^logFC)*100-100 ),position='dodge',size=5,color='black')+
  #geom_boxplot()+
  ylab('Growth rate (relative to OP50Sens), doubling/h')+
  #ylim(-150,300)+
  xlab('Metformin, mM')+
  facet_grid(.~Strains)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_rate_by_metformin_spread.pdf",sep=''),
             width=12,height=5, useDingbats=FALSE)



shift<-25
ggplot(relresults,aes(x=Metformin,y=Prc,fill=Strains))+
  geom_errorbar(aes(ymin=Prc_min,ymax=Prc_max),position="dodge", width=.4)+
  geom_bar(stat="identity", position="dodge", width=0.5)+
  geom_hline(yintercept=0,colour='red',alpha=0.5)+
  geom_text(aes(label=as.character(Stars),y=ifelse(Prc>=0,Prc+shift,Prc-shift)),
            size=5,color='black',angle=90,nudge_x = 0.1,vjust=0.5,hjust=0.5)+
  ylab('Relative growth change in comparison to OP50Sens (by metformin), %')+
  scale_y_continuous(breaks = seq(-1000,1000,by=50))+
  coord_cartesian(ylim=c(-100,300))+
  xlab('Metformin, mM')+
  facet_grid(.~Strains)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Relative_growth_by_metformin_spread_statsSE.pdf",sep=''),
             width=14,height=5, useDingbats=FALSE)






datats<-read.table('Clean/2017-Rosie_Mutants/Data.csv',sep=',',quote = '"',header = TRUE)
#datats<-subset(datats,Data=='600nm_b')

datatm<-melt(datats,id.vars = c('Plate','Well','Data','Media','Media_conc','Metformin_conc_mM','Replicate','Strains'),value.name = 'OD',variable.name = 'Time_s')
datatm$Time_s<-as.numeric(gsub('X','',datatm$Time_s))
datatm$Time_h<-datatm$Time_s/3600
datatm$Strains<-factor(datatm$Strains, levels=strainlist,labels=strainlist)
datatm$Metformin<-as.factor(datatm$Metformin_conc_mM)

tssum<-ddply(datatm,.(Strains,Metformin,Data,Time_h),summarise,Mean=mean(OD),SD=sd(OD))


ggplot(subset(tssum,Data=="600nm_f"),aes(x=Time_h,fill=Metformin))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD))+
  geom_line(aes(y=Mean),colour="black")+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  ylab("OD@600nm")+
  labs(fill="Metformin, mM")+
  xlab("Time, h")+
  facet_grid(.~Strains)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_curves_Ribbons.pdf",sep=''),
             width=14,height=5, useDingbats=FALSE)


ggplot(subset(tssum,Data=="600nm_log"),aes(x=Time_h,fill=Metformin))+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD))+
  geom_line(aes(y=Mean),colour="black")+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  labs(fill="Metformin, mM")+
  ylab("log_OD@600nm")+
  xlab("Time, h")+
  facet_grid(.~Strains)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_log_curves_Ribbons.pdf",sep=''),
             width=14,height=5, useDingbats=FALSE)




