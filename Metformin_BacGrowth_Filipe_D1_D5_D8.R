library(tidyverse)
library(readxl)

library(RColorBrewer)
library(ggrepel)


devtools::install_github("PNorvaisas/PFun")
library(PFun)

remove(adjustments)

theme_set(theme_light())



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



data<-read_xlsx('Filipe_bac_growth/Bacterial_growth_curve_D1_D5_D8.xlsx',sheet= 'Clean') %>%
  mutate_at(c('Day','Strain','Replicate'),as.factor) %>%
  gather(Time_min,OD_raw,matches("[[:digit:]]")) %>%
  mutate_at(c('Time_min','OD_raw'),as.numeric) %>%
  group_by(Day, Strain, Replicate) %>%
  mutate(Time_h=Time_min/60,
         Base=mean(OD_raw[Time_h<1]),
         OD=OD_raw-Base) %>%
  ungroup





ggplot(data,aes(x=Time_h,y=OD))+
  geom_line(aes(group=interaction(Replicate)))+
  xlab('Time, h')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_grid(Strain~Day,labeller = labeller(.rows = label_both, .cols = label_both))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=9,height=6)



data.sum<-data %>%
  group_by(Day,Strain,Time_min,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))

ggplot(data.sum,aes(x=Time_h,y=OD_Mean,color=Strain,fill=Strain))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_grid(Day~.,labeller=labeller(.rows = label_both))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_by_day.pdf",sep=''),
             width=6,height=6)




data.int<-data %>%
  group_by(Day,Strain,Replicate) %>%
  summarise(OD_Int=sum(OD)*5/60) %>%
  mutate(logOD_Int=log2(OD_Int),
         Sample=paste(gsub('-','',Strain),Day,sep='_')) %>%
  gather(Measure,Value,OD_Int,logOD_Int)
  
  

sel.samples<-unique(data.int$Sample)




contrasts<-read.contrasts2('!Contrasts_Filipe_Bac.xlsx')

contrasts$Contrasts.table

contrasts.desc<-contrasts$Contrasts.table %>%
  select(Description:Strain)

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

contrasts.desc


lmdata<-data.int%>%
  group_by(Measure) %>%
  do(hypothesise2(.,formula='Value~0+Sample',contr.matrix))
  

results<-contrasts.desc %>%
  left_join(lmdata) %>%
  adjustments %>%
  mutate(Contrast=factor(Contrast,levels=contrasts.desc$Contrast)) %>%
  select(Measure,Contrast,Description:Strain,everything())



results.m<-results %>%
  gather(Stat,Value,logFC:logFDR)

results.castfull<-results.m %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Measure,CS,Value) %>%
  spread(CS,Value)

results.cast<-results.m %>%
  filter(Stat %in% c('logFC','FDR')) %>%
  arrange(Contrast,desc(Stat)) %>%
  unite(CS,Contrast,Stat) %>%
  select(Measure,CS,Value) %>%
  spread(CS,Value)

head(results.castfull)
head(results.cast)


write.csv(results,paste(odir,'/Growth_results.csv',sep=''),row.names = FALSE)



ODstars_day<-results %>%
  filter(Measure=='OD_Int' & Contrast_type=='Strain') %>%
  select(Day,pStars)

ODstars_day$Strain<-'N2'

head(data.int)


dwidth<-0.5

data.int %>%
  filter(Measure=='OD_Int') %>%
  ggplot(aes(x=Day,y=Value,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('OD integral, OD*h')+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=ODstars_day,aes(x=Day,y=4.5,label=as.character(pStars)))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_OD_int_by_day.pdf",sep=''),
             width=6,height=4)



logODstars_day<-results %>%
  filter(Measure=='logOD_Int' & Contrast_type=='Strain') %>%
  select(Day,pStars)

logODstars_day$Strain<-'N2'


data.int %>%
  filter(Measure=='logOD_Int') %>%
  ggplot(aes(x=Day,y=Value,fill=Strain))+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('log2 OD integral, log2(OD*h)')+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2)+
  geom_text(data=logODstars_day,aes(x=Day,y=2.25,label=as.character(pStars)))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_logOD_int_by_day.pdf",sep=''),
             width=6,height=4)



#By strain

ODstars_strain<-results %>%
  filter(Measure=='OD_Int' & Contrast_type=='Time' & !Contrast %in% c('N2_8-5','Phm2_8-5') )


data.int %>%
  filter(Measure=='OD_Int') %>%
  ggplot(aes(x=Day,y=Value,fill=Strain))+
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





logODstars_strain<-results %>%
  filter(Measure=='logOD_Int' & Contrast_type=='Time' & !Contrast %in% c('N2_8-5','Phm2_8-5') )


data.int %>%
  filter(Measure=='logOD_Int') %>%
  ggplot(aes(x=Day,y=Value,fill=Strain))+
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



