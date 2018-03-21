library(tidyverse)
library(readxl)

library(RColorBrewer)
library(ggrepel)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



theme_set(theme_Publication())


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


#save.image('BGA_Gut.RData')
#load('BGA_Gut.RData')



odir<-'Summary_Gut'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")





data1<-read_xlsx('Filipe_bac_growth/Bacterial growth curve D1,D5,D8.xlsx',sheet='Export') %>%
  filter(Worm=='N2') %>%
  mutate(Worm=NULL,
         File='1')

data2<-read_xlsx('Filipe_bac_growth/Bacterial growth curve D2,D6,D9 gut function.xlsx',sheet='Export') %>%
  mutate(File='2')


colnames(data1)
colnames(data2)

data<-rbind(data1,data2) %>%
  group_by(Day, Strain, Metformin_mM) %>%
  mutate(Replicate=row_number()) %>%
  ungroup %>%
  gather(Time_min,OD_raw,matches("[[:digit:]]")) %>%
  mutate_at(c('Time_min','OD_raw'),as.numeric) %>%
  group_by(File,Strain, Metformin_mM,Replicate) %>%
  mutate(Time_h=Time_min/60,
         Base=mean(OD_raw[Time_h<1]),
         OD=OD_raw-Base) %>%
  ungroup %>%
  mutate_at(c('Day','Strain','Metformin_mM','Replicate'),as.factor) %>%
  filter(File=='1') %>%
  mutate(Strain=recode(Strain,'OP50'='OP50-C'), #,'OP50-MR'='MR'
    Sample=paste(gsub('-','',recode(Strain,'OP50'='C','OP50-MR'='MR')),Day,sep='_'),
    Group=paste(recode(Strain,'OP50-C'='C','OP50-MR'='MR'),recode(Metformin_mM,'0'='C','50'='T'),sep='_'),
    DDay=paste0('D',Day))


unique(data$Replicate)
unique(data$Strain)
unique(data$Group)
#View(data)


data %>%
  group_by(File,Day, Strain, Metformin_mM,Replicate) %>%
  summarise(Count=n())




unique(data$Day)

unique(data$Strain)

unique(data$Metformin_mM)

unique(data$Replicate)




Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin, mM'


ggplot(data,aes(x=Time_h,y=OD,color=Metformin_mM,linetype=File))+
  geom_line(aes(group=interaction(Replicate,Metformin_mM)))+
  xlab('Time, h')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  scale_colour_manual(name = Metlab,values =Metcols)+
  facet_grid(Strain~Day,labeller = labeller(.rows = label_both, .cols = label_both))


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=9,height=6)



data.sum<-data %>%
  group_by(Day,Strain,Metformin_mM,Time_min,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))





ggplot(data.sum,aes(x=Time_h,y=OD_Mean,color=Metformin_mM,fill=Metformin_mM))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_fill_manual(name = Metlab,values =Metcols)+
  xlab('Time, h')+
  ylab('OD')+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_grid(Strain~Day,labeller=label_both)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary.pdf",sep=''),
             width=9,height=6)




data.int<-data %>%
  group_by(Day,Strain,Metformin_mM,Replicate,Sample,Group,DDay) %>%
  summarise(OD_Int=sum(OD)*5/60) %>%
  mutate(logOD_Int=log2(OD_Int)) %>%
  gather(Measure,Value,OD_Int,logOD_Int)



data.int.sum<-data.int %>%
  group_by(Measure,Day,Strain,Metformin_mM) %>%
  summarise(Mean=mean(Value),SD=sd(Value)) %>%
  group_by(Measure) %>%
  mutate(C_norm=Mean-Mean[Day=='1' & Strain=='OP50-C' & Metformin_mM=='0'],
         Prc=2^C_norm*100,
         PrcSD=(2^SD-1)*100)







contrasts<-read.contrasts('!Contrasts_BacGrowth_Gut.xlsx')

contrasts$Contrasts.table

contrasts.desc<-contrasts$Contrasts.table %>%
  select(Description:Strain)

contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

contrasts.desc


lmdata_S<-data.int%>%
  group_by(Measure,Day) %>%
  do(hypothesise(.,formula='Value~0+Group',contr.matrix)) %>%
  getresults(contrasts.desc %>% select(-Day)) %>%
  pluck('results') %>%
  mutate(ID=paste(Day,Contrast,sep='_'),
         Group='Any')

  
lmdata_D<-data.int%>%
  group_by(Measure,Group,Strain,Metformin_mM) %>%
  do(hypothesise(.,formula='Value~0+DDay',contr.matrix)) %>%
  getresults(contrasts.desc %>% select(-Strain,-Metformin_mM)) %>%
  pluck('results') %>%
  mutate(ID=paste(Group,Contrast,sep='_'))
  
  

results<-rbind(lmdata_S,lmdata_D) %>%
  mutate(Contrast=factor(Contrast,levels=contrasts.desc$Contrast)) %>%
  select(Measure,ID,Day,Group,Strain,Metformin_mM,Contrast,Contrast_type,Description,everything())



View(results)

# 
# results.m<-results %>%
#   gather(Stat,Value,logFC:logFDR)
# 
# results.castfull<-results.m %>%
#   arrange(ID,desc(Stat)) %>%
#   unite(CS,ID,Stat) %>%
#   select(Measure,CS,Value) %>%
#   spread(CS,Value)
# 
# results.cast<-results.m %>%
#   filter(Stat %in% c('logFC','FDR')) %>%
#   arrange(ID,desc(Stat)) %>%
#   unite(CS,ID,Stat) %>%
#   select(Measure,CS,Value) %>%
#   spread(CS,Value)

# head(results.castfull)
# head(results.cast)


write.csv(results,paste(odir,'/Growth_results.csv',sep=''),row.names = FALSE)



ODstars_treatment<-results %>%
  filter(Measure=='OD_Int' & Contrast_type=='Treatment') %>%
  select(Strain,Metformin_mM,Day,pStars)%>%
  mutate_all(as.factor)
ODstars_treatment


ODstars_day<-results %>%
  filter(Measure=='OD_Int' & Contrast_type=='Time' & Contrast!='Day_8') %>%
  select(Strain,Metformin_mM,Day,pStars)%>%
  mutate_all(as.factor)


ODstars_day

head(data.int)


dwidth<-0.5

data.int %>%
  filter(Measure=='OD_Int') %>%
  ggplot(aes(x=Day,y=Value,fill=Metformin_mM,color=Metformin_mM))+
  geom_errorbar(stat = 'summary', position = position_dodge(width=dwidth), width = 0.2,color='black')+
  geom_bar(stat = 'summary', fun.y = 'mean',position = position_dodge(width=dwidth),width=0.5) +
  ylab('AUC, OD*h')+
  xlab('Day')+
  labs(fill='Metformin, mM',
       color='Metformin, mM')+
  geom_text(data=ODstars_treatment,aes(x=Day,y=4,label=as.character(pStars)),color='black')+
  geom_text(data=ODstars_day,aes(x=Day,y=ifelse(Metformin_mM=='0',5,4.5),label=as.character(pStars)))+
  facet_grid(~Strain)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary_OD_int_by_day_bar.pdf",sep=''),
             width=8,height=4)





data.int.sum %>%
  filter(Measure=='OD_Int') %>%
  ggplot(aes(x=Day,y=Mean,color=Metformin_mM))+
  stat_summary(aes(group=interaction(Metformin_mM)),fun.y=sum, geom="line")+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width = 0.2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  geom_point(size=2)+
  ylab('AUC, OD*h')+
  xlab('Day')+
  scale_y_continuous(breaks=0:10,limits=c(2,4.2))+
  labs(color='Metformin, mM')+
  geom_text(data=ODstars_day,aes(x=Day,y=ifelse(Metformin_mM=='0',4.2,4.1),label=as.character(pStars)),nudge_x = -0.5)+
  geom_text(data=ODstars_treatment,aes(x=Day,label=as.character(pStars)),y=4,color='black')+
  facet_grid(~Strain)



dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Growth_Summary_OD_int_by_day_line.pdf",sep=''),
             width=5,height=3)


data.int.sum %>%
  filter(Measure=='OD_Int') %>%
  ggplot(aes(x=Day,y=Prc,color=Metformin_mM))+
  stat_summary(aes(group=interaction(Metformin_mM)),fun.y=sum, geom="line")+
  geom_errorbar(aes(ymin=Prc-PrcSD,ymax=Prc+PrcSD),width = 0.2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  geom_point(size=2)+
  ylab('Growth AUC vs OP50-C Day 1 Control, %')+
  xlab('Day')+
  scale_y_continuous(breaks=seq(0,300,by=50),limits=c(50,250))+
  labs(color='Metformin, mM')+
  geom_text(data=ODstars_day,aes(x=Day,y=ifelse(Metformin_mM=='0',240,250),label=as.character(pStars)),nudge_x = -0.5)+
  geom_text(data=ODstars_treatment,aes(x=Day,label=as.character(pStars)),y=150,color='black')+
  facet_grid(~Strain)


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Growth_Summary_Prc_by_day_line.pdf",sep=''),
             width=5,height=3)

