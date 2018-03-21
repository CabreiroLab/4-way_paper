library(tidyverse)
library(broom)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_Publication())


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


#save.image('BGA_Media.RData')
#load('BGA_Media.RData')



odir<-'Summary_Media'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


medias<-c('Bacto peptone','Soy peptone','LB','MRS')

metf<-c("0","25","50","75","100","150")


MeasNames<-data.frame(Measure=c("AUC","logAUC","Dph"),
                      MeasName=c("Growth AUC, OD*h","Growth log AUC, log2(OD*h)","Growth rate, doubling/h"),
                      MeasShort=c("Growth AUC","Growth log AUC","Growth rate") )


data<-read_csv('Media growth assays/Media_data/Summary.csv') %>%
  select(Plate:TReplicate,AUC=Int_600nm_f,logAUC=Int_600nm_log,Dph=a_log,-Data) %>%
  filter(!is.na(Media)) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Media=factor(Media,levels=medias),
         Metformin_mM=factor(Metformin_mM,levels=metf )) %>%
  gather(Measure,Value,AUC:Dph) %>%
  group_by(Media,Measure) %>%
  mutate(Ref_M=mean(Value[Metformin_mM=='0'])) %>%
  group_by(Measure) %>%
  mutate(Norm_M=Value-Ref_M) %>%
  ungroup %>%
  gather(Normalisation,Value,Value,Norm_M) %>%
  left_join(MeasNames)


head(data)

data_ts<-read_csv('Media growth assays/Media_data/Data.csv') %>%
  filter(Data=='600nm_f' & !is.na(Media) ) %>%
  gather(Time_s,OD,contains('.0')) %>%
  select(-Data) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Time_s=as.numeric(Time_s),
         Time_h=Time_s/3600,
         Media=factor(Media,levels=medias),
         Metformin_mM=factor(Metformin_mM,levels=metf ))


head(data_ts)


Metcols <- colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data$Metformin_mM)
Metlab<-'Metformin, mM'

ggplot(data_ts,aes(x=Time_h,y=OD,color=Metformin_mM))+
  geom_line(aes(group=interaction(Replicate,TReplicate,Metformin_mM)))+
  xlab('Time, h')+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_x_continuous(breaks=seq(0,24,by=2))+
  facet_wrap(~Media)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=9,height=6)

data.sum<-data_ts %>%
  group_by(Media,Metformin_mM,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))

ggplot(data.sum,aes(x=Time_h,y=OD_Mean,color=Metformin_mM,fill=Metformin_mM))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_colour_manual(name = Metlab,values = Metcols)+
  scale_fill_manual(name = Metlab,values = Metcols)+
  scale_x_continuous(breaks=seq(0,18,by=2))+
  facet_wrap(~Media)

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Summary.pdf",sep=''),
             useDingbats=FALSE,
             width=6,height=4)





fit<-lm(logAUC~Media*Metformin_mM,data)
summary(fit)



stat<-data %>%
  group_by(Media,Measure,MeasName,MeasShort,Normalisation) %>%
  do(tidy(lm(Value~0+Metformin_mM,data=.))) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  mutate(Metformin_mM=str_extract(term,'[[:digit:]]{1,3}'),
    Metformin_mM=factor(Metformin_mM,levels=metf ),
         PE=logFC+SE,
         NE=logFC-SE,
         Prc=2^logFC*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100,
    pStars=pStars(p.value)) 


stat %>%
  filter(Measure=="logAUC" & Normalisation=="Norm_M") %>%
  ggplot(aes(x=Metformin_mM,y=Prc,color=Metformin_mM))+
  geom_line(aes(group=1))+
  geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25)+
  geom_point()+
  scale_colour_manual(name = Metlab,values = Metcols)+
  ylab("Growth AUC, %")+
  xlab("Metformin, mM")+
  geom_text(aes(label=ifelse(Metformin_mM!="0",as.character(pStars),"" )),color='black',nudge_y = 20)+
  scale_y_continuous(breaks=seq(0,100,by=25),limits=c(0,120))+
  facet_wrap(~Media)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Comparison.pdf",sep=''),
             useDingbats=FALSE,
             width=6,height=4)







Medcols <- c("red","blue4", viridis::viridis(4,option="magma"))
names(Medcols) <- medias
Medlab<-"Media"


PlotComp<-function(data,stat,meas,measname) {
  stars<-stat %>%
    filter(Measure==meas & Normalisation=="Norm_M")
  
  data %>%
    ggplot(aes(x=Metformin_mM,y=Prc,color=Media,fill=Media))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,130))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Media),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Media)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs Control")+
    scale_colour_manual(name = Medlab,values = Medcols)+
    scale_fill_manual(name = Medlab,values = Medcols)+
    ylab(paste0(measname," vs Control, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=20-as.numeric(Media)*4 ),nudge_y = 110, show.legend = FALSE)
}



stat %>%
  filter(Normalisation=="Norm_M" & Measure=="logAUC") %>% 
  PlotComp(.,stat,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) )


Compplots<-stat %>%
  filter(Normalisation=="Norm_M" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) )

Compwrap<-stat %>%
  filter(Normalisation=="Norm_M" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) + facet_wrap(~Media))


map2(paste0(odir,"/Growth_comparison_vs_Control_condensed_",as.character(Compplots$Measure),".pdf"),
     Compplots$plot,
     width=5,height=3, useDingbats=FALSE, ggsave)


map2(paste0(odir,"/Growth_comparison_vs_Control_",as.character(Compwrap$Measure),".pdf"),
     Compwrap$plot,
     width=7,height=5, useDingbats=FALSE, ggsave)


