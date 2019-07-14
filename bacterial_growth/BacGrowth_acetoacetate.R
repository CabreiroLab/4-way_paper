#Figure numbering might have been changed
library(tidyverse)
library(broom)
library(ggrepel)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir<-'Summary_Acetoacetate'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#save.image('BGA_Acetoacetate.RData')
#load('BGA_Acetoacetate.RData')



strainlist<-c('OP50-C')
metf<-c("0","50")
aclevels<-c('0','1','5','10','20')


NormNames<-data.frame(Normalisation=c("Value","Norm_A","Norm_M","Norm_AM"), NormName=c("Absolute","By Acetoacetate=0","By Metformin=0","By Acetoacetate=0 and Metformin=0") )
MeasNames<-data.frame(Measure=c("AUC","logAUC","Dph",'logDph'),
                      MeasName=c("Growth AUC, OD*h","Growth log AUC, log2(OD*h)","Growth rate, doubling/h","Growth rate, log2(doubling/h)"),
                      MeasShort=c("Growth AUC","Growth log AUC","Growth rate","Growth rate log") )



data<-read_csv('Acetoacetate/Data/Summary.csv') %>%
  filter(!is.na(Strain) ) %>%
  select(Plate:TReplicate,logAUC=Int_600nm_log,Dph=a_log,-Data) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist,labels=strainlist),
         Acetoacetate_mM=factor(Acetoacetate_mM,levels=aclevels),
         logDph=log2(Dph)) %>%
  gather(Measure,Value,logAUC,logDph) %>%
  group_by(Metformin_mM,Measure) %>%
  mutate(Ref_A=mean(Value[Acetoacetate_mM=='0']) ) %>%
  group_by(Acetoacetate_mM,Measure) %>%
  mutate(Ref_M=mean(Value[Metformin_mM=='0'])) %>%
  group_by(Measure) %>%
  mutate(Ref_AM=mean(Value[Metformin_mM=='0' & Acetoacetate_mM=='0']),
         Norm_A=Value-Ref_A,
         Norm_M=Value-Ref_M,
         Norm_AM=Value-Ref_AM
         ) %>%
  ungroup %>%
  gather(Normalisation,Value,Value,Norm_A,Norm_M,Norm_AM) %>%
  left_join(NormNames) %>%
  left_join(MeasNames) %>%
  mutate(Value=ifelse(Value %in% c(Inf,-Inf),NA,Value ))


data.sum<-data %>%
  group_by(Measure,Normalisation,Strain,Metformin_mM,Acetoacetate_mM,NormName,MeasName,MeasShort) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value)) %>%
  mutate(PE=Mean+SD,
         NE=Mean-SD,
         Prc=2^Mean*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100)

data.sum2<-multiplex(data.sum,c('Strain','Metformin_mM','Acetoacetate_mM'),dim=2)


data %>%
  select(-c(NormName,MeasName,MeasShort)) %>%
  write_csv(paste0(odir,"/Raw_data_Summary.csv"))


  
data_ts<-read_csv('Acetoacetate/Data/Data.csv') %>%
  filter(!is.na(Strain) ) %>%
  filter(Data=='600nm_f') %>%
  gather(Time_s,OD,contains('.0')) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  select(-Data) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist),
         Time_s=as.numeric(Time_s),
         Time_h=Time_s/3600)


Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin,\nmM'



ggplot(data_ts,aes(x=Time_h,y=OD,color=Metformin_mM))+
  geom_line(aes(group=interaction(Plate,Replicate,TReplicate,Metformin_mM)))+
  xlab('Time, h')+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  facet_wrap(~Acetoacetate_mM,ncol=5)


ggsave(file=paste0(odir,"/Growth_overview.pdf"),
             width=110,height=41,units='mm',scale=2,family="Arial",device=cairo_pdf)




datats.sum<-data_ts %>%
  group_by(Strain,Acetoacetate_mM,Metformin_mM,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))


Accols <- colorRampPalette(c("green4", "cyan"))(5)

names(Accols) <- levels(data_ts$Acetoacetate_mM)
Aclab<-'Acetoacetate, mM'

ggplot(datats.sum,aes(x=Time_h,y=OD_Mean,color=Acetoacetate_mM,fill=Acetoacetate_mM))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_colour_manual(name = Aclab,values = Accols)+
  scale_fill_manual(name = Aclab,values = Accols)+
  scale_x_continuous(breaks=seq(0,18,by=6))+
  facet_grid(Metformin_mM~.)+
  theme(legend.position = "top")


ggsave(file=paste0(odir,"/Growth_Summary_tiny.pdf"),
       width=25,height=41,units="mm",scale=2,family="Arial",device=cairo_pdf)



PlotBox<-function(data,yvar,ytitle,title) {
  data %>%
  ggplot()+
    aes(x=Acetoacetate_mM,y=Value,color=Acetoacetate_mM)+
    geom_hline(aes(yintercept = 0),color='red4',alpha=0.5,linetype='longdash')+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "dodge",alpha=0.2)+
    geom_point()+
    ylab(ytitle)+
    ggtitle(paste0("Normalisation: ",title))+
    xlab("Acetoacetate, mM")+
    scale_colour_manual(name = Aclab,values = Accols)+
    scale_y_continuous(breaks = seq(-10,10,by=1))+
    coord_cartesian(ylim = c(-3,3)) +
    facet_wrap(~Metformin_mM)
}


Boxplots<-data %>%
  group_by(Measure,Normalisation) %>%
  do(plot=PlotBox(.,unique(as.character(.$Measure)),unique(as.character(.$MeasName)), unique(as.character(.$NormName))))
  

map2(paste0(odir,"/Growth_comparisons_",as.character(Boxplots$Measure),"_",as.character(Boxplots$Normalisation),".pdf"),
     Boxplots$plot,
     scale=2,
     width=120,height=60,units="mm", useDingbats=FALSE, ggsave)





stat<-data %>%
  group_by(Strain,Metformin_mM,Measure,MeasName,MeasShort,Normalisation,NormName) %>%
  filter(!is.na(Value))%>%
  do(tidy(lm(Value~0+Acetoacetate_mM,data=.))) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  mutate(Acetoacetate_mM=str_replace(term,'Acetoacetate_mM',''),
         Acetoacetate_mM=factor(Acetoacetate_mM,levels=aclevels),
         PE=logFC+SE,
         NE=logFC-SE,
         Prc=2^logFC*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100,
         pStars=pStars(p.value)) 

stat2<-data %>%
  group_by(Strain,Measure,MeasName,MeasShort,Normalisation,NormName) %>%
  filter(!is.na(Value))%>%
  do(tidy(lm(Value~Metformin_mM*Acetoacetate_mM,data=.))) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  filter(term!='(Intercept)') %>%
  filter(str_detect(term,':') ) %>%
  mutate(term=str_replace_all(term,'Metformin_mM|Acetoacetate_mM',''),
         PE=logFC+SE,
         NE=logFC-SE,
         Prc=2^logFC*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100,
         pStars=pStars(p.value)) %>%
  separate(term,c('Metformin_mM','Acetoacetate_mM'),sep=':') %>%
  mutate_at(c("Metformin_mM",'Acetoacetate_mM'),as.factor)




head(data)
stat.new<-data %>%
  filter(Normalisation=='Value' ) %>%
  group_by(Measure) %>%
  lmtest("Value~Acetoacetate_mM*Metformin_mM")


stat.new %>%
  write_csv(paste0(odir,"/Stats_Summary.csv"))


stat %>%
  filter(Measure=='logAUC' & Normalisation=="Norm_AM" & Metformin_mM=='50')



stat2 %>%
  filter(Measure=='logAUC' & Normalisation=="Norm_A" )




stat %>%
  filter(Measure=='logAUC' & Normalisation=="Norm_AM") %>%
  ggplot(aes(x=Acetoacetate_mM,y=Prc,color=Acetoacetate_mM))+
  geom_hline(aes(yintercept = 100),color='red4',alpha=0.5,linetype='longdash')+
  scale_y_continuous(breaks=seq(0,400,by=25))+
  coord_cartesian(ylim=c(0,130))+
  geom_line(aes(group=Strain))+
  geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25)+
  geom_point()+
  ggtitle("Comparison vs OP50-C in control")+
  scale_colour_manual(name = Aclab,values = Accols)+
  ylab("Growth AUC vs Acetoacetate=0mM, %")+
  xlab("Acetoacetate, mM")+
  geom_text(data=stat %>%
              filter(Measure=='logAUC' & Normalisation=="Norm_A"),
            aes(label=pStars),color='black',y = 125)+
  facet_wrap(~Metformin_mM,ncol=2)



ggsave(file=paste(odir,"/Growth_Comparison_vs_Acetoacetate0.pdf",sep=''),
             useDingbats=FALSE,scale=2,
             width=120,height=60,units='mm')


PlotComp<-function(data,stat,meas,measname) {
  stars<-stat %>%
    filter(Measure==meas,Normalisation=="Norm_A")

  data %>%
    ggplot(aes(x=Acetoacetate_mM,y=Prc,color=Metformin_mM,fill=Metformin_mM))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,110))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Metformin_mM),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Metformin_mM)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs OP50-C")+
    scale_colour_manual(values = Metcols)+
    scale_fill_manual(values = Metcols)+
    labs(color='Metformin,\nmM',fill='Metformin,\nmM')+
    ylab(paste0(measname," vs Control, %") )+
    xlab("Acetoacetate, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=10-as.numeric(Metformin_mM)*4 ),nudge_y = 70, angle=90,show.legend = FALSE,size=5)
  #+theme(legend.position = "top")
}


data.sum %>%
  filter(Normalisation=="Norm_AM" & Measure=="logAUC") %>% 
  PlotComp(.,stat2,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) )

ggsave(file=paste(odir,"/Growth_comparison_vs_Metf0Acet0_condensed_logAUC_SD_tiny.pdf",sep=''),
       useDingbats=FALSE,scale=2,
       width=40,height=41,units='mm')


Compplots<-data.sum %>%
  filter(Normalisation=="Norm_AM" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat2,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) )

Compwrap<-data.sum %>%
  filter(Normalisation=="Norm_AM" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat2,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) + facet_wrap(~Metformin_mM,ncol=5))

 
map2(paste0(odir,"/Growth_comparison_vs_Metf0Acet0_condensed_",as.character(Compplots$Measure),"_SD.pdf"),
     Compplots$plot,
     width=55,height=41,scale=2,units ='mm',family="Arial",device=cairo_pdf, ggsave)


map2(paste0(odir,"/Growth_comparison_vs_Metf0Acet0_",as.character(Compwrap$Measure),"_SD.pdf"),
     Compwrap$plot,
     scale=2,
     width=110,height=40, units ='mm',family="Arial",device=cairo_pdf, ggsave)



errcolor<-'grey80'

data.sum2 %>%
  filter(x_Measure=='logAUC' & x_Normalisation=='Norm_M' &
           y_Measure=='logDph' & y_Normalisation=='Norm_M') %>%
  #filter(Strain=='OP50-C')%>%
  ggplot(aes(x=x_Prc,y=y_Prc))+
  geom_vline(xintercept = 100,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 100,color='gray70',alpha=0.5,linetype='longdash')+
  scale_x_continuous(breaks=seq(0,200,by=20))+
  scale_y_continuous(breaks=seq(0,200,by=20))+
  coord_cartesian(xlim = c(0,120),ylim=c(0,120))+
  geom_errorbar(aes(ymin=y_PrcNE,ymax=y_PrcPE),color=errcolor)+
  geom_errorbarh(aes(xmin=x_PrcNE,xmax=x_PrcPE),color=errcolor)+
  scale_colour_manual(name = Metlab,values = Metcols)+
  #scale_fill_manual(name = Metlab,values = Metcols)+
  geom_line(aes(group=interaction(Metformin_mM),color=Metformin_mM,size=Acetoacetate_mM),alpha=0.5)+
  geom_point(aes(color=Metformin_mM,size=Acetoacetate_mM))+
  labs(size="Acetoacetate,\nmM")+
  xlab('Relative growth AUC vs OP50-C Control, %')+
  ylab('Relative growth rate vs OP50-C Control, %')

#hc

ggsave(file=paste0(odir,"/Growth_AUC_and_rate_comparison.pdf"),
       device=cairo_pdf,family="Arial",
             scale=2,
             width=70,height=41,units='mm')


