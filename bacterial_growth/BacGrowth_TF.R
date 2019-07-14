#Figure numbering might have been changed
library(tidyverse)
library(broom)
library(ggrepel)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


odir<-'Summary_TF'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#save.image('BGA_TF.RData')
#load('BGA_TF.RData')



strainlist<-c('OP50-C','crp','cra','argR','ntrC','arcA','csiR','fur','gcvA','marA','mlc','nac')
metf<-c("0","25","50","75","100","150")


NormNames<-data.frame(Normalisation=c("Value","Norm_C","Norm_M","Norm_CM"), NormName=c("Absolute","By OP50-C","By Metformin=0","By OP50-C and Metformin=0") )
MeasNames<-data.frame(Measure=c("AUC","logAUC","Dph",'logDph'),
                      MeasName=c("Growth AUC, OD*h","Growth log AUC, log2(OD*h)","Growth rate, doubling/h","Growth rate, log2(doubling/h)"),
                      MeasShort=c("Growth AUC","Growth log AUC","Growth rate","Growth rate log") )

mutdata<-read_csv('Resistance_mutants_growth_assays/2017-Rosie_Mutants/Summary.csv') %>%
  filter(Strain %in% c('crp','mlc','OP50-C'))


newdata<-read_csv('TF/Data_new/Summary.csv')

newdatats<-read_csv('TF/Data_new/Data.csv')


data<-read_csv('TF/Data/Summary.csv') %>%
  rbind(newdata) %>%
  rbind(mutdata) %>%
  filter(!is.na(Strain) ) %>%
  select(Plate:TReplicate,logAUC=Int_600nm_log,Dph=a_log,-Data) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist,labels=strainlist),
         logDph=log2(Dph)) %>%
  gather(Measure,Value,logAUC,logDph) %>%
  group_by(Metformin_mM,Measure) %>%
  mutate(Ref_C=mean(Value[Strain=='OP50-C']) ) %>%
  group_by(Strain,Measure) %>%
  mutate(Ref_M=mean(Value[Metformin_mM=='0'])) %>%
  group_by(Measure) %>%
  mutate(Ref_CM=mean(Value[Metformin_mM=='0' & Strain=='OP50-C']),
         Norm_C=Value-Ref_C,
         Norm_M=Value-Ref_M,
         Norm_CM=Value-Ref_CM) %>%
  ungroup %>%
  gather(Normalisation,Value,Value,Norm_C,Norm_M,Norm_CM) %>%
  left_join(NormNames) %>%
  left_join(MeasNames) %>%
  mutate(Value=ifelse(Value %in% c(Inf,-Inf),NA,Value ))


data.sum<-data %>%
  group_by(Measure,Normalisation,Strain,Metformin_mM,NormName,MeasName,MeasShort) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value)) %>%
  mutate(PE=Mean+SD,
         NE=Mean-SD,
         Prc=2^Mean*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100)

data.sum2<-multiplex(data.sum,c('Strain','Metformin_mM'),dim=2)


data %>%
  select(-c(NormName,MeasName,MeasShort)) %>%
  write_csv(paste0(odir,"/Raw_data_Summary.csv"))


mutdatats<-read_csv('Resistance_mutants_growth_assays/2017-Rosie_Mutants/Data.csv') %>%
  filter(Strain %in% c('crp','mlc','OP50-C'))
  
  
data_ts<-read_csv('TF/Data/Data.csv') %>%
  rbind(newdatats) %>%
  rbind(mutdatats) %>%
  filter(!is.na(Strain) ) %>%
  filter(Data=='600nm_f') %>%
  gather(Time_s,OD,contains('.0')) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  select(-Data) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist),
         Time_s=as.numeric(Time_s),
         Time_h=Time_s/3600)



Metcols <- colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin,\nmM'

ggplot(data_ts,aes(x=Time_h,y=OD,color=Metformin_mM))+
  geom_line(aes(group=interaction(Plate,Replicate,TReplicate,Metformin_mM)))+
  xlab('Time, h')+
  scale_colour_manual(name = Metlab,values =Metcols)+
  scale_x_continuous(breaks=seq(0,24,by=6))+
  facet_wrap(~Strain,ncol=6)


ggsave(file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=110,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")




datats.sum<-data_ts %>%
  group_by(Strain,Metformin_mM,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))

Metlab2<-"Metformin, mM"

datats.sum %>%
  #filter(Strain %in% c("OP50-C","crp","cra","argR","ntrC","marA","nac") ) %>%
  ggplot(aes(x=Time_h,y=OD_Mean,color=Metformin_mM,fill=Metformin_mM))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_colour_manual(name = Metlab,values = Metcols)+
  scale_fill_manual(name = Metlab,values = Metcols)+
  scale_x_continuous(breaks=seq(0,18,by=6))+
  facet_wrap(~Strain,ncol=6)

ggsave(file=paste0(odir,"/Growth_Summary_horizontal.pdf"),
       width=110,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")

# ggsave(file=paste0(odir,"/Growth_Summary_horizontal.pdf"),
#        width=90,height=25,units='mm',scale=2,device=cairo_pdf,family="Arial")



PlotBox<-function(data,yvar,ytitle,title) {
  data %>%
  ggplot()+
    aes(x=Metformin_mM,y=Value,color=Metformin_mM)+
    geom_hline(aes(yintercept = 0),color='red4',alpha=0.5,linetype='longdash')+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "dodge",alpha=0.2)+
    geom_point()+
    ylab(ytitle)+
    ggtitle(paste0("Normalisation: ",title))+
    xlab("Metformin, mM")+
    scale_colour_manual(name = Metlab,values = Metcols)+
    scale_y_continuous(breaks = seq(-10,10,by=1))+
    coord_cartesian(ylim = c(-3,3)) +
    facet_wrap(~Strain)
}


Boxplots<-data %>%
  group_by(Measure,Normalisation) %>%
  do(plot=PlotBox(.,unique(as.character(.$Measure)),unique(as.character(.$MeasName)), unique(as.character(.$NormName))))
  

map2(paste0(odir,"/Growth_comparisons_",as.character(Boxplots$Measure),"_",as.character(Boxplots$Normalisation),".pdf"),
     Boxplots$plot,
     width=120,height=60,units="mm",scale=2,device=cairo_pdf,family="Arial", ggsave)





stat<-data %>%
  group_by(Metformin_mM,Measure,MeasName,MeasShort,Normalisation,NormName) %>%
  filter(!is.na(Value))%>%
  do(tidy(lm(Value~0+Strain,data=.))) %>%
  rename(SE=std.error,
         logFC=estimate) %>%
  mutate(Strain=str_replace(term,'Strain',''),
         Strain=factor(Strain,levels=strainlist,labels=strainlist),
         PE=logFC+SE,
         NE=logFC-SE,
         Prc=2^logFC*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100,
         pStars=pStars(p.value)) 


stat2<-data %>%
  filter(Measure=="logAUC" & Normalisation=="Value") %>%
  do(tidy(lm(Value~Strain*Metformin_mM,data=.))) %>%
  filter(str_detect(term,":")) %>%
  rename(logFC=estimate,SE=std.error,t.value=statistic,p=p.value) %>%
  mutate(term=str_replace_all(term,"Strain|Metformin_mM","")) %>%
  separate(term,c("Strain","Metformin_mM"),sep=":") %>%
  mutate(Comparison="Interaction",
         Strain=factor(Strain,levels=strainlist),
         pStars=pStars(p)) %>%
  ungroup


stat %>%
  write_csv(paste0(odir,"/Stats_Summary.csv"))



stat %>%
  filter(Measure=='logAUC' & Normalisation=="Norm_CM") %>%
  ggplot(aes(x=Metformin_mM,y=Prc,color=Metformin_mM))+
  geom_hline(aes(yintercept = 100),color='red4',alpha=0.5,linetype='longdash')+
  scale_y_continuous(breaks=seq(0,400,by=25))+
  coord_cartesian(ylim=c(0,130))+
  geom_line(aes(group=Strain))+
  geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25)+
  geom_point()+
  ggtitle("Comparison vs OP50-C in control",subtitle = "Significance shown for difference vs OP50-C")+
  scale_colour_manual(name = Metlab,values = Metcols)+
  ylab("Growth logAUC, %")+
  xlab("Metformin, mM")+
  geom_text(data=stat %>%
              filter(Measure=='AUC' & Normalisation=="Norm_C"),
            aes(label=pStars),color='black',y = 125)+
  facet_wrap(~Strain,ncol=6)

ggsave(file=paste(odir,"/Growth_Comparison_vs_OP50-C_Control.pdf",sep=''),
             width=120,height=60,units='mm',scale=2,device=cairo_pdf,family="Arial")


Straincols<-ggthemes::tableau_color_pal(palette = "tableau20")(12)
names(Straincols) <- strainlist


Straincols.new<-ggthemes::tableau_color_pal(palette = "tableau20")(13)
names(Straincols.new) <- strainlist.new



PlotComp<-function(data,stars,measname) {

  data %>%
    ggplot(aes(x=Metformin_mM,y=Prc,color=Strain,fill=Strain))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,175))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Strain),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Strain)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs OP50-C")+
    scale_colour_manual(name = Strainlab,values = Straincols.new)+
    scale_fill_manual(name = Strainlab,values = Straincols.new)+
    #labs(color='Strain',fill='Strain')+
    ylab(paste0(measname," vs OP50-C Control, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=50-as.numeric(Strain)*5 ),nudge_y = 135, show.legend = FALSE,size=5)
}




PlotComp2<-function(data,stars,measname) {
  
  data %>%
    ggplot(aes(x=Metformin_mM,y=Prc,color=Strain,fill=Strain))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,160))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Strain),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Strain)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs OP50-C")+
    scale_colour_manual(name = Strainlab,values = Straincols)+
    scale_fill_manual(name = Strainlab,values = Straincols)+
    labs(color='Strain',fill='Strain')+
    ylab(paste0(measname," vs OP50-C Control, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=50-as.numeric(Strain)*5 ),nudge_y = 130 , show.legend = FALSE,size=5) #nudge.y=100
}



data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure=="logAUC") %>%
  PlotComp(filter(stat,Measure=='logAUC', Normalisation=="Norm_C"),'Growth AUC' )


ggsave(file=paste(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_logAUC_fix.pdf",sep=''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


MRdata<-read_csv('Summary_resistance/Data_Summary.csv') %>%
  filter(Strain=='OP50-MR')%>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf))

MRstat<-read_csv('Summary_resistance/Stats_Summary.csv') %>%
  filter(Strain=='OP50-MR') %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf))


strainlist.new<-c('OP50-C','OP50-MR','crp','cra','argR','ntrC','arcA','csiR','fur','gcvA','marA','mlc','nac')


stat<-stat %>%
  ungroup %>%
  bind_rows(MRstat) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist.new))

stat %>%
  write_csv(paste0(odir,"/Stats_Summary_with_MR.csv"))




data.sum<-data.sum %>%
  ungroup %>%
  bind_rows(MRdata) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist.new))


FS6<-c('OP50-C','arcA','csiR','fur','gcvA','marA','mlc','nac')



PlotCompfS5<-function(data,stars,measname) {
  
  data %>%
    ggplot(aes(x=Metformin_mM,y=Prc,color=Strain,fill=Strain))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,150))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Strain),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Strain)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs OP50-C")+
    scale_colour_manual(name = Strainlab,values = Straincols.new)+
    scale_fill_manual(name = Strainlab,values = Straincols.new)+
    #labs(color='Strain',fill='Strain')+
    ylab(paste0(measname," vs OP50-C Control, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=50-as.numeric(Strain)*5 ),nudge_y = 135, show.legend = FALSE,size=5)
}


data.sum %>%
  ungroup %>%
  filter(Normalisation=="Norm_CM" & Measure=="logAUC" & Strain %in% FS6) %>%
  mutate(Strain=factor(Strain,levels=FS6)) %>%
  PlotCompfS5(filter(stat,Measure=='logAUC', Normalisation=="Norm_C" & Strain %in% FS6),'Growth AUC' )


ggsave(file=paste(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_FigS6.pdf",sep=''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


F6<-c('OP50-C','OP50-MR','crp','cra','argR','ntrC')

PlotCompf5<-function(data,stars,measname) {
  
  data %>%
    ggplot(aes(x=Metformin_mM,y=Prc,color=Strain,fill=Strain))+
    scale_y_continuous(breaks=seq(0,200,by=25))+
    coord_cartesian(ylim=c(0,150))+
    geom_hline(aes(yintercept = 100),color='gray75',linetype='longdash')+
    geom_ribbon(aes(group=interaction(Strain),
                    ymin=PrcNE,
                    ymax=PrcPE),
                alpha=0.75,
                color=NA,
                show.legend=FALSE)+
    geom_line(aes(group=interaction(Strain)))+
    #geom_errorbar(aes(ymin=PrcNE,ymax=PrcPE),width=0.25,alpha=0.5)+
    geom_point()+
    #ggtitle("Significance shown for differences vs OP50-C")+
    scale_colour_manual(name = Strainlab,values = Straincols.new)+
    scale_fill_manual(name = Strainlab,values = Straincols.new)+
    #labs(color='Strain',fill='Strain')+
    ylab(paste0(measname," vs OP50-C Control, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=50-as.numeric(Strain)*6 ),nudge_y = 110, show.legend = FALSE,size=5)
}



data.sum %>%
  ungroup %>%
  filter(Normalisation=="Norm_CM" & Measure=="logAUC" & Strain %in% F6) %>%
  mutate(Strain=factor(Strain,levels=F6)) %>%
  PlotCompf5(filter(stat,Measure=='logAUC', Normalisation=="Norm_C" & Strain %in% F6),'Growth AUC' )+
  theme(legend.position = "top")


ggsave(file=paste(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_Fig6.pdf",sep=''),
       width=35,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


FT<-c('OP50-C','crp','cra','argR','ntrC')

data.sum %>%
  ungroup %>%
  filter(Normalisation=="Norm_CM" & Measure=="logAUC" & Strain %in% FT) %>%
  mutate(Strain=factor(Strain,levels=FT)) %>%
  PlotComp(filter(stat,Measure=='logAUC', Normalisation=="Norm_C" & Strain %in% FT),'Growth AUC' )


ggsave(file=paste(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_Thesis.pdf",sep=''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



Compplots.I<-data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat2, unique(as.character(.$MeasShort)) ) )

map2(paste0(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_",as.character(Compplots$Measure),"_Interaction_SD.pdf"),
     Compplots.I$plot,
     width=55,height=41,scale=2,units ='mm',device=cairo_pdf,family="Arial",ggsave)



Compwrap<-data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure!="AUC") %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat %>% filter(Measure==unique(as.character(.$Measure)), Normalisation=="Norm_C"),unique(as.character(.$MeasShort)) ) + facet_wrap(~Strain,ncol=6)) %>%
  mutate(file=)

map2(paste0(odir,"/Growth_comparison_vs_OP50-C_Control_",as.character(Compwrap$Measure),"_SD.pdf"),
     Compwrap$plot,
     width=110,height=41,units ='mm', scale=2,device=cairo_pdf,family="Arial",ggsave)




errcolor<-'grey80'

data.sum2 %>%
  filter(x_Measure=='logAUC' & x_Normalisation=='Norm_CM' &
           y_Measure=='logDph' & y_Normalisation=='Norm_CM') %>%
  #filter(Strain=='OP50-C')%>%
  ggplot(aes(x=x_Prc,y=y_Prc))+
  geom_vline(xintercept = 100,color='gray70',alpha=0.5,linetype='longdash')+
  geom_hline(yintercept = 100,color='gray70',alpha=0.5,linetype='longdash')+
  scale_x_continuous(breaks=seq(0,200,by=20))+
  scale_y_continuous(breaks=seq(0,200,by=20))+
  coord_cartesian(xlim = c(0,120),ylim=c(0,120))+
  geom_errorbar(aes(ymin=y_PrcNE,ymax=y_PrcPE),color=errcolor)+
  geom_errorbarh(aes(xmin=x_PrcNE,xmax=x_PrcPE),color=errcolor)+
  scale_colour_manual(name = Strainlab,values = Straincols)+
  scale_fill_manual(name = Strainlab,values = Straincols)+
  geom_line(aes(group=interaction(Strain),color=Strain,size=Metformin_mM),alpha=0.5)+
  geom_point(aes(color=Strain,size=Metformin_mM))+
  labs(size="Metformin,\nmM")+
  xlab('Relative growth AUC vs OP50-C Control, %')+
  ylab('Relative growth rate vs OP50-C Control, %')

#hc

ggsave(file=paste(odir,"/Growth_AUC_and_rate_comparison.pdf",sep=''),
             width=100,height=80,units='mm',scale=2,device=cairo_pdf,family="Arial")


