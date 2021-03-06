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


odir<-'Summary_PTS'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#save.image('BGA_PTS.RData')
#load('BGA_PTS.RData')


strainlist<-c('OP50-C','OP50-MR','crp','crr','cyaA','mlc','cpdA','ptsG')
metf<-c("0","25","50","75","100","150")


NormNames<-data.frame(Normalisation=c("Value","Norm_C","Norm_M","Norm_CM"), NormName=c("Absolute","By OP50-C","By Metformin=0","By OP50-C and Metformin=0") )
MeasNames<-data.frame(Measure=c("AUC","logAUC","Dph"),
                      MeasName=c("Growth AUC, OD*h","Growth log AUC, log2(OD*h)","Growth rate, doubling/h"),
                      MeasShort=c("Growth AUC","Growth log AUC","Growth rate") )


data<-read_csv('Resistance_mutants_growth_assays/2017-Rosie_Mutants/Summary.csv') %>%
  select(Plate:TReplicate,AUC=Int_600nm_f,logAUC=Int_600nm_log,Dph=a_log,-Data) %>%
  mutate_at(vars(Plate:TReplicate),as.factor) %>%
  mutate(Metformin_mM=factor(Metformin_mM,levels=metf),
         Strain=factor(Strain,levels=strainlist,labels=strainlist)) %>%
  gather(Measure,Value,AUC:Dph) %>%
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
  left_join(MeasNames)


data.sum<-data %>%
  group_by(Measure,Normalisation,Strain,Metformin_mM,NormName,MeasName,MeasShort) %>%
  summarise(Mean=mean(Value),
            SD=sd(Value)) %>%
  mutate(PE=Mean+SD,
         NE=Mean-SD,
         Prc=2^Mean*100,
         PrcNE=2^NE*100,
         PrcPE=2^PE*100)

  
data_ts<-read_csv('Resistance_mutants_growth_assays/2017-Rosie_Mutants/Data.csv') %>%
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
  scale_x_continuous(breaks=seq(0,24,by=2))+
  facet_wrap(~Strain)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_overview.pdf",sep=''),
             width=9,height=6)




datats.sum<-data_ts %>%
  group_by(Strain,Metformin_mM,Time_h) %>%
  summarise(OD_Mean=mean(OD),
            OD_SD=sd(OD))


datats.sum %>%
  filter(!Strain %in% c("OP50-MR","mlc")) %>%
  ggplot(aes(x=Time_h,y=OD_Mean,color=Metformin_mM,fill=Metformin_mM))+
  geom_line()+
  geom_ribbon(aes(ymin=OD_Mean-OD_SD,
                  ymax=OD_Mean+OD_SD),alpha=0.5,color=NA)+
  xlab('Time, h')+
  ylab('OD')+
  scale_colour_manual(name = Metlab,values = Metcols)+
  scale_fill_manual(name = Metlab,values = Metcols)+
  scale_x_continuous(breaks=seq(0,18,by=6))+
  facet_wrap(~Strain,ncol=3)

ggsave(file=paste0(odir,"/Growth_Summary_horizontal.pdf"),
             width=110,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")


ggsave(file=paste0(odir,"/Growth_Summary_horizontal_tiny.pdf"),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")



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
     width=9,height=6, useDingbats=FALSE, ggsave)





stat<-data %>%
  group_by(Metformin_mM,Measure,MeasName,MeasShort,Normalisation,NormName) %>%
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

stat %>%
  write_csv(paste0(odir,"/Stats_Summary.csv"))


stat %>%
  filter(Measure=='AUC' & Normalisation=="Norm_CM") %>%
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
  facet_wrap(~Strain)


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Growth_Comparison_vs_OP50-C_Control.pdf",sep=''),
             useDingbats=FALSE,
             width=9,height=6)

install.packages("viridis")
library(viridis)




Straincols <- c("red","blue4", colorRampPalette(c("orange", "black"))(6))
Straincols <- c("red4","blue4", rainbow(6))
Straincols <- c("red4","blue4", terrain.colors(6))




Straincols <- c("red","blue", viridis::viridis(5,option="inferno"))
names(Straincols) <- setdiff(strainlist,c("OP50-MR","mlc"))
Strainlab<-"Strain"




PlotComp<-function(data,stat,meas,measname) {
  stars<-stat %>%
    filter(Measure==meas,Normalisation=="Norm_C")

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
    # scale_fill_viridis(discrete=TRUE)+
    # scale_color_viridis(discrete=TRUE)+
    # scale_colour_brewer(palette="Set1")+
    # scale_fill_brewer(palette="Set1")+
    ylab(paste0(measname," vs OP50-C, %") )+
    xlab("Metformin, mM")+
    geom_text(data=stars,
              aes(label=pStars,y=40-as.numeric(Strain)*4 ),nudge_y = 120, show.legend = FALSE,size=5)
}


data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure=="logAUC" & !Strain %in% c("OP50-MR","mlc")) %>% 
  PlotComp(.,stat %>% filter(!Strain %in% c("OP50-MR","mlc")) ,unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) )



#Vs OP50-C
Compplots<-data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure!="AUC"  & !Strain %in% c("OP50-MR","mlc")) %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat %>% filter(!Strain %in% c("OP50-MR","mlc")),unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) )

Compwrap<-data.sum %>%
  filter(Normalisation=="Norm_CM" & Measure!="AUC"  & !Strain %in% c("OP50-MR","mlc")) %>%
  group_by(Measure) %>%
  do(plot=PlotComp(.,stat %>% filter(!Strain %in% c("OP50-MR","mlc")),unique(as.character(.$Measure)),unique(as.character(.$MeasShort)) ) + facet_wrap(~Strain))

 


map2(paste0(odir,"/Growth_comparison_vs_OP50-C_Control_condensed_",as.character(Compplots$Measure),"_SD.pdf"),
     Compplots$plot,
     width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial",ggsave)


map2(paste0(odir,"/Growth_comparison_vs_OP50-C_Control_",as.character(Compwrap$Measure),"_SD.pdf"),
     Compwrap$plot,
     width=110,height=82,units='mm',scale=2,device=cairo_pdf,family="Arial", ggsave)
