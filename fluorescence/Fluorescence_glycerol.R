#Figure numbering might have been changed
library(tidyverse)

#devtools::install_github("PNorvaisas/PFun")
library(PFun)

setwd("~/Dropbox/Projects/Metformin_project/Fluorescence microscopy/")


odir<-'Summary_Glycerol'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#load("Fluorescence_Glycerol.RData")
#save.image('Fluorescence_Glycerol.RData')

#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau

Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin,\nmM'

Slevels<-c('OP50-C','OP50-C-Gly','glpK','glpK-Gly')
Slabels<-c('OP50-C','OP50-C+Glyc','dglpK','dglpK+Glyc')
strains<-c('OP50-C','glpK')

data <- read_csv('All_raw_data.csv') %>%
  filter(Type %in% c('Glycerol')) %>%
  mutate_at(c('Type','Gene','Strain','Metformin_mM','Supplement','SGroup','Condition','ID','Replicate'),as.factor) %>%
  mutate(SGroup=factor(SGroup,levels=Slevels,labels=Slabels),
         Strain=factor(Strain,levels=strains),
         Supplement=factor(Supplement,levels=c('None','Glycerol'),labels=c('','+ Glycerol') ))



plotBox<-function(data,yval,ylb){
  data %>%
  ggplot+
    aes_string(x="SGroup",y=yval,color="Metformin_mM")+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_jitter(width=0.25)+
    ylab(ylb)+
    xlab('Bacterial strain + Supplement')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_colour_manual(name = Metlab,values =Metcols)+
    facet_grid(Gene~Type+Replicate,scale="free_y")
}

plotBoxC<-function(data,yval,ylb){
  data %>%
    ggplot+
    aes_string(x="SGroup",y=yval,color="Metformin_mM")+
    stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
    geom_jitter(width=0.25)+
    ylab(ylb)+
    xlab('Bacterial strain + Supplement')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_colour_manual(name = Metlab,values =Metcols)+
    facet_grid(Gene~.,scale="free_y")
}



Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBox(.,"Value",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_raw_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=10,height=12, useDingbats=FALSE, ggsave)


Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBox(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=10,height=12, useDingbats=FALSE, ggsave)
  

Fluorplots<-data %>%
  group_by(Measure) %>%
  do(plot=plotBoxC(.,"Norm",as.character(unique(.$Measure))) )


map2(paste0(odir,"/Fluorescence_norm_joined_",as.character(Fluorplots$Measure),".pdf"),
     Fluorplots$plot,
     width=6,height=12, useDingbats=FALSE, ggsave)


data %>%
  filter(Measure=='Log') %>%
  ggplot(aes(x=interaction(Strain,Supplement,sep = ' '),y=Norm,color=Metformin_mM))+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_jitter(width=0.25)+
  ylab('C. elegans fluorescence (normalised), log2 A.U.')+
  xlab('Bacterial strain + Supplement')+
  labs(color='Metformin, mM')+
  scale_y_continuous(breaks=seq(-10,10,by=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_colour_manual(name = Metlab,values =Metcols)+
  facet_wrap(~Gene,scale="free_y")


dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Fluorescence_comparison.pdf",sep=''),
             width=10,height=6, useDingbats=FALSE)



#
sum.c<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene,ID,Strain,Supplement,Metformin_mM,SGroup) %>%
  summarise(SD=sd(Value,na.rm = TRUE),Mean=mean(Value,na.rm = TRUE)) %>%
  mutate(Index=paste(Gene,ID),
         VarPrc=ifelse(is.na(SD) ,Inf, (2^(SD)-1)*100 ) ) %>%
  arrange(VarPrc) %>%
  data.frame %>%
  mutate(Index=factor(Index, levels=Index,labels=Index))

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=10, useDingbats=FALSE)



ggplot(sum.c,aes(x=Index,y=VarPrc))+
  geom_point()+
  scale_y_continuous(breaks = seq(0,1000,by=25),limits=c(0,100))+
  ylab('In-group variation, %')+
  coord_flip()

dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Ingroup_variation_Percentage.pdf",sep=''),
             width=7,height=10, useDingbats=FALSE)



#Linear modelling

contrasts<-read.contrasts('!Contrasts_fluorescence.xlsx')

contrasts$Contrasts.table
contrasts.desc<-contrasts$Contrasts.table%>%
  select(Description:Metformin_mM) %>%
  mutate(SGroup=factor(SGroup,levels=Slevels,labels=Slabels))


contr.matrix<-contrasts$Contrasts.matrix
contr.matrix

results.all<-data %>%
  filter(Measure=='Log') %>%
  group_by(Gene) %>%
  do(hypothesise(.,"Norm~0+ID",contr.matrix)) %>%
  getresults(contrasts.desc)


results<-results.all$results
results.castfull<-results.all$castfull
results.cast<-results.all$cast
results.multi<-results.all$multi


write.csv(results,paste(odir,'/All_results.csv',sep=''),row.names = FALSE)
write.csv(results.cast,paste(odir,'/All_results_sidebyside.csv',sep=''),row.names = FALSE)
write.csv(results.castfull,paste(odir,'/All_results_sidebyside_full.csv',sep=''),row.names = FALSE)


#Generate table for heatmap

heatsum<-results %>%
  filter(Contrast_type %in% c('Treatment','Interaction' ) ) %>%
  select(Description,Gene,logFC) %>%
  spread(Description,logFC) %>%
  data.frame(check.names = FALSE,check.rows = FALSE)


rownames(heatsum)<-heatsum$Gene
heatsum$Gene<-NULL


d<-dist(as.matrix(heatsum),method = "euclidean")
h<-hclust(d)
ordmet<-rownames(heatsum[h$order,])

if (length(ordmet)!=length(unique(ordmet))){
  print("Non unique metabolites!")
}


max(heatsum)
min(heatsum)

amp<-6

minv<- -amp
maxv<- amp

nstep<-maxv-minv
nstep<-8

clrbrks<-seq(-amp,amp,by=2)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)




results %>%
  #filter(Contrast %in% c('OP50_T','OP50-MR_T','OP50-MR_I')) %>%
  ggplot(aes(x=Description,y=Gene))+
  geom_tile(aes(fill=logFC))+
  #geom_point(aes(size=FDRc,colour=logFC),alpha=0.9)+
  theme_minimal()+
  geom_text(aes(label=as.character(pStars)))+
  #scale_size_discrete(range = c(2,4))+#,breaks=brks
  scale_fill_gradientn(colours = clrscale,
                       breaks=clrbrks,limits=c(-amp,amp))+
  #scale_fill_gradient2(low = "purple", mid = "gray", high = "red", midpoint = 0, breaks = clrbrks)+
  xlab("Comparison")+
  theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Comparison_Heatmap_All.pdf',sep = ''),
             width=6,height=4,useDingbats=FALSE)



descriptions<-c('SM-S'='Treatment effect on OP50','RM-R'='Treatment effect on OP50-MR','SM-S-(RM-R)'='Interaction between treatment and OP50-MR')
translations<-c('F44G3.2'='argk-1','C05D11.7'='atgl-1')

genes<-unique(as.character(results$Gene))

showstats<-results %>%
  filter(Contrast %in% c('OP50_T',"glpK_T","OP50Gly_T","glpKGly_T","OP50Gly_I","gklpKGly_I"))


data %>%
  filter(Measure=='Log') %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(ymin=min(NormAbs),
            ymax=max(NormAbs),
            logymin=log2(ymin),
            logymax=log2(ymax))

blank_data <- data %>%
  filter(Measure=='Log') %>%
  group_by(Gene,Strain,SGroup) %>%
  summarise(NormAbs=2^(min(Norm)+(max(Norm)-min(Norm))*1.3 ) )%>% #ymin+2^( log2( max(NormAbs)/ymin )*1.5 )  #ymin+(max(NormAbs)-ymin)*1.5
  mutate(Metformin_mM="0",
         Metformin_mM=factor(Metformin_mM,levels=c("0","50"),labels=c("0",'50')))

hj<-0.8
vj<-2
nx<--0.5

quartz()
data %>%
  filter(Measure=='Log') %>%
  ggplot+
  aes(x=SGroup,y=NormAbs,color=Metformin_mM)+
  geom_jitter(aes(fill=Metformin_mM),width=0.25,size=1,alpha=0.5)+
  stat_summary(fun.data=MinMeanSDMax, geom="boxplot",position = "identity",alpha=0.5) +
  geom_blank(data = blank_data,aes(x=SGroup,y=NormAbs))+
  scale_y_continuous(trans = 'log2',
                     breaks = scales::trans_breaks('log2', function(x) 2^x),
                     labels = scales::trans_format('log2', scales::math_format(2^.x))) + 
  scale_colour_manual(name = "Metformin, mM",values =Metcols)+
  scale_fill_manual(name = "Metformin, mM",values =Metcols)+
  geom_text(data=showstats %>% filter(Contrast_type=="Treatment"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,nudge_x=nx, vjust=vj,angle=45)+
  geom_text(data=showstats %>% filter(Contrast_type=="Interaction"),aes(label=pStars,y=Inf),show.legend = FALSE,size=5,color="red",nudge_x=nx, vjust=vj+0.5,angle=45,hjust=hj)+
  ylab('Mean fluorescence per worm, a.u.')+
  xlab('Strain & Supplement')+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Gene,scale = "free_y",ncol=3)


ggsave(file=paste(odir,'/Fluorescence_logScale2.pdf',sep = ''),
       width=55,height=41,units='mm',scale=2,device=cairo_pdf,family="Arial")
