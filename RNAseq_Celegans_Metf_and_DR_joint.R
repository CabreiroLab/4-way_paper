library(tidyverse)
library(ggthemes)


#devtools::install_github("PNorvaisas/PFun")
library(PFun)

# theme_Publication <- function(base_size=14) {
#   (theme_foundation(base_size=base_size)
#    + theme(plot.title = element_text(face = "bold",
#                                      size = rel(1.2), hjust = 0.5),
#            text = element_text(),
#            panel.background = element_rect(colour = NA),
#            plot.background = element_rect(colour = NA),
#            #panel.border = element_rect(colour = NA),
#            axis.title = element_text(face = "bold",size = rel(1)),
#            axis.title.y = element_text(angle=90,vjust =2),
#            axis.title.x = element_text(vjust = -0.2),
#            axis.text = element_text(), 
#            axis.line = element_line(colour="black"),
#            axis.ticks = element_line(),
#            panel.grid.major = element_line(colour="#f0f0f0"),
#            panel.grid.minor = element_blank(),
#            legend.key = element_rect(colour = NA),
#            #legend.position = "bottom",
#            #legend.direction = "horizontal",
#            #legend.key.size= unit(0.2, "cm"),
#            #legend.margin = unit(0, "cm"),
#            legend.title = element_text(face="italic"),
#            #plot.margin=unit(c(10,5,5,5),"mm"),
#            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#            strip.text = element_text(face="bold")
#    ))
#   
# }


theme_set(theme_Publication())


setwd("~/Dropbox/Projects/2015-Metformin/RNAseq")


odir<-'Summary_Celegans_Metf_DR'
#odir<-'~/Desktop/RNAseq'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")


#save.image('Metformin_Celegans_RNAseq_Metf_and_DR.RData')
#load('Metformin_Celegans_RNAseq_Metf_and_DR.RData')


#Enrichment data
metf.KEGG.old<-read_csv('Celegans_metformin/Results/All_results_KEGG.csv') %>%
  filter(Threshold=='None') %>%
  mutate(Type='Metf_Old')


metf.KEGG.new<-read_csv('Celegans_metformin/Results_1thrs_newannot/All_results_KEGG.csv') %>%
  filter(Threshold=='None') %>%
  mutate(Type='Metf_New')


DR.KEGG<-read_csv('Celegans_DR_Heintz/Results_no_polyA_trimming_conservative/All_results_KEGG.csv') %>%
  filter(Threshold=='None') %>%
  mutate(Type='DR')




enrbrks<-c(0,-log(0.05,10),2,3,4,5,100)
enrlbls<-c('N.S.','<0.05','<0.01','<0.001','<0.0001','<0.000001')

enrcols<-colorRampPalette(c("gray90","steelblue1","blue4"))(n = length(enrbrks))



enr<-metf.KEGG.old %>%
  rbind(metf.KEGG.new) %>%
  rbind(DR.KEGG) %>%
  mutate(Threshold=NULL) %>%
  rename(Count_Up=Up,Count_Down=Down,
         FDR_Up=P.Up,FDR_Down=P.Down) %>%
  gather(Measure,Value,Count_Up:FDR_Down) %>%
  separate(Measure,c('Measure','Direction')) %>%
  spread(Measure,Value) %>%
  unite(CD,Comparison,Direction,remove = FALSE) %>%
  mutate(logFDR=abs(-log10(FDR)),
         FDRbin=cut(logFDR,breaks=enrbrks,labels=enrlbls,right=FALSE),
         Direction=factor(Direction,levels=c('Up','Down'))) %>%
  select(Type,Direction,Comparison,CD,ID,Description,everything())


clustorder<-function(data,rows,variable,value,dst.method='euclidean',cl.method='ward.D2',fill=0,reverse=FALSE,scalesel='none'){
  heatsum<-data %>%
    ungroup %>%
    select_(rows,variable,value) %>%
    spread_(variable,value) %>%
    data.frame
  
  heatsum[is.na(heatsum)]<-fill
  
  rownames(heatsum)<-heatsum[,rows]
  heatsum[,rows]<-NULL
  
  if (scalesel=='none') {
    heatmat<-as.matrix(heatsum)
  } else if (scalesel=='col') {
    heatmat<-scale(as.matrix(heatsum))
  } else if (scalesel=='row') {
    heatmat<-t(scale(t(as.matrix(heatsum))))
  } else {
    stop(paste0('Unknown scaling!: ',scalsel), call. = FALSE)
  }
  
  d<-dist(heatmat,method = dst.method)
  h<-hclust(d,method=cl.method)
  ordered<-rownames(heatsum[h$order,])
  if (reverse){
    ordered<-rev(ordered)
  }
  return(ordered)
}




enr %>%
  filter(Type=='DR' & !ID %in% DRexclude$ID) %>%
  clustorder(.,'Description','CD','logFDR')



heatsum<-datasel %>%
  select_('Description','CD','logFDR') %>%
  spread_('CD','logFDR') %>%
  data.frame



enr %>%
  filter(is.na(FDRbin))


unique(enr$Comparison)

sel.conts<-list('Metf_Old'=c('SM-S','RM-R',"SM-S-(RM-R)"),
             'Metf_New'=c('SM-S','RM-R',"SM-S-(RM-R)"),
             'DR'=c('DR_3','Age_N2','Age_eat2','DR_15','DR_27'))

DRexclude<-enr %>%
  filter(str_detect(Description,'Selenocompound|Glycosaminoglycan')) %>%
  group_by(Description,ID) %>%
  summarise()



enr %>%
  filter(Type=='DR' & !ID %in% DRexclude$ID) %>%
  filter(Comparison %in% sel.conts[[unique(as.character(Type))]]) %>%
  mutate(Comparison=factor(Comparison,levels=sel.conts[[unique(as.character(Type))]])) %>%
  group_by(Description) %>%
  mutate(Select=sum(FDR<0.05)>1) %>%
  ungroup %>%
  filter(Select)%>%
  mutate(Description=factor(Description,levels=clustorder(.,'Description','CD','logFDR',reverse=TRUE))) %>%
  ggplot(aes(x=Direction,y=Description))+
  geom_tile(aes(fill=FDRbin))+
  scale_fill_manual(values =enrcols)+
  xlab("Direction of change")+
  ylab('KEGG pathway')+
  labs(fill='FDR')+
  theme(axis.ticks=element_blank(),
          panel.border=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(colour = NA),
          axis.line.x = element_line(colour = NA),
          axis.line.y = element_line(colour = NA),
          strip.text = element_text(colour = 'black', face='bold',size=10),
          axis.text.x= element_text(face='bold', colour='black', size=10, angle = 90, hjust = 1))+
  facet_grid(~Comparison)


dev.copy2pdf(device=cairo_pdf,file=paste(odir,'/Enrichment_heatmap_DR_KEGG_aging.pdf',sep = ''),
             width=8,height=6,useDingbats=FALSE)



# 
# unique(enr$Type)
# 
# enr %>%
#   filter(Type=='Metf_Old' & !ID %in% DRexclude$ID) %>%
#   #filter(Type=='DR' & !ID %in% DRexclude$ID) %>%
#   filter(Comparison %in% sel.conts[[unique(as.character(Type))]]) %>%
#   mutate(Comparison=factor(Comparison,levels=sel.conts[[unique(as.character(Type))]])) %>%
#   group_by(Description) %>%
#   mutate(Select=sum(FDR<0.05)>1) %>%
#   ungroup %>%
#   filter(Select)%>%
#   mutate(Description=factor(Description,levels=clustorder(.,'Description','CD','logFDR',reverse=TRUE))) %>%
#   ggplot(aes(x=Direction,y=Description))+
#   geom_tile(aes(fill=FDRbin))+
#   scale_fill_manual(values =enrcols)+
#   xlab("Direction of change")+
#   ylab('KEGG pathway')+
#   labs(fill='FDR')+
#   theme(axis.ticks=element_blank(),
#         panel.border=element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         axis.line = element_line(colour = NA),
#         axis.line.x = element_line(colour = NA),
#         axis.line.y = element_line(colour = NA),
#         strip.text = element_text(colour = 'black', face='bold',size=10),
#         axis.text.x= element_text(face='bold', colour='black', size=10, angle = 90, hjust = 1))+
#   facet_grid(~Comparison)
# 
# 
# 
# 
# metf.GO<-read_csv('Celegans_metformin/Results/All_results_GO.csv') %>%
#   filter(Comparison %in% c('RM-R','SM-S','SM-S-(RM-R)') & Threshold=='None') %>%
#   rename(Count_Up=Up,Count_Down=Down,
#          FDR_Up=P.Up,FDR_Down=P.Down) %>%
#   gather(Measure,Value,Count_Up:FDR_Down) %>%
#   separate(Measure,c('Measure','Direction')) %>%
#   spread(Measure,Value) %>%
#   unite(CD,Comparison,Direction,remove = FALSE) %>%
#   mutate(logFDR=abs(-log10(FDR)),
#          FDRbin=cut(logFDR,breaks=enrbrks,labels=enrlbls,right=FALSE),
#          Direction=factor(Direction,levels=c('Up','Down'))) %>%
#   select(Direction,Comparison,CD,ID,Description,everything()) %>%
#   select(Ont,Description,ID,CD,FDR) %>%
#   spread(CD,FDR)


