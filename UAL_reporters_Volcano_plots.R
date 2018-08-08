library(tidyverse)
library(PFun)
library(broom)
library(viridis) #Nice colors
library(readxl)
library(ggrepel)


library(PFun)


#New default theme
theme_set(theme_PN(base_size = 12))
scale_colour_discrete <- ggthemes::scale_colour_tableau
scale_fill_discrete <- ggthemes::scale_fill_tableau






setwd("~/Dropbox/Projects/2015-Metformin/GFP_reporters/Summary/")


#Volcano plots
cnames<-c('SM-S'='OP50-C treatment','RM-R'='OP50-MR treatment','R-S'='Strain difference','SM-S-(RM-R)'='Longevity effect')

data <- read_csv('All_results.csv') %>%
  group_by(Contrast) %>%
  arrange(desc(logFDR)) %>%
  mutate(Name=ifelse(row_number()<=10 & FDR<0.05,Promoter,NA ))


unique(data$Comparison)

data %>%
  filter(Contrast %in% c('PM-PC','PGM-PGC_adj','PGC-PC')) %>%
  ggplot(aes(x=logFC,y=logFDR))+
  geom_hline(yintercept = -log10(0.05),color='red4',alpha=0.5)+
  geom_point(size=1,color='gray50' )+
  geom_text_repel(aes(label=Name),color='red2')+
  ylab('-log10(FDR)')+
  facet_wrap(~Description,ncol = 3)

ggsave(file='UAL_Volcano_plot_clean.pdf',
       width=120,height=60,units='mm',scale=2,device=cairo_pdf,family="Arial")




