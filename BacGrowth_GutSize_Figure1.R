library(tidyverse)
library(readxl)

library(RColorBrewer)
library(ggrepel)


devtools::install_github("PNorvaisas/PFun")
library(PFun)


theme_set(theme_Publication())


cwd<-"~/Dropbox/Projects/Metformin_project/Bacterial Growth Assays/"
setwd(cwd)


#save.image('BGA_Gut_size.RData')
#load('BGA_Gut_size.RData')



odir<-'Summary_Gut_size'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



data<-read_xlsx('Gut_size_PPT/Chart_data.xlsx',sheet='Data') %>%
  mutate_at(c('Strain','Metformin_mM','Day'),as.factor)

stars<-read_xlsx('Gut_size_PPT/Chart_data.xlsx',sheet='Stars') %>%
  mutate_at(c('Contrast','Strain','Metformin_mM','Day'),as.factor)

stars_day<-stars %>%
  filter(Contrast=='Day')

stars_treat<-stars %>%
  filter(Contrast=='Treatment')



Metcols <- c("#FF0000","#32006F")#colorRampPalette(c("red", "blue4"))(6)
names(Metcols) <- levels(data_ts$Metformin_mM)
Metlab<-'Metformin, mM'


data %>%
  ggplot(aes(x=Day,y=Size,color=Metformin_mM))+
  stat_summary(aes(group=interaction(Metformin_mM)), geom="line")+
  geom_errorbar(aes(ymin=Size-SE,ymax=Size+SE),width = 0.2)+
  geom_point(size=2)+
  scale_colour_manual(name = Metlab,values =Metcols)+
  labs(y = expression(paste("Gut size (normalised), ",mu, "m"^2)),
       x='Day',
       color='Metformin, mM')+
  #scale_y_continuous(breaks=0:10,limits=c(2,4.2))+
  geom_text(data=stars_day,aes(x=Day,y=ifelse(Metformin_mM=='0',32,35),label=as.character(Stars)),nudge_x = -0.5,legend=FALSE)+
  geom_text(data=stars_treat,aes(x=Day,label=as.character(Stars)),y=27,color='black',legend=FALSE)+
  facet_grid(~Strain)


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Gut_size.pdf",sep=''),
             width=5,height=3)

