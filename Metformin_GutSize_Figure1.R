library(tidyverse)
library(readxl)

library(RColorBrewer)
library(ggrepel)


devtools::install_github("PNorvaisas/PFun")
library(PFun)



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
           #legend.title = element_text(face="italic"),
           #plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

theme_set(theme_Publication())


cwd<-"~/Dropbox/Projects/Metformin_project/Figures/Figure_1/Subfigures/1e/"
setwd(cwd)


#save.image('Gut_size.RData')
#load('BGA_BacGrowth_Gut.RData')



odir<-'.'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



data<-read_xlsx('Chart_data.xlsx',sheet='Data') %>%
  mutate_at(c('Strain','Metformin_mM','Day'),as.factor)

stars<-read_xlsx('Chart_data.xlsx',sheet='Stars') %>%
  mutate_at(c('Contrast','Strain','Metformin_mM','Day'),as.factor)

stars_day<-stars %>%
  filter(Contrast=='Day')

stars_treat<-stars %>%
  filter(Contrast=='Treatment')

data %>%
  ggplot(aes(x=Day,y=Size,color=Metformin_mM))+
  stat_summary(aes(group=interaction(Metformin_mM)), geom="line")+
  geom_errorbar(aes(ymin=Size-SE,ymax=Size+SE),width = 0.2)+
  geom_point(size=2)+
  labs(y = expression(paste("Gut size (normalised), ",mu, "m"^2)),
       x='Day',
       color='Metformin, mM')+
  #scale_y_continuous(breaks=0:10,limits=c(2,4.2))+
  geom_text(data=stars_day,aes(x=Day,y=ifelse(Metformin_mM=='0',32,35),label=as.character(Stars)),nudge_x = -0.5)+
  geom_text(data=stars_treat,aes(x=Day,label=as.character(Stars)),y=27,color='black')+
  facet_grid(~Strain)


dev.copy2pdf(device=cairo_pdf,
             useDingbats=FALSE,
             file=paste(odir,"/Gut_size.pdf",sep=''),
             width=5,height=3)

