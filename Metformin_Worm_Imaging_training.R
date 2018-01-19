library(gplots)
library(ggplot2)
library(xlsx)
library(plyr)
library(reshape2)
library(ggbiplot)
library(ggrepel)
library(car)
library(heatmap3)
library(gtools)

library(RColorBrewer)
library(grid)
library(gridExtra)

library(plot3D)

# library(rgl)
# library(pca3d)

library(multcomp)
library(contrast)
library(scales)


library(ggrepel)


MinMeanSDMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


theme_set(theme_light())

theme_update(panel.background = element_rect(colour = "black"),
             axis.text = element_text(colour = "black"))

setwd("~/Dropbox/Projects/2015-Metformin/Biolog_Met_NGM/Celegans/")



theme_set(theme_light())



odir<-'Summary'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")



theme_set(theme_light())


thres.man<-read.csv('Thresholds_6-28_1625_from_manual.csv',sep=',')
thres.ran<-read.csv('Thresholds_6-28_1652_random.csv',sep=',')
thres<-rbind(thres.man,thres.ran)


peaks.man<-read.csv('Thresholds_data_7-3_1314_manual_peaks.csv',sep='\t')
peaks.ran<-read.csv('Thresholds_data_7-3_1317_random_peaks.csv',sep='\t')
peaks.all<-rbind(peaks.man,peaks.ran)

data.man<-read.csv('Thresholds_data_6-28_1633_from_manual.csv',sep='\t')
data.ran<-read.csv('Thresholds_data_6-28_1654_random.csv',sep='\t')
data<-rbind(data.man,data.ran)

data.m<-melt(data,id.vars = c('Replicate','Plate','Index','Layer'),
             variable.name = 'Brightness',value.name = 'Count')
data.m$Brightness<-as.numeric(gsub('X','',data.m$Brightness))


data.hv<-dcast(data.m,Replicate+Plate+Index+Brightness~Layer,value.var = 'Count',mean)


comb<-merge(data.hv,thres,by=c('Replicate','Plate','Index'))

findpeak<-function(x,y){
  peak<-x[y==max(y)]
  return(peak)
}


findpeak(comb$Brightness,comb$h)
peaks<-ddply(comb,.(Replicate,Plate,Index,Threshold),summarise,Peak_h=findpeak(Brightness,h),Peak_v=findpeak(Brightness,v) )


data.old<-read.csv('Thresholds_all.csv',sep=',')
data.old$ID<-paste(data.old$Replicate,data.old$Plate,data.old$Index,sep='_')



data.new<-read.csv('Thresholds_7-14_1737_manual_part4.csv',sep=',')
data.new$ID<-paste(data.new$Replicate,data.new$Plate,data.new$Index,sep='_')
data.new$Replicate<-as.factor(data.new$Replicate)

intersect(data.old$ID,data.new$ID)



data.oldf<-subset(data.old,!ID %in% intersect(data.old$ID,data.new$ID))



data.all<-rbind(data.oldf,data.new)
write.csv(data.all,'Thresholds_all_3.csv',row.names = FALSE)







data.nclean<-subset(data.new,Threshold!=0 & Mu>0)
data.rep56<-subset(data.new,Replicate %in% c(5,6) & Threshold!=0)






ggplot(data.nclean,aes(x=Mu,y=Threshold,color=Replicate))+
  geom_abline(slope=0.995318,intercept=0.020434)+
  geom_abline(slope=0.995318,intercept=0.05,color='red')+
  geom_point()





model<-lm(Threshold~Peak_h,peaks)
result<-summary(model)
result
result$r.squared


ggplot(peaks,aes(x=Peak_h,y=Threshold,shape=as.factor(Replicate)))+
  geom_point()+
  geom_abline(aes(intercept=0.023005,slope=1.001005))
dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hue_threshold_model.pdf",sep=''),
             width=9,height=6)




testing<-merge(peaks.all,thres,by=c('Replicate','Plate','Index'))




model<-lm(Threshold~Mu+SD,data.rep56)
result<-summary(model)
result
result$r.squared



ggplot(testing,aes(x=Mu,y=Threshold,color=SD,shape=as.factor(Replicate)))+
  geom_point()+
  geom_abline(aes(intercept=0.020434,slope=0.995318  ))



dev.copy2pdf(device=cairo_pdf,
             file=paste(odir,"/Hue_threshold_model_with_SD.pdf",sep=''),
             width=9,height=6)







dim(comb)


peakfind<-function (x, thresh = 0) {
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
  if (!missing(thresh)) {
    pks[x[pks - 1] - x[pks] > thresh]
  }
  else pks
}


cleanh<-subset()

comb.hv<-dcast(comb,Replicate+Plate+Index+Brightness~Layer,value.var = 'Count')
dim(comb.hv)



max(comb.hv$h)


plate<-'PM1'
rep<-1

datasel<-subset(comb,Plate==plate & Replicate==rep)

ggplot(datasel,aes(x=Brightness,y=Count,color=Layer))+
  geom_line()+
  facet_wrap(~Index)


datasel.hv<-subset(comb.hv,Plate==plate & Replicate==rep)
ggplot(datasel.hv,aes(x=Brightness,y=h))+
  geom_line()+
  facet_wrap(~Index)










#Find peaks and their widths
source("https://bioconductor.org/biocLite.R")
biocLite('alsace')
data(tea)
new.lambdas <- seq(260, 500, by = 2)
tea <- lapply(tea.raw, preprocess, dim2 = new.lambdas)
tea.split <- splitTimeWindow(tea, c(12, 14), overlap = 10)

Xl <- tea.split[[2]]
Xl.opa <- opa(Xl, 4)

Xl.als <- doALS(Xl, Xl.opa)

tpoints <- getTime(Xl.als)
plot(tpoints, Xl.als$CList[[2]][,2], type = "l", col = "gray")
pk.pos <- findpeaks(Xl.als$CList[[2]][,2], span = 11)
abline(v = tpoints[pk.pos], col = 4)

pks <- fitpeaks(Xl.als$CList[[2]][,2], pk.pos)
apply(pks, 1,
      function(pkmodel) {
        lines(tpoints,
              dnorm(1:length(tpoints), pkmodel["rt"], pkmodel["sd"]) *
                pkmodel["area"],
              col = 2)
        invisible()
      })
## reasonably close fit, apart from the small peak in the middle...


