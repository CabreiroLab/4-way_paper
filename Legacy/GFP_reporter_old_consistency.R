library('ggplot2')
library('gplots')
library('plyr')
#Vennerable installation: install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library('Vennerable')
library('GOplot')

difff<-function(x)abs(diff(x))

folders<-c('13/13_Control','13/13_Exp','17/17_Exp','1/1_Control','1/1_Exp')
it<-1
for (f in folders) {
  print(unlist(strsplit(f,'/'))[1])
  f1<-read.csv(paste(f,'/Control_Fluorescence_norm_Shift.csv',sep=''),sep='\t',quote = '"',header = TRUE)
  f2<-read.csv(paste(f,'/Experiment_Fluorescence_norm_Shift.csv',sep=''),sep='\t',quote = '"',header = TRUE)
  f1$Replicate<-1
  f2$Replicate<-2
  f12<-merge(f1,f2,all.x=TRUE,all.y=TRUE)
  f12$Plate<-unlist(strsplit(f,'/'))[1]
  if (length(grep("_Exp",f))>0){
    f12$Type<-"With Metformin"
  } else {
    f12$Type<-"Control"
  }
  if (it==1){
    duplicates_temp<-f12
  } else {
    duplicates_temp<-merge(f12,duplicates_temp,all.x=TRUE,all.y=TRUE)
  }
  
  it<-it+1
}

sdcalc<-function(frame,x,y,st=1,sto1=2,sto2=3.5,step=0.1,window=0.25,fill=13,smooth=0.5){
  #step=0.1
  #window=0.25
  Mean=numeric()
  SD=numeric()
  Mean<-c(st)
  SD<-c(sd(subset(frame,x>st & x<st+window)[,y]))
  for (i in seq(st,sto1-window/2,step)){
    Mean<-c(Mean,mean(c(i,i+window)))
    SD<-c(SD,sd(subset(frame,x>i & x<i+window)[,y]))
  }
  if (sto2>sto1){
    SDconst<-sd(subset(frame,x>sto1 & x<sto2)[,y])
    Mean<-c(Mean,seq(sto1,sto2,length.out=fill))
    SD<-c(SD,rep(SDconst,fill))
  }
  print(Mean)
  print(SD)
  SD<-loess(SD ~ Mean,span=smooth)
  sigma<-SD$fitted
  stat<-data.frame(Mean,sigma)
  stat$sigma2<-stat$sigma*2
  stat$sigma2.5<-stat$sigma*2.5
  stat$sigma3<-stat$sigma*3
  return(stat)
}

Ztest<-function(Means,Diffs,SD){
  SDs<-approx(SD$Mean,SD$sigma,xout=Means,rule=2:2)
  Z<-abs(Diffs)/SDs$y
  return(list(Z=Z,SD=SDs$y))
}


datacalc<-function(frame,valuesrange){
  frame$Max<-apply(frame[,valuesrange],1,max)
  frame$AUC<-apply(frame[,valuesrange],1,sum)
  frame$AUC<-(frame$AUC)/12
  return(frame)
}

dupprep<-function(frame){
  sided<-merge(subset(frame,Replicate==1),subset(frame,Replicate==2),by=c('Type','Plate','Gene'),all.x=TRUE,all.y=TRUE)
  frame_temp<-ddply(frame, .(Type,Gene,Replicate,Plate), summarise, Max=mean(Max),AUC=mean(AUC), AUCnorm=mean(AUCnorm))
  frame_stat<-ddply(frame_temp, .(Type,Gene,Plate), summarise, Max_max=max(Max),Mean_max=mean(Max),Dif_max=diff(Max),Max_AUC=max(AUC), Mean_AUC=mean(AUC),Dif_AUC=diff(AUC),Max_AUCnorm=max(AUCnorm),Mean_AUCnorm=mean(AUCnorm),Dif_AUCnorm=diff(AUCnorm))
  return(list(main=frame,sided=sided,stat=frame_stat))
}

dataprep<-function(frame){
  frame_temp<-ddply(frame, .(Type,Gene), summarise, Max=mean(Max),AUC=mean(AUC),AUCnorm=mean(AUCnorm))
  sided<-merge(frame_temp[frame_temp$Type=='Control',],frame_temp[frame_temp$Type=='With Metformin',],by='Gene',suffixes = c("_C","_M"))
  sided$Max_max<-apply(sided[,c('Max_M','Max_C')],1,max)
  sided$Mean_max<-apply(sided[,c('Max_M','Max_C')],1,mean)
  sided$Diff_max<-apply(sided[,c('Max_C','Max_M')],1,diff)
  sided$Max_AUC<-apply(sided[,c('AUC_M','AUC_C')],1,max)
  sided$Mean_AUC<-apply(sided[,c('AUC_M','AUC_C')],1,mean)
  sided$Diff_AUC<-apply(sided[,c('AUC_C','AUC_M')],1,diff)
  sided$Max_AUCnorm<-apply(sided[,c('AUCnorm_M','AUCnorm_C')],1,max)
  sided$Mean_AUCnorm<-apply(sided[,c('AUCnorm_M','AUCnorm_C')],1,mean)
  sided$Diff_AUCnorm<-apply(sided[,c('AUCnorm_C','AUCnorm_M')],1,diff)
  return(sided)
}

duplicates<-duplicates_temp
duplicates<-subset(duplicates,Plate!=1)
duplicates[, 2:242 ]<-log10(sapply(duplicates[1:nrow(duplicates), 2:242 ], as.numeric))
duplicates[duplicates=='NaN']<--Inf
duplicates0<-duplicates
duplicates0[duplicates0<0]<-0
duplicates[duplicates<1]<-1

duplicates<-datacalc(duplicates,2:242)
duplicates0<-datacalc(duplicates0,2:242)


duplicates$AUCnorm<-(duplicates$AUC-241/12)/duplicates$Max
duplicates0$AUCnorm<-(duplicates0$AUC)/duplicates$Max


dup0<-dupprep(duplicates0)
dup<-dupprep(duplicates)


# 
# 
# repl_comp<-ggplot(dup$sided,aes(x=Max.x,y=Max.y,color=Type))+geom_point()+geom_rug(size=0.1)+geom_smooth(aes(group=1),method = "lm")+geom_abline(intercept=0,slope=1,color='red',alpha=0.2)
# repl_comp<-repl_comp+xlab('Maximum promoter activation in replicate 1, a.u.')+labs(color='Type')+ggtitle('Comparison of duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Maximum promoter activation in replicate 2, a.u.')
# repl_comp
# dev.copy2pdf(device=cairo_pdf,file="Stats/replcomp.pdf",width=11.69,height=8.27)
# 
# repl_comp0<-ggplot(dup0$sided,aes(x=Max.x,y=Max.y,color=Type))+geom_point()+geom_rug(size=0.1)+geom_smooth(aes(group=1),method = "lm")+geom_abline(intercept=0,slope=1,color='red',alpha=0.2)
# repl_comp0<-repl_comp0+xlab('Maximum promoter activation in replicate 1, a.u.')+labs(color='Type')+ggtitle('Comparison of duplicate experiments for plates 13, 17 (Cutoff=0)')+ylab('Maximum promoter activation in replicate 2, a.u.')
# repl_comp0
# dev.copy2pdf(device=cairo_pdf,file="Stats/replcomp0.pdf",width=11.69,height=8.27)
# 
# 
# 
# 
# 
# #ggplot(subset(dup_stat,Mean_max>1.5 & Mean_max<1.75),aes(x=Dif_max,fill=Type))+geom_histogram(binwidth=0.05,position='identity',alpha=0.5)#geom_density(alpha=0.5)
# repl<-ggplot(dup$stat,aes(x=Mean_max,y=Dif_max,color=Type))+geom_point()+geom_abline(intercept=-2,slope=2,color='red',alpha=0.2)+geom_abline(intercept=2,slope=-2,color='red',alpha=0.2)+ylim(-0.75,0.75)
# repl<-repl+xlab('Average promoter activation, a.u.')+labs(color='Type')+ggtitle('Differences between duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Absolute difference between replicates, a.u.')
# repll<-repl+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# repl2D<-repl+geom_density2d()
# repl
# dev.copy2pdf(device=cairo_pdf,file="Stats/repl.pdf",width=11.69,height=8.27)
# repl2D
# dev.copy2pdf(device=cairo_pdf,file="Stats/repl2D.pdf",width=11.69,height=8.27)
# repll
# dev.copy2pdf(device=cairo_pdf,file="Stats/repll.pdf",width=11.69,height=8.27)
# 
# 
# 
# repl0<-ggplot(dup0$stat,aes(x=Mean_max,y=Dif_max,color=Type))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)+ylim(-1.5,1.5)
# repl0<-repl0+xlab('Average promoter activation, a.u.')+labs(color='Type')+ggtitle('Differences between duplicate experiments for plates 13, 17 (Cutoff=0)')+ylab('Absolute difference between replicates, a.u.')
# repll0<-repl0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# repl02D<-repl0+geom_density2d()
# repl02D
# dev.copy2pdf(device=cairo_pdf,file="Stats/repl02D.pdf",width=11.69,height=8.27)
# repll0
# dev.copy2pdf(device=cairo_pdf,file="Stats/repll0.pdf",width=11.69,height=8.27)
# repl0
# dev.copy2pdf(device=cairo_pdf,file="Stats/repl0.pdf",width=11.69,height=8.27)
# 
# 
# 
# Diffsig<-sdcalc(dup$stat,dup$stat$Mean_max,'Dif_max',st=1,sto1=1.75,sto2=3.5,step=0.1,window=0.25)
# Diffsigy=Diffsig$sigma2
# Diffl<-'2 sigma\nthreshold'
# replcut<-repl+geom_line(data=Diffsig,aes(x=Mean,y=Diffsigy,color=Diffl))+geom_line(data=Diffsig,aes(x=Mean,y=Diffsigy*-1,color=Diffl))+xlim(1,3)
# replcut
# dev.copy2pdf(device=cairo_pdf,file="Stats/replcut.pdf",width=11.69,height=8.27)
# 
# Diffsig0<-sdcalc(dup0$stat,dup0$stat$Mean_max,'Dif_max',st=0.5,sto1=1.75,sto2=3.5,step=0.1,window=0.25)
# Diffsig0y=Diffsig0$sigma2
# repl0cut<-repl0+geom_line(data=Diffsig0,aes(x=Mean,y=Diffsig0y,color=Diffl))+geom_line(data=Diffsig0,aes(x=Mean,y=Diffsig0y*-1,color=Diffl))+xlim(0,3)
# repl0cut
# dev.copy2pdf(device=cairo_pdf,file="Stats/repl0cut.pdf",width=11.69,height=8.27)
# 
# 
# replp<-ggplot(dup$stat,aes(x=Mean_max,y=Dif_max,color=Plate))+geom_point()+geom_abline(intercept=-2,slope=2,color='red',alpha=0.2)+geom_abline(intercept=2,slope=-2,color='red',alpha=0.2)+ylim(-1,1)
# replp<-replp+xlab('Average promoter activation, a.u.')+labs(color='Plate')+ggtitle('Differences between duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Absolute difference between replicates, a.u.')
# replp2D<-replp+geom_density2d()
# replp
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp.pdf",width=11.69,height=8.27)
# replp2D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp2D.pdf",width=11.69,height=8.27)
# 
# 
# replp0<-ggplot(dup0$stat,aes(x=Mean_max,y=Dif_max,color=Plate))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)+ylim(-1,1)
# replp0<-replp0+xlab('Average promoter activation, a.u.')+labs(color='Plate')+ggtitle('Differences between duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Absolute difference between replicates, a.u.')
# replp02D<-replp0+geom_density2d()
# replp0
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp0.pdf",width=11.69,height=8.27)
# replp02D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp02D.pdf",width=11.69,height=8.27)
# 
# 
# replp0<-ggplot(dup0$stat,aes(x=Mean_max,y=Dif_max,color=Plate))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)+ylim(-1,1)
# replp0<-replp0+xlab('Average promoter activation, a.u.')+labs(color='Plate')+ggtitle('Differences between duplicate experiments for plates 13, 17 (Cutoff=0)')+ylab('Absolute difference between replicates, a.u.')
# replp02D<-replp0+geom_density2d()
# replp0
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp0.pdf",width=11.69,height=8.27)
# replp02D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replp02D.pdf",width=11.69,height=8.27)
# 
# 
# replAUC<-ggplot(dup$stat,aes(x=Mean_AUC,y=Dif_AUC,color=Type))+geom_point()+geom_abline(intercept=-40,slope=2,color='red',alpha=0.2)+geom_abline(intercept=40,slope=-2,color='red',alpha=0.2)
# replAUC<-replAUC+xlab('Average promoter activity (AUC), a.u./h')+labs(color='Type')+ggtitle('Differences in promoter activity (AUC) between duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Absolute difference in activity (AUC), a.u./h')
# replAUC2D<-replAUC+geom_density2d()
# repllAUC<-replAUC+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# replAUC
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUC.pdf",width=11.69,height=8.27)
# replAUC2D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUC2D.pdf",width=11.69,height=8.27)
# repllAUC
# dev.copy2pdf(device=cairo_pdf,file="Stats/repllAUC.pdf",width=11.69,height=8.27)
# 
# AUCsig<-sdcalc(dup$stat,dup$stat$Mean_AUC,'Dif_AUC',st=20,sto1=30,sto2=45,step=1,window=5,smooth=0.5,fill=23)
# AUCsigy=AUCsig$sigma2
# replAUCcut<-replAUC+geom_line(data=AUCsig,aes(x=Mean,y=AUCsigy,color='2 sigma\nthreshold'))+geom_line(data=AUCsig,aes(x=Mean,y=AUCsigy*-1,color='2 sigma\nthreshold'))+xlim(20,50)
# replAUCcut
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCcut.pdf",width=11.69,height=8.27)
# 
# replAUC0<-ggplot(dup0$stat,aes(x=Mean_AUC,y=Dif_AUC,color=Type))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
# replAUC0<-replAUC0+xlab('Average promoter activity (AUC), a.u./h')+labs(color='Type')+ggtitle('Differences in promoter activity (AUC) between duplicate experiments for plates 13, 17 (Cutoff=0)')+ylab('Absolute difference in activity (AUC), a.u./h')
# replAUC02D<-replAUC0+geom_density2d()
# repllAUC0<-replAUC0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# replAUC0
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUC0.pdf",width=11.69,height=8.27)
# replAUC02D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUC02D.pdf",width=11.69,height=8.27)
# repllAUC0
# dev.copy2pdf(device=cairo_pdf,file="Stats/repllAUC.pdf",width=11.69,height=8.27)
# 
# AUCsig0<-sdcalc(dup0$stat,dup0$stat$Mean_AUC,'Dif_AUC',st=0,sto1=30,sto2=45,step=1,window=5,smooth=0.5,fill=23)
# AUCsig0y=AUCsig0$sigma2
# replAUC0cut<-replAUC0+geom_line(data=AUCsig0,aes(x=Mean,y=AUCsig0y,color='2 sigma\nthreshold'))+geom_line(data=AUCsig0,aes(x=Mean,y=AUCsig0y*-1,color='2 sigma\nthreshold'))
# replAUC0cut
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUC0cut.pdf",width=11.69,height=8.27)
# 
# 
# 
# 
# 
# 
# replAUCnorm<-ggplot(dup$stat,aes(x=Mean_AUCnorm,y=Dif_AUCnorm,color=Type))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
# replAUCnorm<-replAUCnorm+xlab('Normalised average promoter activity (AUC), a.u./h')+labs(color='Type')+ggtitle('Differences in normalised promoter activity (AUC) between duplicate experiments for plates 13, 17 (Cutoff=1)')+ylab('Normalised difference in activity (AUC), a.u./h')
# replAUCnorm2D<-replAUCnorm+geom_density2d()
# repllAUCnorm<-replAUCnorm+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# replAUCnorm
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnorm.pdf",width=11.69,height=8.27)
# replAUCnorm2D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnorm2D.pdf",width=11.69,height=8.27)
# repllAUCnorm
# dev.copy2pdf(device=cairo_pdf,file="Stats/repllAUCnorm.pdf",width=11.69,height=8.27)
# 
# 
# AUCnormsig<-sdcalc(dup$stat,dup$stat$Mean_AUCnorm,'Dif_AUCnorm',st=0,sto1=7.5,sto2=15,step=0.5,window=1,smooth=0.5)
# AUCnormsigy=AUCnormsig$sigma2
# replAUCnormcut<-replAUCnorm+geom_line(data=AUCnormsig,aes(x=Mean,y=AUCnormsigy,color='2 sigma\nthreshold'))+geom_line(data=AUCnormsig,aes(x=Mean,y=AUCnormsigy*-1,color='2 sigma\nthreshold'))+xlim(0,10)
# replAUCnormcut
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnormcut.pdf",width=11.69,height=8.27)
# 
# 
# replAUCnorm0<-ggplot(dup0$stat,aes(x=Mean_AUCnorm,y=Dif_AUCnorm,color=Type))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
# replAUCnorm0<-replAUCnorm0+xlab('Normalised average promoter activity (AUC), a.u./h')+labs(color='Type')+ggtitle('Differences in normalised promoter activity (AUC) between duplicate experiments for plates 13, 17 (Cutoff=0)')+ylab('Normalised difference in activity (AUC), a.u./h')
# replAUCnorm02D<-replAUCnorm0+geom_density2d()
# repllAUCnorm0<-replAUCnorm0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=5)
# replAUCnorm0
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnorm0.pdf",width=11.69,height=8.27)
# replAUCnorm02D
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnorm02D.pdf",width=11.69,height=8.27)
# repllAUCnorm0
# dev.copy2pdf(device=cairo_pdf,file="Stats/repllAUCnorm0.pdf",width=11.69,height=8.27)
# 
# 
# AUCnormsig0<-sdcalc(dup0$stat,dup0$stat$Mean_AUCnorm,'Dif_AUCnorm',st=0,sto1=15,sto2=40,step=0.5,window=2.5,smooth=0.5)
# AUCnormsig0y=AUCnormsig0$sigma2
# replAUCnorm0cut<-replAUCnorm0+geom_line(data=AUCnormsig0,aes(x=Mean,y=AUCnormsig0y,color='2 sigma\nthreshold'))+geom_line(data=AUCnormsig0,aes(x=Mean,y=AUCnormsig0y*-1,color='2 sigma\nthreshold'))+xlim(0,17.5)
# replAUCnorm0cut
# dev.copy2pdf(device=cairo_pdf,file="Stats/replAUCnorm0cut.pdf",width=11.69,height=8.27)
# 
# 


Cfluor<-read.csv('../Allnew/Control_Fluorescence_norm_Shift.csv',sep='\t',quote = '"',header = TRUE)
Efluor<-read.csv('../Allnew/Experiment_Fluorescence_norm_Shift.csv',sep='\t',quote = '"',header = TRUE)
Cfluor$Type<-'Control'
Efluor$Type<-'With Metformin'

if ("Index" %in% colnames(Cfluor)) {
	valuesrange<-3:243
} else {
	valuesrange<-2:242
}

Fluor_temp<-merge(Cfluor,Efluor,all.x=TRUE,all.y=TRUE)

Allfluor<-Fluor_temp
#Allfluor[, valuesrange ]<-log10(sapply(Allfluor[1:nrow(Allfluor), valuesrange ], as.numeric))
Allfluor[Allfluor=='NaN']<--Inf
Allfluor0<-Allfluor
Allfluor0[Allfluor0<0]<-0
Allfluor[Allfluor<1]<-1


Allfluor<-datacalc(Allfluor,valuesrange)
Allfluor0<-datacalc(Allfluor0,valuesrange)

Allfluor$AUCnorm<-(Allfluor$AUC-241/12)/Allfluor$Max
Allfluor0$AUCnorm<-Allfluor0$AUC/Allfluor$Max

Allfluor0$Cutoff<-0
Allfluor$Cutoff<-1

Reffluor<-merge(Allfluor,Allfluor0,all.x=TRUE,all.y=TRUE)
Reffluor<-Reffluor[Reffluor$Gene %in% c('lacZ','wrbA','serA'),]

Refstat<-ddply(Reffluor,.(Type,Gene,Cutoff),summarize, Mean_max=mean(Max), SD_max=sd(Max), Mean_AUC=mean(AUC), SD_AUC=sd(AUC), Mean_AUCnorm=mean(AUCnorm), SD_AUCnorm=sd(AUCnorm))



RefPeak<-ggplot(data=subset(Reffluor,Cutoff=0), aes(x=factor(Gene, levels=c('lacZ','wrbA','serA')),y=Max,fill=Type))+
  geom_boxplot(notch=TRUE)+
  #ylim(0,10)+
  xlab('Promoter')+
  labs(fill='Type',alpha='Cutoff')+
  ylab('Promoter activation level (PNF), a.u.')+
  ggtitle('c',subtitle = 'Before log2 transformation')+
  theme(legend.position="none")
RefPeak
dev.copy2pdf(device=cairo_pdf,file="RefPeak.pdf",width=4,height=4)

# RefAUC<-ggplot(data=Reffluor, aes(x=factor(Gene, levels=c('lacZ','wrbA','serA')),y=AUC,fill=Type,alpha=factor(Cutoff)))+
#   geom_boxplot(notch=TRUE)+
#   xlab('Promoter')+
#   labs(fill='Type',alpha='Cutoff')+
#   ylab('Promoter activity (AUC), a.u.')+
#   ggtitle('Reference promoter activity (AUC) variation through the plates with different cutoffs')
# RefAUC
# dev.copy2pdf(device=cairo_pdf,file="RefAUC_log.pdf",width=4,height=4)

RefAUCnorm<-ggplot(data=subset(Reffluor,Cutoff==4), aes(x=factor(Gene, levels=c('lacZ','wrbA','serA')),y=AUCnorm,fill=Type))+
  geom_boxplot(notch=TRUE)+
  xlab('Promoter')+
  labs(fill='Type',alpha='Cutoff')+
  ylab('Promoter activation timing (AUCnorm), a.u.')+ggtitle('d')
RefAUCnorm
dev.copy2pdf(device=cairo_pdf,file="RefAUCnorm.pdf",width=6,height=4)

RefPeakSD<-ggplot(data=subset(Refstat,Cutoff==4),aes(y=SD_max,x=Mean_max,shape=Gene,color=Type))+
  geom_point(size=5)+
  xlab('Average promoter activation, a.u.')+
  ylab('SD')+
  labs(shape='Promoter',color='Type')+
  ggtitle('d')
RefPeakSD

RefPeakSD0<-ggplot(data=subset(Refstat,Cutoff==0),aes(y=SD_max,x=Mean_max,shape=Gene,color=Type))+
  geom_smooth(aes(group=1),method = "lm")+
  geom_point(size=5)+#+xlim(0,50)+ylim(0,8.5)+
  xlab('Average promoter activation, a.u.')+
  ylab('SD')+
  labs(shape='Promoter',color='Type')+
  ggtitle('d',subtitle='Before log2 transformation')
RefPeakSD0
dev.copy2pdf(device=cairo_pdf,file="Lincor.pdf",width=6,height=4)


All_sided<-dataprep(Allfluor)
All_sided0<-dataprep(Allfluor0)
All_sided$Cutoff<-1
All_sided0$Cutoff<-0


All_sided$Z_max<-Ztest(All_sided$Mean_max,All_sided$Diff_max,Diffsig)$Z
All_sided$Z_AUC<-Ztest(All_sided$Mean_AUC,All_sided$Diff_AUC,AUCsig)$Z
All_sided$Z_AUCnorm<-Ztest(All_sided$Mean_AUCnorm,All_sided$Diff_AUCnorm,AUCnormsig)$Z

All_sided0$Z_max<-Ztest(All_sided0$Mean_max,All_sided0$Diff_max,Diffsig0)$Z
All_sided0$Z_AUC<-Ztest(All_sided0$Mean_AUC,All_sided0$Diff_AUC,AUCsig0)$Z
All_sided0$Z_AUCnorm<-Ztest(All_sided0$Mean_AUCnorm,All_sided0$Diff_AUCnorm,AUCnormsig0)$Z


All<-merge(All_sided,All_sided0,by='Gene',suffixes = c("1","0"))
All_listed<-merge(All_sided,All_sided0,all.x=TRUE,all.y=TRUE)


#All$Z_max1<-Ztest(All$Mean_max_1,All$Diff_max_1,Diffsig)$Z
#All$Z_max0<-Ztest(All$Mean_max_0,All$Diff_max_0,Diffsig)$Z
#All$Z_AUC1<-Ztest(All$Mean_AUC_1,All$Diff_AUC_1,AUCsig)$Z
#All$Z_AUC0<-Ztest(All$Mean_AUC_0,All$Diff_AUC_0,AUCsig0)$Z
#All$Z_AUCnorm1<-Ztest(All$Mean_AUCnorm_1,All$Diff_AUCnorm_1,AUCnormsig)$Z
#All$Z_AUCnorm0<-Ztest(All$Mean_AUCnorm_0,All$Diff_AUCnorm_0,AUCnormsig0)$Z

#g<-approx(Diff_cut$Mean,Diff_cut$sigma2,xout=c(3.5),rule=2:2)

#Avoid using automatic summary
#All_stat<-ddply(All_temp, .(Gene), summarise, Max_max=max(Max), Diff_max=diff(Max))
All_diff0c<-ggplot(All,aes(x=Max_C0,y=Max_M0))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_diff0c<-All_diff0c+xlab('Maximum promoter activation in Control, a.u.')+ggtitle('Comparison of promoter activation with and without metformin (Cutoff=0)')+ylab('Maximum promoter activation with Metformin, a.u.')
All_diff0c
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diff0c.pdf",width=11.69,height=8.27)

All_AUC0c<-ggplot(All,aes(x=AUC_C0,y=AUC_M0))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_AUC0c<-All_AUC0c+xlab('Promoter activity (AUC) in Control, a.u.')+ggtitle('Comparison of promoter activity (AUC) with and without metformin (Cutoff=0)')+ylab('Promoter activity (AUC) with Metformin, a.u.')
All_AUC0c
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUC0c.pdf",width=11.69,height=8.27)

All_AUCnorm0c<-ggplot(All,aes(x=AUCnorm_C0,y=AUCnorm_M0))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_AUCnorm0c<-All_AUCnorm0c+xlab('Normalised promoter activity (AUC) in Control, a.u.')+ggtitle('Comparison of normalised promoter activity (AUC) with and without metformin (Cutoff=0)')+ylab('Normalised promoter activity (AUC) with Metformin, a.u.')
All_AUCnorm0c
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnorm0c.pdf",width=11.69,height=8.27)

All_diffc<-ggplot(All,aes(x=Max_C1,y=Max_M1))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_diffc<-All_diffc+xlab('Maximum promoter activation in Control, a.u.')+ggtitle('Comparison of promoter activation with and without metformin (Cutoff=1)')+ylab('Maximum promoter activation with Metformin, a.u.')
All_diffc
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diffc.pdf",width=11.69,height=8.2)

All_AUCc<-ggplot(All,aes(x=AUC_C1,y=AUC_M1))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_AUCc<-All_AUCc+xlab('Promoter activity (AUC) in Control, a.u.')+ggtitle('Comparison of promoter activity (AUC) with and without metformin (Cutoff=1)')+ylab('Promoter activity (AUC) with Metformin, a.u.')
All_AUCc
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCc.pdf",width=11.69,height=8.2)

All_AUCnormc<-ggplot(All,aes(x=AUCnorm_C1,y=AUCnorm_M1))+geom_point()+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)+geom_smooth(color='navy')+geom_abline(intercept=0,slope=1,color='red',alpha=0.5)
All_AUCnormc<-All_AUCnormc+xlab('Normalised promoter activity (AUC) in Control, a.u.')+ggtitle('Comparison of normalised promoter activity (AUC) with and without metformin (Cutoff=1)')+ylab('Normalised promoter activity (AUC) with Metformin, a.u.')
All_AUCnormc
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnormc.pdf",width=11.69,height=8.2)

Diffsigy=Diffsig$sigma2
All_diff<-ggplot(All,aes(x=Mean_max1,y=Diff_max1))+geom_point()+geom_abline(intercept=-2,slope=2,color='red',alpha=0.2)+geom_abline(intercept=2,slope=-2,color='red',alpha=0.2)
All_diff<-All_diff+geom_line(data=Diffsig,aes(x=Mean,y=Diffsigy,color='2 sigma\nthreshold'),size=2)+geom_line(data=Diffsig,aes(x=Mean,y=Diffsigy*-1,color='2 sigma\nthreshold'),size=2)+labs(color='Significance')+xlim(1,3.5)
All_diff<-All_diff+geom_smooth(color='navy',size=1.5)+xlab('Maximum promoter activation, a.u.')+ggtitle('Differences in promoter activation with and without metformin (Cutoff=1)')+ylab('Absolute difference, a.u.')
All_diff<-All_diff+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_diff2D<-All_diff+geom_density2d()
All_diff
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diff.pdf",width=11.69,height=8.2)
All_diff2D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diff2D.pdf",width=11.69,height=8.2)



All_diff0<-ggplot(All,aes(x=Mean_max0,y=Diff_max0))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
All_diff0<-All_diff0+geom_line(data=Diffsig0,aes(x=Mean,y=Diffsig0y,color='2 sigma\nthreshold'),size=2)+geom_line(data=Diffsig0,aes(x=Mean,y=Diffsig0y*-1,color='2 sigma\nthreshold'),size=2)+labs(color='Significance')
All_diff0<-All_diff0+geom_smooth(color='navy',size=1.5)+xlab('Maximum promoter activation, a.u.')+ggtitle('Differences in promoter activation with and without metformin (Cutoff=0)')+ylab('Absolute difference, a.u.')
All_diff0<-All_diff0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_diff02D<-All_diff0+geom_density2d()
All_diff0
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diff0.pdf",width=11.69,height=8.2)
All_diff02D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_diff02D.pdf",width=11.69,height=8.2)



All_AUC<-ggplot(All,aes(x=Mean_AUC1,y=Diff_AUC1))+geom_point()+geom_abline(intercept=-40,slope=2,color='red',alpha=0.2)+geom_abline(intercept=40,slope=-2,color='red',alpha=0.2)
All_AUC<-All_AUC+geom_smooth(color='navy')+xlab('Average promoter activity (AUC), a.u./h')+ggtitle('Differences in promoter activity (AUC) with and without metformin (Cutoff=1)')+ylab('Absolute difference (AUC), a.u./h')
All_AUC<-All_AUC+geom_line(data=AUCsig,aes(x=Mean,y=AUCsigy,color='2 sigma\nthreshold'),size=2,alpha=0.5)+geom_line(data=AUCsig,aes(x=Mean,y=AUCsigy*-1,color='2 sigma\nthreshold'),size=2,alpha=0.5)+labs(color='Significance')+xlim(20,42.5)
All_AUC<-All_AUC+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_AUC2D<-All_AUC+geom_density2d()
All_AUC
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUC.pdf",width=11.69,height=8.2)
All_AUC2D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUC2D.pdf",width=11.69,height=8.2)

All_AUC0<-ggplot(All,aes(x=Mean_AUC0,y=Diff_AUC0))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
All_AUC0<-All_AUC0+geom_smooth(color='navy')+xlab('Average promoter activity (AUC), a.u./h')+ggtitle('Differences in promoter activity (AUC) with and without metformin (Cutoff=0)')+ylab('Absolute difference (AUC), a.u./h')
All_AUC0<-All_AUC0+geom_line(data=AUCsig0,aes(x=Mean,y=AUCsig0y,color='2 sigma\nthreshold'),size=2,alpha=0.5)+geom_line(data=AUCsig0,aes(x=Mean,y=AUCsig0y*-1,color='2 sigma\nthreshold'),size=2,alpha=0.5)+labs(color='Significance')+xlim(0,42.5)
All_AUC0<-All_AUC0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_AUC02D<-All_AUC0+geom_density2d()
All_AUC0
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUC0.pdf",width=11.69,height=8.2)
All_AUC02D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUC02D.pdf",width=11.69,height=8.2)

All_AUCnorm<-ggplot(All,aes(x=Mean_AUCnorm1,y=Diff_AUCnorm1))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
All_AUCnorm<-All_AUCnorm+geom_smooth(color='navy')+xlab('Average normalised promoter activity (AUC), a.u./h')+ggtitle('Differences in normalised promoter activity (AUC) with and without metformin (Cutoff=1)')+ylab('Normalised difference (AUC), a.u./h')
All_AUCnorm<-All_AUCnorm+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_AUCnorm<-All_AUCnorm+geom_line(data=AUCnormsig,aes(x=Mean,y=AUCnormsigy,color='2 sigma\nthreshold'),size=2,alpha=0.5)+geom_line(data=AUCnormsig,aes(x=Mean,y=AUCnormsigy*-1,color='2 sigma\nthreshold'),size=2,alpha=0.5)+labs(color='Significance')+xlim(0,12.5)
All_AUCnorm2D<-All_AUCnorm+geom_density2d()
All_AUCnorm
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnorm.pdf",width=11.69,height=8.2)
All_AUCnorm2D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnorm2D.pdf",width=11.69,height=8.2)

All_AUCnorm0<-ggplot(All,aes(x=Mean_AUCnorm0,y=Diff_AUCnorm0))+geom_point()+geom_abline(intercept=0,slope=2,color='red',alpha=0.2)+geom_abline(intercept=0,slope=-2,color='red',alpha=0.2)
All_AUCnorm0<-All_AUCnorm0+geom_smooth(color='navy')+xlab('Average normalised promoter activity (AUC), a.u./h')+ggtitle('Differences in normalised promoter activity (AUC) with and without metformin (Cutoff=0)')+ylab('Normalised difference (AUC), a.u./h')
All_AUCnorm0<-All_AUCnorm0+geom_text(aes(label=Gene),hjust=0, vjust=0,size=3)
All_AUCnorm0<-All_AUCnorm0+geom_line(data=AUCnormsig0,aes(x=Mean,y=AUCnormsig0y,color='2 sigma\nthreshold'),size=2,alpha=0.5)+geom_line(data=AUCnormsig0,aes(x=Mean,y=AUCnormsig0y*-1,color='2 sigma\nthreshold'),size=2,alpha=0.5)+labs(color='Significance')+xlim(0,20)
All_AUCnorm02D<-All_AUCnorm0+geom_density2d()
All_AUCnorm0
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnorm0.pdf",width=11.69,height=8.2)
All_AUCnorm02D
dev.copy2pdf(device=cairo_pdf,file="Stats/All_AUCnorm02D.pdf",width=11.69,height=8.2)

Z_Peak1<-ggplot(All,aes(x=Diff_max1,y=Z_max1,color=ifelse(Diff_max0>0,'Upregulated','Downregulated')))+geom_point()+geom_text(aes(label=ifelse(Z_max1>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)#+ylim(-2,4)
Z_Peak1<-Z_Peak1+ggtitle('Z-scored differences in cumulative promoter activation between control and treatment (Cutoff=1)')+ylab('Z-score')+xlab('Difference in peak promoter activation, a.u.')+labs(color='With metformin:')
Z_Peak1
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_Peak1.pdf",width=11.69,height=8.2)

Z_Peak0<-ggplot(All,aes(x=Diff_max0,y=Z_max0,color=ifelse(Diff_max0>0,'Upregulated','Downregulated')))+geom_point(alpha=0.9)+geom_text(aes(label=ifelse(Z_max0>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)
Z_Peak0<-Z_Peak0+ggtitle('Z-scored differences in promoter activation between control and treatment (Cutoff=0)')+ylab('Z-score')+xlab('Difference in peak promoter activation, a.u.')+labs(color='With metformin:')
Z_Peak0
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_Peak0.pdf",width=11.69,height=8.2)

Z_AUC1<-ggplot(All,aes(x=Diff_AUC1,y=Z_AUC1,color=ifelse(Diff_AUC1>0,'Upregulated','Downregulated')))+geom_point()+geom_text(aes(label=ifelse(Z_AUC1>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)
Z_AUC1<-Z_AUC1+ggtitle('Z-scored differences in promoter activation between control and treatment (Cutoff=1)')+ylab('Z-score')+xlab('Difference in promoter activity (AUC), a.u.')+labs(color='With metformin:')
Z_AUC1
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_AUC1.pdf",width=11.69,height=8.2)

Z_AUC0<-ggplot(All,aes(x=Diff_AUC0,y=Z_AUC0,color=ifelse(Diff_AUC0>0,'Upregulated','Downregulated')))+geom_point(alpha=0.9)+geom_text(aes(label=ifelse(Z_AUC0>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)
Z_AUC0<-Z_AUC0+ggtitle('Z-scored differences in promoter activation between control and treatment (Cutoff=0)')+ylab('Z-score')+xlab('Difference in peak promoter activity (AUC), a.u.')+labs(color='With metformin:')
Z_AUC0
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_AUC0.pdf",width=11.69,height=8.2)

Z_AUCnorm1<-ggplot(All,aes(x=Diff_AUCnorm1,y=Z_AUCnorm1,color=ifelse(Diff_AUCnorm1>0,'Upregulated','Downregulated')))+geom_point()+geom_text(aes(label=ifelse(Z_AUCnorm1>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)
Z_AUCnorm1<-Z_AUCnorm1+ggtitle('Z-scored differences in promoter activation between control and treatment (Cutoff=1)')+ylab('Z-score')+xlab('Difference in normalised promoter activity (AUC), a.u.')+labs(color='With metformin:')
Z_AUCnorm1
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_AUCnorm1.pdf",width=11.69,height=8.2)

Z_AUCnorm0<-ggplot(All,aes(x=Diff_AUCnorm0,y=Z_AUCnorm0,color=ifelse(Diff_AUCnorm0>0,'Upregulated','Downregulated')))+geom_point(alpha=0.9)+geom_text(aes(label=ifelse(Z_AUCnorm0>2,as.character(Gene),'')),hjust=0,just=0,size=3)+geom_hline(yintercept=2,color='red',alpha=0.5)
Z_AUCnorm0<-Z_AUCnorm0+ggtitle('Z-scored differences in promoter activation between control and treatment (Cutoff=0)')+ylab('Z-score')+xlab('Difference in normalised promoter activity (AUC), a.u.')+labs(color='With metformin:')
Z_AUCnorm0
dev.copy2pdf(device=cairo_pdf,file="Stats/Z_AUCnorm0.pdf",width=11.69,height=8.2)

write.csv(All,file='../Results/All_results.csv')
write.csv(All_listed,file='../Results/All_listed_results.csv')

#subset(All, Z_max1>2.5 & !Z_max0>2.5)[,'Gene']

Thres<-2

GPeak1<-subset(All,Z_max1>Thres)
GPeak0<-subset(All,Z_max0>Thres)
GAUC1<-subset(All,Z_AUC1>Thres)
GAUC0<-subset(All,Z_AUC0>Thres)
GAUCnorm1<-subset(All,Z_AUCnorm1>Thres)
GAUCnorm0<-subset(All,Z_AUCnorm0>Thres)

Lists1<-list('Peak'=GPeak1$Gene,'AUC'=GAUC1$Gene,'AUCnorm'=GAUCnorm1$Gene)
Lists0<-list('Peak'=GPeak0$Gene,'AUC'=GAUC0$Gene,'AUCnorm'=GAUCnorm0$Gene)

Peakc<-list('Peak0'=GPeak0$Gene,'Peak1'=GPeak1$Gene)
AUCc<-list('AUC0'=GAUC0$Gene,'AUC1'=GAUC1$Gene)
AUCnormc<-list('AUCnorm0'=GAUCnorm0$Gene,'AUCnorm1'=GAUCnorm1$Gene)

Peakcc<-list('Peak0_UP'=subset(GPeak0,Diff_max0>0)$Gene,'Peak0_DOWN'=subset(GPeak0,Diff_max0<0)$Gene,'Peak1_UP'=subset(GPeak1,Diff_max1>0)$Gene,'Peak1_DOWN'=subset(GPeak1,Diff_max1<0)$Gene)
AUCcc<-list('AUC0_UP'=subset(GAUC0,Diff_AUC0>0)$Gene,'AUC0_DOWN'=subset(GAUC0,Diff_AUC0<0)$Gene,'AUC1_UP'=subset(GAUC1,Diff_AUC1>0)$Gene,'AUC1_DOWN'=subset(GAUC1,Diff_AUC1<0)$Gene)
AUCnormcc<-list('AUCnorm0_UP'=subset(GAUCnorm0,Diff_AUCnorm0>0)$Gene,'AUCnorm0_DOWN'=subset(GAUCnorm0,Diff_AUCnorm0<0)$Gene,'AUCnorm1_UP'=subset(GAUCnorm1,Diff_AUCnorm1>0)$Gene,'AUCnorm1_DOWN'=subset(GAUCnorm1,Diff_AUCnorm1<0)$Gene)

plot(Venn(Peakcc),type='ellipses', doWeights = FALSE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file="Stats/Peak_comparison.pdf",width=11.69,height=8.2)
plot(Venn(AUCcc),type='ellipses', doWeights = FALSE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file="Stats/AUC_comparison.pdf",width=11.69,height=8.2)
plot(Venn(AUCnormcc),type='ellipses', doWeights = FALSE)#,type='ellipses'
dev.copy2pdf(device=cairo_pdf,file="Stats/AUCnorm_comparison.pdf",width=11.69,height=8.2)


plot(Venn(Lists1), doWeights = FALSE)#,type='ellipses'


final<-list("Peak_UP"=unique(c(as.character(Peakcc$Peak0_UP), as.character(Peakcc$Peak1_UP))),
               "Peak_DOWN"=unique(c(as.character(Peakcc$Peak0_DOWN),as.character(Peakcc$Peak1_DOWN))),
               "AUC_UP"=unique(c(as.character(AUCcc$AUC0_UP),as.character(AUCcc$AUC1_UP))),
               "AUC_DOWN"=unique(c(as.character(AUCcc$AUC0_DOWN),as.character(AUCcc$AUC1_DOWN))),
                "AUCnorm_UP"=as.character(AUCnormcc$AUCnorm1_UP),"AUCnorm_DOWN"=as.character(AUCnormcc$AUCnorm1_DOWN)
               )

final<-c(final,list("UP"=unique(c(as.character(final$Peak_UP), as.character(final$AUC_UP), as.character(final$AUCnorm_UP))),
                    "DOWN"=unique(c(as.character(final$Peak_DOWN), as.character(final$AUC_DOWN), as.character(final$AUCnorm_DOWN))),
                    "UP_consensus"=intersect(final$Peak_UP,intersect(final$AUC_UP,final$AUCnorm_UP)),
                    "DOWN_consensus"=intersect(final$Peak_DOWN,intersect(final$AUC_DOWN,final$AUCnorm_DOWN)),
                    "AUCnorm_DOWN_exclusive"=setdiff(final$AUCnorm_DOWN,union(final$Peak_DOWN,final$AUC_DOWN)),
                    "AUC_DOWN_exclusive"=setdiff(final$AUC_DOWN,union(final$Peak_DOWN,final$AUCnorm_DOWN))
                    ))
        




plot(Venn(final[c('Peak_UP','AUC_UP','AUCnorm_UP')]), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file="Stats/Upregulated_Venn.pdf",width=11.69,height=8.2)
plot(Venn(final[c('Peak_DOWN','AUC_DOWN','AUCnorm_DOWN')]), doWeights = FALSE)
dev.copy2pdf(device=cairo_pdf,file="Stats/Downregulated_Venn.pdf",width=11.69,height=8.2)

lwrite<-function(x,flnm) {
  for (l in names(x)) {
      write(paste(l) ,file=flnm, append=TRUE, ncolumns=1000)
      write(unlist(x[l]),file=flnm, append=TRUE, ncolumns=1000)
  }
  
}

lwrite(final, "Stats/Gene_lists.txt")
