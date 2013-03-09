
library(ggplot2)
library(scales)
library(flowCore)
library(gridExtra)
library(reshape)
library(plyr)
library(grImport)
library(Hmisc)
library(corrplot)
source("/Users/jonathan/Documents/flowset2ggplot.r")



# generates Figure 1

#Figure 1a Clonal Workflow
color<-TRUE


if(color) {
	PostScriptTrace("Clonal_workflow_color.ps")
	workflow<-readPicture("Clonal_workflow_color.ps.xml")
	
	
}else {
	PostScriptTrace("Clonal_workflow_grayscale.ps")
	workflow<-readPicture("Clonal_workflow_grayscale.ps.xml")
	}	

F1a<-pictureGrob(workflow);

#Figure 1b Bulk Gating
LGbulk_GFP<-read.table("/Users/jonathan/Documents/JF_SD_paperFigs/Data/LGbulk.txt")
NJ_GFP<-read.table("/Users/jonathan/Documents/JF_SD_paperFigs/Data/NJ.txt")
LGBulk<-cbind(name=rep("LGM2 Polyclonal",length(LGbulk_GFP)),rename(LGbulk_GFP,c(V1="GFP")))
NJ<-cbind(name=rep("Autofluorescence",length(NJ_GFP)),rename(NJ_GFP,c(V1="GFP")))

#Shifts RFU top 10^0-10^4 rather than 10^-1-10^3, doesn't really matter
#Bulk_Auto<-as.data.frame(rbind(LGBulk,NJ))
#Bulk_Auto$GFP<-10^(log10(Bulk_Auto$GFP+1));

F1b<-ggplot(Bulk_Auto,aes(x=GFP,fill=name))+geom_density(alpha=0.4)
if(color) F1b<-F1b+scale_fill_manual(values=c("green","grey")) else F1b<-F1b+scale_fill_grey()
F1b<-F1b+scale_x_log10(name="GFP RFU",breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)))
F1b<-F1b+theme_bw(16)+theme(legend.position="bottom")+guides(fill=guide_legend(title=NULL))
F1b<-F1b+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())

gateMin<-2
gateMax<-95
gateVert<-0.9
F1b<-F1b+geom_segment(aes(x=gateMin,y=gateVert,xend=gateMax,yend=gateVert))+geom_segment(aes(x=gateMin-0.05,xend=gateMin-0.05,y=gateVert-0.1,yend=gateVert+0.1))+geom_segment(aes(x=gateMax+0.05,xend=gateMax+0.05,y=gateVert-0.1,yend=gateVert+0.1))

F1b<-F1b+annotate("text", x = 12.5, y = 1, label = "Unbiased Gate")

#Figure 1C sorted clones


  rfuTransform<-function(transformId,r=256,min=0){
 	t=new("transform",.Data=function(x) 10^((x/r)+min))
 	 t@transformationId = transformId
  t
}

channelTrans<-rfuTransform("256 per decade",min=-1) #trasform channel numbers to 10^-1-10^3


plates<-list.dirs("/Users/jonathan/Documents/JF_SD_paperFigs/Data/LGM2_allClones_10kreads/")
LGM2.Clones<-read.flowSet(path=plates[[2]],dataset=2,transformation=FALSE)

for (i in 3:length(plates)) {
	LGM2.Clones<-rbind2(LGM2.Clones,read.flowSet(path=plates[[i]],dataset=2,transformation=FALSE))
}

colnames(LGM2.Clones)[1:3]<-c("FSC","SSC","GFP")
LGM2.Clones.gated<-Subset(LGM2.Clones,norm2Filter("SSC","FSC"))
LGM2.Clones.filtered=Subset(LGM2.Clones.gated,boundaryFilter(x=c("GFP"),side="both"))

LGM2.Clones.RFU<-transform(LGM2.Clones.filtered,`GFP`=channelTrans(`GFP`))
LGM2.Clones.df<-flowset2ggplot(LGM2.Clones.RFU[,c(1:3)])

cloneMeans<-rename(ddply(LGM2.Clones.df, c("name"), function(df)mean(df$GFP)),c(V1="GFP"))
cloneMeans<-cloneMeans[with(cloneMeans, order(GFP)), ]

desired_order<-cloneMeans$name
LGM2.Clones.df$name <- factor( as.character(LGM2.Clones.df$name), levels=desired_order )
LGM2.Clones.df <- LGM2.Clones.df[order(LGM2.Clones.df$name),]
F1c<-ggplot(LGM2.Clones.df,aes(y=name,x=GFP))+ stat_density(aes(fill=..density..), geom="tile", position="identity")
F1c<-F1c+scale_x_log10(name="GFP RFU",breaks=trans_breaks("log10",function(x) 10^x),labels=trans_format("log10",math_format(10^.x)))
F1c<-F1c
F1c<-F1c+theme_bw(18)+theme(axis.text.y = element_blank())+scale_y_discrete(breaks=NULL)+
ylab("227 LGM2 Clones")+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())


	pdf("Figure1.pdf",width=10,height=10,colormodel="cmyk")
	grid.arrange(F1a,arrangeGrob(F1b,F1c,ncol=2),ncol=1)
	dev.off()






#generate Figure 2


fano<-function(x){
	return(var(x)/mean(x))
}

cv<-function(x){
	return(sqrt(var(x))/mean(x))
}


LGM2.Stats<-read.table("/Users/jonathan/Documents/JF_SD_paperFigs/Data/LGM2_NSA_withErrors_11202012.csv",sep=",",header=TRUE)
FISHnames<-as.matrix(LGM2.Stats$Clone)

FISHnames<-apply(FISHnames,1,function(x){strsplit(as.character(x),"_")[[1]][[2]]})
LGM2.Stats$Clone<-FISHnames

LGM2.SelectedClones.df<-subset(LGM2.Clones.df,name %in% FISHnames)
cloneMoments.all<-rename(ddply(LGM2.Clones.df, c("name"), function(df)c(mean(df$GFP),var(df$GFP),cv(df$GFP))),c(V1="Mean",V2="Var",V3="CV"))
cloneMoments.FISH<-rename(ddply(LGM2.SelectedClones.df, c("name"), function(df)c(mean(df$GFP),var(df$GFP),cv(df$GFP))),c(V1="Mean",V2="Var",V3="CV"))
pointSize<-1.3
poissonVar<-function(x){x}
scalingVarAll<-function(x){2*x-0.92}
scalingVarSubset<-function(x){2*x-0.73}
telegraphScaling<-function(x){sqrt(x)}
telegraphMean<-apply(as.matrix(cloneMoments.all$Var),1,telegraphScaling)
cloneMoments.all<-cbind(cloneMoments.all,telegraphMean)




F2a<-ggplot(cloneMoments.all,aes(x=log10(Mean),y=log10(Var)))+geom_point(size=pointSize)+geom_smooth(method=lm)+ylab("Log10(GFP Variance)")+xlab("Log10(<GFP>)")+stat_function(fun=poissonVar,colour="black",linetype="dashed")+stat_function(fun=scalingVarAll,colour="green",linetype="dashed")+theme_grey(14)

F2b<-ggplot(cloneMoments.all,aes(x=log10(Mean),y=log10(CV)))+geom_point(size=pointSize)+geom_smooth(method=lm)+ylab("Log10(GFP CV)")+xlab("Log10(<GFP>)")+theme_grey(14)

F2c<-ggplot(cloneMoments.FISH,aes(x=log10(Mean),y=log10(Var)))+geom_point(size=pointSize)+geom_smooth(method=lm)+ylab("Log10(GFP Variance)")+xlab("Log10(<GFP>)")+stat_function(fun=poissonVar,colour="black",linetype="dashed")+stat_function(fun=scalingVarSubset,colour="green",linetype="dashed")+theme_grey(14)

F2d<-ggplot(cloneMoments.FISH,aes(x=log10(Mean),y=log10(CV)))+geom_point(size=pointSize)+geom_smooth(method=lm)+ylab("Log10(GFP CV)")+xlab("Log10(<GFP>)")+theme_grey(14)


pdf("Figure2_allClones.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(F2a,F2b,F2c,F2d,ncol=2,nrow=2)
dev.off();




#Figure 3 subpanel

BC6data<-read.csv("/Users/jonathan/Documents/JF_SD_paperFigs/Data/LGM2_GFP_FISH_csv/BC6.csv")
scatter<-ggplot(BC6data,aes(x=RNA,y=GFP))+geom_point()
RNAhist<-ggplot(BC6data)+geom_histogram(aes(RNA,colour=red))

pdf("Figure3x.pdf",width=10,height=10,colormodel="cmyk")              
grid.arrange(RNAhist,scatter,nrow=2)
dev.off()             




#Figure 4

cloneMoments.FISH$name<-as.character(droplevels(cloneMoments.FISH)$name)
cloneMoments.FISH<-arrange(cloneMoments.FISH,name)
LGM2.Stats<-arrange(LGM2.Stats,Clone)

F4b<-ggplot(LGM2.Stats,aes(x=log10(RNA.Mean),y=log10(RNA.Var)))+geom_point()+geom_smooth(method=lm)+ylab("Log10(RNA Variance)")+xlab("Log10(<RNA>)")+theme_grey(16)+geom_errorbarh(aes(xmin=log10(RNA.Mean.Min),xmax=log10(RNA.Mean.Max),height=0.02))+geom_errorbar(aes(ymin=log10(RNA.Var.Min),ymax=log10(RNA.Var.Max),width=0.02))

F4c<-ggplot(LGM2.Stats,aes(x=log10(RNA.Mean),y=log10(RNA.CV)))+geom_point()+geom_smooth(method=lm)+ylab("Log10(RNA CV)")+xlab("Log10(<RNA>)")+theme_grey(16)+geom_errorbarh(aes(xmin=log10(RNA.Mean.Min),xmax=log10(RNA.Mean.Max),width=0.02))+geom_errorbarh(aes(xmin=log10(RNA.Mean.Min),xmax=log10(RNA.Mean.Max),height=0.02))+geom_errorbar(aes(ymin=log10(RNA.CV.Min),ymax=log10(RNA.CV.Max)), height=0.02)




F4d<-ggplot(LGM2.Stats,aes(x=log10(RNA.Mean),y=log10(GFP.Mean)))+geom_point()+geom_smooth(method=lm)+ylab("Log10(<GFP>)")+xlab("Log10(<RNA>)")+theme_grey(16)+geom_errorbarh(aes(xmin=log10(RNA.Mean.Min),xmax=log10(RNA.Mean.Max),height=0.02))

F4e<-ggplot(LGM2.Stats,aes(x=log10(RNA.Var),y=log10(GFP.Var)))+geom_point()+geom_smooth(method=lm)+ylab("Log10(GFP Variance)")+xlab("Log10(RNA Variance)")+theme_grey(16)+geom_errorbarh(aes(xmin=log10(RNA.Var.Min),xmax=log10(RNA.Var.Max),height=0.02))


pdf("Figure4bcde.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(F4b,F4c,F4d,F4e,nrow=2,ncol=2)
dev.off();

LGM2.Stats.outliersRemoved<-subset(LGM2.Stats,!Clone %in% c("ID3","CA2"))
LGM2.Stats.noID3<-subset(LGM2.Stats,!Clone %in% c("ID3"))
LGM2.Stats.noCA2<-subset(LGM2.Stats,!Clone %in% c("CA2"))
LGM2.Stats.outliersRemovedB<-subset(LGM2.Stats,!Clone %in% c("BC6","CA2"))

##Figure 5
 ktd<-0.34
 ktdMin<-0.26
 ktdMax<-0.43

F5a<-ggplot(LGM2.Stats.noID3,aes(x=B,y=RNA.Mean))+geom_point()+geom_smooth(method=lm)+ylab("<RNA>")+xlab("Burst Size")+theme_grey(18)+geom_errorbar(aes(ymin=RNA.Mean.Min,ymax=RNA.Mean.Max,width=0.5))+geom_errorbarh(aes(xmin=B.min,xmax=B.max,height=0.5))

F5b<-ggplot(LGM2.Stats.noID3,aes(y=RNA.Mean,x=ka/ktd))+geom_point()+geom_smooth(method=lm)+ylab("<RNA>")+xlab("Normalized On rate")+theme_grey(18)+geom_errorbar(aes(ymin=RNA.Mean.Min,ymax=RNA.Mean.Max,width=0.05))+geom_errorbarh(aes(xmin=ka.Min/ktd,xmax=ka.Max/ktd,height=0.5))

F5c<-ggplot(LGM2.Stats.noID3,aes(y=RNA.CV,x=B))+geom_point()+geom_smooth(method=lm)+ylab("RNA CV")+xlab("Burst Size")+theme_grey(18)+geom_errorbar(aes(ymin=RNA.CV.Min,ymax=RNA.CV.Max,width=0.5))+geom_errorbarh(aes(xmin=B.min,xmax=B.max,height=0.02))

F5d<-ggplot(LGM2.Stats.noID3,aes(y=RNA.CV,x=ka/ktd))+geom_point()+geom_smooth(method=lm)+ylab("RNA CV")+xlab("Normalized On Rate")+theme_grey(18)+geom_errorbar(aes(ymin=RNA.CV.Min,ymax=RNA.CV.Max,width=0.2))+geom_errorbarh(aes(xmin=ka.Min/ktd,xmax=ka.Max/ktd,height=0.02))



pdf("Figure5.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(F5a,F5b,F5c,F5d,ncol=2,nrow=2)
dev.off();

##Figure 7


F7a<-ggplot(LGM2.Stats.noCA2,aes(x=Nuc1,y=ka/ktd))+geom_point()+geom_smooth(method=lm)+xlab("Nuc-1 Chromatin Inaccessibility")+ylab("Normalized On Rate")+theme_grey(18)+geom_errorbarh(aes(xmin=Nuc1-Nuc1.Error,xmax=Nuc1+Nuc1.Error,height=0.05))+geom_errorbar(aes(ymin=ka.Min/ktd,ymax=ka.Max/ktd,width=0.02))

F7b<-ggplot(LGM2.Stats.NSA,aes(x=log(HSS),y=kr/ktd))+geom_point()+geom_smooth(method=lm)+xlab("Log(HSS Chromatin Inaccessibility)")+ylab("Normalized Off Rate")+theme_grey(16)+ylim(-10,40)+labs(title="Exponential Fit")

F7c<-ggplot(LGM2.Stats.NSA,aes(x=log(HSS),y=kr/ktd))+geom_point()+geom_smooth()+xlab("Log(HSS Chromatin Inaccessibility)")+ylab("Normalized Off Rate")+theme_grey(16)+ylim(-10,40)+labs(title="LOESS")

pdf("Figure7.pdf",width=6,height=3,colormodel="cmyk")
grid.arrange(F7a,F7b,F7c,nrow=1)
dev.off()


#Supp Figure- FSC vs. GFP

LGM2.Clones.Channels.df<-flowset2ggplot(LGM2.Clones.filtered[,c(1:3)])
LGM2.FSCvGFP<-ddply(LGM2.Clones.Channels.df, c("name"), function(df)(c(Pearson.r=rcorr(df$GFP,df$FSC,type=("pearson"))$r[[1,2]],Spearman.r=rcorr(df$GFP,df$FSC,type=("spearman"))$r[[1,2]],Pearson.p=rcorr(df$GFP,df$FSC,type=("pearson"))$P[[1,2]],Spearman.p=rcorr(df$GFP,df$FSC,type=("pearson"))$P[[1,2]],r.squared=summary(lm(df$GFP~df$FSC))$r.squared)))
FSCGFP_pearsonR2r<-ggplot(LGM2.FSCvGFP,aes(x=Pearson.r,y=r.squared))+geom_point(aes(colour=name))+geom_rug()+guides(colour=FALSE)
FSCGFP_pearsonR2p<-ggplot(LGM2.FSCvGFP,aes(x=Pearson.p,y=r.squared))+geom_point(aes(colour=name))+geom_rug()+guides(colour=FALSE)
FSCGFP_spearmanR2r<-ggplot(LGM2.FSCvGFP,aes(x=Spearman.r,y=r.squared))+geom_point(aes(colour=name))+geom_rug()+guides(colour=FALSE)
FSCGFP_spearmanR2p<-ggplot(LGM2.FSCvGFP,aes(x=Spearman.p,y=r.squared))+geom_point(aes(colour=name))+geom_rug()+guides(colour=FALSE)
pdf("AllClones_FSCvGFP_uncorrelated.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(FSCGFP_pearsonR2r,FSCGFP_spearmanR2r,FSCGFP_pearsonR2p,FSCGFP_spearmanR2p, nrow=2,ncol=2)
dev.off()


pdf("SelectedClones_FSCvGFP.pdf",width=10,height=10,colormodel="cmyk")
ggplot(LGM2.SelectedClones.df,aes(y=log10(GFP),x=log10(FSC)))+geom_point(shape=1)+facet_wrap(~name,ncol=4)+geom_smooth(method=lm,aes(colour="red"))
dev.off()

pdf("FISHclones_gating.pdf",width=10,height=10,colormodel="cmyk")
print(xyplot(`SSC` ~ `FSC`,LGM2.Clones[which(sampleNames(LGM2.Clones) %in% FISHnames)],filter=norm2Filter("SSC","FSC")))
dev.off()

#Supp Figure-RNA Degradation

RNAdeg<-read.csv("/Users/jonathan/Documents/JF_SD_paperFigs/Data/RNA_degradation.csv")
dfc <- summarySE(RNAdeg, measurevar="GFP", groupvars=c("Time","Condition"))
LGM2degmodel<-lm(log(RNAdeg[16:24,]$GFP)~RNAdeg[16:24,]$Time)
BactModel<-lm(log(RNAdeg[25:36,]$GFP)~RNAdeg[25:36,]$Time)
LGM2cntrlModel<-lm(log(RNAdeg[1:12,]$GFP)~RNAdeg[1:12,]$Time)
BActcntrlModel<-lm(log(RNAdeg[37:51,]$GFP)~RNAdeg[37:51,]$Time)



LGM2deg<-   ggplot(RNAdeg[1:24,], aes(x=Time, y=log(GFP), colour=Condition))+geom_point()+stat_smooth(data=RNAdeg[16:24,],method=lm)+stat_smooth(data=RNAdeg[1:12,],method=lm)+theme_grey(16)+theme(legend.position="bottom")+ylab("ln(LGM2 RNA) Arbitrary Units")
BActdeg<-    ggplot(RNAdeg[25:51,], aes(x=Time, y=log(GFP), colour=Condition))+geom_point()+stat_smooth(data=RNAdeg[25:36,],method=lm)+stat_smooth(data=RNAdeg[37:51,],method=lm)+theme_grey(16)+theme(legend.position="bottom")+ylab("ln(B-Actin RNA) Arbitrary Units")
pdf("RNA_degradation.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(LGM2deg,BActdeg,ncol=1)
dev.off()

#Supp Figure-Spearman Matrix

LGM2.Stats.reduced<-data.frame(LGM2.Stats.outliersRemoved$RNA.Mean,LGM2.Stats.outliersRemoved$RNA.Var,LGM2.Stats.outliersRemoved$RNA.CV,LGM2.Stats.outliersRemoved$GFP.Mean,LGM2.Stats.outliersRemoved$GFP.Var,LGM2.Stats.outliersRemoved$GFP.CV,LGM2.Stats.outliersRemoved$ka,LGM2.Stats.outliersRemoved$B,LGM2.Stats.outliersRemoved$Nuc0,LGM2.Stats.outliersRemoved$HSS,LGM2.Stats.outliersRemoved$Nuc1)
colnames(LGM2.Stats.reduced)<-c("RNA.Mean","RNA.Var","RNA.CV","GFP.Mean","GFP.Var","GFP.CV","ka","B","Nuc0","HSS","Nuc1")

pdf("Spearman_visualization_p025.pdf",width=10,height=10,colormodel="cmyk")
corrplot.mixed(rcorr(as.matrix(LGM2.Stats.reduced),type=c("spearman"))$r,lower=c("number"),upper=c("circle"),p.mat=rcorr(as.matrix(LGM2.Stats.reduced),type=c("spearman"))$P,tl.pos=c("lt"),sig.level=0.025,insig=c("pch"),diag=c("n"))
dev.off()


#supp Figure Cell Size

files<-list.files(path="../Data/CellSize/")
cellSizes<-read.csv(paste0("../Data/CellSize/",files[1]))
cellSizes<-cbind(Clone=as.matrix(rep(strsplit(files[1],".csv")[1]),dim(cellSizes)[1]),cellSizes)

for (i in 2:length(files)){
	cellSize<-read.csv(paste0("../Data/CellSize/",files[i]))
	cellSizes<-rbind(cbind(Clone=as.matrix(rep(strsplit(files[i],".csv")[1]),dim(cellSize)[1]),cellSize),cellSizes)
	
}

cellSizes$Clone<-as.character(cellSize$Clone)
CA2<-subset(cellSizes,Clone %in% c("CA2"))
IC4<-subset(cellSizes,Clone %in% c("IC4"))
IB4<-subset(cellSizes,Clone %in% c("IB4"))
AD1<-subset(cellSizes,Clone %in% c("AD1"))
EC5<-subset(cellSizes,Clone %in% c("EC5"))


CA2plot<-ggplot(CA2,aes(x=CellSize))+geom_histogram(aes(y=..density..))+stat_function(fun=dnorm, args=list(mean=mean(CA2$CellSize),sd=sd(CA2$CellSize)),aes(colour="red"))
AD1plot<-ggplot(AD1,aes(x=CellSize))+geom_histogram(aes(y=..density..))+stat_function(fun=dnorm, args=list(mean=mean(AD1$CellSize),sd=sd(AD1$CellSize)),aes(colour="red"))
IC4plot<-ggplot(IC4,aes(x=CellSize))+geom_histogram(aes(y=..density..))+stat_function(fun=dnorm, args=list(mean=mean(IC4$CellSize),sd=sd(IC4$CellSize)),aes(colour="red"))
IB4plot<-ggplot(IB4,aes(x=CellSize))+geom_histogram(aes(y=..density..))+stat_function(fun=dnorm, args=list(mean=mean(IB4$CellSize),sd=sd(IB4$CellSize)),aes(colour="red"))
EC5plot<-ggplot(EC5,aes(x=CellSize))+geom_histogram(aes(y=..density..))+stat_function(fun=dnorm, args=list(mean=mean(EC5$CellSize),sd=sd(EC5$CellSize)),aes(colour="red"))

cellSize.boxplot<-ggplot(cellSizes,aes(factor(Clone),CellSize))+geom_boxplot()

pdf("CellSize_wNorm.pdf",width=10,height=10,colormodel="cmyk")
grid.arrange(cellSize.boxplot,arrangeGrob(AD1plot,CA2plot,IC4plot,EC5plot,nrow=2,ncol=2),nrow=2)
dev.off()






#Supp Figure FISH intensity

FISHIntensity<-rename(read.csv("/Users/jonathan/Documents/JF_SD_paperFigs/Data/FISHIntensity.csv",header=FALSE),c(V1="Intensity"))

pdf("FISH_IntensitywNorm.pdf",width=5,height=5,colormodel="cmyk")
ggplot(FISHIntensity,aes(x=Intensity))+geom_histogram(aes(y=..density..),binwidth=3)+stat_function(fun=dnorm, args=list(mean=mean(FISHIntensity$Intensity),sd=sd(FISHIntensity$Intensity)),aes(colour="red"))+theme(legend.position="none")+xlab("Intensity Arbitrary Units")
dev.off()


#Supp Figure ka v. B

pdf("ka_B_regression.pdf",width=5,height=5,colormodel="cmyk")
ggplot(LGM2.Stats.noID3,aes(y=B,x=ka/ktd))+geom_point()+geom_smooth(method=lm)+ylab("Burst Size")+xlab("Normalized On rate")+theme_grey(16)+geom_errorbar(aes(ymin=B.min,ymax=B.max,width=0.2))+geom_errorbarh(aes(xmin=ka.Min/ktd,xmax=ka.Max/ktd,height=0.2))
dev.off()

#Supp  Figure Nuc0,HS vs. ka

Nuc0ka<-ggplot(LGM2.Stats.noCA2,aes(x=Nuc0,y=ka/ktd))+geom_point()+geom_smooth(method=lm)+xlab("Nuc0 Chromatin Inaccessibility")+ylab("Normalized On Rate")+theme_grey(16)+geom_errorbarh(aes(xmin=Nuc0-Nuc0.Error,xmax=Nuc0+Nuc0.Error,height=0.05))+geom_errorbar(aes(ymin=ka.Min/ktd,ymax=ka.Max/ktd,width=0.02))
HSSka<-ggplot(LGM2.Stats.noCA2,aes(x=HSS,y=ka/ktd))+geom_point()+geom_smooth(method=lm)+xlab("HSS Chromatin Inaccessibility")+ylab("Normalized On Rate")+theme_grey(16)+geom_errorbarh(aes(xmin=HSS-HSS.Error,xmax=HSS+HSS.Error,height=0.05))+geom_errorbar(aes(ymin=ka.Min/ktd,ymax=ka.Max/ktd,width=0.02))
pdf("Nuc0HSS_ka.pdf",width=8,height=4,colormodel="cmyk")
grid.arrange(Nuc0ka,HSSka,ncol=2)
dev.off()


#Regression Models...this should be refactored by creating a mixed data frame and then feeding vector of regression pairs to
#some function that generates all this in a nice compact fashion.... this is a copy paste hell

MuVarModel<-lm(log10(cloneMoments.all$Var)~log10(cloneMoments.all$Mean))
MuVarTgraphModel<-lm(log10(cloneMoments.all$Var)~log10(cloneMoments.all$telegraphMean))
MuVarCorr<-rcorr(log10(cloneMoments.all$Var),log10(cloneMoments.all$Mean),type=c("spearman"))

MuCVModel<-lm(log10(cloneMoments.all$CV)~log10(cloneMoments.all$Mean))
MuCVCorr<-rcorr(log10(cloneMoments.all$CV),log10(cloneMoments.all$Mean),type=c("spearman"))

SubsetMuVarModel<-lm(log10(cloneMoments.FISH$Var)~log10(cloneMoments.FISH$Mean))
SubsetMuVarCorr<-rcorr(log10(cloneMoments.FISH$Var),log10(cloneMoments.FISH$Mean),type=c("spearman"))

SubsetMuCVModel<-lm(log10(cloneMoments.FISH$CV)~log10(cloneMoments.FISH$Mean))
SubsetMuCVCorr<-rcorr(log10(cloneMoments.FISH$CV),log10(cloneMoments.FISH$Mean),type=c("spearman"))

RNAMuCVModel<-lm(log10(LGM2.Stats$RNA.CV)~log10(LGM2.Stats$RNA.Mean))
RNAMuCVCorr<-rcorr(log10(LGM2.Stats$RNA.CV),log10(LGM2.Stats$RNA.Mean),type=c("spearman"))

RNAMuVarModel<-lm(log10(LGM2.Stats$RNA.Var)~log10(LGM2.Stats$RNA.Mean))
RNAMuVarCorr<-rcorr(log10(LGM2.Stats$RNA.Var),log10(LGM2.Stats$RNA.Mean),type=c("spearman"))

RNAMuGFPMuModel<-lm(log10(LGM2.Stats$GFP.Mean)~log10(LGM2.Stats$RNA.Mean))
RNAMuGFPMuCorr<-rcorr(log10(LGM2.Stats$GFP.Mean),log10(LGM2.Stats$RNA.Mean),type=c("spearman"))

RNAVarGFPVarModel<-lm(log10(LGM2.Stats$GFP.Var)~log10(LGM2.Stats$RNA.Var))
RNAVarGFPVarCorr<-rcorr(log10(LGM2.Stats$GFP.Var),log10(LGM2.Stats$RNA.Var),type=c("spearman"))
MuBModel<-lm(LGM2.Stats.noID3$RNA.Mean~LGM2.Stats.noID3$B)
MuBCorr<-rcorr(LGM2.Stats.noID3$RNA.Mean,LGM2.Stats.noID3$B,type=c("spearman"))

MukaModel<-lm(LGM2.Stats.noID3$RNA.Mean~LGM2.Stats.noID3$ka)
MukaCorr<-rcorr(LGM2.Stats.noID3$RNA.Mean,LGM2.Stats.noID3$ka,type=c("spearman"))

CVBModel<-lm(LGM2.Stats.noID3$RNA.CV~LGM2.Stats.noID3$B)
CVBCorr<-rcorr(LGM2.Stats.noID3$RNA.CV,LGM2.Stats.noID3$B,type=c("spearman"))
CVkaModel<-lm(LGM2.Stats.noID3$RNA.CV~LGM2.Stats.noID3$ka)
CVkaCorr<-rcorr(LGM2.Stats.noID3$RNA.CV,LGM2.Stats.noID3$ka,type=c("spearman"))
Nuc1kaModel<-lm(LGM2.Stats.noCA2$ka~LGM2.Stats.noCA2$Nuc1)
Nuc1kaCorr<-rcorr(LGM2.Stats.noCA2$Nuc1,LGM2.Stats.noCA2$ka,type=c("spearman"))
Nuc1CVModel<-lm(LGM2.Stats.noCA2$RNA.CV~LGM2.Stats.noCA2$Nuc1)
Nuc1CVCorr<-rcorr(LGM2.Stats.noCA2$RNA.CV,LGM2.Stats.noCA2$Nuc1,type=c("spearman"))
HSSkaModel<-lm(LGM2.Stats.noCA2$ka~LGM2.Stats.noCA2$HSS)
HSSkaCorr<-rcorr(LGM2.Stats.noCA2$HSS,LGM2.Stats.noCA2$ka,type=c("spearman"))
Nuc0kaModel<-lm(LGM2.Stats.noCA2$ka~LGM2.Stats.noCA2$Nuc0)
Nuc0kaCorr<-rcorr(LGM2.Stats.noCA2$Nuc0,LGM2.Stats.noCA2$ka,type=c("spearman"))


kaBModel<-lm(LGM2.Stats.noID3$B~LGM2.Stats.noID3$ka)
kaBCorr<-rcorr(LGM2.Stats.noID3$B,LGM2.Stats.noID3$ka,type=c("spearman"))

models<-list(MuVarModel,MuCVModel,SubsetMuVarModel,SubsetMuCVModel,RNAMuCVModel,RNAMuVarModel,RNAMuGFPMuModel,RNAVarGFPVarModel,MuBModel,MukaModel,CVBModel,CVkaModel,Nuc1kaModel,Nuc1CVModel,kaBModel)

corrs<-list(MuVarCorr,MuCVCorr,SubsetMuVarCorr,SubsetMuCVCorr,RNAMuCVCorr,RNAMuVarCorr,RNAMuGFPMuCorr,RNAVarGFPVarCorr,MuBCorr,MukaCorr,CVBCorr,CVkaCorr,Nuc1kaCorr,Nuc1CVCorr,kaBCorr)



modelStats<-lapply(models,function(x){data.frame(explanatory=as.character(x$call$formula)[[3]],response=as.character(x$call$formula)[[2]],intercept=x$coefficients[[1]],slope=x$coefficients[[2]],slope95CI=x$coefficients[[2]]-confint(x)[[2,1]],r.squared=summary(x)$r.squared,regression.p=summary(x)$coefficients[[2,4]])})


corrStats<-lapply(corrs,function(x){data.frame(Spearman.r=x$r[[1,2]],corelation.p=x$P[[1,2]])})
corrStats.df=corrStats[[1]];

modelStats.df=cbind(modelStats[[1]],corrStats[[1]])
for(i in 2:length(modelStats)){
	modelStats.df=rbind(modelStats.df,cbind(modelStats[[i]],corrStats[[i]]))
}
