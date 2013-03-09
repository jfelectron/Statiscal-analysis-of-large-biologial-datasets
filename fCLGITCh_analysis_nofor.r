 library(flowCore)
 library(flowQ)
 library(flowViz)
 library(flowStats)
 library(GMD)
 library(psych)

 
 #import data
flowData<-read.flowSet(path="/Users/jonathan/Desktop/100412 fCLGITCh Clones 2nd read GFP-mCherry/Clones/",transformation=FALSE)
 #replace column names to change from PMT names to desired ones
 colnames(flowData) <-c("FSC","SSC","GFP","Tat-mCherry","Time")
 print(xyplot(`SSC` ~ `FSC`,flowData[1:25],filter=norm2Filter("SSC","FSC")))
 dev.new()
 #apply gate to data
 gatedData<-Subset(flowData,norm2Filter("SSC","FSC"))

#remove events in 1st and last channels 
filteredData=Subset(gatedData,boundaryFilter(x=c("GFP","Tat-mCherry"),side="both"))
 #generate plot of GFP v. mCherry for first 12 clones
 print(xyplot(`GFP`~`Tat-mCherry`,filteredData[1:25]))
 dev.new()
 #generate side by side histograms of the two channesl
 print(densityplot(~`GFP`+`Tat-mCherry`,breaks=256,filteredData[1:25]))
 
 mCherryBins=matrix(data=0,nrow=length(filteredData),ncol=1024)
 GFPBins=matrix(data=0,nrow=length(filteredData),ncol=1024)
 
 
for (i in 1:dim(mCherryBins)[1]) {
 	mCh<-exprs(filteredData[[i]])[,"Tat-mCherry"]
 	GFP<-exprs(filteredData[[i]])[,"GFP"]
 	 
	for(j in 1:length(mCh)){
	 	mCherryBins[i,mCh[j]]<-mCherryBins[i,mCh[j]]+1
	 	GFPBins[i,GFP[j]]<-GFPBins[i,GFP[j]]+1
	 	}
	 }
mCherryData<-list(mCherryBins[1,])
GFPData<-list(GFPBins[1,])
for(k in 2:dim(mCherryBins)[1]){
	mCherryData[k]<-list(mCherryBins[k,])
	GFPData[k]<-list(GFPBins[k,])
	}
mDistances<-gmdm(mCherryData,keyword(filteredData,"GUID"),sliding=FALSE)
gDistances<-gmdm(GFPData,keyword(filteredData,"GUID"),sliding=FALSE)
mDist<-gmdm2dist(mDistances)
gDist<-gmdm2dist(gDistances)
mhClust<-hclust(mDist,method="ward")
ghClust<-hclust(gDist,method="ward")
mCSSmulti <- css.hclust(mDist,mhClust)
gCSSmulti <- css.hclust(gDist,ghClust)
mElbow <- elbow.batch(mCSSmulti,ev.thres=0.90,inc.thres=0.05)
gElbow<-elbow.batch(gCSSmulti,ev.thres=0.90,inc.thres=0.05)
gCutree <- cutree(ghClust,k=30)
phenoData(filteredData)$gCluster<-as.vector(gCutree)
clusters<-split(filteredData,pData(filteredData)$gClust)
pdf("fCLGITCh_Ward_Clusters_Density.pdf",width=7,height=7)
for(i in 1:length(clusters)){
print(densityplot(~`GFP`+`Tat-mCherry`,main=paste("Cluster",i,sep=" "),groups=gCluster,clusters[[i]]))
}
dev.off()

pdf("fCLGITCh_Ward_Clusters_SmoothedScatter",width=7,height=7)
for(i in 1:length(clusters)){
print(xyplot(`GFP`~`Tat-mCherry`,main=paste("Cluster",i,sep=" "),clusters[[i]]))
}
dev.off()
 
 
 
for(i in 1:length(clusters)){
	for(j in 1:length(clusters[[i]])) {
		write.table(exprs(clusters[[i]][[j]]), file=paste("Cluster",i,keyword(clusters[[i]][[j]],"GUID"),sep="_"),sep=",",row.names=FALSE,col.names=FALSE)
}

}

  bimodality<-function(x){skew(x)^2/(kurtosi(x)+1)}
  clusterMeans<-lapply(clusters,function(x){fsApply(x[,3:4],each_col,mean)})
  clusterVars<-lapply(clusters,function(x){fsApply(x[,3:4],each_col,var)})
  clusterSkews<-lapply(clusters,function(x){fsApply(x[,3:4],each_col,skew)})
  clusterBimodality<-lapply(clusters,function(x){fsApply(x[,3:4],each_col,bimodality)})
  
 
 
 plot(as.vector(lapply(clusterMeans,function(x){x[,1]})),as.vector(lapply(clusterBimodality,function(x){x[,1]})),xlab="<GFP>",ylab="GFP Bimodality coefficient")