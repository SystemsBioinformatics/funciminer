

data<-read.table("path/ko_taxfull_abund.csv")

library(apcluster)
apres1<-apcluster(corSimMat(),data,q=0.9)
aggres2 <- aggExCluster(x=apres1)
#plot histogram to idetifiy best number of clusters
plot(aggres2)
apres1<-apclusterK(corSimMat(),data,K=6)

heatmap(apres1,col=terrain.colors(12),margins=c(5, 10,0),dendScale=.3,main="",cexRow=0,cexCol=1,Colv=NA)
