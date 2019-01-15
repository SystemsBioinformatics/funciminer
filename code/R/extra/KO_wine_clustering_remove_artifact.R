data<-read.table("path/ko_tax_abund.csv")

library(apcluster)
test_data<-data[,-grep("*_Vitis",colnames(data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Vitis",cexRow=0,cexCol=1.0,Colv=NA,legend="col")
test_data<-data[,-grep("*_Sacch",colnames(data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Sacch",cexRow=0,cexCol=1.0,Colv=NA,legend="col")



test_data<-data[,-grep("*_Vitis",colnames(data))]
test_data<-test_data[,-grep("*_Sacch",colnames(test_data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Vitis and Sacch",cexRow=0,cexCol=1.0,Colv=NA,legend="col")




test_data<-data[,-grep("*_Oenoc",colnames(data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Oenoc",cexRow=0,cexCol=1.0,Colv=NA,legend="col")
test_data<-data[,-grep("*_Lacto",colnames(data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Lacto",cexRow=0,cexCol=1.0,Colv=NA,legend="col")

test_data<-data[,-grep("*_Oenoc",colnames(data))]
test_data<-test_data[,-grep("*_Lacto",colnames(test_data))]
aprer_data<-apclusterK(corSimMat(),test_data,K=6)
heatmap(aprer_data,col=terrain.colors(12),dendScale=.3,main="without Oenoc and Lacto",cexRow=0,cexCol=1.0,Colv=NA,legend="col")
