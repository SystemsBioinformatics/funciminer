######################################find lbp count ko per genus 
#load ko_tax_table and lbp ko =lbp_ko
data<-read.table("path/ko_taxfull_abund.csv")

read <- read.table("path/lbp_ko.txt",sep="\t", header=FALSE, fill = TRUE)
lbp_ko<-grep("^K",as.matrix(read),value=TRUE)
lbp_ko<-unique(lbp_ko)

unique_genera<-unique(gsub(".*_","",colnames(data)))


data_lac<-data[,which(gsub("_.*","",colnames(data))%in%lbp_ko==T)]
tab<-c()
toremove<-c()
for (i in 1:length(unique_genera)){
  datatemp<-data[,which(gsub(".*_","",colnames(data))==unique_genera[i])]
  if(is.data.frame(datatemp)==T){
    t<-apply(datatemp, 1, function(x) length(which(gsub("_.*","",colnames(datatemp)[which(x>0)])%in%lbp_ko==T)))
    tab<-cbind(tab,t)
  }else{
    toremove<-c(toremove,i)
  }
}

unique_genera<-unique_genera[-toremove]
colnames(tab)<-unique_genera


t<-apply(tab, 2, function(x) if (length(which(x>300)>0)){T}else{F})
tabf<-tab[,which(t==T)]

library(apcluster)
apres1<-apcluster(corSimMat(),tabf,q=0.9)
aggres2 <- aggExCluster(x=apres1)
plot(aggres2)
apres1<-apclusterK(corSimMat(),t(tabf),K=15)
library(lattice)
levelplot (as.matrix(tabf)[,unlist(apres1@clusters)], xlab="",ylab="",main="" ,col.regions = colorRampPalette(c("lightblue","blue","yellow","orange","red", "black")), aspect = "fill",scale=list(y=list(cex=.8),x=list(cex=1,rot=45)))
