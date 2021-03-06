read_data<-function(path){
  data <- read.csv(path)
  rownames(data)<-data[,1]
  data<-data[,-1]
  data<-t(data)
  data
}

wine_guys_pathways<-read_data("path/pw_count_ko_numbers_oe_sce_lp_lb.csv")

rownames(wine_guys_pathways)<-c("L.plantarum","O.oeni_airen","O.oeni_bobal","O.oeni_tempranillo","genus.sacchar.","L.brevis","number_kos")

data<-wine_guys_pathways
data<-data[,!is.na(as.character(data[nrow(data),]))]

pathways_norm<-apply(rbind(data[-nrow(data),],data[nrow(data),]), MARGIN = 2, function(x) x/x[length(x)]*100)
pathways_norm<-pathways_norm[-nrow(pathways_norm),]
#remove yeast
pathways_norm<-pathways_norm[-5,]
#
 t<-apply(pathways_norm,2,function(x) all(as.integer(names(table(x)))<=2))
 datar<-pathways_norm[,!t]


library(KEGGREST)
templist<-list()
# high l abund
datat<-datar[c(1,5,2,3,4),]
colnames(datat)<-gsub("_.*","",colnames(datat))

for (j in 1:length(colnames(datat))){
  templist[[j]]<-keggGet(colnames(datat)[j])
}
pathwaynames<-lapply(templist, function(x) x[[1]]$NAME)
colnames(datat)<-paste(unlist(pathwaynames),gsub(".*_","",colnames(datat)),sep="-")

#filter high sd
sd_data<-apply(datat, MARGIN = 2, function(x) sd(x))
plot(sd_data)
data_r<-datat[, which(sd_data>mean(sd_data))]


library(apcluster)
apres1<-apcluster(corSimMat(),t(data_r),q=0.9)
aggres2 <- aggExCluster(x=apres1)
plot(aggres2)
apres1<-apclusterK(corSimMat(),t(data_r),K=7)
library(lattice)
levelplot (as.matrix(data_r)[,unlist(apres1@clusters)], xlab="",ylab="",main="" ,col.regions = colorRampPalette(c("lightblue","blue","yellow","orange","red", "black")), aspect = "fill",scale=list(y=list(cex=.8),x=list(cex=1,rot=45)))

