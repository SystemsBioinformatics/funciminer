
#create overview heatmap for genera 
data<-read.table("path/tax_of_kos_abund_lbpnorm.csv")

sd_data<-apply(data, MARGIN = 2, function(x) sd(x))
plot(sd_data)
data_r<-data[, which(sd_data>mean(sd_data))]

library(apcluster)
apres<-apcluster(corSimMat(),t(data_r),q=0.9)
aggres <- aggExCluster(x=apres)
plot(aggres)
apres<-apclusterK(corSimMat(),t(data_r),K=6)


library(lattice)
levelplot (as.matrix(data_r)[,unlist(apres@clusters)], xlab="Samples",ylab="Genera",main="" ,col.regions = colorRampPalette(c("lightblue","blue","green", "yellow","orange","red", "black")), aspect = "fill",scale=list(x=list(cex=0.7,rot=45)))


meta_data <- read.csv("path/meta_data_ko.csv", sesp="")

######## diverstiy calculation
Variety<-meta_data$variety
inx<-1:length(Variety)
abundance.matrix<-data
Richness <- rowSums(abundance.matrix>0)
Shannon <- diversity(abundance.matrix)
datagg<-as.data.frame(cbind(Richness,Shannon,inx))
Samples<-rownames(datagg)
datagg<-cbind(datagg,Variety,Samples)

ggplot(datagg, aes(x=Samples, y=Richness,color=Variety) ) +
  geom_point(shape=19,cex=3)+scale_color_manual(values=c("yellow", "red", "purple"))+theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1))+scale_x_discrete(limits= datagg$Samples)+
  theme(legend.text=element_text(size=15),legend.title=element_text(size=20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size=6))

ggplot(datagg, aes(x=Samples, y=Shannon,color=Variety) ) +
  geom_point(shape=17,cex=3)+scale_color_manual(values=c("yellow", "red", "purple"))+theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1))+scale_x_discrete(limits= datagg$Samples)+
  theme(legend.text=element_text(size=15),legend.title=element_text(size=20),axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size=6))
