R_input_ph_exp <- read.csv("~/Desktop/Skin_experiment/R_input_skin_exp_mean.csv")
x<-R_input_ph_exp
rownames(x)<-x[,1]
x<-x[,-1]
x<-as.matrix(x)
x[,1:3]<-log10(x[,1:3])
library(ggplot2)
dx<-as.data.frame(x)
dx$Var<-as.factor(dx$Var)
dxl<-split(dx,dx$Var)
library(reshape)
plotlist<-list()
plotlist1<-list()
names<-rownames(x)[seq(2,length(rownames(x)),2)]
for (i in 1:length(dxl)){
  mdata <- melt(dxl[[i]], id=c("Var","skin")) 
  colnames(mdata)[4]<-"malic"
  colnames(mdata)[3]<-"Day_of_sampling"
  mdata$skin[mdata$skin==1]="yes"
  mdata$skin[mdata$skin==0]="no"
  plotlist1[[i]]<-mdata
  p<-ggplot(data=mdata, aes(y=malic, x=Day_of_sampling, group=skin )) +
    geom_line(aes(color=skin),size=1.3)+
    geom_point(aes(color=skin),size=3)+theme(legend.position="none")+ ylim(0,max(plotlist1[[i]]$malic))+theme(axis.title.x=element_blank(),
                                                                                     axis.title.y=element_blank() )+ggtitle(names[i])
  plotlist[[i]]<-p
  
}
                
library(gridExtra)
grid.arrange(plotlist[[1]],plotlist[[3]],plotlist[[4]],plotlist[[5]],plotlist[[6]] , plotlist[[7]] , plotlist[[8]] , plotlist[[9]],plotlist[[10]], nrow = 3,ncol=3)
grid.arrange(plotlist[[2]],plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]] , plotlist[[15]], nrow = 2,ncol=3)

