outputfinalnormb <- read.csv("path/7PVA_Bobal_idba/merged_all_bins/table_final_abund.csv")
temp<-colnames(outputfinalnormb)[-ncol(outputfinalnormb)]
outputfinalnormb<-outputfinalnormb[,-1]
colnames(outputfinalnormb)<-temp
outputfinalnormtemp<- read.table("path/THLR_Tempranillo_maxbin/tank1_res_maxbin/table_final_abund.csv" ,header=TRUE, sep="\t")
outputfinalnormtemp<-outputfinalnormtemp[,-1]
outputfinalnorma <- read.csv("path/VECX_Airen_idba/tank1plus_maxbin/table_final_abund.csv",sep = "\t")

lbp_bin_abund<-unlist(c(outputfinalnormb[5,],outputfinalnormtemp[13,],outputfinalnorma[14,]))
lbb_bin_abund<-unlist(c(rep(0,75-length(outputfinalnorma[3,])),outputfinalnorma[3,]))

all_lab<-rbind(lbp_bin_abund,lbb_bin_abund)

data<-all_lab

#rownames(data)<-c("Presense","Absense")
datatest<-cbind(data[,1:8],c(NA,NA),data[,9:17],c(NA,NA),data[,18:24],c(NA,NA),data[,25:34],c(NA,NA),data[,35:42],c(NA,NA),data[,43:51],c(NA,NA),data[,52:60],c(NA,NA),data[,61:68],c(NA,NA),data[,69:75])

t1<-adjustcolor( "#33a02c", alpha.f = 1)
t2<-adjustcolor( "#b2df8a", alpha.f = 1)
bp1<-barplot(datatest,main = "",xaxt= "n",border=c("#33a02c","#b2df8a"),col=c(t1,t2),axes = T)

grid(nx=NA, ny=NULL)
vecb<-c("B_Inoc.",rep("",7),"","B_Control1",rep("",8),"","B_Control2",rep("",6),"","T_Inoc.",rep("",9),"","T_Control1",rep("",7),"","T_Control2",rep("",8),"","A_Inoc.",rep("",8),"","A_Control1",rep("",7),"","A_Control2",rep("",6))
ylim=c(0,100)
axis( 1, at=bp1, labels=F, las= 1,lwd.ticks=1,col.ticks ="black",cex.axis=1)
axis( 1, at=bp, labels=vecb, las= 2,lwd.ticks=2,col.ticks ="black",cex.axis=1)
#
legend(ncol(datatest)+10,0.6, col=c("#33a02c","#b2df8a") , pch = 20,legend = c("Lb.P abund.","Lb.B abund."),cex=1.5,pt.cex =3)


##############################################KO pres/abs###########################################################################

##################################################################################################################################
#read KO binary table
x_ko<-read.table("path/ko_wine_binary.tsv")
#read L.plantarum KOs
read <- read.table("path/lbp_ko.txt",sep="\t", header=FALSE, fill = TRUE)
lbp_ko<-grep("^K",as.matrix(read),value=TRUE)
lbp_ko<-unique(lbp_ko)

temp_test<-apply(x_ko, 1, function(x) lbp_ko%in%names(which(x==1)))
temp_test_table<-apply(temp_test, 2, function(x) table(x))


data<-temp_test_table[c(2,1),]
rownames(data)<-c("Presense","Absense")
datatestko<-cbind(data[,1:8],c(NA,NA),data[,9:17],c(NA,NA),data[,18:24],c(NA,NA),data[,25:34],c(NA,NA),data[,35:42],c(NA,NA),data[,43:51],c(NA,NA),data[,52:60],c(NA,NA),data[,61:68],c(NA,NA),data[,69:75])

t1<-adjustcolor( "#1f78b4", alpha.f = 1)
t2<-adjustcolor( "#a6cee3", alpha.f = 1)
bp1<-barplot(datatestko[c(1,2),],xaxt= "n",border=c("#1f78b4","#a6cee3"),col=c(t1,t2),axes = T)

grid(nx=NA, ny=NULL)
vecb<-c("B_Inoc.",rep("",7),"","B_Control1",rep("",8),"","B_Control2",rep("",6),"","T_Inoc.",rep("",9),"","T_Control1",rep("",7),"","T_Control2",rep("",8),"","A_Inoc.",rep("",8),"","A_Control1",rep("",7),"","A_Control2",rep("",6))
ylim=c(0,100)
axis( 1, at=bp1, labels=F, las= 1,lwd.ticks=1,col.ticks ="black",cex.axis=1)
axis( 1, at=bp, labels=vecb, las= 2,lwd.ticks=2,col.ticks ="black",cex.axis=1)

legend(1,1, col=c("#1f78b4","#a6cee3") , pch = 20,legend = c("Presense KOs","Absense KOs"),cex=1.5,pt.cex = 3)

