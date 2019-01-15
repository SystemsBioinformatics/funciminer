
##############################################binning visulization################################################################

library("robCompositions", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("Rtsne", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
library("Biostrings")
library(rgl)
library(evd)
library(scatterplot3d)
library(apcluster)
library(gridExtra)
library(ShortRead)

setwd("path/VECX_Airen_idba/tank1plus_maxbin")

inp <- list.files(pattern = "*.fasta")
counter<-1
bin<-c()
names_all<-c()
All<-DNAStringSet()
names<-list()
data<-list()
counterb<-1
for (input in inp){
  data[[counterb]]<-readFasta(input)
  bin<-c(bin,rep(counter,length(data[[counterb]])))
  counter<-counter+1
  All<-c(All,data[[counterb]]@sread)
  names_all<-c(names_all,as.character(data[[counterb]]@id))
  names[[counterb]]<-as.character(data[[counterb]]@id)
  counterb<-counterb+1
}

fivemers<-oligonucleotideFrequency(All,5,step = 1)
fivemers<-fivemers+1
fivemers_clr<-cenLR(fivemers)
fivemers_clr_m <- unique(as.matrix(fivemers_clr$x.clr))
set.seed(42) # Set a seed if you want reproducible results

tsne_All_2d_U <- Rtsne(fivemers_clr_m,dims = 2, theta =0.7, perplexity =30) # Run TSNE
tsne_All_3d_U <-Rtsne(fivemers_clr_m,dims = 3, theta =0.7, perplexity =30) # Run TSNE

#load output of megan following multi-metagenome pipeline. c
finaltax_c<-findtax("path/tank1plusf.orfs.hmm.blast-ex.txt",6)
finaltax_s<-findtax("path/tank1plusf.orfs.hmm.blast-ex.txt",9)


#####################3d 
summar <- read.delim("path/bins.summary")
completenes<-summar$Completeness
completenes<-as.numeric(gsub("%","",as.character(completenes)))

centers3<-c()
for(i in 1:length(unique(bin)))
{
  
  group<-tsne_All_3d_U$Y[which(bin==i),]
  
  cent<-t(as.matrix(colMeans(group)))
  
  centers3<-rbind(centers3,cent)
  
}
indord<-order(completenes,decreasing = TRUE)
centers3_ord<-centers3[indord,]
completenes_ord<-completenes[indord]
finaltax_ord<-finaltax_c[indord]
finaltax_ords<-finaltax_s[indord]
n<-length(table(finaltax_c))
names_c<-names(table(finaltax_c))

col_new<-palette(c(rgb(216,197,243, maxColorValue=255),
                   rgb(142,184,240, maxColorValue=255),
                   rgb(244,177,137, maxColorValue=255),
                   rgb(90,208,218, maxColorValue=255),
                   rgb(242,162,186, maxColorValue=255),
                   rgb(165,227,173, maxColorValue=255),
                   rgb(201,181,247, maxColorValue=255),
                   rgb(132,192,140, maxColorValue=255),
                   rgb(219,171,224, maxColorValue=255),
                   rgb(209,239,179, maxColorValue=255),
                   rgb(173,181,240, maxColorValue=255),
                   rgb(173,189,126, maxColorValue=255),
                   rgb(231,165,204, maxColorValue=255),
                   rgb(121,221,203, maxColorValue=255),
                   rgb(233,167,155, maxColorValue=255),
                   rgb(124,204,238, maxColorValue=255),
                   rgb(239,227,174, maxColorValue=255),
                   rgb(216,197,243, maxColorValue=255),
                   rgb(140,197,165, maxColorValue=255),
                   rgb(228,186,217, maxColorValue=255),
                   rgb(176,240,210, maxColorValue=255)
                   # rgb(237,185,194, maxColorValue=255),
                   # rgb(169,236,233, maxColorValue=255),
                   # rgb(223,184,155, maxColorValue=255),
                   # rgb(184,207,242, maxColorValue=255),
                   # rgb(204,210,161, maxColorValue=255),
                   # rgb(168,176,214, maxColorValue=255),
                   # rgb(204,186,143, maxColorValue=255),
                   # rgb(139,198,194, maxColorValue=255),
                   # rgb(174,211,171, maxColorValue=255)
))
# col_new<-palette(c( rgb(216,197,243, maxColorValue=255),
#                     rgb(142,184,240, maxColorValue=255),
#                     rgb(244,177,137, maxColorValue=255),
#                     rgb(90,208,218, maxColorValue=255),
#                     rgb(242,162,186, maxColorValue=255),
#                     rgb(165,227,173, maxColorValue=255),
#                     rgb(201,181,247, maxColorValue=255)
                    
#))
centers<-centers3
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
t<-scatterplot3d(centers3_ord,pch=20,color=col_new[factor(finaltax_ord)],cex.symbols = 6,angle = 200,xlab = "x",ylab = "y", zlab = "z", main = "Binning overview - tSNE")
text(t$xyz.convert(centers3_ord), labels = 1:nrow(centers3_ord),
     cex= 0.9, col = "black",font=2)
legend("topleft",legend = unique(finaltax_ord),pch=20,bty = "n",col = col_new,pt.cex=3,cex=1.65,title = "Until class level")

colors<-rep("azure2",length(unique(bin)))
n <- length(which(completenes>70))
ramp <- colorRamp(c("#00007F", "blue", "#007FFF", "cyan",
                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
temp_col<-rgb( ramp(seq(0, 1, length = n)), max = 255)
########
centers<-centers3
colf<-col_new[factor(finaltax_ord)]
col_2<-palette(c(   rgb(165,227,173, maxColorValue=255),
                   rgb(90,208,218, maxColorValue=255),
                   rgb(142,184,241, maxColorValue=255),
                   rgb(135,206,250, maxColorValue=255),
                   rgb(229,204,255, maxColorValue=255),
                   rgb(65,105,225, maxColorValue=255),
                   rgb(50,205,50, maxColorValue=255),
                   rgb(244,177,137, maxColorValue=255),
                   rgb(239,227,174, maxColorValue=255),
                   rgb(216,197,243, maxColorValue=255),
                   rgb(140,197,165, maxColorValue=255),
                   rgb(228,186,217, maxColorValue=255)
                   
))
colf[1:12]<-col_2
colf[13:50]<-adjustcolor("gray100", alpha=0)
t<-scatterplot3d(centers3_ord,pch=20,color=colf,cex.symbols = 6,angle = 200,xlab = "x",ylab = "y", zlab = "z",main = "Binning best - tSNE")
text(t$xyz.convert(centers3_ord[1:12,]), labels = 1:nrow(centers3_ord[1:12,]),
     cex= 0.9, col = "black",font=2)
legend("topleft",legend = finaltax_ords[1:12],pch=20,bty = "n",col = colf[1:12],pt.cex=3,cex=1.5,title = "Until species level")


#function to assign taxonomy to bins 5 for class, 9 for species
findtax <-function(pathfile,tax_level){
  nom <- read.delim(pathfile, header=FALSE)
  tax<-as.character(nom[,1])
  tax_temp<-gsub('_([^_]*)$', "", tax)
  nomsplit<-strsplit(as.character(nom$V2),";")
  nenom<-c()
  counter<-1
  for(i in 1:length(nomsplit)){
    if(length(nomsplit[[i]])<=3){
      nomsplit[[i]]<-NA
      tax_temp[i]<-NA
    }else{
      nenom[counter]<-nomsplit[[i]][tax_level]
      counter<-counter+1
    }
    
  }
  tax_temp<-tax_temp[ !is.na( tax_temp ) ]
  tax<-cbind(tax_temp,nenom)
  
  
  abad_name_list<-list()
  abad_tax_list<-list()
  for (i in 1:length(names)){
    abad_id2<-match(names[[i]], as.character(tax[,1]))
    abad_tax_list[[i]]<-tax[abad_id2,]
  }
  
  num_unq<-c()
  final_tax<-c()
  for (i in 1:length(names)){
    temp<-table(abad_tax_list[[i]][,2])[table(abad_tax_list[[i]][,2])>0]
    num_unq<-c(num_unq,length(temp))
    if(length(temp)>0){
      final_tax<-c(final_tax,names(temp[which.max(temp)]))
    }else{
      final_tax<-c(final_tax,"unclassified")
    }
    
  }
  
  finaltax<-final_tax
}
