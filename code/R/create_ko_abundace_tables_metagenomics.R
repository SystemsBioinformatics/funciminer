
#functions
norm_lacto<-function(table_tax){
  ind_l<-which(colnames(table_tax)=="Lactobacillus")
  table_tax<-t(apply(table_tax, MARGIN = 1, function(x) x/sum(x)))
  table_tax_nol<-t(apply(table_tax, MARGIN = 1, function(x)  x[-ind_l]/(1-x["Lactobacillus"])))
  table<-cbind(table_tax_nol,table_tax[,"Lactobacillus"])
  colnames(table)[ncol(table)]<-"Lactobacillus"
  table
}

norm_lacto_ko<-function(table_tax){
  ind_l<-grep("_Lacto",colnames(table_tax))
  table_tax<-t(apply(table_tax, MARGIN = 1, function(x) x/sum(x)))
  table_tax_nol<-t(apply(table_tax, MARGIN = 1, function(x)  x[-ind_l]/(1-mean(x[ind_l][x[ind_l]!=0]))))
  table<-cbind(table_tax_nol,table_tax[,ind_l])
  table
}

names_ordered<-c("B5_0" ,    "B5_1"  ,   "B5_3"   ,  "B5_4"  ,   "B5_5"   ,  "B5_6"   ,  "B5_7" , "Bott_B5"  ,  "B7_0"    , "B7_1"  ,  "B7_2"  ,  "B7_3" ,    "B7_4"  ,  
                 "B7_5_"  ,  "B7_6"   ,  "B7_7",  "Bott_B7"   ,  "B8_1"   ,  "B8_2"  ,   "B8_3"  ,   "B8_4"   ,  "B8_5"   ,  "B8_6"  ,   "B8_7"   ,  
                 "T1_0" ,    "T1_1" ,  "T1_2"  , "T1_3"  ,  "T1_4"    ,   "T1_5"   ,  "T1_6"  ,   "T1_7"  , "T1_8",   "Bott_T1" ,       "T3_0"   , 
                 "T3_1"  ,   "T3_2"   ,  "T3_3"    , "T3_5"  ,   "T3_6"  ,   "T3_8"  ,"Bott_T3"  , "T4_0"   ,  "T4_1"   ,  "T4_2"   ,  "T4_3"   ,  "T4_4"  ,   "T4_5"    ,
                 "T4_6"   ,  "T4_7"   ,  "T4_8", "9_0","A9_0","A9_1","A9_2",
                 "A9_3","A9_4","A9_5","A9_6"  ,"Bott_A9","A11_0","A11_1","A11_2","A11_3",
                 "A11_4","A11_5","A11_6","Bott_A11" ,"A12_0" ,"A12_1","A12_2","A12_3","A12_4","A12_5","A12_6")  



isEmpty <- function(x) {
  return(length(x)==0)
}
access_kegg_ko<-function(ko_s){
  
  if(ncol(ko_s)>10){
    temp_all<-list()
    names<-colnames(ko_s)
    names<-names[-index]
    for (j in 1:length(names)){
      print(j)
      temp_all[[j]]<-keggGet(names[j])
    }
    print("enter")
    temp<-sapply(temp_all, function(x) c(x[[1]]$ENTRY,x[[1]]$NAME,if(isEmpty(x[[1]]$DEFINITION)){"NA"}else{paste(x[[1]]$DEFINITION,collapse = " / ")},if(isEmpty(x[[1]]$MODULE)){"NA"}else{paste(x[[1]]$MODULE,collapse = " / ")},if(isEmpty(x[[1]]$PATHWAY)){"NA"}else{paste(x[[1]]$PATHWAY,collapse = " / ")},if(isEmpty(x[[1]]$BRITE)){"NA"}else{paste(x[[1]]$BRITE,collapse = " / ")}),simplify = F)
  }else{
    print("enter")
    temp_all<-keggGet(colnames(ko_s))
    temp<-sapply(temp_all, function(x) c(x$ENTRY,x$NAME,if(isEmpty(x$DEFINITION)){"NA"}else{paste(x$DEFINITION,collapse = " / ")},if(isEmpty(x$MODULE)){"NA"}else{paste(x$MODULE,collapse = " / ")},if(isEmpty(x[[1]]$PATHWAY)){"NA"}else{paste(x[[1]]$PATHWAY,collapse = " / ")},if(isEmpty(x$BRITE)){"NA"}else{paste(x$BRITE,collapse = " / ")}),simplify = F)
    
    
  }
  temp.table<-t(data.frame(temp))
  rownames( temp.table)<-temp.table[,1]
  colnames( temp.table)<-c("ENTRY","NAME","DEFINITION","MODULE","PATHWAY","BRITE")
  temp.table
}

access_kegg_module<-function(module_s){
  if(ncol(module_s)>10){
    temp_all<-list()
    for (j in 1:ncol(module_s)){
      temp_all[[j]]<-keggGet(colnames(module_s)[j])
    }
    temp_a<-sapply(temp_all, function(x) c(x[[1]]$ENTRY,x[[1]]$NAME,if(isEmpty(x[[1]]$DEFINITION)){"NA"}else{paste(x[[1]]$DEFINITION,collapse = " / ")},if(isEmpty(x[[1]]$CLASS)){"NA"}else{paste(x[[1]]$CLASS,collapse = " / ")},if(isEmpty(x[[1]]$PATHWAY)){"NA"}else{paste(x[[1]]$PATHWAY,collapse = " / ")}),simplify = F)
  }else{
    temp_all<-keggGet(colnames(module_s))
    temp_a<-sapply(temp_all, function(x) c(x$ENTRY,x$NAME,x$DEFINITION,if(isEmpty(x$CLASS)){"NA"}else{paste(x$CLASS,collapse = " ")},if(isEmpty(x$PATHWAY)){"NA"}else{paste(x$PATHWAY,collapse = " ")}),simplify = F)
  }
  
  lp_modules.table<-t(data.frame(temp_a))
  rownames( lp_modules.table)<-lp_modules.table[,1]
  colnames( lp_modules.table)<-c("ENTRY","NAME","DEFINITION","CLASS","PATHWAY")
}


clean_wine_wnames<-function(x){
  x[grep("*_Bobal_Finished_Inoc",x)]<-"Bott_B5"
  x[grep("*_Bobal_Finished_CTRL",x)]<-"Bott_B7"
  x[grep("*_Airen_Finished_Inoc",x)]<-"Bott_A9"
  x[grep("*_Airen_Finished_CTRL",x)]<-"Bott_A11"
  x[grep("*_Tempranillo_Finished_Inoc",x)]<-"Bott_T1"
  x[grep("*_Tempranillo_Finished_CTRL",x)]<-"Bott_T3"
  
  x<-gsub(pattern ="THLR_" ,replacement = "",x = x)
  x<-gsub(pattern ="X7PVA_" ,replacement = "",x = x)
  x<-gsub(pattern ="7PVA_" ,replacement = "",x = x)
  x<-gsub(pattern ="VECX_" ,replacement = "",x = x)
  x<-gsub(pattern ="*_KO.txt$" ,replacement = "",x = x) 
  x<-gsub(pattern ="*_sort.csv$" ,replacement = "",x = x)
  x<-gsub(pattern ="Sample_TOG_" ,replacement = "",x = x)
  x<-gsub(pattern ="*_R1_val_1.fq$" ,replacement = "",x = x)
  x<-gsub(pattern ="*.top$" ,replacement = "",x = x)
  x<-gsub(pattern ="*_KO$" ,replacement = "",x = x)
  x<-gsub(pattern ="*.csv$" ,replacement = "",x = x)
  x<-gsub(pattern ="*.faa$" ,replacement = "",x = x)
  x<-gsub(pattern ="*.fa$" ,replacement = "",x = x)
  
  x
}

lb_p_image<-function(x,name){
  col_names<-c(rep("firebrick1",6),rep("firebrick",13),rep("khaki",21))
  par( mar = par( "mar" ) + c( 4, 4, 0, 8 ) )
  image( x, xaxt= "n", yaxt= "n" ,main=name,col=c("white","black"))
  axis( 2, at=seq(0,1,length.out=ncol( x ) ), labels= colnames( x ), las= 2 ,font=2,cex.axis=1.2)
  vecb<-c("Inoc.",rep("",7),"Control1",rep("",8),"Control2",rep("",6))
  axis( 1, at=seq(0,1,length.out=nrow( x ) )[1:24], labels=vecb, las= 2,lwd.ticks=10,col.ticks ="red",cex.axis=1)
  vect<-c("Inoc.",rep("",9),"Control1",rep("",7),"Control2",rep("",8))
  axis( 1, at=seq(0,1,length.out=nrow( x ) )[25:51], labels=vect, las= 2,lwd.ticks=10,col.ticks ="red3",cex.axis=1)
  veca<-c("Inoc.",rep("",8),"Control1",rep("",7),"Control2",rep("",6))
  axis( 1, at=seq(0,1,length.out=nrow( x ) )[52:75], labels=veca, las= 2,lwd.ticks=10,col.ticks ="khaki",cex.axis=1)
  par(xpd=TRUE)
  
  legend(x=1.05,y=1.05,col = c("white","black"),pch = 20,legend = c("Absence","Presence"),bg="gainsboro",cex=1.4,pt.cex = 2)
  legend(x=1.05,y=0.6,col = c("red","red3","khaki"),pch = 20,legend = c("Bobal","Tempranillo","Airen"),cex=1.4,pt.cex = 2)
  
}

#path to the files with read count per varity for all samples
totalread_persample_ord<-function(){
  bobal_wine_QQcon_f <- read.csv("~/path/bobal_wine_QQcon_f.csv", sep="")
  bcounts<-bobal_wine_QQcon_f[,8]
  bnames_counts<-bobal_wine_QQcon_f[,3]
  
  Tempranillo_wine_QQconf <- read.csv("~/path/Tempranillo_wine_QQconf.csv", sep="")
  tcounts<-Tempranillo_wine_QQconf[,8]
  tnames_counts<-Tempranillo_wine_QQconf[,3]
  
  Airen_wine_QQcon_f <- read.csv("~/path/Airen_wine_QQcon_f.csv", sep="")
  acounts<-Airen_wine_QQcon_f[,8]
  anames_counts<-Airen_wine_QQcon_f[,3]
  #combine total
  all_count<-c(bcounts,tcounts,acounts)
  all_names_counts<-c(as.character(bnames_counts),as.character(tnames_counts),as.character(anames_counts))
  #clean
  all_names_counts<-clean_wine_wnames(all_names_counts)
  all_count_ord<-(all_count[match(names_ordered,all_names_counts)])
  
  all_count_ord
}

#path to shotgun sample mapped to L.plantarum genome
lb_p_abund_norm<-function(){
  
  output<-"/path/lp_iso/"
  folders<-c("7PVA_Bobal","THLR_Tempranillo","VECX_Airen")
  map_list_iso<-list()
  map_names_iso<-c()
  counter<-1
  for(j in folders){
    setwd(paste(output,j,sep="/"))
    flist <- list.files(pattern = "*.csv$")
    
    for(i in flist){
      map_names_iso[[counter]]<-i
      map_list_iso[[counter]]<-sum(read.delim(paste(getwd(),i,sep="/"), header=FALSE)[,2])
      
      counter<-counter+1
      
    }
  }
  #normilize in respect to L.plantarum abundance, total read count per samples.
  all_count_ord<-totalread_persample_ord()
  all_count_ord<-all_count_ord/2
  map_names_iso<-clean_wine_wnames(map_names_iso)
  map_list_iso_ord<-(unlist(map_list_iso)[match(names_ordered,map_names_iso)])
  y_norm<-round((unlist(map_list_iso_ord)/all_count_ord)*1000000)
  y_norm
}


######################################################################################################


#create table of KO from species with abundance


###########################################################################################
library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("data.table", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
library("plyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#path to shotgun KO annotations
path<-"/path/"
make_ko_table<-function(path){
  folders<-c("7PVA_Bobal_maxbin_input","THLR_Tempranillo_maxbin_input","VECX_Airen_maxbin_input")
  list_top<-list()
  count<-1
  names_list_top<-c()
  for(j in folders){
    setwd(paste(path,j,sep=""))
    flist <- list.files(pattern = "*fa_folder$")
    for(i in flist){
      setwd(paste(path,j,i,sep="/"))
      ffile <- list.files(pattern = "*_KO.txt$")
      ab <-read.table(ffile,sep="\t", header=FALSE, fill = TRUE)
      # ab <- as.character(ab[grep("K[[:digit:]]{5}",ab),])
      if(ncol(ab)==2){
        list_top[[count]]<-levels(ab$V2)
      }else{
        list_top[[count]]<-grep("K[[:digit:]]{5}",as.matrix(ab),value=TRUE)
      }
      names_list_top<-c(names_list_top,ffile)
      print(count)
      count<-count+1
      
    }
  }
  data<-list(list_top,names_list_top)
}
#path to shotgun abundancne mapping to contigs/scaffolds
make_abundko_table<-function(path, type){
  folders<-c("7PVA_Bobal_maxbin_input","THLR_Tempranillo_maxbin_input","VECX_Airen_maxbin_input")
  list_top<-list()
  list_top1<-list()
  count<-1
  names_list_top<-c()
  KO_all<-c()
  for(j in folders){
    setwd(paste(path,j,sep=""))
    flist <- list.files(pattern = "*fa_folder$")
    for(i in flist){
      setwd(paste(path,j,i,sep="/"))
      ffile <- list.files(pattern = "*.top$")
      top <- read.table(ffile,sep="\t", header=FALSE, fill = TRUE)
      colnames(top)<-c("scaffold","KO","tax2","tax3","genus","geneid","score")
      top[,1]<-sub("user:","",top[,1])
      ffile <- list.files(pattern = "*.csv$")
      ab <- read.delim(ffile, header=FALSE)
      colnames(ab)<-c("scaffold","abundance")
      
      join<-full_join(top,ab, by="scaffold")
      joinf<-filter(join,KO !="")
      
      if (type==1){
        joinf[,2]<-paste(joinf[,2],joinf[,5],sep = "_")
        joinf_r<-joinf[,c(2,8,7)]
        joinf_rr<-filter(joinf_r,abundance !=0)
        joinf_rrr<-ddply(joinf_rr, "KO", numcolwise(sum))
        list_top[[count]]<-joinf_rrr
        KO_all<-c(KO_all,as.character(joinf_rrr[,1]))
        names_list_top<-c(names_list_top,ffile)
        print(count)
        count<-count+1
      }else if (type==2){
        joinf_r_withtax<-joinf[,c(2,8)]
        joinf_rr_withtax<-filter(joinf_r_withtax,abundance !=0)
        joinf_rrr_withtax<-ddply(joinf_rr_withtax, "KO", numcolwise(sum))
        list_top[[count]]<-joinf_rrr_withtax
        KO_all<-c(KO_all,as.character(joinf_rrr_withtax[,1]))
        names_list_top<-c(names_list_top,ffile)
        print(count)
        count<-count+1
        
      }else if (type==3){
        joinf_tax<-joinf[,c(5,8,7)]
        joinf_taxr<-filter(joinf_tax,abundance !=0)
        joinf_taxr<-filter(joinf_taxr,score > 100)
        joinf_taxrr<-ddply(joinf_taxr, "genus", numcolwise(sum))
        list_top[[count]]<-joinf_taxrr
        KO_all<-c(KO_all,as.character(joinf_taxrr[,1]))
        names_list_top<-c(names_list_top,ffile)
        print(count)
        count<-count+1
      }
      
      
      
      
    }
  }
  data<-list(list_top,KO_all,names_list_top)
}
###################################################
# create matrix of KO, KO + taxonomy, KO + abundance +taxonomy, taxonomy + abundance
##################################################
data_ko_tax<-make_abundko_table(path,type=1)
# data_ko<-make_abundko_table(path,type=2)
# data_tax<-make_abundko_table(path,type=3)
# data_tax1<-make_abundko_table(path,type=3)
# data_ko_binary<-make_ko_table(path)

list_top<-data_ko_tax[[1]]
KO_all<-data_ko_tax[[2]]
names_list_top<-data_ko_tax[[3]]

###################score filtering 
names_list_top<-data_ko_binary[[2]]
list_top<-data_ko_binary[[1]]
KO_all_u<-unique(unlist(list_top))
KO_all_u<-KO_all_u[-1]

###
KO_all_u<-unique(KO_all)
finaltable<-c()
for (i in 1:length(names_list_top)){
  temp<-rep(0,length(KO_all_u))
  ind<-match(KO_all_u, list_top[[i]],nomatch = 0)
  indx<-which(ind!=0)
  temp[indx]<-1
  finaltable<-rbind(finaltable,temp)
}
rownames(finaltable)<-clean_wine_wnames(names_list_top)
colnames(finaltable)<-KO_all_u
finaltable_ord<-(finaltable[match(names_ordered,rownames(finaltable)),])

#########abund

####create table with abundance for the threee cases, ko_tax, ko, tax . ordered
create_table<-function(list){
  KO_all<-list[[2]]
  list_top<-list[[1]]
  names_list_top<-list[[3]]
  
  KO_all_u<-unique(KO_all)
  finaltable<-c()
  for (i in 1:length(names_list_top)){
    temp<-rep(0,length(KO_all_u))
    ind<-match(KO_all_u, list_top[[i]]$genus,nomatch = 0)
    indx<-which(ind!=0)
    temp[indx]<-list_top[[i]]$abundance[ind[ind!=0]]
    finaltable<-rbind(finaltable,temp)
  }
  rownames(finaltable)<-clean_wine_wnames(names_list_top)
  colnames(finaltable)<-KO_all_u
  finaltable_ord<-(finaltable[match(names_ordered,rownames(finaltable)),])
  
  total_reads<-totalread_persample_ord()/2
  finaltable_ord_norm<-finaltable_ord/total_reads *min(total_reads)
  finaltable_ord_norm
}

table_ko_tax<-create_table(data_ko_tax)
table_ko<-create_table(data_ko)

table_tax<-create_table(data_tax)
table_tax1<-create_table(data_tax1)
#remove if all are 0
ind_rem<-apply(table_tax,MARGIN = 2, function(x) all(x==0))
table_ko_tax_red<-table_tax[,!ind_rem]
#remove if appear only once
ind_nozero<-apply(table_ko_tax_red,MARGIN = 2, function(x) length(which(x!=0)))
finaltable<-table_ko_tax_red[,-which(ind_nozero==1)] 
finaltable_ord_norm_reduces_removed<-table_ko_tax_red[,which(ind_nozero==1)] 

finaltable_n<-apply(finaltable, 1, function(x) x/sum(x))
write.table(finaltable_n,"/path/tax_of_kos_abund_snorm.csv", sep="\t", col.names = TRUE, row.names = TRUE)
finaltable_norm<-norm_lacto(finaltable)
write.table(finaltable_norm,"/path/tax_of_kos_abund_lbpnorm.csv", sep="\t", col.names = TRUE, row.names = TRUE)

write.table(finaltable,"/path/tax_of_kos_abund.csv", sep="\t", col.names = TRUE, row.names = TRUE)
test<-read.table("/path/ko_tax_abund.csv")

