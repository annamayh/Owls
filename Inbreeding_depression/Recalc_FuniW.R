.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library(hierfstat)
library(gaston)
library(parallel)
library(tidyverse)


## function to get Funi weighted
get.funiwn<-function(dos){
  
  if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
  p<-dos@p
  het<-2*p*(1-p)
  res<-apply(gaston::as.matrix(dos),1,function(x) {
    nas<-which(is.na(x)); 
    if(length(nas)>0){
      xs<-x[-nas];
      ps<-p[-nas];
      hets<-het[-nas];
      sum(xs^2-(1+2*ps)*xs+2*ps^2)/sum(hets);
    }
    else sum(x^2-(1+2*p)*x+2*p^2)/sum(het);
  }
  )
  return(list(Funi=unlist(res),het=sum(het)))
}


SuperScaffolds = read.table("./input_dfs/GoodQuality_ss_minus22.txt", sep = ',')[[1]]

#create variables to populate
Funiw.t<-rep(0,3085)
shet<-0
fun_per_ss=list()


for(ss in SuperScaffolds){
  
  print(paste0("Starting super scaffold ", ss))
  
  vcf_file_name=paste0("/scratch/ahewett1/3k_VCFs/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_Super-Scaffold_",ss,"_GenPOSplus10.vcf.gz")
 
  bed=read.vcf(vcf_file_name, convert.chr = FALSE)
  fun_per_ss[[ss]]=get.funiwn(bed)#cacl Funi per chr
  
  Funiw.t<-Funiw.t+fun_per_ss[[ss]]$Funi*fun_per_ss[[ss]]$het#cacl total funi cpmbined with all chr 
  shet<-shet+fun_per_ss[[ss]]$het
  
}

Funiw.t<-Funiw.t/shet

Funi_recalc=as.data.frame(Funiw.t)%>%
  rownames_to_column(var="RingId")


write.table(Funi_recalc,
            file="./outputs/FuniW_goodQ_ss_only.txt",
            row.names = F, col.names = FALSE, quote = F, 
            sep = ",",na = "NA")

