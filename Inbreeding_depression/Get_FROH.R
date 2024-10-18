##  this is a modified script from taken from Eleonore Lavanchy to calculate FROH using the package RZooRoH with different generations back  

library(RZooRoH,lib = "/users/ahewett1/R")

#read sample IDs Libnames from Elu's files
sampleIDs = as.vector(read.table("./elavanc1/ID_in3K/data/inputFILES/All3085_NewNamesCORRECTED.list")[,1])

#read ss list
SuperScaffolds = as.vector(read.table("./elavanc1/ID_in3K/data/inputFILES/SS_AUTOSAUMES_3K.list")[,1])

#### F ####
#Create DF with just sample IDs that we'll complete with each F interval
F_HBD_genint = as.data.frame(sampleIDs)


#Read the 53 ss R ENV files and extract F
for(ss in SuperScaffolds){
  
  print(paste0("Starting SS ", ss))
  
  #Load the big R file from Elu per scaffold
  load(list.files(path="./elavanc1/ID_in3K/data/4_ROHs", pattern=paste0("EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_",ss,".RData"), full.names = T))
  
  #Estimate FROH 
  # T = double the number of generations
  Fgenint = cumhbd(zres=loc_mod, T = 32)
  
  #merge with sample IDs
  Fgen = as.data.frame(cbind(loc_mod@sampleids,Fgenint))
  colnames(Fgen)[1] = "sampleIDs"
  
  #merge the full dataframe
  F_HBD_genint = merge(F_HBD_genint, Fgen, by = "sampleIDs")
  #Set new column name
  colnames(F_HBD_genint)[ncol(F_HBD_genint)] = paste0("FHBD_ss_",ss)
  
}

#Read the file with nb of SNPs per interval (for F weigthed mean)
SNPsNB = as.vector(read.table("./ahewett/ID_owls/scripts/4_gen_arch/Get_FROH/SNPsnb_per_ss.list", h = F)[,2])

#pass FHBD columns of F_HBD_genint  numeric
for(col in 2:ncol(F_HBD_genint)){F_HBD_genint[,col] = as.numeric(as.character(F_HBD_genint[,col]))}

#Create new vectors with weighted means (per intervals)
F_HBD = as.data.frame(cbind(INDVs=F_HBD_genint[,1], FHBDgen=as.numeric(apply(F_HBD_genint[,2:ncol(F_HBD_genint)],1 ,function(x){mean(x, weigths = SNPsNB)}))))

#pass FHBD as numeric
F_HBD[,2] = as.numeric(as.character(F_HBD[,2]))

colnames(F_HBD) = c("RingId","FHBD_16gens")


#save the FHBD "final" table
write.table(F_HBD, "./ahewett/ID_owls/outputs/0_FROH_ests/3k_FHBD_16gens.txt", sep = "\t", quote = F, col.names = T, row.names = F)