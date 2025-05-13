.libPaths(c("/work/FAC/FBM/DEE/jgoudet/barn_owl/ahewett/R", .libPaths())) #specify library cluster path

library("RZooRoH")
library("tidyverse")


# Only calculating genomewide FROH using ss that are of good quality 
SuperScaffolds = read.table("./input_dfs/GoodQuality_ss_minus22.txt", sep = ',')[[1]]

#SuperScaffolds=SuperScaffolds[c(24, 28, 35, 33)]

list_Fs=list()
list_snps=list()


for(ss in SuperScaffolds){
  
  #removing loaded in files at the start of the loop so that R doesnt crash with two copies of things
  if (exisits("Mod")){
  rm(list = c("data_Rohs", "loc_mod", "loc_table", "Mod"))}
  
  print(paste0("Starting ", ss))
  load(paste0("/scratch/ahewett1/RzooROH_Rsess/EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_Super-Scaffold_",ss,".RData"))
  #Estimate F of 512 generations (1/2 T)
  list_Fs[[ss]] = cumhbd(zres=loc_mod, T = 1024)
  #get number of SNPs per scaffold
  list_snps[[ss]]=data_Rohs@nsnps
  
}

#unlist number of snps per ss
snpN=list_snps%>%
  unlist()

# get ring ids (note: this just uses the last loaded scaffold)
RingId = loc_mod@sampleids

# unlist and create dataframe of all FROH per ss and calculates genomewide FROH weighted by number of snps
ss_Fs_df=list_Fs%>%
  bind_cols() %>% # unlists all ss FROHs 
  rowwise() %>% # treating each line seperatly 
  mutate(FROH = weighted.mean(c_across(everything()), w = snpN)) %>% #calculating average weighted by number of snps
  ungroup()%>%
  add_column(RingId, .before = 1) # adding Ringid
  


write.table(ss_Fs_df,
            file="./outputs/FHBD_512gen_goodQ_ss_only.txt",
            row.names = F, col.names = FALSE, quote = F, 
            sep = ",",na = "NA")


