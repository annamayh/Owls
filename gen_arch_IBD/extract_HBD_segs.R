library(RZooRoH)

setwd("/Users/ahewett1/Documents")



load("Inbreeding_depression_owls/ROH_regions/EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_Super-Scaffold_12.RData")

loc_mod@hbdp[[3]][,1:15]

data_Rohs@bp[1] ## bp of snps .. 1st SNP is 0.5Mb in to the chr?
data_Rohs@bp[56547] ## bp of snps last snp is at ~25Mb 


id_ibc=list()


    y1 <- probhbd(zres = loc_mod, ## created from zoorun
                  zooin = data_Rohs, ## zoo data
                  id = 1, 
                  chrom = 1, ## as it has been run per chromosome, chr number always =1
                  startPos = data_Rohs@bp[1], ## default is 0
                  endPos = data_Rohs@bp[10])
    
    mean(y1) ## mean pr of being HBD for specific section 
