
library(RZooRoH,lib = "/users/ahewett1/R")

args <- commandArgs(trailingOnly = TRUE)
ss=args[1]

load(paste0("./elavanc1/ID_in3K/data/4_ROHs/EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_",ss,".RData"))



ids_HBD_chr=list()
number_ids=data_Rohs@nind


for (id in 1:number_ids){
  
          y1 <- probhbd(zres = loc_mod, ## created from zoorun
                        zooin = data_Rohs, ## zoo data
                        id = id, ## specific id
                        chrom = 1, ## as it has been run per chromosome, chr number always =1
                        #startPos = data_Rohs@bp[window], ## snp bp at start of window 
                       # endPos = data_Rohs@bp[window+overlap]
                        ) ## by snp window size
          ## T = base pop for HBD segs .. defualt all HBD segs 

          ids_HBD_chr[[id]]=mean(y1) # get mean pr HBD for window
      

      if (id %% 500 == 0) {
        print(paste0("Finished ID ", id))
      }
}


m_HBD_segs=do.call(rbind, ids_HBD_chr) ## 

row.names(m_HBD_segs)=data_Rohs@sample_ids
colnames(m_HBD_segs) <- c(paste0("pr_HBD_", ss , "_wind_"))


saveRDS(m_HBD_segs,file=paste0("./ahewett/ID_owls/outputs/4_gen_arch/HBD_per_CHR/HBD_perCHR_perID_ss",ss,".RDS"))

