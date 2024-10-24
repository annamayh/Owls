library(RZooRoH,lib = "/users/ahewett1/R")

args <- commandArgs(trailingOnly = TRUE)

ss=args[1]

load(paste0("./elavanc1/ID_in3K/data/4_ROHs/EntireRsession_All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_GenPOSplus10_Model13HBDclasses_ss_",ss,".RData"))


snp_window_size=2500
ids_HBD_chr=list()
number_ids=data_Rohs@nind


for (id in 1:number_ids){
  
      id_ibc_seg_list <- list()
      seed=0
      for (window in seq(from = 1, to = which.max(data_Rohs@bp) , by = snp_window_size)){ # from first snp position to last snp position in chr 
          
          seed=seed+1
          overlap=snp_window_size-1 ## no overalap 
          y1 <- probhbd(zres = loc_mod, ## created from zoorun
                        zooin = data_Rohs, ## zoo data
                        id = id, ## specific id
                        chrom = 1, ## as it has been run per chromosome, chr number always =1
                        startPos = data_Rohs@bp[window], ## snp bp at start of window 
                        endPos = data_Rohs@bp[window+overlap]#, ## by snp window size
                        #T = 40  # 20 generations back but max is 16
                        )

          seg_ibc=mean(y1) # get mean pr HBD for window
          id_ibc_seg_list[[seed]] <-  seg_ibc ## mean pr of being HBD for specific section 

      }
      
      ids_HBD_chr[[id]] <- do.call(cbind, id_ibc_seg_list) ## adding ids in list
      
      if (id %% 500 == 0) {
        print(paste0("Finished ID ", id))
      }
}


m_HBD_segs=do.call(rbind, ids_HBD_chr) ## 

row.names(m_HBD_segs)=data_Rohs@sample_ids
colnames(m_HBD_segs) <- c(paste0("pr_HBD_", ss , "_wind_", 1:seed))


saveRDS(m_HBD_segs,file=paste0("./ahewett/ID_owls/outputs/4_gen_arch/HBD_per_window_2500_16gens/16_gens_HBD_perID_",snp_window_size,"-wind_",ss,".RDS"))

