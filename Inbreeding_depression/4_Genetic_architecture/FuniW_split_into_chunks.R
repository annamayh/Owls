library(janitor)
library(dplyr)
library(tibble)
library(gaston)
library(stringr)
setwd("/Users/ahewett1/Documents")

## script to split chromosome inbreeding coeffs into chunks so can run a more efficient array job
output_path="Inbreeding_depression_owls/ROH_regions/Funi_segs_out/"
funi_segs <- readRDS("Inbreeding_depression_owls/ROH_regions/Funi_segs_out/Funiw_3kowls_perScaffold_2500SNPs_windows.RDS")

good_sc=read.table("Barn_owl_general/ss_goodQ.txt")
scaffold_names=good_sc$V1%>%
  str_sub(start=7,end = -1)

filtered_funi_segs <- funi_segs[grepl(paste(scaffold_names, collapse = "|"), names(funi_segs))]

names_vcf=read.vcf("Inbreeding_depression_owls/ROH_regions/Funi_segs_out/All3085_NewNamesCORRECTED_AUTOSAUMES_RP502SNPs_ss_Super-Scaffold_32_GenPOSplus10.vcf.gz")
RingId=names_vcf@ped$id


for (i in 1:length(filtered_funi_segs)){
  
  Funi_segs_list_chr <- filtered_funi_segs[i]%>%
    as.data.frame()%>%
    add_column(RingId, .before = 1)%>%
    select_if(~ !any(is.na(.)))# remove the columns with Nas at the end of chromosomes
  #change col names to simpler thing
  chr_num=colnames(Funi_segs_list_chr)[2]%>%
    str_sub(start=66,end = -23)
  

  colnames(Funi_segs_list_chr)[2:ncol(Funi_segs_list_chr)] <- str_extract(colnames(Funi_segs_list_chr)[2:ncol(Funi_segs_list_chr)], "Super\\.Scaffold\\_\\d+") %>%
    paste0("_chunk", 2:(ncol(Funi_segs_list_chr))-1)  # Append "chunk" and the column number
  
  split_df_by_cols_save <- function(df, chunk_size) {
    n_cols <- ncol(df)  # Number of columns in the dataframe
    col_indices <- split(seq(2, n_cols), ceiling(seq(2, n_cols) / chunk_size))  # Split column indices, excluding first column
    
    # Loop through each chunk, adding the first column and saving as .RDS
    for (i in seq_along(col_indices)) {
      chunk <- df %>% 
        select(1, all_of(col_indices[[i]]))  # Include the first column
      saveRDS(chunk, file = paste0(output_path,"localFuniW_ss",chr_num,"_chunk_", i, ".RDS") )   # Save chunk to .RDS
    }
    
  }
  
  chunk_size = 30  # Number of columns per chunk
  split_df_by_cols_save(Funi_segs_list_chr, chunk_size)
  
}
