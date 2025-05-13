library(janitor)
library(dplyr)
library(tibble)

setwd("/Users/ahewett1/Documents")

## script to split chromosome inbreeding coeffs into chunks so can run a more efficient array job
output_path="Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_2500snpwind_chunks/"

hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_per_window_2500/HBD_perID_2500*Super-Scaffold_*.RDS"), readRDS)

for (i in 1:length(hbd_segs_list)){
  
hbd_segs_list_chr <- hbd_segs_list[i]%>%
  as.data.frame()%>%
  rownames_to_column(var='RingId')%>%
  select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes

chr_num=colnames(hbd_segs_list_chr)[2]%>%
  str_sub(start=23,end = -8)
  

split_df_by_cols_save <- function(df, chunk_size) {
  n_cols <- ncol(df)  # Number of columns in the dataframe
  col_indices <- split(seq(2, n_cols), ceiling(seq(2, n_cols) / chunk_size))  # Split column indices, excluding first column
  
  # Loop through each chunk, adding the first column and saving as .RDS
  for (i in seq_along(col_indices)) {
    chunk <- df %>% 
      select(1, all_of(col_indices[[i]]))  # Include the first column
    saveRDS(chunk, file = paste0(output_path,"prHBD_ss",chr_num,"_chunk_", i, ".RDS") )   # Save chunk to .RDS
  }
  
}

chunk_size = 30  # Number of columns per chunk
split_df_by_cols_save(hbd_segs_list_chr, chunk_size)

}
