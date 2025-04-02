library(tidyverse)
library(data.table)

setwd("/Users/ahewett1/Documents")

hbd_segs_list <- lapply(Sys.glob("Inbreeding_depression_owls/ROH_regions/HBD_segs_out/HBD_perCHR*.RDS"), readRDS)

hbd_segs_matrix <- map(hbd_segs_list, ~ as.data.frame(.x))%>% # convert to tibble 
  map(~ rownames_to_column(.x, var="RingId"))%>% # make column of ring id
  reduce(full_join, by="RingId")%>% # reducing list into df with cols of ids
  select_if(~ !any(is.na(.))) # remove the columns with Nas at the end of chromosomes


bill_df_read=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)%>%
  mutate(age_acc=184*exp(-0.99*0.932^age_days)) ## account for age using gompertz growth
 
#######################################################################################################
######        now only for chr that have good coverage #######
##########################################################################

good_sc=read.table("Barn_owl_general/ss_goodQ.txt")
scaffold_names=good_sc$V1

scaffolds_to_keep <- paste0(scaffold_names, collapse = "|")

## linking HBD segs to phenotype df and creating a matrix 
hbd_segs_matrix_rm=hbd_segs_matrix%>%
  select(RingId, matches(scaffolds_to_keep)) ## removes about 300 windows 

names(hbd_segs_matrix_rm) <- gsub(x = names(hbd_segs_matrix_rm), pattern = "-", replacement = "_")  
names(hbd_segs_matrix_rm)[2:ncol(hbd_segs_matrix_rm)] <- substr(names(hbd_segs_matrix_rm)[2:ncol(hbd_segs_matrix_rm)], 4, 
                                                                nchar(names(hbd_segs_matrix_rm)[2:ncol(hbd_segs_matrix_rm)]) - 6)

## joining the hbd matrix to the bill df as we have mutliple records per id
froh_mat=hbd_segs_matrix_rm%>%
    right_join(bill_df_read, by = "RingId") #%>% # to account for the repeated records
    #select(matches("pr_hbd")) %>%## removing ring id to convert to df (BUT BE CAREFUL IT IS EXACTLY THE SAME AS DF)
    #data.matrix()

#froh_mat[1:3, 1:3]
## calculating sum and sum div of hbd segs
sum_hbd_segs=hbd_segs_matrix_rm%>%
    mutate(hbd_sum = rowSums(.[2:ncol(hbd_segs_matrix_rm)], na.rm = TRUE))%>% ## remove nas from last col of each chr
    mutate(hbd_sum_div=hbd_sum/(ncol(hbd_segs_matrix_rm)-1))%>%
    select(RingId, hbd_sum, hbd_sum_div)
  
bill_df=froh_mat%>%
    left_join(sum_hbd_segs)%>%
    as.data.frame()

nrow(bill_df)
ncol(bill_df)
  
#HBD segs as matrix for multi-memebership
nrow(froh_mat) # check froh matrix has same number of rows as original data
ncol(froh_mat)# check has expected number of columns
#bill_df$froh_mat_linked=froh_mat


any(is.na(froh_mat))

saveRDS(bill_df, file="Inbreeding_depression_owls/ROH_regions/df_with_HBD_mat/CHR_bill_df_OnlyGoodSS.RDS")
