library(tidyverse)


bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_brmsGWAS_FINAL/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_brmsGWAS_FINAL/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_brmsGWAS_FINAL/*.RDS.RDS"



gwas_output_list <- lapply(Sys.glob(paste0(bill_files)), readRDS)


unlisted_gwas=do.call(rbind,gwas_output_list)%>%
  as.data.frame()%>%
  mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
  arrange(super_scaffold)%>%
  group_by(super_scaffold)%>%
  mutate(id = cur_group_id())

bonf=0.05/nrow(unlisted_gwas)

unlisted_gwas$sig <- ifelse(unlisted_gwas$p_val < bonf, "yes", "no") 

mean(unlisted_gwas$estimate)

blw_zero=length(which(sign(unlisted_gwas$estimate)%in%'-1'))## number of times estimate is below 0
abv_zero=length(which(sign(unlisted_gwas$estimate)%in%'1'))## number of times estimate is above 0

blw_zero/(blw_zero+abv_zero)


data_tarsus=read.table("Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",sep=",", header=T)
## % change = old - new / old * 100

normal=mean(data_tarsus$LeftTarsus)
inb=mean(data_tarsus$LeftTarsus)-20

((normal-inb)/normal)*100






data_bill=read.table("Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",sep=",", header=T)
## % change = old - new / old * 100

normal=mean(data_bill$BillLength)
inb=mean(data_bill$BillLength)-4.3

((normal-inb)/normal)*100


