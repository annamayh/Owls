library(tidyverse)
setwd("/Users/ahewett1/Documents")


## Inline numerical values ###

## Results 3.1 ##

## LMM ##

## BILL##
id_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_Funi_unscaled.RDS")
id_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.1.ID_bill_GRM_FROH_unscaled.RDS")
## TARSUS ##
id_tarsus_funi=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.4.ID_tarsus_GRM_Funi.RDS")
id_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/1_unscaled_ID_wGRM/1.4.ID_tarsus_GRM_FROH.RDS")

summary(id_tarsus_froh)



get_perc_change_inbred_etc_LMM=function(f, model, fcoeff, age_function){

    standard=(fixef(model, pars = "Intercept")[,1]+
      fixef(model, pars = "rank")[,1]+
      age_function)#taken from age_acc varible 
    
    inbred=standard+(f*fixef(model, pars = paste0(fcoeff))[,1])
    
    perc_change = ((inbred-standard)/inbred) *100
    perc_change
    
    print(c(perc_change, inbred, standard, inbred-standard))

}


## remember bill needs to be converted to mm by * 0.1
get_perc_change_inbred_etc_LMM(0.25, id_bill_froh, "FHBD512gen", (184*exp(-0.99*0.932^30)))
get_perc_change_inbred_etc_LMM(0.25, id_bill_funi, "FuniWE", (184*exp(-0.99*0.932^30)))


get_perc_change_inbred_etc_LMM(0.25, id_tarsus_froh, "FHBD512gen",1) ## check and re-run without mean centering??









### NLM ###
## Bill and tarsus show sig diff growth rates ##

f=0.25

growth_bill_funi=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_totalfixedb.RDS")
#growth_bill_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.1.bill_gr_FROH_totalfixedb.RDS")

summary(growth_bill_funi)


asym=fixef(growth_bill_funi, pars = "asym_Intercept")[,1]
b=fixef(growth_bill_funi, pars = "b_Intercept")[,1]
c=fixef(growth_bill_funi, pars = "c_Intercept")[,1]

standard_age1= (asym*0.1)*exp(-b*(c)^0)  # (* 0.1 to get bill length in mm)
standard_age2= (asym*0.1)*exp(-b*(c)^30)

perc_change_growth_standard=((standard_age2-standard_age1)/standard_age1)*100

#asym_inbred=asym+(f*fixef(growth_bill_funi, pars = "asym_FuniWE")[,1])
c_inbred=c+(f*fixef(growth_bill_funi, pars = "c_FuniWE")[,1])

inbred_age1=(asym*0.1)*exp(-b*(c_inbred)^0)
inbred_age2=(asym*0.1)*exp(-b*(c_inbred)^30)

perc_change_growth_inbred=((inbred_age2-inbred_age1)/inbred_age1)*100

perc_change_growth_standard ## standard indiv grows 162% between ages 0 and 30 days 
perc_change_growth_inbred ## inbred indiv grows 156% between ages 0 and 30 days 



## tarsus 


growth_tarsus_funi=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr.RDS")
growth_tarsus_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.4_tarsus_gr_FROH_fixedb.RDS")


summary(growth_tarsus_froh)

asym1=fixef(growth_tarsus_froh, pars = "asym1_Intercept")[,1]
asym=fixef(growth_tarsus_froh, pars = "asym_Intercept")[,1]
b=fixef(growth_tarsus_froh, pars = "b_Intercept")[,1]
c=fixef(growth_tarsus_froh, pars = "c_Intercept")[,1]

standard_age1= asym1+(asym*exp(-b*(c)^0))
standard_age2= asym1+(asym*exp(-b*(c)^20))

perc_change_growth_standard=((standard_age2-standard_age1)/standard_age1)*100

c_inbred=c+(f*fixef(growth_tarsus_froh, pars = "c_FHBD512gen")[,1])

inbred_age1=asym1+(asym*exp(-b*(c_inbred)^0))
inbred_age2=asym1+(asym*exp(-b*(c_inbred)^20))

perc_change_growth_inbred=((inbred_age2-inbred_age1)/inbred_age1)*100

perc_change_growth_standard ## standard indiv grows 326% between ages 0 and 20 days 
perc_change_growth_inbred ## inbred indiv grows 310% between ages 0 and 30 days 



### mass ###

# sig asymptote for funi 


growth_mass_funi=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_fixedb.RDS")
#growth_mass_froh=readRDS("Inbreeding_depression_owls/Model_outputs/3_growth/3.2_mass_gr_FROH_totalfixedb.RDS")

summary(growth_mass_funi)

asym=fixef(growth_mass_funi, pars = "asym_Intercept")[,1]
b=fixef(growth_mass_funi, pars = "b_Intercept")[,1]
c=fixef(growth_mass_funi, pars = "c_Intercept")[,1]


asym_inbred=asym+(f*fixef(growth_mass_funi, pars = "asym_FuniWE")[,1])


standard_mass_final= asym*exp(-b*(c)^40)  # 
standard_mass_final
inbred_mass_final= asym_inbred*exp(-b*(c)^40)  # 
inbred_mass_final

inbred_mass_final-standard_mass_final ## 13g lighter 

perc_change_mass = ((inbred_mass_final-standard_mass_final)/inbred_mass_final) *100
perc_change_mass

####################
### Results 3.2 ###
####################

####################
### Results 3.3 ###
####################

## Gen Arch ##


get_estimates=function(input_files){
  
  gwas_output_list <- lapply(Sys.glob(paste0(input_files)), readRDS)
  
  unlisted_gwas=do.call(rbind,gwas_output_list)%>%
    as.data.frame()%>%
    mutate(super_scaffold = as.numeric(stringr::str_split(Window, '_', simplify = TRUE)[, 5]))%>%
    arrange(super_scaffold)%>%
    group_by(super_scaffold)%>%
    mutate(id = cur_group_id())
  
  bonf=0.05/nrow(unlisted_gwas)
  
  unlisted_gwas$sig <- ifelse(unlisted_gwas$p_val < bonf, "yes","no")      

  
  unlisted_gwas
  
}


bill_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.1_bill_brmsGWAS_FINAL/*.RDS.RDS"
mass_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.2_mass_brmsGWAS_FINAL/*.RDS.RDS"
tarsus_files="Inbreeding_depression_owls/Model_outputs/4_Gen_arch/4.4_tarsus_brmsGWAS_FINAL/*.RDS.RDS"


bill_gen_arch=get_estimates(bill_files)


ggplot(bill_gen_arch, aes(x=estimate)) +
  geom_density(fill=)+
  theme_classic()+
  geom_vline(xintercept = 0)






