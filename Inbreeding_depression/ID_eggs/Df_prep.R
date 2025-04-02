## prepping dfs for first-step ##


library(tidyverse)
library(lme4)
library(brms)
library(WriteXLS)

setwd("/Users/ahewett1/Documents")


##Fgrm and FROH from elu
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)



egg_info=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Egg.csv", header = T)%>%
  select(BornClutchId, RaisedClutchId, Rank, WidthMM, LengthMM, WeightGrams, ObservationDate)%>%
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA), too_few = "align_start")%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  separate_wider_delim(obs_date, delim = "-", names = c("Season","Month",NA))

hist(egg_info$WidthMM, breaks = 50) ## checking dist of the egg widths etc 
hist(egg_info$LengthMM, breaks = 50)


## labelling eggs 
eggs=egg_info%>%
  group_by(BornClutchId)%>%
  drop_na(Rank, WidthMM, LengthMM)%>%
  distinct(Rank, .keep_all = T)%>% #only keeping one record per egg (or rank)
  ungroup()%>%
  filter(WidthMM>100&LengthMM>100)%>% #remove eggs less that 100mm 
  filter(LengthMM>WidthMM)# remove eggs where lenght is longer than width (probably recorded wrong)


egg_errors=egg_info%>%
  filter(WidthMM<100|LengthMM<100|LengthMM<WidthMM)

# WriteXLS::WriteXLS(egg_errors,
#             ExcelFileName = "teaching/FS_ID_eggs/Egg_errors.xls")
# 


clutch_info=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Clutch.csv", header = T)%>%
  select(Id, MaleRing, FemaleRing, Season)%>%
  mutate(Season=as.factor(Season))
head(clutch_info)


owl_birb=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, GrowthStageWhenRinged, HatchDate, RingDate)%>%
  mutate(HatchingDate=as.Date(HatchDate, format = "%d/%m/%Y"))%>%
  separate_wider_delim(HatchingDate, delim = "-", names = c("MotherBirthYear",NA,NA))%>%
  mutate(MotherBirthYear=as.numeric(MotherBirthYear))

check=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv", header = T)

mother_age_clutch=clutch_info%>%
  mutate(Season=as.numeric(as.character(Season)))%>%
  filter(FemaleRing!="")%>%
  left_join(owl_birb, by = c("FemaleRing"="RingId"))%>%
  mutate(mother_age=Season-MotherBirthYear)%>%
  mutate(naive_mother=case_when(
    mother_age==1 ~ "Naive",
    mother_age>1 ~ "Adult", 
    is.na(mother_age) ~ "Unknown"
  ))%>%
  select(Id, MaleRing, FemaleRing, Season, naive_mother)%>%
  mutate(Season=as.factor(Season))



## 3 main dfs to work with:
# inbreeding coeffs 
# egg info
# clutch info to link eggs with parents

# write.table(eggs,
#             file = "teaching/FS_ID_eggs/Egg_phenotypes.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 
# write.table(inb,
#             file = "teaching/FS_ID_eggs/Inbreeding_coefficients.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 
# write.table(clutch_info,
#             file = "teaching/FS_ID_eggs/Clutch_info.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")


clutch_full_parents=left_join(eggs, mother_age_clutch, by = c("BornClutchId"="Id", 'Season'))


mother_InbDep=clutch_full_parents%>%
  right_join(inb, by = c('FemaleRing'='RingId'))%>%
  rename(mother_FHBD512gen= FHBD512gen, mother_FuniWE= FuniWE)%>%
  mutate(Rank=as.numeric(Rank), 
         Month=as.numeric(Month))

mean(mother_InbDep$LengthMM, na.rm = T)


father_InbDep=clutch_full_parents%>%
  right_join(inb, by = c('MaleRing'='RingId'))%>%
  rename(father_FHBD512gen = FHBD512gen, father_FuniWE = FuniWE)
  

prior_length=c(prior(student_t(3, 394,17), class = "Intercept"), ## 
               prior(student_t(3,0,17), class = "sd"),
               prior(student_t(3,0,17), class = "sigma"))

get_prior(LengthMM ~ 1 + mother_FuniWE + Rank + (1 | Season) + (1 | BornClutchId) + (1 | Month), 
          data = mother_InbDep)
  
mod_length=brm(LengthMM ~ 1 + mother_FuniWE + Rank + naive_mother + (1 | Season) + (1 | BornClutchId) + (1 | Month), 
    data = mother_InbDep, 
    cores = 4, 
    chains = 4, 
    prior = prior_length,
    control=list(adapt_delta=0.95),
    iter = 10000)

## yes for funi, no for froh

summary(mod_length)
plot(mod_length)


prior_wid=c(prior(student_t(3, 306,10), class = "Intercept"), ## 
            prior(student_t(3,0,10), class = "sd"),
            prior(student_t(3,0,10), class = "sigma"))

mod_wid=brm(WidthMM ~ 1 + mother_FuniWE + Rank + (1 | Season) + (1 | BornClutchId) + (1 | Month), 
               data = mother_InbDep, 
               cores = 4, 
               chains = 4, 
               prior = prior_wid,
               control=list(adapt_delta=0.95),
               iter = 10000)

summary(mod_wid)
plot(mod_wid)




n_distinct(eggs$BornClutchId)
