## prepping dfs for first-step ##


library(tidyverse)

setwd("/Users/ahewett1/Documents")


##Fgrm and FROH from elu
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)



egg_info=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Egg.csv", header = T)%>%
  select(BornClutchId, RaisedClutchId, WidthMM, LengthMM, WeightGrams,CrossFosteredClutchId, ObservationDate)%>%
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA), too_few = "align_start")%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  separate_wider_delim(obs_date, delim = "-", names = c("Season",NA,NA))

## labelling eggs 
eggs=egg_info%>%
  group_by(BornClutchId)%>%
  mutate(egg_Id=1:n())%>%
  ungroup()

clutch_info=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Clutch.csv", header = T)%>%
  select(Id, MaleRing, FemaleRing, Season)
head(clutch_info)

## 3 main dfs to work with:
# inbreeding coeffs 
# egg info
# clutch info to link eggs with parents

write.table(eggs,
            file = "teaching/FS_ID_eggs/Egg_phenotypes.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(inb,
            file = "teaching/FS_ID_eggs/Inbreeding_coefficients.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(clutch_info,
            file = "teaching/FS_ID_eggs/Clutch_info.txt",
            row.names = F, quote = F, sep = ",",na = "NA")







n_distinct(eggs$BornClutchId)
