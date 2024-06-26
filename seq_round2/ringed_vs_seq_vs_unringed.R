
library(tidyverse)


setwd("/Users/ahewett1/Documents")


bird=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(-RingType, -Species, -RingAdministrator,-CH1903X,-CH1903Y,
         -CH1903Altitude,-RingedBy,-GeneticSex,-PhenotypeSex,-RFRingId,
         -RFRingNumber,-RFRingDate,-RFRingDateErrorMargin,-Remarks)%>%
  unique()

sequenced_ids=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)%>%
  rename(RingId=IID)%>%
  add_column(seq_round_1="yes")

#6348 born after 2010 and before 2020
born_after_2010=bird%>%
  filter(as.Date(HatchDate, format = "%d/%m/%Y") > as.Date("01/01/2010", format = "%d/%m/%Y"))#%>%
  #filter(as.Date(HatchDate, format = "%d/%m/%Y") < as.Date("01/01/2020", format = "%d/%m/%Y"))


length(which(all_blood_after2010$RingId=="")) ## 601 unringed


lab=read.csv("Barn_owl_general/Data_base_dfs/Lab_20240523105823_LabIndividualTissue.csv")%>%
  select(IndividualIdentifier,OtherIdentifiers, SampleType,LabIndividualTissueId,TissueIdentifier,LabIndividualId)%>%
  rename(RingId=IndividualIdentifier)

lab_indiv=read.csv("Barn_owl_general/Data_base_dfs/Lab_20240523105823_LabIndividual.csv")%>%
  select(LabIndividualId,IndividualIdentifier,OtherIdentifiers)


per_yr=born_after_2010%>%
  mutate(hatch_date = as.Date(HatchDate, format = "%d/%m/%Y"))%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))%>%
  select(RingId, year)%>%
  left_join(sequenced_ids)%>%
  mutate(seq_round_1 = replace_na(seq_round_1, 'no'))



seq_per_yr=per_yr%>%
  group_by(year)%>%
  summarise(total=n(), 
            Ringed_and_seq=sum((seq_round_1 == 'yes'&RingId != ''), na.rm = TRUE), 
            Ringed_unseq=sum((seq_round_1 == 'no'&RingId != ''), na.rm = TRUE), 
            Unringed_unseq=sum(RingId == '', na.rm = TRUE))%>%
  mutate(all=Ringed_and_seq+Ringed_unseq+Unringed_unseq)

seq_per_yr$all == seq_per_yr$total

plotting=seq_per_yr%>%
  pivot_longer(cols=c(Ringed_and_seq,Ringed_unseq,Unringed_unseq),values_to = 'number_ids')%>%
  mutate(name=fct_relevel(name, 'Unringed_unseq', 'Ringed_unseq', 'Ringed_and_seq'))

plotting%>%
  ggplot(aes(x=year,y=number_ids,fill=name)) + 
    geom_bar(position = "stack", stat="identity")+
    theme_bw()+
  scale_fill_brewer(palette = 'Set2', direction = -1)




all_blood_after2010=born_after_2010%>%
  left_join(lab)%>%
  filter(!(SampleType%in%c("swab","Swab", "Feather", "",NA)))%>%
  distinct(RingId, .keep_all = T)%>% #
  filter(is.na(seq_round_1)) ## Removing ring ids that have already been sequenced
  

## ids with blood samples that have been ringed
proposed_ringed_ids=all_blood_after2010%>%
  mutate(hatch_date = as.Date(HatchDate, format = "%d/%m/%Y"))%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))


proposed_ids_per_yr=proposed_ringed_ids%>%
  group_by(year)%>%
  summarise(Ringed_unseq_with_blood=sum((RingId != ''), na.rm = TRUE))%>%
  full_join(seq_per_yr, by = 'year')%>%
  mutate(Ringed_unseq_DNA_not_good=Ringed_unseq-Ringed_unseq_with_blood)%>% ## getting number of ids with useless blood samples (like feathers etc)
  mutate(Ringed_unseq_with_blood = if_else(is.na(Ringed_unseq_with_blood),Ringed_unseq, Ringed_unseq_with_blood))





plotting_2=proposed_ids_per_yr%>%
  pivot_longer(cols=c(Ringed_and_seq, Ringed_unseq_with_blood, Ringed_unseq_DNA_not_good,Unringed_unseq),values_to = 'number_ids')%>%
  mutate(name=fct_relevel(name, 'Ringed_unseq_DNA_not_good','Unringed_unseq', 'Ringed_unseq_with_blood','Ringed_and_seq'))

  ggplot(plotting_2,aes(x=year,y=number_ids,fill=name)) + 
  geom_bar(position = "stack", stat="identity")+
  theme_bw()+
  scale_fill_brewer(palette = 'Set2', direction = -1)

proposed_ids_per_yr_b42019= proposed_ids_per_yr%>%
    filter(year<2021)

sum(proposed_ids_per_yr_b42019$Ringed_unseq_with_blood)## 1364
sum(proposed_ids_per_yr_b42019$Unringed_unseq) #484

## how many have a recorded death in proposed??
## (before 2021)
death_recorded=proposed_ringed_ids%>%
  filter(DeathDate!='')%>%
  filter(as.Date(HatchDate, format = "%d/%m/%Y") < as.Date("01/01/2021", format = "%d/%m/%Y"))
nrow(death_recorded) #n=239 

## and born after 2020

born_after_2020=bird%>%
  filter(as.Date(HatchDate, format = "%d/%m/%Y") > as.Date("01/01/2020", format = "%d/%m/%Y"))

