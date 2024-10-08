library(tidyverse)
library(WriteXLS)

setwd("/Users/ahewett1/Documents")
##########################
## priority 2 individuals ##
##########################

# including:

# 1. Confirmed dead ids that died within 90 days of hatching (n=224, ext = 220)
# 2. Missing (almost certainly) dead ids (n=119, ext = 118)
# 3. Un-ringed individuals with known death dates (n=73, ext = 55)
# 6. Siblings of previously sequenced ids that had high inbreeding coefficients (n=75, ext = 74)


clutch=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Clutch.csv")%>%
  select(Id, Season, Nestbox, MaleRing, FemaleRing, CrossFostered)

bird=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(-RingType, -Species, -RingAdministrator,-CH1903X,-CH1903Y,
         -CH1903Altitude,-RingedBy,-GeneticSex,-PhenotypeSex,-RFRingId,
         -RFRingNumber,-RFRingDate,-RFRingDateErrorMargin,-Remarks)

sequenced_ids=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)%>%
  rename(RingId=IID)%>%
  add_column(seq_round_1="yes")

big_boy_lab=read.csv("Barn_owl_general/Data_base_dfs/Ta_Switzerland_Tissue_DNA _BIG_FILE_11.01.24.csv")

samples_2021=read.csv("Barn_owl_general/Data_base_dfs/Ta_2021_and_more_pour_DB_rect_12.05.22.csv")%>%
  select(RingId, NB, Date.extr, emplacement.ADN, Boite.ADN, place.ADN, DNA.Tube.n.,conc.ng.ul)%>%
  filter(RingId!='')%>%
  drop_na(RingId)%>%
  filter(Date.extr!="-")
  
head(samples_2023)

samples_2022=read.csv("Barn_owl_general/Data_base_dfs/Ta_2022.csv")%>%
  select(Ring.Number.NB, NB.ring, Date.extr, DNA.Freezer.Rack.Box.Localisation, 
         DNA_Box_Name, DNA_Tube_.Place, 
         DNA_Tube_Nb, conc.ng.ul)%>%
  rename(RingId=Ring.Number.NB)

samples_2023=read.csv("Barn_owl_general/Data_base_dfs/Pour_Clara_Ta_2023.csv")#%>%
  select(RingId,extr.date,DNA.Freezer.Rack.Box.Localisation, DNA_Box_Name, 
         DNA_Tube_.Place, DNA.n., Hidex.C.ng.ul..26.01.24)%>%
  filter(RingId!='')%>%
  drop_na(RingId)%>%
  rename(conc.ng.ul=Hidex.C.ng.ul..26.01.24)

  

p1_ids=read.table("New_owls_to_seq/P1_Ids_to_seq_DNA_extracted.txt",
                   sep = ",",na = "NA", header = T)

## getting individuals with known death and birth dates 
dead_exact_after2021=bird%>%
  filter(!is.na(HatchDate) & !is.na(DeathDate) & DeathDate != "") %>%
  mutate(
    hatch_date = as.Date(HatchDate, format = "%d/%m/%Y"), #mutating dates to correct formatting
    death_date = as.Date(DeathDate, format = "%d/%m/%Y"),
    age_at_death = death_date - hatch_date)%>%
  filter(age_at_death<=90, #died before 90 days
         RingId != "", # has ring id
         hatch_date > as.Date("2021-01-01")) # hatched after / in 2010




p2=big_boy_lab%>%
  rename(RingId=Sample_Name_1)%>%
  right_join(dead_exact_after2021, relationship = "many-to-many")%>%
  select(RingId, BornClutchId, Date.extraction, DNA.Freezer.Rack.Box.Localisation, 
         DNA_Box_Name, DNA_Tube_.Place, DNA_Tube_Nb, DNA.C.ng.ul.)%>%
  filter(Date.extraction!=''&Date.extraction!='-'&DNA.Freezer.Rack.Box.Localisation!=' ')


n_distinct(p2$RingId) # 186 ids that died 


## and ~1 sibling per clutch that survived, to distinguish year effects and clutch effects.
sibs_2021=p2%>%
  select(BornClutchId)%>% 
  distinct(BornClutchId)%>% ## 100 different clutches
  left_join(bird)%>%
  filter(RingId!='')%>%
  anti_join(p2, by='RingId')%>%
  group_by(BornClutchId)%>%
  sample_n(1)%>%
  select(RingId)

sibs_2021_plus_lab=big_boy_lab%>%
  rename(RingId=Sample_Name_1)%>%
  right_join(sibs_2021, relationship = "many-to-many")%>%
  select(RingId, BornClutchId, Date.extraction, DNA.Freezer.Rack.Box.Localisation, 
         DNA_Box_Name, DNA_Tube_.Place, DNA_Tube_Nb, DNA.C.ng.ul.)%>%
  filter(Date.extraction!=''&Date.extraction!='-'&DNA.Freezer.Rack.Box.Localisation!="") # 90 siblings 


p1_ids_ids=p1_ids$RingId
#combine both dead ids and sibs to one df 

p2_all_dead_and_sibs=sibs_2021_plus_lab%>%
  rbind(p2)%>%
  group_by(RingId)%>%
  mutate(DNA.C.ng.ul.=as.numeric(DNA.C.ng.ul.))%>%
  slice(which.max(DNA.C.ng.ul.))%>%
  ungroup()%>%
  filter(DNA.C.ng.ul.>1)%>%
  select(-BornClutchId)%>%
  anti_join(p1_ids_ids)

n_distinct(p2_all_dead_and_sibs$RingId) # 270 ids (just under 3 plates)

final_p2=p2_all_dead_and_sibs[!p2_all_dead_and_sibs$RingId %in% p2_all_dead_and_sibs, ] #just a quick check they werent done before 
final_p2=p2_all_dead_and_sibs[!p2_all_dead_and_sibs$RingId %in% p1_ids_ids, ] #just a quick check they werent done before 



WriteXLS(p2_all_dead_and_sibs, "New_owls_to_seq/P2_dead_2021_onwards_and_sibs.xlsx")




########## 






