
library(tidyverse)

setwd("/Users/ahewett1/Documents")

clutch=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Clutch.csv")%>%
  select(Id, Season, Nestbox, MaleRing, FemaleRing, CrossFostered)

sequenced_ids=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)%>%
  rename(RingId=IID)%>%
  add_column(seq_round_one="yes")

  

mum_sequenced_ids=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)%>%
  rename(FemaleRing=IID)%>%
  add_column(mother_seq="yes")


dad_sequenced_ids=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)%>%
  rename(MaleRing=IID)%>%
  add_column(dad_seq="yes")


## number of clutches where both parents have been sequenced.
list_of_parents_seqd=clutch%>%
  left_join(mum_sequenced_ids, by = 'FemaleRing')%>%
  left_join(dad_sequenced_ids, by = 'MaleRing')%>%
  filter(dad_seq=='yes'&mother_seq=='yes')
#863


non_cf_parents_seqd=list_of_parents_seqd%>%
  filter(!CrossFostered==1|is.na(CrossFostered))

  
bird=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(-RingType, -Species, -RingAdministrator,-CH1903X,-CH1903Y,
         -CH1903Altitude,-RingedBy,-GeneticSex,-PhenotypeSex,-RFRingId,
         -RFRingNumber,-RFRingDate,-RFRingDateErrorMargin,-Remarks)


# check number of ids that have been sequenced in the clutches where both parents have been seq

chicks_seq=non_cf_parents_seqd%>%
  rename(BornClutchId=Id)%>%
  left_join(bird)%>%
  drop_na(RingId)%>%
  inner_join(sequenced_ids, by = 'RingId')


n_distinct(chicks_seq$RingId) ## 710 chicks in 223 clutches 
n_distinct(chicks_seq$BornClutchId)
n_distinct(chicks_seq$FemaleRing)
n_distinct(chicks_seq$MaleRing)




chicks_to_seq=non_cf_parents_seqd%>%
  rename(BornClutchId=Id)%>%
  left_join(bird)%>%
  drop_na(RingId)%>%
  anti_join(sequenced_ids, by = 'RingId', copy = T)

