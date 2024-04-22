load("sequioa/pedigree_FULL_recontr_out26.2.RData")


ped_corrected2=read.table("sequioa/pedigreeCORRECTED.tab", header = T)%>%
  mutate(dadid=case_when(
    momid=="M026267" ~"M026267", 
    dadid=="M026658" ~ NA,
    dadid=="M031195" ~ "M041195",
    TRUE ~ dadid
  ))%>%
  mutate(momid=case_when( ## found during first parentage assigment
    dadid=="M026658" ~ "M026658",## M026267   recorded as female but actually a male 
    momid=="M027000" ~ "M022696",       ## this id was ringed twice 
    
    momid=="M026267" ~ NA,       ## M026658 recorded as male but genetic sex is female
    TRUE ~ momid 
    
  ))%>%
  mutate(id=case_when(
    id=="M011867" ~ "889913",
    TRUE ~ id))
    




sequenced_ids_read=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)

sequenced_ids=as.character(sequenced_ids_read[,1])

full_ped_recontr=seq_reconstr_out$Pedigree

full_compare_parentage=PedCompare(
  Ped1 = ped_corrected2,
  Ped2 = full_ped_recontr,
  DumPrefix = c("XF", "XM"),
  SNPd = sequenced_ids
  
)


consensus_ped=full_compare_parentage$ConsensusPed

corrected_consensus_ped=consensus_ped%>%
  select(id, dam.c, sire.c)%>%
  rename(RingId=id, mumID=dam.c, dadID=sire.c)

mis=full_compare_parentage$Mismatch


owl_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, HatchDate)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA), too_few = "align_start")%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))%>%
  select(-hatch_date)%>%
  unique()

  owl_bird%>%
    filter(RingId=="M026076")

 
  corrected_consensus_ped%>%
    filter(RingId=="M026161")
  
  dum_match=full_compare_parentage$DummyMatch
  
  dum_match%>%
    filter(id.1=="M026903")
  
  ped_corrected2%>%
    filter(id=="M026903")
  
  
    
  # ID  MOTHER  FATHER cohort Mum.cohort Dad.cohort
  # 7025 M032164 M032156 M032167     22         22         NA
  # 7027 M032166 M032156 M032167     22         22         NA
  # 7373 M032168 M032156 M032167     22         22         NA
  # 7513 M032163 M032156 M032167     22         22         NA
  # 7515 M032165 M032156 M032167     22         22         NA
  # 730  M026076 M011734 M026903     17         12         17
  # 738  M026161 M011734 M026903     17         12         17
  
  
  
  
  
  ## hartch date errors
  owl_hatch=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, HatchDate, RingDate, GrowthStageWhenRinged)%>%
    separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA), too_few = "align_start")%>%
    mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
    separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("birth_year",NA,NA))%>%
    
    separate_wider_delim(RingDate, delim = " ", names = c("ring_date",NA), too_few = "align_start")%>%
    mutate(ring_date=as.Date(ring_date, format = "%d/%m/%Y"))%>%
    separate_wider_delim(ring_date, delim = "-", cols_remove=F,names = c("ring_year",NA,NA))%>%
    mutate(ring_year=as.numeric(ring_year),birth_year=as.numeric(birth_year) )%>%
    mutate(discrepancy=ring_year-birth_year)


  
  
  
ped_birth_yr=ped_corrected2%>%
    rename(RingId=id)%>%
    left_join(owl_bird, relationship = "many-to-many")%>%
    mutate(year=case_when(
      RingId=="M032156" ~ "2014", 
      TRUE ~ year
      ))


owl_hatch%>%
  filter(RingId=="M032156")
