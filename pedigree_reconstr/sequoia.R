library(tidyverse)
library(sequoia)

setwd("/Users/ahewett1/Documents/")
## read in filtered 

input=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  column_to_rownames(var="IID")

GenoM=as.matrix(input[6:460])
geno_check=CheckGeno(GenoM, DumPrefix = c("XF", "XM"))

## life history data from database
life=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(RingId, HatchDate, DeathDate, PhenotypeSex, RingDate)%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("birth_date",NA),too_few = "align_end")%>%
  separate_wider_delim(DeathDate, delim = " ", names = c("death_date",NA),too_few = "align_end")%>%
  separate_wider_delim(RingDate, delim = " ", names = c("ring_date",NA),too_few = "align_end")%>%
  mutate(birth_date=as.Date(birth_date, format = "%d/%m/%Y"),
         death_date=as.Date(death_date, format = "%d/%m/%Y"),
         ring_date=as.Date(ring_date, format = "%d/%m/%Y"))%>%
  separate_wider_delim(birth_date, delim = "-", cols_remove=T,names = c("BirthYear",NA,NA))%>%
  separate_wider_delim(death_date, delim = "-", cols_remove=T,names = c("DeathYear",NA,NA))%>%
  separate_wider_delim(ring_date, delim = "-", cols_remove=T,names = c("RingYear",NA,NA))%>%
  drop_na(RingId)%>%
  filter(!RingId=="")%>%
  unique()

life$phen_sex=if_else(life$PhenotypeSex=="Male", 1, 
                      if_else(life$PhenotypeSex=="Female",2, 
                              NA))

sex=read.table("Barn_owl_general/GeneticSex_3K_RP548.txt", header=T)

## some mismatched sexes found from first parentage assignment 
sexes=life%>%
  left_join(sex)

mismatched_sexes <- sexes %>%
  filter(!is.na(phen_sex) & !is.na(GeneticSex) & phen_sex != GeneticSex)
print(mismatched_sexes)


all_life_hist=life%>%
  left_join(sex, by='RingId')%>%
  mutate(BirthYear=as.integer(BirthYear))%>%
  mutate(sex_combined=case_when(
    GeneticSex==1 ~ 2, ## in sequoia female = 1 and male = 2
    GeneticSex==2 ~ 1, 
    (is.na(GeneticSex)&phen_sex==1) ~ 2,
    (is.na(GeneticSex)&phen_sex==2) ~ 1,
    (GeneticSex==1&phen_sex==2) ~ 2,
    (GeneticSex==2&phen_sex==1) ~ 1))%>%
  mutate(year_last=case_when(
    !is.na(DeathYear) ~ DeathYear))%>%
  add_column(min=NA)%>%
  as.data.frame()

all_life_hist=all_life_hist[,c(1,8,2,10,5,9)]

all_life_hist$BirthYear=as.integer(all_life_hist$BirthYear)
all_life_hist$min=as.integer(all_life_hist$min)
all_life_hist$RingYear=as.integer(all_life_hist$RingYear)
all_life_hist$year_last=as.integer(all_life_hist$year_last)

head(all_life_hist)


################################################
### RUN PARENTAGE ASSIGNMENT ##############################
################################################
## first just the parentage to check any mismatches
seq_parentage_out=sequoia(GenoM, 
                  DummyPrefix = c("XF", "XM"),
                  LifeHistData=all_life_hist,
                  Module='par', 
                  args.AP=list(MaxAgeParent=15)
                  )

#################################################
parentage_ped=seq_parentage_out$Pedigree

## corrected pedigree found during first parentage assigment
ped_corrected=read.table("sequioa/pedigreeCORRECTED.tab", header = T)%>%
  mutate(dadid=case_when(
    momid=="M026267" ~"M026267", 
    momid=="M038112" ~ "M038112",
    dadid=="M026658" ~ NA,
    dadid=="M031195" ~ "M041195",
    TRUE ~ dadid
  ))%>%
  mutate(momid=case_when( ## found during first parentage assigment
    dadid=="M026658" ~ "M026658",## M026267 , M038112  recorded as female but actually a male 
    momid=="M026267" ~ NA,       ## M026658 recorded as male but genetic sex is female
    momid=="M038112" ~ NA,
    TRUE ~ momid 
    
  ))
           
          
ped_corrected=ped_corrected[,c(1,3,2)]




sequenced_ids_read=read.table("sequioa/3.inputfile_for_sequoia.raw", header=T)%>%
  select(IID)
sequenced_ids=as.character(sequenced_ids_read[,1])


compare_parentage=PedCompare(
  Ped1 = ped_corrected,
  Ped2 = parentage_ped,
  DumPrefix = c("XF", "XM"),
  SNPd = sequenced_ids
  
)


mismatch_parents=compare_parentage$Mismatch
print(mismatch_parents)
# > print(mismatch_parents)
# id   dam.1  sire.1   dam.2  sire.2    id.r dam.r sire.r id.dam.cat id.sire.cat dam.class sire.class
# 332  M010989 M027000 M010971 M022696    <NA> M010989    NA     NA         GG          GD  Mismatch     P1only
# 410  M011867 M011634 M011635  871649  870404 M011867    NA     NA         GG          GG  Mismatch   Mismatch
# 613  M026190 M027000 M010971 M022696    <NA> M026190    NA     NA         GG          GD  Mismatch     P1only
# 622  M026231 M027000 M010971 M022696    <NA> M026231    NA     NA         GG          GD  Mismatch     P1only
# 623  M026232 M027000 M010971 M022696    <NA> M026232    NA     NA         GG          GD  Mismatch     P1only
# 2377 M040592 M028977 M043668    <NA> M040587 M040592    NA     NA         GG          GG    P1only   Mismatch
# 2406 M040633 M031605 M031195 M031605 M041195 M040633    NA     NA         GG          GG     Match   Mismatch
# 2407 M040634 M031605 M031195 M031605 M041195 M040634    NA     NA         GG          GG     Match   Mismatch
# 2408 M040636 M031605 M031195    <NA> M041195 M040636    NA     NA         GG          GG    P1only   Mismatch
# 2905 M043834 M039430 M027966  881868    <NA> M043834    NA     NA         GG          GG  Mismatch     P1only




# CalcOHLLR(Pedigree=   ped_corrected,
#           GenoM=GenoM,
#           CalcLLR = FALSE,
#           LifeHistData=all_life_hist,
#           DumPrefix = c("XF", "XM")
#           )




########################################################################################
##### RUN FULL PEDIGREE RECONSTRUCTION ################################################
#######################################################################################


## try with vector of error rate 

error=c(0.01, 0.01, 0.01)

seq_reconstr_out=sequoia(GenoM, 
                          DummyPrefix = c("XF", "XM"),
                          LifeHistData=all_life_hist, 
                          Module = 'ped',
                          Err = error, 
                         args.AP=list(MaxAgeParent=15)
                         
                         
)

## created 220 dummy ids

#save(seq_reconstr_out, seq_parentage_out, file = "sequioa/pedigree_FULL_recontr_out26.2.RData")
load("sequioa/pedigree_FULL_recontr_out26.2.RData")

full_ped_recontr=seq_reconstr_out$Pedigree

compare_pedigrees=PedCompare(
  Ped1 = ped_corrected,
  Ped2 = full_ped_recontr,
  DumPrefix = c("XF", "XM"),
  SNPd = sequenced_ids
  
)

## pedigree diagnostics 




full_mismatch_parents=compare_pedigrees$Mismatch
print(full_mismatch_parents)


dummys=compare_pedigrees$DummyMatch%>%
  filter(id.1=="nomatch")

##127 dummy ids created with no matching id\

consensus_pedf=compare_pedigrees$ConsensusPed%>%
  filter(grepl("XF", dam.c))

consensus_pedm=compare_pedigrees$ConsensusPed%>%
  filter(grepl("XM", sire.c))

consensus_pedi=compare_pedigrees$ConsensusPed%>%
  filter(grepl("XM|XF", id))


##365 parental gaps filled


life2=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(RingId, HatchDate, BornClutchId, HatchClutchId, RaisedClutchId, RingDate, SiteId)

##mother of this id was created
life2%>%
  filter(RingId=="M026522")

df_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Clutch.csv")%>%
  select(Id, Season, Nestbox, MaleRing, FemaleRing, CrossFostered)

df_bird%>%
  filter(Id=="1308")
