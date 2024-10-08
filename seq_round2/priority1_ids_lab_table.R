library(tidyverse)
library(WriteXLS)
##########################
## priority 1 individuals ##
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


#################
### ~~ 1 ~~ #####
#################

## getting individuals with known death and birth dates 
age_dead_exact=bird%>%
  filter(!is.na(HatchDate) & !is.na(DeathDate) & DeathDate != "") %>%
  mutate(
    hatch_date = as.Date(HatchDate, format = "%d/%m/%Y"), #mutating dates to correct formatting
    death_date = as.Date(DeathDate, format = "%d/%m/%Y"),
    age_at_death = death_date - hatch_date)

dead_wRing_knownDd_blood=age_dead_exact%>%
  filter(age_at_death<=90, #died before 90 days
         RingId != "", # has ring id
         hatch_date > as.Date("2010-01-01"))%>% # hatched after / in 2010
  anti_join(sequenced_ids, by = "RingId")%>% # remove ids sequenced last time
  left_join(lab,  by = "RingId")%>% #joing with lab info
  filter(SampleType%in%c("blood","Blood cell"))%>% #check which ones have blood samples
  distinct(RingId, .keep_all = TRUE) #

################
### ~~ 2 ~~ ### 
###############
missings=read.csv("Barn_owl_general/missing_ids/Confirmed_missingIDs_with_blood_and_rings.csv")%>%
  anti_join(sequenced_ids, by = "RingId")%>% # remove ids previously sequenced (n=1)
  anti_join(dead_wRing_knownDd_blood, by = 'RingId')

cfs_missings=read.csv("Barn_owl_general/missing_ids/Confirmed_cross_fost_missingIDs_with_blood.csv")%>%
  anti_join(sequenced_ids, by = "RingId")%>% # remove ids previously sequenced (n=1)
  anti_join(dead_wRing_knownDd_blood, by = 'RingId')

total_missings=rbind(cfs_missings,missings)%>%
  unique()

#############
## ~~ 6 ~~ ##
#############

inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)
#top 90% inbred individuals from 1st round of sequencing
top_ibcs=inb%>%
  filter(FuniWE > quantile(FuniWE, 0.9))%>% # top % ibcs
  left_join(bird)%>% ## joining with bird info to get DOB
  filter(as.Date(HatchDate, format = "%d/%m/%Y") > as.Date("2010-01-01")) # hatched after / in 2010

#n_distinct(top_ibcs$BornClutchId) ##  number of different clutches that have highest ibcs 

#inbred clutches
inbred_clutches=top_ibcs%>%
  select(BornClutchId)%>%
  distinct(BornClutchId)%>%
  na.omit()
# inbred sibling within those inbred clutches
inbred_sibs=bird%>%
  right_join(inbred_clutches)%>% #
  anti_join(sequenced_ids, by = "RingId")%>% # remove ids previously sequenced
  anti_join(dead_wRing_knownDd_blood, by = 'RingId')%>% # remove ids from stage 1.1
  anti_join(total_missings, by = 'RingId')%>% # remove ids from stage 1.2
  left_join(lab,  by = "RingId")%>% # join with lab info
  filter(SampleType%in%c("blood","Blood cell")) #filter for only ids with blood samples




######################################################################
#### linking to lab data ####

missing=total_missings$RingId
dead_ids=dead_wRing_knownDd_blood$RingId
inbred_ids=inbred_sibs$RingId

all_Ringids_plate_1=c(missing, dead_ids, inbred_ids)%>%
  unique()%>%
  as.data.frame()%>%
  rename(RingId='.')


DNA_extr_1=big_boy_lab%>%
  select(Sample_Name_1, Tissue.Type, Tissue.Freezer.Rack.Box.Localisation, Tissue_Container_Name, 
         Tissue_Container_Place, Tissue.Number, 
         Date.extraction, DNA.Freezer.Rack.Box.Localisation, 
         DNA_Box_Name, DNA_Tube_.Place, DNA_Tube_Nb, DNA.C.ng.ul.)%>%
  rename(RingId=Sample_Name_1)%>%
  filter(RingId!='')%>%
  right_join(all_Ringids_plate_1, by = 'RingId')%>%
  filter(Tissue.Type%in%c("blood","Blood", "Blood cell", "Blood & CORT","blood spot" ))%>%
  filter(Date.extraction!=''| is.na(Date.extraction))

n_distinct(DNA_extr_1$RingId)

write.table(DNA_extr_1,
            file = "New_owls_to_seq/P1_Ids_to_seq_DNA_extracted.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

WriteXLS(DNA_extr_1, "New_owls_to_seq/P1_Ids_to_seq_DNA_extracted.xlsx")





################################################################################


###############
### ~~ 3 ~~ ###
###############
All_after2010_noRing=bird%>%
  filter(as.Date(HatchDate, format = "%d/%m/%Y") > as.Date("2010-01-01"))%>% # hatched after / in 2010
  filter(as.Date(HatchDate, format = "%d/%m/%Y") < as.Date("2021-01-01"))%>% #hatched before 2020 (to keep with 3k)
  filter(RingId=="")

clutch_born=clutch%>%
  rename(BornClutchId=Id)


dead_noRing=All_after2010_noRing%>%
  filter(!(DeathDate==""))%>%
  left_join(clutch_born, by = 'BornClutchId')%>%
  select(BirdId, Mark, BornClutchId, Season, Nestbox.y, CrossFostered, HatchDate, DeathDate)

NBs=big_boy_lab%>%
  filter(str_starts(Sample_Name_1, 'NB')| str_starts(Sample_Name_2, 'NB'))%>%
  filter(!str_starts(Sample_Name_2, 'M') &  !str_starts(Sample_Name_3, 'M') & !str_starts(Sample_Name_1,'M'))%>%
  filter(Tissue.Type%in%c("blood","Blood", "Blood cell", "Blood & CORT","blood spot" ))%>%
  
# First sample same seq: NB, Season, N, nestbox, Mark
# and match to sample name 2
dead_noRing_knownDd_NB_S_N_nest_M=dead_noRing%>%
  drop_na(BornClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_2 = paste0('NB', Season, 'N', nestbox_edit, Mark))

sample_name2_match_NB_S_N_nest_M=dead_noRing_knownDd_NB_S_N_nest_M%>%
  inner_join(NBs, by = 'Sample_Name_2', relationship = "many-to-many")

# NB, nestbox, Mark
# and match to sample name 1 & 2
dead_noRing_knownDd_NB_nest_mark=dead_noRing%>%
  drop_na(BornClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_2 = paste0('NB', nestbox_edit,Mark))

sample_name2_match_NB_nest_mark=dead_noRing_knownDd_NB_nest_mark%>%
  inner_join(NBs, by = 'Sample_Name_2',relationship = "many-to-many")


## sample name 1
dead_noRing_knownDd_NB_nest_mark=dead_noRing%>%
  drop_na(BornClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_1 = paste0('NB', nestbox_edit,Mark))

sample_name1_match_NB_nest_mark=dead_noRing_knownDd_NB_nest_mark%>%
  inner_join(NBs, by = 'Sample_Name_1',relationship = "many-to-many")


# N, nestbox, Mark
# and match to sample name 1
dead_noRing_knownDd_N_nest_M=dead_noRing%>%
  # filter(CrossFostered==0 | is.na(CrossFostered))%>%
  drop_na(BornClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_1 = paste0('N', nestbox_edit, Mark))

sample_name1_match_N_nest_M=dead_noRing_knownDd_N_nest_M%>%
  inner_join(NBs, by = 'Sample_Name_1',relationship = "many-to-many")


all_unringed_linked_blood=rbind(sample_name2_match_NB_S_N_nest_M,sample_name1_match_NB_nest_mark)%>%
  rbind(sample_name2_match_NB_nest_mark)%>%
  rbind(sample_name1_match_N_nest_M)%>%
  unique()%>%
  rename(clutch=BornClutchId)

##########################################################################
## and CFd ids ##
clutch_raised=clutch%>%
  rename(RaisedClutchId=Id)

dead_noRing_CF=All_after2010_noRing%>%
  filter(!(DeathDate==""))%>%
  left_join(clutch_raised, by = 'RaisedClutchId')%>%
  filter(CrossFostered==1)%>%
  select(BirdId, Mark, RaisedClutchId, Season, Nestbox.y, CrossFostered, HatchDate, DeathDate)

# First sample same seq: NB, Season, N, nestbox, Mark
# and match to sample name 2
CF_dead_noRing_knownDd_NB_S_N_nest_M=dead_noRing_CF%>%
  drop_na(RaisedClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_2 = paste0('NB', Season, 'N', nestbox_edit, Mark))

CF_sample_name2_match_NB_S_N_nest_M=CF_dead_noRing_knownDd_NB_S_N_nest_M%>%
  inner_join(NBs, by = 'Sample_Name_2', relationship = "many-to-many")

# NB, nestbox, Mark
# and match to sample name 1 & 2
CF_dead_noRing_knownDd_NB_nest_mark=dead_noRing_CF%>%
  drop_na(RaisedClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_2 = paste0('NB', nestbox_edit,Mark))

CF_sample_name2_match_NB_nest_mark=CF_dead_noRing_knownDd_NB_nest_mark%>%
  inner_join(NBs, by = 'Sample_Name_2',relationship = "many-to-many")


## sample name 1
CF_dead_noRing_knownDd_NB_nest_mark=dead_noRing_CF%>%
  drop_na(RaisedClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_1 = paste0('NB', nestbox_edit,Mark))

CF_sample_name1_match_NB_nest_mark=CF_dead_noRing_knownDd_NB_nest_mark%>%
  inner_join(NBs, by = 'Sample_Name_1',relationship = "many-to-many")


# N, nestbox, Mark
# and match to sample name 1
CF_dead_noRing_knownDd_N_nest_M=dead_noRing_CF%>%
  # filter(CrossFostered==0 | is.na(CrossFostered))%>%
  drop_na(RaisedClutchId)%>% #remove 2 ids with no clutch info
  mutate(nestbox_edit=gsub("[A-Za-z]", "", Nestbox.y))%>%
  mutate(Sample_Name_1 = paste0('N', nestbox_edit, Mark))

CF_sample_name1_match_N_nest_M=CF_dead_noRing_knownDd_N_nest_M%>%
  inner_join(NBs, by = 'Sample_Name_1',relationship = "many-to-many")


CF_all_unringed_linked_blood=rbind(CF_sample_name2_match_NB_S_N_nest_M,CF_sample_name1_match_NB_nest_mark)%>%
  rbind(CF_sample_name2_match_NB_nest_mark)%>%
  rbind(CF_sample_name1_match_N_nest_M)%>%
  unique()%>%
  rename(clutch=RaisedClutchId)

all_unringed_ids=rbind(all_unringed_linked_blood,CF_all_unringed_linked_blood)%>%
  unique()%>%
  select(Sample_Name_1, Sample_Name_2, BirdId, Tissue.Type, Tissue.Freezer.Rack.Box.Localisation, 
         Tissue_Container_Name, Tissue_Container_Place, Tissue.Number, 
         Date.extraction, DNA.Freezer.Rack.Box.Localisation, 
         DNA_Box_Name, DNA_Tube_.Place, DNA_Tube_Nb, DNA.C.ng.ul.)%>%
  filter(Date.extraction!='')

n_distinct(all_unringed_ids$BirdId)

write.table(all_unringed_ids,
            file = "New_owls_to_seq/P1_UNRINGED_Ids_to_seq_DNA_extracted.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

WriteXLS(all_unringed_ids, "New_owls_to_seq/P1_UNRINGED_Ids_to_seq_DNA_extracted.xlsx")
