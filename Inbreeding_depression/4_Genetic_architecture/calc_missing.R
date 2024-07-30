library(tidyverse)


clutch=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Clutch.csv")%>%
  select(Id, Season, CrossFostered)%>%
  rename(BornClutchId=Id)

mesurements=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv")%>%
  select(RingId, LeftWing, Mass, ObservationDate)%>%
    mutate(ObservationDate = as.Date(ObservationDate, format = "%d/%m/%Y"))
  
bird=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv")%>%
  select(RingId, HatchDate, DeathDate, BornClutchId, RaisedClutchId, LatestBirdObservationDate)%>%
  mutate( ## getting dates into date format
    HatchDate = as.Date(HatchDate, format = "%d/%m/%Y"),
    DeathDate = as.Date(DeathDate, format = "%d/%m/%Y"),
    LatestBirdObservationDate = as.Date(LatestBirdObservationDate, format = "%d/%m/%Y"))


lab=read.csv("Barn_owl_general/Data_base_dfs/Lab_20240523105823_LabIndividualTissue.csv")%>%
  select(IndividualIdentifier,OtherIdentifiers, SampleType,LabIndividualTissueId,
         TissueIdentifier,LabIndividualId, BoxCode, BoxName, TubeIdentifier)%>%
  rename(RingId=IndividualIdentifier)


all_records_clutch=mesurements%>%
  right_join(bird, by = 'RingId', relationship = "many-to-many")%>%
  right_join(clutch, by='BornClutchId', relationship = "many-to-many")%>%
  filter(CrossFostered==0|is.na(CrossFostered))%>%
  filter(RingId!=''&!is.na(RingId),
         HatchDate > as.Date("2010-01-01"))%>% # hatched after / in 2010
  drop_na(Mass)


## only for non-cross fostered ids for now cos its a bit of a headache otherwise
missings <- all_records_clutch %>%
  filter(ObservationDate <= as.Date(paste(Season+1, 1, 1, sep = "-"))) %>% ## filter for only observations taken in the same season born
  mutate(age_at_latest=as.integer(LatestBirdObservationDate - HatchDate))%>%
  group_by(BornClutchId) %>%
  ## Checking if there's any observations within the same clutch later than the last obs of the individual i.e. a clutch is visited again but there are no records of the id here
  # conditional on whether the observation date is within 35 days of the latest bird obs (in case observations of others in clutch could have occured after id fledged)
  mutate(obs_after = as.integer(
    sapply(LatestBirdObservationDate, function(x) any((ObservationDate > x)& (ObservationDate <= x + 35)))
  ))%>%
    ## define missing ids if the last obs is before other obs taken in the clutch, the chick is within 45 days old and the age of the latest observation was <=35 days 
  # so couldnt possibly have fledged but records stop for that id (i.e. likely because it is dead)
  mutate(missing = case_when(
    obs_after==1 & (LatestBirdObservationDate<HatchDate+45) & (age_at_latest<=35) ~1
  ))%>%
  ungroup()%>%
  filter(missing==1)


missings_with_blood=missings%>%
left_join(lab,  by = "RingId", relationship = "many-to-many")%>% #joing with lab info
  filter(SampleType%in%c("blood","Blood cell"))%>%
  distinct(RingId, .keep_all = TRUE) %>% #
  filter(is.na(DeathDate))



n_distinct(missings_with_blood$RingId)


write.csv(missings_with_blood,"Barn_owl_general/missing_ids/missingIDs_with_blood.csv")




