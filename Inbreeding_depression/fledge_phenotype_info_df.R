## making df for : ###
### modelling inbreeding depression in tarsus length for fledelings ###
## using glmm, and animal models (ped RM and GRM)
setwd("C:/Users/s1881212/Downloads")

##Fgrm and FROH from elo
inb=read.table("All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)
## reading in phenotype info
owl_mes=read.csv("BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
  select(RingId,LeftTarsus, Mass, LeftWing, BillLength, Observer,ObservationDate,EstimatedGrowthStage)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  filter(EstimatedGrowthStage=="Fledgling"|EstimatedGrowthStage=="Nestling or Fledgling") ##only taking nestlings/fledgelings i.e. first year of life 

owl_sex=read.table("GeneticSex_3K_RP548.txt", header = T) # sex info 
## info on clutch 
owl_bird=read.csv("BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, BornClutchId, RaisedClutchId, SiteId, HatchDate)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA))%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  unite("clutch_merge", BornClutchId:RaisedClutchId, remove = F)# merging raised and born clutch to one variable 
  

##merge all df together by RingID
fledge_merge=inb%>%
  inner_join(owl_mes)%>%
  inner_join(owl_sex)%>%
  inner_join(owl_bird,relationship = "many-to-many")%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))

fledge_merge$date_diff <- as.Date(as.character(fledge_merge$obs_date))-
  as.Date(as.character(fledge_merge$hatch_date)) ## getting rough age of id when observation was taken based on hatch date

## also getting info on rank i.e. what order they were born in in the nest
ranked=fledge_merge%>%
  select(RingId,hatch_date,RaisedClutchId)%>%
  distinct(RingId, .keep_all = TRUE)%>%
  group_by(RaisedClutchId)%>%
  mutate(rank = rank(hatch_date, ties.method= "min")) %>%
  arrange(hatch_date)

fledge_mes_rep=fledge_merge%>%
  inner_join(ranked,relationship = "many-to-many")

###  
fledge_pheno_df=fledge_merge%>%
  mutate(age_days=as.numeric(date_diff))%>%
  mutate(age_days_sq=age_days^2)%>%
  filter(age_days<90)%>%
  as.data.frame()

n_distinct(fledge_pheno_df$RingId) ## 2299 diff ids and 14867 obs but some have missing info for pheno traits

## ~~ Tarsus df ~~ ###
fledge_tarsus_df=fledge_pheno_df%>%
  select(-Mass, -LeftWing, -BillLength)%>% ##remove phenotype info we may not have 
  na.omit() #remove NAs

n_distinct(fledge_tarsus_df$RingId) ##2262 ids and 6490 records


## ~~ Mass df ~~ ####
fledge_mass_df=fledge_pheno_df%>%
  select(-LeftTarsus, -LeftWing, -BillLength)%>% ##remove phenotype info we may not have 
  na.omit() #remove NAs

n_distinct(fledge_mass_df$RingId) ##2294 ids with mass records, 9864 records

## ~~ Bill length  ~~ ####
fledge_bill_df=fledge_pheno_df%>%
  select(-LeftTarsus, -LeftWing, -Mass)%>% ##remove phenotype info we may not have 
  na.omit() #remove NAs

n_distinct(fledge_bill_df$RingId) ##2262 ids with records, 6128 records


## ~~ wing length ~~  ####
### more records for wing length than tarsus
fledge_wing_df=fledge_pheno_df%>%
  select(-LeftTarsus, -BillLength, -Mass)%>% ##remove phenotype info we may not have 
  na.omit() #remove NAs

n_distinct(fledge_wing_df$RingId) ##2296 ids with records, 9410 records


