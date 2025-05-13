## making df for : ###
### modelling inbreeding depression in tarsus length  ###
## using glmm, and animal models (ped RM and GRM)
## as we take into account clutch info most ids are ones that we have fledeling records from 
library(tidyverse)

setwd("/Users/ahewett1/Documents")


##Fgrm and FROH from elo
inb_FROH=read.table("Inbreeding_depression_owls/pheno_df/FHBD_512gen_goodQ_ss_only.txt",  sep=',')%>%
  select('V1', 'V39')%>% #only get genome-wide FROH (other cols are per chr FROH)
  rename(RingId=V1, FROH=V39) 

inb_FuniW=read.table("Inbreeding_depression_owls/pheno_df/FuniW_goodQ_ss_only.txt", sep=',')%>%
  rename(RingId=V1, FuniW=V2) 

#table(owl_mes$EstimatedGrowthStage)

## reading in phenotype info
owl_mes=read.csv("Barn_owl_general/Data_base_dfs//BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
  select(RingId,LeftTarsus, Mass, LeftWing, BillLength, HookLength, Observer,ObservationDate,EstimatedGrowthStage)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  mutate(stage=case_when(
    grepl("Adult", EstimatedGrowthStage) ~ "adult",
    EstimatedGrowthStage=="Fledgling"~"juvenile", 
    EstimatedGrowthStage=="Nestling or Fledgling"~"juvenile", 
    ))

#table(owl_mes$EstimatedGrowthStage)



owl_sex_g=read.table("Inbreeding_depression_owls/GeneticSex_3K_RP548.txt", header = T) # sex info 
owl_sex_p=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, PhenotypeSex)

owl_sex_p$phen_sex=if_else(owl_sex_p$PhenotypeSex=="Male", 1, 
                           if_else(owl_sex_p$PhenotypeSex=="Female",2, 
                                   NA))

sex=owl_sex_g%>%
  left_join(owl_sex_p, by='RingId')%>%
  mutate(sex=case_when(
    GeneticSex==1 ~ 1, ## male = 1, female =2
    GeneticSex==2 ~ 2, 
    (is.na(GeneticSex)&phen_sex==1) ~ 1,
    (is.na(GeneticSex)&phen_sex==2) ~ 2,
    (GeneticSex==1&phen_sex==2) ~ 1,
    (GeneticSex==2&phen_sex==1) ~ 2))%>%
  as.data.frame() %>%
  unique()%>%
  select(RingId, sex)


# owl_site=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Site.csv", header = T)%>%
#   select(SiteId, CH1903X, CH1903Y)%>%
#   mutate(SiteId=as.character(SiteId))
# 
# check=owl_bird%>%
#    filter(HatchDateErrorMargin>180)%>%
#   drop_na(BornClutchId)%>%
#   inner_join()
#   

## info on clutch 
owl_bird=read.csv("Barn_owl_general/Data_base_dfs/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, BornClutchId, RaisedClutchId, SiteId, Nestbox, HatchDate, HatchDateErrorMargin)%>%
 # filter(HatchDateErrorMargin<180|is.na(HatchDateErrorMargin))%>% ## for now this isnt needed as ids with genetic info are all adults (so wouldnt matter what their age is in this analysis adults are just adults)
  select(-HatchDateErrorMargin)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA))%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  unite("clutch_merge", BornClutchId:RaisedClutchId, remove = F)%>%# merging raised and born clutch to one variable 
  mutate(SiteId=as.character(SiteId))%>%
  mutate(nestboxID=case_when( ## creating new variable for nestbox id for when sites have 2 nestboxes
    Nestbox=="" ~ SiteId, 
    Nestbox!="" ~ Nestbox, 
  ))%>%
  left_join(owl_site)


##merge all df together by RingID
fledge_merge=inb_FROH%>%
  left_join(inb_FuniW)%>%
  inner_join(owl_mes)%>%
  inner_join(sex)%>%
  inner_join(owl_bird,relationship = "many-to-many")%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year","month",NA))%>%
  mutate(month=as.numeric(month))

fledge_merge$date_diff <- as.Date(as.character(fledge_merge$obs_date))-
  as.Date(as.character(fledge_merge$hatch_date)) ## getting rough age of id when observation was taken based on hatch date

## also getting info on rank i.e. what order they were born in in the nest
ranked=fledge_merge%>%
  select(RingId,hatch_date,RaisedClutchId)%>%
  distinct(RingId, .keep_all = TRUE)%>%
  group_by(RaisedClutchId)%>%
  mutate(rank = rank(hatch_date, ties.method= "min")) %>%
  arrange(hatch_date)

#combine with ranked info
fledge_mes_rep=fledge_merge%>%
  inner_join(ranked,relationship = "many-to-many") # because many repeated measures for the same id

###  
all_pheno_df=fledge_mes_rep%>%
  mutate(age_days=as.numeric(date_diff))%>%
  group_by(RingId)%>%
  mutate(min_mes=min(age_days))%>%
  filter(!age_days<0)%>% ## for ids where measurement is before hatch date (9 ids)
  as.data.frame()%>%
  mutate(julian_hatchdate=format(hatch_date, "%j"))%>%
  mutate(gr_stage=case_when(
    age_days<=180 ~ "Juvenile", 
    age_days>180 ~ "Adult"
  ))%>%
  mutate(dev_stage=case_when(
    age_days<=10 ~ "Day_0-10", 
    age_days>10&age_days<=20 ~ "Day_10-20",
    age_days>20&age_days<=40 ~ "Day_20-40", 
    age_days>40&age_days<=60  ~ "Day_40-60",
    age_days>60 ~ "Day_60plus")
  )%>%
  filter(min_mes>0)# filter 




n_distinct(all_pheno_df$RingId) ## 2300 diff ids and total of 17073 obs but some have missing info for pheno traits

### REMOVING THE ~250 IDS THAT WERE SEQ ON MED COV AND HAVE BIAS UPWARDS ESTIMATES
low_cov=read.table("Inbreeding_depression_owls/pheno_df/RingIdsCov5.txt",sep=",", header=T)

## ~~ Mass df ~~ ####
mass_df_all=all_pheno_df%>%
  select(RingId, Mass, FROH, FuniW, age_days,sex,rank, clutch_merge, Observer, year, nestboxID, gr_stage, julian_hatchdate, month, dev_stage, min_mes)%>% ##remove phenotype info we may not have 
  drop_na(Mass,FuniW) %>%#remove NAs
  inner_join(low_cov)%>%
  unique() %>%## duplicates from repeated tracking ids on 
  filter(Mass<1000) %>%#remove 1 incorrect mass record 
  mutate(mass_scale=Mass/100)%>%
  mutate(age_acc=362*exp(-4.01*0.893^age_days)) %>% ## account for age using gompertz growth
  mutate(RingId_pe=RingId) # add permanent environment for repeated measures 

#plot(mass_df_all$age_days[mass_df_all$age_days < 90], mass_df_all$Mass[mass_df_all$age_days < 90])

n_distinct(mass_df_all$RingId) ##1996 ids with mass records,  

table(mass_df_all$gr_stage)


## ~~ Bill length  ~~ ####
bill_df_all=all_pheno_df%>%
  select(RingId, BillLength, FROH, FuniW, age_days,sex,rank, clutch_merge, Observer, year, nestboxID, gr_stage,julian_hatchdate, month, dev_stage, min_mes)%>% ##remove phenotype info we may not have 
  na.omit(BillLength,FuniW) %>%#remove NAs
  inner_join(low_cov)%>%
  filter(!(BillLength<50 & age_days >5))%>% #filtering out 1 ids that are probably incorrect 
  unique() %>%## duplicates from repeated tracking ids on 
  mutate(bill_scale=BillLength/100) %>% ## rescaling for ease of model convergence
  mutate(age_acc=184*exp(-0.99*0.932^age_days)) %>%## account for age using gompertz growth
  mutate(RingId_pe=RingId) # add permanent environment for repeated measures 

plot(bill_df_all$age_days[bill_df_all$age_days < 90], bill_df_all$BillLength[bill_df_all$age_days < 90])



n_distinct(bill_df_all$RingId) ##1996 ids with records, 

table(bill_df_all$gr_stage)


## ~~ Tarsus length ~~ ###
tarsus_df_all=all_pheno_df%>%
  select(RingId, LeftTarsus, FROH, FuniW, age_days,sex,rank, clutch_merge, Observer, year, nestboxID, gr_stage, julian_hatchdate, month, dev_stage,min_mes)%>% ##keep only pheno info interested in
  na.omit(FuniW,LeftTarsus) %>% #remove NAs for important bits
  inner_join(low_cov)%>%
  filter(!(LeftTarsus<200 & age_days >20))%>% #filtering out 6 ids that are probably incorrect i.e. when their tarsus is very small for their age
  unique()%>% ## duplicates from repeated tracking ids on \
  mutate(tarsus_scale=LeftTarsus/100)%>% # scaling for easier model convergence 
  mutate(age_acc= 55 + (662*exp(-2.341*0.89^age_days))) %>% ## account for age unsing gompertz growth
  mutate(RingId_pe=RingId) # add permanent environment for repeated measures 


plot(tarsus_df_all$age_days[tarsus_df_all$age_days < 90], tarsus_df_all$LeftTarsus[tarsus_df_all$age_days < 90])

n_distinct(tarsus_df_all$RingId) ##
## adding extra info from other life stages doesnt increase the number of ids but does inc the number of records by ~ 500 for each  
## most of the ids are ones we have info from clutch etc. 


## save all dfs

write.table(tarsus_df_all,
            file = "Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(bill_df_all,
            file = "Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(mass_df_all,
            file = "Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")







####3
## none have measurements at this observation date ... so perhaps its the obs date that is wrong 
# errors=fledge_mes_rep%>%
#   mutate(age_days=as.numeric(date_diff))%>%
#   group_by(RingId)%>%
#   mutate(min_mes=min(age_days))%>%
#   filter(age_days<0)%>% ## for ids where measurement is before hatch date (9 ids)
#   as.data.frame()
# 

