library(AGHmatrix)
library(tidyverse)
library(brms)
library(pedigree)
library(pedigreemm)


setwd("/Users/ahewett1/Documents")


## reading in phenotype info
owl_mes=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
  select(RingId,LeftTarsus, Mass, LeftWing, BillLength, Observer,ObservationDate,EstimatedGrowthStage, PhenotypeSex)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))
n_distinct(owl_mes$RingId) ## 

## info on clutch 
owl_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, BornClutchId, RaisedClutchId, SiteId, HatchDate)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA))%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  unite("clutch_merge", BornClutchId:RaisedClutchId, remove = F)# merging raised and born clutch to one variable 

## getting info on rank within the clutch (i.e. what order they were born)
rank=owl_bird%>%
  select(RingId,hatch_date,RaisedClutchId)%>%
  distinct(RingId, .keep_all = TRUE)%>%
  group_by(RaisedClutchId)%>%
  mutate(rank = rank(hatch_date, ties.method= "min")) %>%
  arrange(hatch_date)

##merge clutch and pehntpype info
clutch_pheno=owl_mes%>%
  inner_join(owl_bird,relationship = "many-to-many")%>%
  separate_wider_delim(hatch_date, delim = "-", cols_remove=F,names = c("year",NA,NA))

clutch_pheno$date_diff <- as.Date(as.character(clutch_pheno$obs_date))-
  as.Date(as.character(clutch_pheno$hatch_date)) ## getting rough age of id when observation was taken based on hatch date


###  
all_pheno_df=clutch_pheno%>%
  inner_join(rank,relationship = "many-to-many")%>% # because many repeated measures for the same id
  mutate(age_days=as.numeric(date_diff))%>%
  group_by(RingId)%>%
  mutate(min_mes=min(age_days))%>%
  filter(!age_days<0)%>% ## for ids where measurement is before hatch date (9 ids)
  as.data.frame()

n_distinct(all_pheno_df$RingId) ## 9624 diff ids and total of 17073 obs but some have missing info for pheno traits

## ~~ Tarsus df ~~ ###
tarsus_df=all_pheno_df%>%
  select(RingId, LeftTarsus, age_days,rank, clutch_merge, Observer, year, min_mes, PhenotypeSex)%>% ##keep only pheno info interested in
  na.omit(LeftTarsus) %>% #remove NAs for important bits
  unique() ## duplicates from repeated tracking ids on 

n_distinct(tarsus_df_all$RingId) ##7588 ids and 23168 records
## adding extra info from other life stages doesnt increase the number of ids but does inc the number of records by ~ 500 for each  
## most of the ids are ones we have info from clutch etc. 



## pedigree
ped_corrected=read.table("sequioa/pedigreeCORRECTED.tab", header = T)%>%
  select(-sex)%>%
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

ped_corrected[is.na(ped_corrected)] <- 0

a_mat=Amatrix(ped_corrected, ploidy=2)

ped_ordered=orderPed(ped_corrected)
ord<-orderPed(ped_corrected)
ped_ordered=ped_corrected[order(ord),]

Fped=calcInbreeding(ped_ordered)%>%
  as.data.frame()

sum(Fped$.>0)

Fped2=inbreeding(ped_corrected)

## setting up model
tarsus_df[,'clutch_merge']=as.factor(tarsus_df[,'clutch_merge'])
tarsus_df[,'GeneticSex']=as.factor(tarsus_df[,'GeneticSex'])
tarsus_df[,'RingId']=as.factor(tarsus_df[,'RingId'])
tarsus_df[,'year']=as.factor(tarsus_df[,'year'])
tarsus_df[,'Observer']=as.factor(tarsus_df[,'Observer'])
tarsus_df[,'rank']=as.numeric(tarsus_df[,'rank'])


## first 3 priors are from get_priors() function in brms 
prior=c(prior(student_t(3, 694,40), class = "Intercept"), ## for fixed effect of intercept setting mean as the rough weight
        prior(student_t(3,0,40), class = "sd"),
        prior(student_t(3,0,40), class = "sigma"),
        prior(cauchy(0, 2), class = "sd", group="RingId_pe"))## for random effect of pe expect very low variance centered on 0


tarsus_df[,'RingId_pe']=tarsus_df[,'RingId'] ##add permanent env variable to get h2 estimate


mod_tarsus_all.3.ARM <- brm(LeftTarsus ~  1 + FuniWE+GeneticSex+rank+SSgompertz(age_days, 717, 2.27, 0.89)+
                              (1|gr(RingId, cov=Amat))+(1|RingId_pe)+(1|Observer)+(1|clutch_merge)+(1|year),
                            prior = prior,
                            data = tarsus_df,
                           # control=list(adapt_delta=0.98),
                            data2 = list(Amat = a_mat),
                            chains = 4,
                            cores=4,
                            iter = 10000,
                            warmup = 2500)


