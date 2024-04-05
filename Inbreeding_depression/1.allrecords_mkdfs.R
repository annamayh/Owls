## making df for : ###
### modelling inbreeding depression in tarsus length  ###
## using glmm, and animal models (ped RM and GRM)
## as we take into account clutch info most ids are ones that we have fledeling records from 

setwd("/Users/ahewett1/Documents")


##Fgrm and FROH from elo
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)

## reading in phenotype info
owl_mes=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
  select(RingId,LeftTarsus, Mass, LeftWing, BillLength, Observer,ObservationDate,EstimatedGrowthStage)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))

owl_sex=read.table("Inbreeding_depression_owls/GeneticSex_3K_RP548.txt", header = T) # sex info 
## info on clutch 
owl_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
  select(RingId, BornClutchId, RaisedClutchId, SiteId, Nestbox, HatchDate)%>%
  na.omit()%>%
  separate_wider_delim(HatchDate, delim = " ", names = c("hatch_date",NA))%>%
  mutate(hatch_date=as.Date(hatch_date, format = "%d/%m/%Y"))%>%
  unite("clutch_merge", BornClutchId:RaisedClutchId, remove = F)%>%# merging raised and born clutch to one variable 
  mutate(SiteId=as.character(SiteId))%>%
  mutate(nestboxID=case_when( ## creating new variable for nestbox id for when sites have 2 nestboxes
    Nestbox=="" ~ SiteId, 
    Nestbox!="" ~ Nestbox, 
  ))
  

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

#combine with ranked info
fledge_mes_rep=fledge_merge%>%
  inner_join(ranked,relationship = "many-to-many") # because many repeated measures for the same id

###  
all_pheno_df=fledge_mes_rep%>%
  mutate(age_days=as.numeric(date_diff))%>%
  group_by(RingId)%>%
  mutate(min_mes=min(age_days))%>%
  filter(!age_days<0)%>% ## for ids where measurement is before hatch date (9 ids)
  as.data.frame()



n_distinct(all_pheno_df$RingId) ## 2300 diff ids and total of 17073 obs but some have missing info for pheno traits

## ~~ Tarsus df ~~ ###
tarsus_df_all=all_pheno_df%>%
  select(RingId, LeftTarsus, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year, min_mes, nestboxID)%>% ##keep only pheno info interested in
  na.omit(FuniWE,LeftTarsus) %>% #remove NAs for important bits
  unique() ## duplicates from repeated tracking ids on 

n_distinct(tarsus_df_all$RingId) ##2299 ids and 7102 records
## adding extra info from other life stages doesnt increase the number of ids but does inc the number of records by ~ 500 for each  
## most of the ids are ones we have info from clutch etc. 


## ~~ Mass df ~~ ####
mass_df_all=all_pheno_df%>%
  select(RingId, Mass, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year,min_mes, nestboxID)%>% ##remove phenotype info we may not have 
  na.omit(Mass,FuniWE) %>%#remove NAs
  unique() ## duplicates from repeated tracking ids on 

n_distinct(mass_df_all$RingId) ##2300 ids with mass records, 11150 records


## ~~ Bill length  ~~ ####
bill_df_all=all_pheno_df%>%
  select(RingId, BillLength, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year, min_mes, nestboxID)%>% ##remove phenotype info we may not have 
  na.omit(BillLength,FuniWE) %>%#remove NAs
  unique() ## duplicates from repeated tracking ids on 


n_distinct(bill_df_all$RingId) ##2298 ids with records, 6746


# ~~ wing length ~~  ####
## more records for wing length than tarsus so maybe good to look at both??
wing_df_all=all_pheno_df%>%
  select(RingId, LeftWing, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year, min_mes, nestboxID)%>% ##remove phenotype info we may not have 
  na.omit(LeftWing,FuniWE) #remove NAs

n_distinct(wing_df_all$RingId) ##2296 ids with records, 10747 records



write.table(bill_df_all,
            file = "Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(mass_df_all,
            file = "Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(tarsus_df_all,
            file = "Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")

write.table(wing_df_all,
            file = "Inbreeding_depression_owls/pheno_df/wing_all_pheno_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA")



## checking growth curve and extracting parameters for growth 
##

fm1 <- nls(LeftTarsus ~ SSgompertz(age_days, asym, b2, b3),
           data = tarsus_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>30)) ## cannot estimate growth curve using all data 

summary(fm1)

plot(tarsus_df_all$age_days, tarsus_df_all$LeftTarsus, xlim = c(0,90))
curve(717*exp(-2.27*0.89^x), from = 0, to=400, add = TRUE, col="red", lwd=2)




fm2 <- nls(Mass ~ SSgompertz(age_days, asym, b2, b3),
           data = mass_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>10))
summary(fm2)

# Formula: Mass ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.691e+02  7.583e-01   486.7   <2e-16 ***
#   b2   4.417e+00  1.181e-01    37.4   <2e-16 ***
#   b3   8.905e-01  1.358e-03   655.7   <2e-16 ***

plot(mass_df_all$age_days, mass_df_all$Mass, xlim = c(0,90), ylim = c(0,600))
curve(352*exp(-4.29*0.88^x), from = 0, to=90, add = TRUE, col="red", lwd=2)



fm3 <- nls(BillLength ~ SSgompertz(age_days, asym, b2, b3),
           data = bill_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>45))
summary(fm3)

# Formula: BillLength ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 1.835e+02  2.282e-01  804.08   <2e-16 ***
#   b2   1.010e+00  1.427e-02   70.82   <2e-16 ***
#   b3   9.310e-01  8.083e-04 1151.68   <2e-16 ***

plot(bill_df_all$age_days, bill_df_all$BillLength, xlim = c(0,90))
curve(184*exp(-1.01*0.93^x), from = 0, to=90, add = TRUE, col="red", lwd=2)





fm4 <- nls(LeftWing ~ SSgompertz(age_days, asym, b2, b3),
           data = wing_df_all
           )
summary(fm4)

# Formula: LeftWing ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.175e+02  7.035e-01   451.4   <2e-16 ***
#   b2   3.815e+00  1.507e-02   253.1   <2e-16 ***
#   b3   9.512e-01  2.041e-04  4660.8   <2e-16 ***

plot(wing_df_all$age_days, wing_df_all$LeftWing, xlim = c(0,90))
curve(317*exp(-4.19*0.95^x), from = 0, to=90, add = TRUE, col="red", lwd=2)



errors=fledge_mes_rep%>%
  mutate(age_days=as.numeric(date_diff))%>%
  group_by(RingId)%>%
  mutate(min_mes=min(age_days))%>%
  filter(age_days<0)%>% ## for ids where measurement is before hatch date (9 ids)
  as.data.frame()
