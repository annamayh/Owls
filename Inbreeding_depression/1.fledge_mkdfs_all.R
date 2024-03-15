## making df for : ###
### modelling inbreeding depression in tarsus length for fledelings ###
## using glmm, and animal models (ped RM and GRM)

setwd("/Users/ahewett1/Documents")


##Fgrm and FROH from elo
inb=read.table("Inbreeding_depression_owls/All3085_FuniWE_FHBD512g.txt", header = T)%>%
  rename(RingId=INDVs)

## reading in phenotype info
owl_mes=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_BirdMeasurement.csv",header = T)%>%
  select(RingId,LeftTarsus, Mass, LeftWing, BillLength, Observer,ObservationDate,EstimatedGrowthStage)%>% ##only selecting the paramters we need
  separate_wider_delim(ObservationDate, delim = " ", names = c("obs_date",NA))%>%
  mutate(obs_date=as.Date(obs_date, format = "%d/%m/%Y"))%>%
  filter(EstimatedGrowthStage=="Fledgling"|EstimatedGrowthStage=="Nestling or Fledgling")##only taking nestlings/fledgelings i.e. first year of life 
  

owl_sex=read.table("Inbreeding_depression_owls/GeneticSex_3K_RP548.txt", header = T) # sex info 
## info on clutch 
owl_bird=read.csv("Barn_owl_general/BarnOwls_Legacy_20231010153920/BarnOwls_Legacy_20231010153920_Bird.csv", header = T)%>%
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

#combine with ranked info
fledge_mes_rep=fledge_merge%>%
  inner_join(ranked,relationship = "many-to-many") # because many repeated measures for the same id

###  
fledge_pheno_df=fledge_mes_rep%>%
  mutate(age_days=as.numeric(date_diff))%>%
  mutate(age_days_sq=age_days^2)%>%
  filter(age_days<90)%>% ## removes ~30 ids that are below 90 days old
  as.data.frame()

n_distinct(fledge_pheno_df$RingId) ## 2299 diff ids and total of 14867 obs but some have missing info for pheno traits

## ~~ Tarsus df ~~ ###
fledge_tarsus_df=fledge_pheno_df%>%
 # mutate(tarsus_scale=LeftTarsus/100)%>%
  select(RingId, LeftTarsus, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year)%>% ##keep only pheno info interested in
  na.omit(LeftTarsus,FuniWE) %>% #remove NAs for important bits
  unique() ## duplicates from repeated tracking ids on 

n_distinct(fledge_tarsus_df$RingId) ##2262 ids and 6345 records


## ~~ Mass df ~~ ####
fledge_mass_df=fledge_pheno_df%>%
  select(RingId, Mass, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year)%>% ##remove phenotype info we may not have 
 # mutate(mass_scale=Mass/100)%>%
  na.omit(Mass,FuniWE) %>%#remove NAs
  unique() ## duplicates from repeated tracking ids on 

n_distinct(fledge_mass_df$RingId) ##2294 ids with mass records, 9675 records


## ~~ Bill length  ~~ ####
fledge_bill_df=fledge_pheno_df%>%
  select(RingId, BillLength, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year)%>% ##remove phenotype info we may not have 
  na.omit(BillLength,FuniWE) %>%#remove NAs
  unique() ## duplicates from repeated tracking ids on 


n_distinct(fledge_bill_df$RingId) ##2262 ids with records, 5983


# ~~ wing length ~~  ####
## more records for wing length than tarsus so maybe good to look at both??
fledge_wing_df=fledge_pheno_df%>%
  select(RingId, LeftWing, FHBD512gen, FuniWE, age_days,GeneticSex,rank, clutch_merge, Observer, year)%>% ##remove phenotype info we may not have 
  na.omit(LeftWing,FuniWE) #remove NAs

n_distinct(fledge_wing_df$RingId) ##2296 ids with records, 9410 records



# write.table(fledge_bill_df, 
#             file = "Inbreeding_depression_owls/pheno_df/bill_fledge_pheno_df.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 
# write.table(fledge_mass_df, 
#             file = "Inbreeding_depression_owls/pheno_df/mass_fledge_pheno_df.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 
# write.table(fledge_tarsus_df, 
#             file = "Inbreeding_depression_owls/pheno_df/tarsus_fledge_pheno_df.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 
# write.table(fledge_wing_df, 
#             file = "Inbreeding_depression_owls/pheno_df/wing_fledge_pheno_df.txt",
#             row.names = F, quote = F, sep = ",",na = "NA")
# 


## checking growth curve for all 

fm1 <- nls(LeftTarsus ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_tarsus_df)
summary(fm1)

plot(fledge_tarsus_df$age_days, fledge_tarsus_df$LeftTarsus, xlim = c(0,90))
curve(720*exp(-2.26*0.89^x), from = 0, to=90, add = TRUE, col="red", lwd=2)
#curve(720-(720*exp(-0.08*x)), from = 0, to=90, add = TRUE, col="blue", lwd=2)




fm2 <- nls(Mass ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_mass_df)
summary(fm2)

# Formula: Mass ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.691e+02  7.583e-01   486.7   <2e-16 ***
#   b2   4.417e+00  1.181e-01    37.4   <2e-16 ***
#   b3   8.905e-01  1.358e-03   655.7   <2e-16 ***

plot(fledge_mass_df$age_days, fledge_mass_df$Mass, xlim = c(0,90))
curve(370*exp(-4.42*0.89^x), from = 0, to=90, add = TRUE, col="red", lwd=2)



fm3 <- nls(BillLength ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_bill_df)
summary(fm3)

# ormula: BillLength ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 1.798e+02  2.501e-01  719.02   <2e-16 ***
#   b2   1.132e+00  1.821e-02   62.17   <2e-16 ***
#   b3   9.203e-01  9.991e-04  921.16   <2e-16 ***

plot(fledge_bill_df$age_days, fledge_bill_df$BillLength, xlim = c(0,90))
curve(180*exp(-1.13*0.92^x), from = 0, to=90, add = TRUE, col="red", lwd=2)





fm4 <- nls(LeftWing ~ SSgompertz(age_days, asym, b2, b3),
           data = fledge_wing_df)
summary(fm4)

# Formula: LeftWing ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.175e+02  7.035e-01   451.4   <2e-16 ***
#   b2   3.815e+00  1.507e-02   253.1   <2e-16 ***
#   b3   9.512e-01  2.041e-04  4660.8   <2e-16 ***

plot(fledge_wing_df$age_days, fledge_wing_df$LeftWing, xlim = c(0,90))
curve(318*exp(-3.82*0.95^x), from = 0, to=90, add = TRUE, col="red", lwd=2)

