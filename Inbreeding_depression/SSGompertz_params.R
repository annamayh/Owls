library(tidyverse)

setwd("/Users/ahewett1/Documents")

bill_df_all = read.table(
            file = "Inbreeding_depression_owls/pheno_df/bill_all_pheno_df.txt",
           sep=",", header=T)

mass_df_all = read.table(
            file = "Inbreeding_depression_owls/pheno_df/mass_all_pheno_df.txt",
           sep=",", header=T)

tarsus_df_all = read.table(
            file = "Inbreeding_depression_owls/pheno_df/tarsus_all_pheno_df.txt",
           sep=",", header=T)

wing_df_all = read.table(
            file = "Inbreeding_depression_owls/pheno_df/wing_all_pheno_df.txt",
           sep=",", header=T)


####3
## checking growth curve and extracting parameters for growth 
##


##### MASS ######

fm2 <- nls(Mass ~ SSgompertz(age_days, asym, b2, b3),
           data = mass_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>15))
summary(fm2)

# Formula: Mass ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 3.620e+02  1.278e+00  283.24   <2e-16 ***
#   b2   4.013e+00  1.390e-01   28.88   <2e-16 ***
#   b3   8.932e-01  2.087e-03  427.94   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot(mass_df_all$age_days, mass_df_all$Mass, xlim = c(0,90), ylim = c(0,600))
curve(362*exp(-4.01*0.893^x), from = 0, to=90, add = TRUE, col="red", lwd=2)

# 
#####.  BILL ######## 
# 
fm3 <- nls(BillLength ~ SSgompertz(age_days, asym, b2, b3),
           data = bill_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>45))
summary(fm3)
# 
# Formula: BillLength ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 1.836e+02  2.341e-01  784.36   <2e-16 ***
#   b2   9.986e-01  1.442e-02   69.23   <2e-16 ***
#   b3   9.317e-01  8.260e-04 1127.90   <2e-16 ***
#   ---
# 
plot(bill_df_all$age_days, bill_df_all$BillLength, xlim = c(0,90))
curve(184*exp(-0.99*0.932^x), from = 0, to=90, add = TRUE, col="red", lwd=2)


plot(bill_df_all$age_days, bill_df_all$bill_scale, xlim = c(0,90))
curve(1.84*exp(-0.99*0.932^x), from = 0, to=90, add = TRUE, col="red", lwd=2)



### WING #######

fm4 <- nls(LeftWing ~ SSgompertz(age_days, asym, b2, b3),
           data = wing_df_all%>%
             filter(age_days>=14)
)
summary(fm4)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 2.988e+02  2.292e-01  1303.9   <2e-16 ***
#   b2   4.186e+00  1.430e-02   292.6   <2e-16 ***
#   b3   9.454e-01  1.206e-04  7836.0   <2e-16 ***
#   ---


plot(wing_df_all$age_days, wing_df_all$LeftWing, xlim = c(0,100))
curve(298*exp(-4.22*0.945^x), from = 0, to=100, add = TRUE, col="red", lwd=2)


# plot(wing_df_all$age_days, wing_df_all$wing_scale, xlim = c(0,100))
# curve(2.98*exp(-4.19*0.945^x), from = 0, to=100, add = TRUE, col="red", lwd=2)
# 

##need to take c down to 3 dp 



### TARSUS #####
### tarsus is a less accurate measure than wing ###

fm1 <- nls(LeftTarsus ~ SSgompertz(age_days, asym, b2, b3),
           data = tarsus_df_all%>%
             group_by(RingId)%>%
             filter(!min_mes>30)) ## cannot estimate growth curve using all data

summary(fm1)

# Formula: LeftTarsus ~ SSgompertz(age_days, asym, b2, b3)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# asym 7.172e+02  6.038e-01 1187.74   <2e-16 ***
#   b2   2.341e+00  2.689e-02   87.06   <2e-16 ***
#   b3   8.877e-01  7.084e-04 1252.97   <2e-16 ***
# 
plot(tarsus_df_all$age_days, tarsus_df_all$LeftTarsus, xlim = c(0,90))
curve(55 + (662*exp(-2.341*0.89^x)), from = 0, to=400, add = TRUE, col="red", lwd=2)
curve(717*exp(-2.34*0.888^x), from = 0, to=400, add = TRUE, col="blue", lwd=2)


tarsus_nls=nls(
  LeftTarsus ~ asym1 + (asym * exp(-b2 * b3^x)), 
  data = tarsus_df_all,
  start = list(asym1 = 10, asym = 500, b2 = 0, b3 = 0)
)


prior_tarsus_nls<- c(
  prior(normal(50, 50), nlpar = "asym1",  class="b", coef="Intercept"),## priors for intercept expectations
  prior(normal(600, 50), nlpar = "asym",  class="b", coef="Intercept"),## priors for intercept expectations
  prior(normal(1, 5), nlpar = "b",  class="b", coef="Intercept"), ## 
  prior(normal(0, 2), nlpar = "c",  class="b", coef="Intercept") ## 
)

tarsus_nls_brms=brm(
  ## model 
  bf(LeftTarsus ~ asym1 + asym * exp(-b*(c)^age_days),
     asym1 ~ 1,
     asym ~ 1,
     b ~ 1 ,
     c ~ 1, 
     nl=TRUE), 
  data = tarsus_df_all, 
  prior = prior_tarsus_nls,
  cores = 4, 
  chains = 4, 
  init = 0
  
  )


summary(tarsus_nls_brms)


plot(tarsus_df_all$age_days, tarsus_df_all$LeftTarsus, xlim = c(0,90))
curve(227+ 445*exp(-5.08*0.4^x), from = 0, to=400, add = TRUE, col="red", lwd=2)


# 
# 
# 
# 
