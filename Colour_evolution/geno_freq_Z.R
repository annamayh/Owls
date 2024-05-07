library(SNPRelate)

setwd("/Users/ahewett1/Documents")

# load the data with many phenotypes, not all are useful
NewPhenoColor_Order3K = read.table("Gene_drop/NewPhenoColor_Order3K.txt", header = TRUE)
# load the table with the genotypes, you will meed the package "SNPRelate"
gen_mat <- snpgdsGetGeno('Gene_drop/All3102_NewNames_RP502SNPs.variantsList.gds')

TabSummary = data.frame(NewPhenoColor_Order3K$ID, NewPhenoColor_Order3K$MeanColor, NewPhenoColor_Order3K$SurfaceBlackSpotsCm, NewPhenoColor_Order3K$SexWGS, NewPhenoColor_Order3K$StageBin, gen_mat[,1], gen_mat[,2], gen_mat[,3])
colnames(TabSummary) = c("RingId", "MeanColor", "SurfaceBlackSpotsCm", "SexWGS", "StageBin", "GenoCorin", "GenoMC1R", "GenoZ")


Z_females=TabSummary%>%
  select(RingId, SexWGS, GenoZ)%>%
  filter(SexWGS==2)%>%
  mutate(GenoZ=as.factor(GenoZ))


Z_males=TabSummary%>%
  select(RingId, SexWGS, GenoZ)%>%
  filter(SexWGS==1)


ped_birth_yr2=ped_birth_yr%>%
  unique()%>%
  mutate(year=as.numeric(year))%>%
  filter(year>=1997)%>%
  filter(year<2018)%>%
  mutate(cohort=(year-min(year))+1)%>%
  filter(RingId!="889913") %>% ##sort this out later
  select(RingId, year, cohort)%>%
  mutate(year=as.factor(year))
  


z_fem_yr=Z_females%>%
  left_join(ped_birth_yr2)%>%
  na.omit()%>%
  group_by(year, GenoZ)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = GenoZ, values_from = n)%>%
  rename(A1="2", A2="0")%>%
  mutate(a1_geno_freq=A1/(A1+A2))%>%
  mutate(a2_geno_freq=A2/(A1+A2))


plot(z_fem_yr$year, z_fem_yr$allele_freq)



z_male_yr=Z_males%>%
  left_join(ped_birth_yr2)%>%
  na.omit()%>%
  group_by(year, GenoZ)%>%
  summarise(n=n())%>%
  pivot_wider(names_from = GenoZ, values_from = n)%>%
  rename(A11="2", A12="1", A22="0")%>%
  mutate_at(c("A11", "A12", "A22"), ~replace(., is.na(.), 0))%>%
  mutate(a11_geno_freq=A11/(A11+A12+A22))%>%
  mutate(a12_geno_freq=A12/(A11+A12+A22))%>%
  mutate(a22_geno_freq=A22/(A11+A12+A22))



all_yrs=z_male_yr%>%
  right_join(z_fem_yr)%>%
  mutate(p=(1/3*((2*a11_geno_freq)+a12_geno_freq+a1_geno_freq)))


ggplot(all_yrs, aes(x=year, y=p, group=1))+
  geom_point()+
  geom_line()
