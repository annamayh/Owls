library(tidyverse)
library(pedigree)
library(pedigreemm)


ped <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
                dam  = c(NA,NA,2,NA,3,2), label= 1:6)
inbreeding(ped)


owl_ped_df=read.table("sequioa/coreccted_consensus_pedigree.txt", sep=',', header = T)

owl_ped_df=owl_ped_df[-277,]
owl_ped_df=owl_ped_df[-2329,]


ped_order = orderPed(owl_ped_df)



ped = owl_ped_df[order(ped_order),]

str(ped)

ped$Fped=calcInbreeding(ped)#%>%


ped_corrected3=read.table("sequioa/pedigreeCORRECTED.tab", header = T)

ped_corrected3[ped_corrected3$momid=='M038112',]
