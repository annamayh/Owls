library(tidyverse)
library(WriteXLS)

#seq=read.csv("New_owls_to_seq/Seq_Samples_new_box_positions.csv", header = T)
setwd("/Users/ahewett1/")

seq=read.csv("Downloads/Seq_Samples_new_box_positions_27.09.24_afternoon.csv", header = T)%>%
  rename(Sample_Name=Sample._Name)

p1=read.table("Documents/New_owls_to_seq/P1_Ids_to_seq_DNA_extracted.txt", 
              sep = ",",na = "NA", header= T)%>%
  rename(Sample_Name=RingId)%>%
  rename(TAKEN_FROM.DNA_Tube_Nb=DNA_Tube_Nb)%>%
  rename(TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation=DNA.Freezer.Rack.Box.Localisation)


p1_unr=read.table("Documents/New_owls_to_seq/P1_UNRINGED_Ids_to_seq_DNA_extracted.txt", 
                    sep = ",",na = "NA",header= T)%>%
  rename(TAKEN_FROM.DNA_Tube_Nb=DNA_Tube_Nb)%>%
  rename(TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation=DNA.Freezer.Rack.Box.Localisation)%>%
  rename(Sample_Name=Sample_Name_1)

p1_unr2=read.table("Documents/New_owls_to_seq/P1_UNRINGED_Ids_to_seq_DNA_extracted.txt", 
                  sep = ",",na = "NA",header= T)%>%
  rename(TAKEN_FROM.DNA_Tube_Nb=DNA_Tube_Nb)%>%
  rename(TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation=DNA.Freezer.Rack.Box.Localisation)%>%
  rename(Sample_Name=Sample_Name_2)


head(final_boxes)

final_boxes=seq%>%
  left_join(p1)%>%#by = c('Sample_Name','TAKEN_FROM.DNA_Tube_Nb','TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation' ))%>%
  left_join(p1_unr)%>% #by = c('Sample_Name','TAKEN_FROM.DNA_Tube_Nb','TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation' ))%>%
  left_join(p1_unr2)%>%
  select(today, Ana.File, Sample_Name, NEW_Box_and_position, Final.box, Final.box.position,
         TAKEN_FROM.DNA_Tube_Nb, DNA.C.ng.ul., Date.extraction,  
         TAKEN_FROM..DNA.Freezer.Rack.Box.Localisation,TAKEN_FROM..DNA_Tube_Place, DNA.Vol..50ul, Comments)


WriteXLS(final_boxes, "Documents/New_owls_to_seq/Final_boxes_dna_conc.xlsx")



n_distinct(seq$Sample._Name)


unique=seq%>%
  filter(duplicated(Sample._Name))


dups=seq%>%
  group_by(Sample._Name)%>%
  filter(n()>1)%>%
  ungroup()

n_distinct(unique$Sample._Name)
 ## 38 replicates