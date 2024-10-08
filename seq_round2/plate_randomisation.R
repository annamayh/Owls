## script for creating plates 

library(tidyverse)
library(WriteXLS)

boxes=read.csv("New_owls_to_seq/Seq_Samples_new_box_positions.csv", header = T)%>%
  select(Sample._Name,NEW_Box_and_position,)%>%
  filter(Sample._Name!="")



n_unique_ids=n_distinct(boxes$Sample._Name)

individuals=unique(boxes$Sample._Name)

wells=96
n_slots=wells-3
n_plates=ceiling(n_unique_ids / n_plate_slots) # ceiling becasue we cant have half a plate

control_wells <- c("A1", "H12", "G12")

# Row letters and column numbers
rows <- LETTERS[1:8]  # A to H
columns <- 1:12  # 1 to 12


all_wells <- as.vector(outer(rows, columns, paste0))  # "A1", "A2", ..., "H12"

# Exclude control wells from the assignable wells
available_wells <- setdiff(all_wells, control_wells)

assignments <- data.frame(individuals = character(), Plate = integer(), Well = character(), stringsAsFactors = FALSE)

# Assign individuals to plates
individual_idx <- 1  # Start index for individuals

for (plate in 1:n_plates) {
  # Number of individuals to assign to this plate (max 95 or remaining individuals)
  n_plate_individuals <- min(n_slots, length(individuals) - individual_idx + 1)
  
  # Randomly select wells for the individuals on this plate
  set.seed(42)  # Optional: For reproducibility
  selected_wells <- sample(available_wells, n_plate_individuals)
  
  # Add assignments for this plate
  plate_assignments <- data.frame(
    Sample._Name = individuals[individual_idx:(individual_idx + n_plate_individuals - 1)],
    Plate = plate,
    Well = selected_wells,
    stringsAsFactors = FALSE
  )
  
  # Append control wells to the plate assignments
  control_assignments <- data.frame(
    Sample._Name = c("Control_A1", "Control_H12", "Control_G12"),
    Plate = plate,
    Well = control_wells,
    stringsAsFactors = FALSE
  )
  
  # Combine individuals and controls
  plate_assignments <- rbind(plate_assignments, control_assignments)
  
  # Append the assignments to the final data frame
  assignments <- rbind(assignments, plate_assignments)
  
  # Update the individual index
  individual_idx <- individual_idx + n_plate_individuals
  
  # Break if all individuals have been assigned
  if (individual_idx > length(individuals)) break
}



## merge to original 

plates_df=boxes%>%
  right_join(assignments, by='Sample._Name')

WriteXLS(plates_df, "New_owls_to_seq/P1_plate_assignments.xlsx")
