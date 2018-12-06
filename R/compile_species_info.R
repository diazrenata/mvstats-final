library(dplyr)
# Function to match seed mass data from Valone (2017) with plant species codes

compile_species_info <- function() {
  
  # Load species table from Portal and Valone data
  plant_species_table <- read.csv('data/plant_species_data.csv', stringsAsFactors = F)
  
  seed_size_table <- read.csv('data/valone-seed-size.csv', stringsAsFactors = F)
  
  # Parse species names in Valone and match to species names in Portal table
  seed_size_genuses <- strsplit(seed_size_table$Species, " ")
  
  seed_sizes <- seed_size_table
  seed_sizes$genus <- NA
  seed_sizes$species <- NA
  
  for (i in 1:length(seed_size_genuses)){
    seed_sizes$genus[i] <- unlist(seed_size_genuses[[i]])[1]
    seed_sizes$species[i] <- unlist(seed_size_genuses[[i]])[2]
  }
  
  colnames(seed_sizes) <- c("FullName", "seed_mass", 'type', "genus", "species")
  seed_sizes <- seed_sizes %>%
    select(genus, species, seed_mass)
  
  # Join seed mass table with species codes 
  plant_species_table <- plant_species_table %>%
    select(speciescode, genus, species, community) %>%
    left_join(seed_sizes, by = c('genus', 'species')) %>%
    mutate(speciescode = gsub(" ", ".", speciescode))
  
  write.csv(plant_species_table, 'data/compiled_plant_species_data.csv', row.names = F)
  
}
