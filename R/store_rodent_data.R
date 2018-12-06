library(dplyr)
store_rodent_data <- function(){
  # Import rodent data tables
  primary_tables <- portalr::load_data(path = 'repo')
  
  rodents <- primary_tables[[1]]
  species <- primary_tables[[2]]
  trapping <- primary_tables[[3]]
  newmoons <- primary_tables[[4]]
  plots <- primary_tables[[5]]
  
  # Get plot treatments for each year (treatments have changed over time)
  plots_history <- left_join(plots, trapping, 
                     by = c('year', 'month', 'plot')) 
  plots_history <- plots_history %>%
    select(year, month, plot, treatment, period, sampled)
  
  # Filter to granivores identified to species
  # Combine rodent data & plot treatment tables
  focal_spp <- species %>%
    filter(granivore == 1, unidentified == 0) %>%
    select(species)
  
  these_rodents <- rodents %>%
    filter(species %in% focal_spp$species) %>%
    select(month, year, period, plot, species) %>%
    group_by(month, year, period, plot, species) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    left_join(plots_history, by = c('year', 'month', 'period', 'plot'))
  
  
  write.csv(these_rodents, 'data/rodents-raw.csv', row.names = F)

}