library(dplyr)
store_plant_data <- function(){
 # Get plot treatments per year (treatments have changed over time):
    # Download primary Portal data tables
primary_tables <- portalr::load_data(path = '/Users/renatadiaz/Documents/GitHub/weecology')

    # Get plot treatments for summer & winter months months
plot_treatments_summer <- primary_tables[[5]] %>%
  dplyr::filter(month %in% 4:9) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()

plot_treatments_winter <- primary_tables[[5]] %>%
  dplyr::filter(month %in% c(1:3, 10:12)) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()


# Download all plant data from the Portal repo
all_plants <- portalr::plant_abundance(path = '/Users/renatadiaz/Documents/GitHub/weecology', level = 'Plot', type = 'Annuals', plots = 'Longterm', unknowns = F, 
                          correct_sp = T, shape = 'flat', na_drop = T, effort = T)

# Filter to summer plants, and combine plant data tables with treatment tables
summer <- all_plants %>% 
  dplyr::filter(season == 'summer') %>%
  dplyr::left_join(plot_treatments_summer, by = c('year', 'plot'))

# Do the same for winter plants
winter <- all_plants %>% 
  dplyr::filter(season == 'winter') %>%
  dplyr::left_join(plot_treatments_winter, by = c('year', 'plot'))

# Save
write.csv(summer, 'data/summer-plants-raw.csv', row.names = F)
write.csv(winter, 'data/winter-plants-raw.csv', row.names = F)

}

# Function to download table of plant species names & species codes
# (this table has other info, but that is the part I use here)
store_plant_species_info <- function() {
  species_data <- download.file(url = 'https://raw.githubusercontent.com/weecology/PortalData/master/Plants/Portal_plant_species.csv',
                                destfile = 'data/plant_species_data.csv')
}