#' @importFrom magrittr "%>%"
#' @title store portal plant dataframes
#'
#' Import Portal rodent data using portalr functions.
#'
#' @export

library(dplyr)
store_plant_data <- function(){
primary_tables <- portalr::load_data(path = 'repo')

plot_treatments_summer <- primary_tables[[5]] %>%
  dplyr::filter(month %in% 4:9) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()

plot_treatments_winter <- primary_tables[[5]] %>%
  dplyr::filter(month %in% c(1:3, 10:12)) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()


all_plants <- portalr::plant_abundance(path = 'repo', level = 'Plot', type = 'All', plots = 'Longterm', unknowns = F, 
                          correct_sp = T, shape = 'flat', na_drop = T, effort = T)

summer <- all_plants %>% 
  dplyr::filter(season == 'summer') %>%
  dplyr::left_join(plot_treatments_summer, by = c('year', 'plot'))

winter <- all_plants %>% 
  dplyr::filter(season == 'winter') %>%
  dplyr::left_join(plot_treatments_winter, by = c('year', 'plot'))

write.csv(summer, 'data/summer-plants-raw.csv', row.names = F)
write.csv(winter, 'data/winter-plants-raw.csv', row.names = F)

}

store_plant_species_info <- function() {
  species_data <- download.file(url = 'https://raw.githubusercontent.com/weecology/PortalData/master/Plants/Portal_plant_species.csv',
                                destfile = 'data/plant_species_data.csv')
}