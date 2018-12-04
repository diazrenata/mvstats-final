#' @importFrom magrittr "%>%"
#' @title store portal plant dataframes
#'
#' Import Portal rodent data using portalr functions.
#'
#' @export

library(dplyr)
store_rodent_data <- function(){
  primary_tables <- portalr::load_data(path = 'repo')
  
  rodents <- primary_tables[[1]]
  species <- primary_tables[[2]]
  trapping <- primary_tables[[3]]
  newmoons <- primary_tables[[4]]
  plots <- primary_tables[[5]]
  
  plots_history <- left_join(plots, trapping, 
                     by = c('year', 'month', 'plot')) 
  plots_history <- plots_history %>%
    select(year, month, plot, treatment, period, sampled)
  
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