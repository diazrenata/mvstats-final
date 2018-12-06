# Store rodent data, filtered to years with plant data, and standardized to reflect trapping effort
# @param rodent_data rodent data table saved by store_rodent_data()
# @param plant_data table of all plant data saved by store_plant_data(),
# used to identify which years to use

store_adjusted_rodent_data <- function(rodent_data, plant_data){ 
 # Filter rodents to years with plant data
   rodent_data <- rodent_data %>%
    select(-month) %>%
    filter(year %in% plant_data$year)
  
   # Get trapping effort (per period per treatment)
  rodent_data_effort <- rodent_data %>%
    select(year, period, plot, treatment, sampled) %>%
    distinct() %>%
    group_by(year, period, treatment) %>%
    summarize(effort = sum(sampled)) %>%
    ungroup()
  
  # Standardize rodent data according to trapping effort
  # Summarize abundances as yearly sums
  rodent_data_adjusted <- rodent_data %>%
    select(-sampled) %>%
    group_by(year, period, treatment, species) %>%
    summarize(n = sum(n)) %>% 
    ungroup() %>%
    left_join(rodent_data_effort, 
              by = c("year", "period", "treatment")) %>%
    mutate(adjusted_n =as.integer(ceiling((n / effort) * 8))) %>%
    select(year, species, treatment, adjusted_n) %>%
    group_by(year, species, treatment) %>%
    summarize(adjusted_total = sum(adjusted_n)) %>%
    ungroup() %>%
    tidyr::spread('species', 'adjusted_total', fill = 0) %>%
    filter(treatment %in% c('control', 'exclosure'))
  
  write.csv(rodent_data_adjusted,'data/rodents-adjusted.csv', row.names = F)
  }
  
# Function to store plant data, adjusted to reflect census effort.
# @param plant_data plant data table as saved by store_plant_data(). 
store_adjusted_plant_data <- function(plant_data) {

  # Get sample effort per treatment per season per year
  plant_data_effort <- plant_data %>%
    select(year, plot, season, treatment, quads) %>%
    distinct() %>%
    group_by(year, season, treatment) %>%
    summarize(total_quads = sum(quads)) %>%
    ungroup()
  
  # Standardize species abundances
  plant_data_adjusted <- plant_data %>%
    group_by(year, species, treatment, season) %>%
    summarize(total_n = sum(abundance)) %>%
    ungroup() %>%
    left_join(plant_data_effort, by = c('year', 'season', 'treatment')) %>%
    mutate(adjusted_n = as.integer(ceiling((total_n / total_quads) * 64))) %>%
    select(year, species, treatment, season, adjusted_n) %>%
    tidyr::spread(species, adjusted_n, fill = 0) %>%
    filter(treatment %in% c('exclosure', 'control'))
  
  write.csv(plant_data_adjusted, 'data/plants-adjusted.csv', row.names = F)
  }
  