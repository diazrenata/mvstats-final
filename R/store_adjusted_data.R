
store_adjusted_rodent_data <- function(rodent_data, plant_data){ 
  rodent_data <- rodent_data %>%
    select(-month) %>%
    filter(year %in% plant_data$year)
  
  rodent_data_effort <- rodent_data %>%
    select(year, period, plot, treatment, sampled) %>%
    distinct() %>%
    group_by(year, period, treatment) %>%
    summarize(effort = sum(sampled)) %>%
    ungroup()
  
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
  
store_adjusted_plant_data <- function(plant_data) {

  plant_data_effort <- plant_data %>%
    select(year, plot, season, treatment, quads) %>%
    distinct() %>%
    group_by(year, season, treatment) %>%
    summarize(total_quads = sum(quads)) %>%
    ungroup()
  
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
  