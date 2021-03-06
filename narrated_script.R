## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(vegan)


## ----download data, eval = F---------------------------------------------
## source('R/store_rodent_data.R')
## source('R/store_plant_data.R')
## 
## store_rodent_data()
## store_plant_data()
## store_plant_species_info()

## ----adjust data, eval = F-----------------------------------------------
## rodents <- read.csv('data/rodents-raw.csv', stringsAsFactors = F)
## winter_plants <-read.csv('data/winter-plants-raw.csv', stringsAsFactors = F)
## summer_plants <- read.csv('data/summer-plants-raw.csv', stringsAsFactors = F)
## all_plants <- rbind(winter_plants, summer_plants)
## 
## source('R/store_adjusted_data.R')
## store_adjusted_rodent_data(rodents, all_plants)
## store_adjusted_plant_data(all_plants)
## 
## rm(list=ls())

## ----summer pcoa, echo = T-----------------------------------------------

# Load plant data

all_plants <- read.csv('data/plants-adjusted.csv',
                       stringsAsFactors = F)

# Filter to summer plants on control plots
summer_plants <- all_plants %>%
  filter(season == 'summer', treatment == 'control') %>%
  select(-season, -treatment)

# Remove species not present in the summer
summer_plants <- summer_plants[ , which(colSums(summer_plants) > 0)]

# Wisconsin transform
summer_plants_wis <- vegan::wisconsin(summer_plants[,2:ncol(summer_plants)])

# Bray-Curtis dissimilarity matrix
summer_dist_mat <- vegdist(summer_plants_wis, 'bray')

# PCoA
summer_pcoa <- cmdscale(summer_dist_mat, k = nrow(summer_plants_wis) - 1, eig = T)

# Generate proportion of variance table 
eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants_wis)-1]
propVar <- eigenvalues/sum(eigenvalues)
cumVar <- cumsum(propVar)
Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)


Summer_PCoA_Table[1:15,]


# Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))

# Calculate and save species scores (to look at later)
species_pc <- wascores(summer_pcoa$points[, 1:3], summer_plants_wis)
write.csv(species_pc, 'models/summer_species_scores.csv', row.names = T)


## ----save summer axes, echo = F------------------------------------------
# Save summer axes
summer_vals <- cbind(summer_plants$year, summer_pcoa$points[,1:3])

season_names <- rep("SummerPCoAxis_", 3) %>%
  paste0(1:3)
season_names <- c('year', season_names)
colnames(summer_vals) <- season_names


write.csv(summer_vals, 'models/summer_pcoa_axes.csv', row.names = F)


## ----winter pcoa, echo = T-----------------------------------------------
# Get winter control plants from whole plant dataset
winter_plants <- all_plants %>%
  filter(season == 'winter', treatment == 'control') %>%
  select(-season, -treatment)

# Remove species not present in the winter
winter_plants <- winter_plants[ , which(colSums(winter_plants) > 0)]

# Wisconsin transform
winter_plants_wis <- vegan::wisconsin(winter_plants[,2:ncol(winter_plants)])

# Bray-Curtis matrix
winter_dist_mat <- vegdist(winter_plants_wis, 'bray')

# PCoA
winter_pcoa <- cmdscale(winter_dist_mat, k = nrow(winter_plants_wis) - 1, eig = T)

# Proportion of variance table 
eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants_wis)-1]
propVar <- eigenvalues/sum(eigenvalues)
cumVar <- cumsum(propVar)
winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)

winter_PCoA_Table[1:15,]


# Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))


# Calculate and save winter species scores
species_pc <- wascores(winter_pcoa$points[, 1:3], winter_plants_wis)

write.csv(species_pc, 'models/winter_species_scores.csv', row.names = T)

# Save PCoA tables for paper

write.csv(as.data.frame(Summer_PCoA_Table), 'models/summer_table.csv')
write.csv(as.data.frame(winter_PCoA_Table), 'models/winter_table.csv')


## ----explore axes--------------------------------------------------------
# Load all axes

summer_axes <- read.csv('models/summer_pcoa_axes.csv', stringsAsFactors = F)
winter_axes <- read.csv('models/winter_pcoa_axes.csv', stringsAsFactors = F)
summer_scores <- read.csv('models/summer_species_scores.csv',
                          stringsAsFactors = F)
winter_scores <- read.csv('models/winter_species_scores.csv',
                          stringsAsFactors = F)

colnames(winter_scores) <- c('speciescode', 'waxis1', 
                             'waxis2', 'waxis3')
colnames(summer_scores) <- c('speciescode', 'saxis1', 
                             'saxis2', 'saxis3')

# Compile species seed mass info & load

# # Do not need to run - file already exists
# source('R/compile_species_info.R') 
# compile_species_info()

plant_info <- read.csv('data/compiled_plant_species_data.csv', stringsAsFactors = F)

plant_info <- select(plant_info, speciescode, community, seed_mass)

winter_scores <- left_join(winter_scores, plant_info, by = 'speciescode')

summer_scores <- left_join(summer_scores, plant_info, by = 'speciescode')


## ----seed size-----------------------------------------------------------

# Source function to run a regression of species score against seed mass
source('R/get_seed_mass_regression.R')

# Get p-values for all regressions
regression_ps <- vector(length = 6)

for(i in 1:3) {
  this_regression <- get_seed_mass_regression(summer_scores, i + 1)
  regression_ps[i] <-  this_regression[[1]][2,4]
}


for(i in 1:3) {
  this_regression <- get_seed_mass_regression(winter_scores, i + 1)
  regression_ps[i+3] <-  this_regression[[1]][2,4]
}


## ----plot seed size, echo = F--------------------------------------------
par(mfrow = c(2, 3))
plot(summer_scores$seed_mass, summer_scores$saxis1)
abline(lm(summer_scores$saxis1~summer_scores$seed_mass), col="red")
text(3, -0.3, label = paste0('p = ', substr(as.character(regression_ps[1]), 0, 5)))

plot(summer_scores$seed_mass, summer_scores$saxis2)
abline(lm(summer_scores$saxis2~summer_scores$seed_mass), col="red")
text(3, -0.3, label = paste0('p = ', substr(as.character(regression_ps[2]), 0, 5)))

plot(summer_scores$seed_mass, summer_scores$saxis3)
abline(lm(summer_scores$saxis3~summer_scores$seed_mass), col="red")
text(3, -0.3, label = paste0('p = ', substr(as.character(regression_ps[3]), 0, 5)))

plot(winter_scores$seed_mass, winter_scores$waxis1)
abline(lm(winter_scores$waxis1~winter_scores$seed_mass), col="red")
text(5, 0.2, label = paste0('p = ', substr(as.character(regression_ps[4]), 0, 5)))

plot(winter_scores$seed_mass, winter_scores$waxis2)
abline(lm(winter_scores$waxis2~winter_scores$seed_mass), col="red")
text(5, 0.2, label = paste0('p = ', substr(as.character(regression_ps[5]), 0, 5)))

plot(winter_scores$seed_mass, winter_scores$waxis3)
abline(lm(winter_scores$waxis3~winter_scores$seed_mass), col="red")
text(5, -0.3, label = paste0('p = ', substr(as.character(regression_ps[6]), 0, 5)))


## ----rank scores---------------------------------------------------------

summer_scores_ranks <- summer_scores %>%
  select(speciescode, saxis1, saxis2, saxis3) %>%
  mutate(abs_score = abs(saxis1)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis1_rank = row_number()) %>%
  mutate(abs_score = abs(saxis2)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis2_rank = row_number()) %>%
  mutate(abs_score = abs(saxis3)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis3_rank = row_number()) %>%
  select(-abs_score)
winter_scores_ranks <- winter_scores %>%
  select(speciescode, waxis1, waxis2, waxis3) %>%
  mutate(abs_score = abs(waxis1)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis1_rank = row_number()) %>%
  mutate(abs_score = abs(waxis2)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis2_rank = row_number()) %>%
  mutate(abs_score = abs(waxis3)) %>%
  arrange(desc(abs_score)) %>%
  mutate(axis3_rank = row_number()) %>%
  select(-abs_score)

filter(summer_scores_ranks, speciescode == 'erod.cicu')
filter(winter_scores_ranks, speciescode == 'erod.cicu')


## ----save winter axes, echo = F------------------------------------------
# Save winter axes
winter_vals <- cbind(winter_plants$year, winter_pcoa$points[,1:3])

season_names <- rep("WinterPCoAxis_", 3) %>%
  paste0(1:3)
season_names <- c('year', season_names)
colnames(winter_vals) <- season_names


write.csv(winter_vals, 'models/winter_pcoa_axes.csv', row.names = F)


## ----load data, echo = T-------------------------------------------------
# Load rodent data, winter axes, and summer axes
rodents <- read.csv('data/rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('models/summer_pcoa_axes.csv', 
                        stringsAsFactors = F)

winter_axes <- read.csv('models/winter_pcoa_axes.csv',
                        stringsAsFactors = F)

# Join axes & year to make a predictor values table
pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

# Filter rodents to control plots & year with plant data
rodents <- filter(rodents, year %in% pred_vals$year, 
                  treatment == 'control') %>%
  select(-year, -treatment)

# Hellinger transformation on rodent data
rodents_hel <- decostand(rodents, 'hellinger')


## ----partial RDA---------------------------------------------------------
# Partial RDA with all axes, conditioned on year
rodents_prda <- rda(rodents_hel ~ WinterPCoAxis_1 + 
                      WinterPCoAxis_2 + 
                      WinterPCoAxis_3 +
                      SummerPCoAxis_1 + 
                      SummerPCoAxis_2 + 
                      SummerPCoAxis_3 +
                      Condition(year), pred_vals)

# Evaluate R2 and significance of global model (via permutation test)
R2 <- RsquareAdj(rodents_prda)$r.squared
R2adj <- RsquareAdj(rodents_prda)$adj.r.squared

R2
R2adj


anova(rodents_prda, step = 1000)
anova(rodents_prda, by = "axis", step = 1000)


# Find the most parsimonious model:
set.seed(1)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_prda), 
                           R2scope = F, direction = "forward", pstep = 1000)

# Most parsimonious is Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      SummerPCoAxis_2

# Evaluate the R2 and significance of the most parsimonious model

rod_prda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2 + Condition(year), pred_vals)

pR2p <- RsquareAdj(rod_prda_pars)$r.squared
pR2adjp <- RsquareAdj(rod_prda_pars)$adj.r.squared

pR2p
pR2adjp


anova(rod_prda_pars, step = 1000)
anova(rod_prda_pars, by = "axis", step = 1000)


## ----partial variance partitioning, echo = T-----------------------------
# Run variance partitioning on significant axes of the most parsimonious RDA
rod_ppart <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~  WinterPCoAxis_3, ~ year, data = pred_vals)
rod_ppart

plot(rod_ppart, digits = 2, Xnames = c('Winter1', 'Winter3', 'year'))
title(main="Variance partitioning - Partial RDA")



## ----regular rda, echo = T-----------------------------------------------
# RDA on all variables (no condition)
rodents_rda <- rda(rodents_hel ~ ., pred_vals)

# R2 and significance test for global model
R2 <- RsquareAdj(rodents_rda)$r.squared
R2adj <- RsquareAdj(rodents_rda)$adj.r.squared
R2
R2adj
anova(rodents_rda, step = 1000)
anova(rodents_rda, by = "axis", step = 1000)

# Find the most parsimonious model
set.seed(1)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_rda), 
                           R2scope = F, direction = "forward", pstep = 1000)

# Most parsimonious is rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +  year 

rod_rda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + year, pred_vals)

# R2 and significance for most parsimonious model
R2p <- RsquareAdj(rod_rda_pars)$r.squared
R2adjp <- RsquareAdj(rod_rda_pars)$adj.r.squared
R2p
R2adjp
anova(rod_rda_pars, step = 1000)
anova(rod_rda_pars, by = "axis", step = 1000)

# The first 3 axes are significant (winter 1, winter 3, winter 2)


## ----regular variance partitioning, echo = T-----------------------------
rod_part <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~ WinterPCoAxis_3, ~  WinterPCoAxis_2, data = pred_vals)
rod_part

plot(rod_part, digits = 2, Xnames = c('Winter1', 'Winter3','Winter2'))
title(main="Variance partitioning - unconditioned RDA")


## ----plot winter pcoa1 v year--------------------------------------------


par(mfrow = c(1, 2))
rodents_props <- rodents/rowSums(rodents)

plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
title(main="Winter axis 1")
abline(v = 1990, col = 'red')

plot(pred_vals$year, rodents_props$DS, col = 'blue')
title(main = 'D. spectabilis % abundance')
abline(v = 1990, col = 'red')



# Important species in axis 1
# Save the ranked list of highest absolute value-scoring species for 
# winter axis 1. 
impt_species <- winter_scores_ranks %>%
  select(speciescode, axis1_rank, waxis1) %>%
  arrange(axis1_rank)

write.csv(impt_species, 'models/winter_axis_1.csv')

