
### Years - rodents comparison

```{r plot winter pcoa1 v year}


plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
abline(v = 1990, col = 'red')

colnames(rodents)

rodents_props <- rodents/rowSums(rodents)


plot(pred_vals$year, rodents_props$DS, col = 'blue')
points(pred_vals$year, rodents_props$PE, col = 'green')
points(pred_vals$year, rodents_props$RM, col = 'purple')
points(pred_vals$year, rodents_props$DM, col = 'pink')
abline(v = 1990, col = 'red')

```


### Explore the PCoA axes

```{r explore axes} 
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

## Does seed mass predict scores?

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

# Seed size predicts the score for summer axis 1, but not the other axes. 

#### Rank scores

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

par(mfrow = c(1,3))
plot(summer_scores_ranks$axis1_rank[which(summer_scores_ranks$axis1_rank<5)], summer_scores_ranks$saxis1[which(summer_scores_ranks$axis1_rank<5)], type = "n", xlim = c(0, 5.5))
text(summer_scores_ranks$axis1_rank[which(summer_scores_ranks$axis1_rank<5)], summer_scores_ranks$saxis1[which(summer_scores_ranks$axis1_rank<5)], labels = summer_scores_ranks$speciescode[which(summer_scores_ranks$axis1_rank<5)])

## Ordiplots

ordiplot(summer_scores[, c(2, 3)], type = "n", cex = 1, main = "Summer plant PCoA")

abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

text(summer_scores[, c(2, 3)], summer_scores$speciescode, cex = 0.7)



ordiplot(winter_scores[, c(2, 3)], type = "n", cex = 1, main = "Winter plant PCoA")

abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

text(winter_scores[, c(2, 3)], winter_scores$speciescode, cex = 0.7)


### Axes through time

all_axes <- inner_join(summer_axes, winter_axes, by = 'year')

axismin = min(all_axes[,2:7])
axismax = max(all_axes[,2:7])

axis_cols <- viridis::viridis(6)
cpt_years <- c(1984, 1990, 1999, 2010)
plot(all_axes$year, all_axes$SummerPCoAxis_1, type ='n', ylim = c(axismin, axismax))
for(i in 1:6) {
  lines(all_axes$year, all_axes[,i+1], col = axis_cols[i])
}
for(i in 1:4) {
  abline(v = cpt_years[i], col = 'red')
}
legend(x = 'bottomright', legend = colnames(all_axes)[2:7],fill = axis_cols,
       cex = 0.5)


plot(all_axes$year, all_axes$WinterPCoAxis_1, type = 'l', col = axis_cols[4])
abline(v = cpt_years[2], col = 'red')


rm(list=ls())

```
