Multivar stats final project
================
Renata Diaz
11/23/2018

Get data
--------

-   Download rodent and plant data from the Portal Project repo:

``` r
source('R/store_rodent_data.R')
source('R/store_plant_data.R')

store_rodent_data()
```

    ## Loading in data version 1.66.0

``` r
store_plant_data()
```

    ## Loading in data version 1.66.0
    ## Loading in data version 1.66.0

``` r
store_plant_species_info()
```

Adjust data to compensate for irregular trapping and census effort.

### Plant data

-   Extracted summer & winter plant censuses for all years
-   Standardized plant abundances according to sampling effort
-   Kept seasons separate

### Rodent data

-   Rodent censuses on control plots for all years, restricted to granivores
-   Standardized according to sampling effort (per census period)
-   Summed across all months in each calendar year

``` r
rodents <- read.csv('data/rodents-raw.csv', stringsAsFactors = F)
winter_plants <-read.csv('data/winter-plants-raw.csv', stringsAsFactors = F)
summer_plants <- read.csv('data/summer-plants-raw.csv', stringsAsFactors = F)
all_plants <- rbind(winter_plants, summer_plants)

source('R/store_adjusted_data.R')
store_adjusted_rodent_data(rodents, all_plants)
store_adjusted_plant_data(all_plants)

rm(list=ls())
```

PCoA
----

-   Transformed raw abundance values using the Wisconsin transformation, and then created a Bray-Curtis dissimilarity matrix for each community.
-   Ran PCoA on the dissimilarity matrices, separately for each community.

### Summer PCoA

``` r
all_plants <- read.csv('data/plants-adjusted.csv',
                          stringsAsFactors = F)
summer_plants <- all_plants %>%
  filter(season == 'summer', treatment == 'control') %>%
  select(-season, -treatment)

summer_plants_wis <- vegan::wisconsin(summer_plants[,2:ncol(summer_plants)])

summer_dist_mat <- vegdist(summer_plants_wis, 'bray')
  
summer_pcoa <- cmdscale(summer_dist_mat, k = nrow(summer_plants_wis) - 1, eig = T)
```

    ## Warning in cmdscale(summer_dist_mat, k = nrow(summer_plants_wis) - 1, eig =
    ## T): only 29 of the first 33 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  Summer_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.0448708 0.17950879 0.1795088
    ##  [2,]   1.4511379 0.12738800 0.3068968
    ##  [3,]   0.8779779 0.07707320 0.3839700
    ##  [4,]   0.8292255 0.07279348 0.4567635
    ##  [5,]   0.7310681 0.06417675 0.5209402
    ##  [6,]   0.5734561 0.05034079 0.5712810
    ##  [7,]   0.5363604 0.04708435 0.6183654
    ##  [8,]   0.4492451 0.03943694 0.6578023
    ##  [9,]   0.4385589 0.03849885 0.6963012
    ## [10,]   0.4055310 0.03559950 0.7319007
    ## [11,]   0.3566596 0.03130933 0.7632100
    ## [12,]   0.3444544 0.03023789 0.7934479
    ## [13,]   0.3235685 0.02840443 0.8218523
    ## [14,]   0.2619292 0.02299343 0.8448457
    ## [15,]   0.2480443 0.02177455 0.8666203

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](narrative_files/figure-markdown_github/summer%20pcoa-1.png)

``` r
#ordiplot(scores(summer_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Summer plant PCoA")
## species scores not available
#abline(h = 0, lty = 3)
#abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(summer_pcoa$points[, 1:3], summer_plants_wis)
#text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
```

There seems to be an inflection point in the scree plot around axis 3.

Moving forward, keeping the first 3 axes as predictor variables for the rodent community.

### Winter PCoA

``` r
winter_plants <- all_plants %>%
  filter(season == 'winter', treatment == 'control') %>%
  select(-season, -treatment)

winter_plants_wis <- vegan::wisconsin(winter_plants[,2:ncol(winter_plants)])


winter_dist_mat <- vegdist(winter_plants_wis, 'bray')
  
winter_pcoa <- cmdscale(winter_dist_mat, k = nrow(winter_plants_wis) - 1, eig = T)
```

    ## Warning in cmdscale(winter_dist_mat, k = nrow(winter_plants_wis) - 1, eig =
    ## T): only 32 of the first 34 eigenvalues are > 0

``` r
# Proportion of variance table 
  eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  winter_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.1552658 0.17530090 0.1753009
    ##  [2,]   1.5097990 0.12280115 0.2981021
    ##  [3,]   1.3006726 0.10579163 0.4038937
    ##  [4,]   0.8819023 0.07173049 0.4756242
    ##  [5,]   0.7334878 0.05965903 0.5352832
    ##  [6,]   0.6588629 0.05358933 0.5888725
    ##  [7,]   0.5772049 0.04694759 0.6358201
    ##  [8,]   0.5610690 0.04563517 0.6814553
    ##  [9,]   0.4313966 0.03508812 0.7165434
    ## [10,]   0.3932881 0.03198852 0.7485319
    ## [11,]   0.3706911 0.03015056 0.7786825
    ## [12,]   0.3313629 0.02695176 0.8056343
    ## [13,]   0.2935723 0.02387802 0.8295123
    ## [14,]   0.2930938 0.02383910 0.8533514
    ## [15,]   0.2405956 0.01956911 0.8729205

``` r
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
```

![](narrative_files/figure-markdown_github/winter%20pcoa-1.png)

``` r
#ordiplot(scores(winter_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Winter plant PCoA")
## species scores not available
#abline(h = 0, lty = 3)
#abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(winter_pcoa$points[, 1:3], winter_plants_wis)
# text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
```

There seems to be an inflection point in the scree plot around axis 2 or 3. Since stopping at 2 would only capture 35% of variation, going to go for 3.

### Explore the PCoA axes

``` r
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

# Compile species info & load

source('R/compile_species_info.R')
compile_species_info()

plant_info <- read.csv('data/compiled_plant_species_data.csv', stringsAsFactors = F)

plant_info <- select(plant_info, speciescode, community, seed_mass)

winter_scores <- left_join(winter_scores, plant_info, by = 'speciescode')

summer_scores <- left_join(summer_scores, plant_info, by = 'speciescode')

## Does seed mass predict scores

par(mfrow = c(2, 3))
plot(summer_scores$seed_mass, summer_scores$saxis1)
abline(lm(summer_scores$saxis1~summer_scores$seed_mass), col="red")
plot(summer_scores$seed_mass, summer_scores$saxis2)
abline(lm(summer_scores$saxis2~summer_scores$seed_mass), col="red")
plot(summer_scores$seed_mass, summer_scores$saxis3)
abline(lm(summer_scores$saxis3~summer_scores$seed_mass), col="red")

plot(winter_scores$seed_mass, winter_scores$waxis1)
abline(lm(winter_scores$waxis1~winter_scores$seed_mass), col="red")
plot(winter_scores$seed_mass, winter_scores$waxis2)
abline(lm(winter_scores$waxis2~winter_scores$seed_mass), col="red")
plot(winter_scores$seed_mass, winter_scores$waxis3)
abline(lm(winter_scores$waxis3~winter_scores$seed_mass), col="red")
```

![](narrative_files/figure-markdown_github/explore%20axes-1.png)

``` r
source('R/get_seed_mass_regression.R')

for(i in 1:3) {
this_regression <- get_seed_mass_regression(summer_scores, i + 1)
  print(this_regression)
}
```

    ## [[1]]
    ##                Estimate Std. Error   t value  Pr(>|t|)
    ## (Intercept) -0.12942385 0.11106864 -1.165260 0.2820807
    ## sizes        0.09469495 0.06142757  1.541571 0.1670801
    ## 
    ## [[2]]
    ## [1] 0.2534481
    ## 
    ## [[1]]
    ##                Estimate Std. Error   t value  Pr(>|t|)
    ## (Intercept)  0.09181182 0.08669055  1.059075 0.3247363
    ## sizes       -0.08398597 0.04794504 -1.751714 0.1232844
    ## 
    ## [[2]]
    ## [1] 0.3047624
    ## 
    ## [[1]]
    ##                Estimate Std. Error   t value  Pr(>|t|)
    ## (Intercept) -0.07315974 0.07033771 -1.040121 0.3328700
    ## sizes        0.04211157 0.03890094  1.082533 0.3148903
    ## 
    ## [[2]]
    ## [1] 0.1434038

``` r
for(i in 1:3) {
this_regression <- get_seed_mass_regression(winter_scores, i + 1)
  print(this_regression)
}
```

    ## [[1]]
    ##               Estimate Std. Error   t value     Pr(>|t|)
    ## (Intercept) -0.1930458 0.03529711 -5.469167 9.785684e-06
    ## sizes        0.0238644 0.01651142  1.445327 1.603122e-01
    ## 
    ## [[2]]
    ## [1] 0.07436978
    ## 
    ## [[1]]
    ##                Estimate Std. Error    t value  Pr(>|t|)
    ## (Intercept) -0.03510730 0.03590929 -0.9776662 0.3372492
    ## sizes        0.02508994 0.01679779  1.4936456 0.1473059
    ## 
    ## [[2]]
    ## [1] 0.07902585
    ## 
    ## [[1]]
    ##                Estimate Std. Error    t value  Pr(>|t|)
    ## (Intercept)  0.02651648 0.03483549  0.7611916 0.4533918
    ## sizes       -0.01956670 0.01629548 -1.2007438 0.2406761
    ## 
    ## [[2]]
    ## [1] 0.05253979

``` r
# Insomuch as we have it, seed size data doesn't predict axis scores

# for ordiplots, add colors to scores

community_colors <- summer_scores %>%
  select(community) %>%
  distinct()

community_colors$colour <- viridis::plasma(n = 7)

winter_scores <- left_join(winter_scores, community_colors, by = 'community')

summer_scores <- left_join(summer_scores, community_colors, by = 'community')

par(mfrow = c(1,1))

## Ordiplots

ordiplot(summer_scores[, c(2, 3)], type = "n", cex = 1, main = "Summer plant PCoA")
```

    ## species scores not available

``` r
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

 text(summer_scores[, c(2, 3)], summer_scores$speciescode, cex = 0.7, col = summer_scores$colour)
 legend(x = 'bottomright', legend = community_colors$community,fill = community_colors$colour,
        cex = 0.5)
```

![](narrative_files/figure-markdown_github/explore%20axes-2.png)

``` r
ordiplot(winter_scores[, c(2, 3)], type = "n", cex = 1, main = "Winter plant PCoA")
```

    ## species scores not available

``` r
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

 text(winter_scores[, c(2, 3)], winter_scores$speciescode, cex = 0.7, col = winter_scores$colour)
 legend(x = 'bottomright', legend = community_colors$community,fill = community_colors$colour,
        cex = 0.5)
```

![](narrative_files/figure-markdown_github/explore%20axes-3.png)

``` r
 ### Axes through time
 
all_scores <- inner_join(summer_axes, winter_axes, by = 'year')

axismin = min(all_scores[,2:7])
axismax = max(all_scores[,2:7])

axis_cols <- viridis::viridis(6)
cpt_years <- c(1984, 1990, 1999, 2010)
plot(all_scores$year, all_scores$SummerPCoAxis_1, type ='n', ylim = c(axismin, axismax))
for(i in 1:6) {
  lines(all_scores$year, all_scores[,i+1], col = axis_cols[i])
}
for(i in 1:4) {
  abline(v = cpt_years[i], col = 'red')
}
legend(x = 'bottomright', legend = colnames(all_scores)[2:7],fill = axis_cols,
        cex = 0.5)
```

![](narrative_files/figure-markdown_github/explore%20axes-4.png)

``` r
plot(all_scores$year, all_scores$WinterPCoAxis_1, type = 'l', col = axis_cols[4])
abline(v = cpt_years[2], col = 'red')
```

![](narrative_files/figure-markdown_github/explore%20axes-5.png)
