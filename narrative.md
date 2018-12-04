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

``` r
rm(list=ls())
```

Partial RDA
-----------

Partial redundancy analysis, using combined winter and summer axes, conditioned on year, to predict the rodent community.

Restricted to years with both a winter & summer census (n = 27).

Used rodent data, summarized yearly, transformed via Hellinger transformation.

``` r
rodents <- read.csv('data/rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('models/summer_pcoa_axes.csv', 
                        stringsAsFactors = F)

winter_axes <- read.csv('models/winter_pcoa_axes.csv',
                        stringsAsFactors = F)


pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

rodents <- filter(rodents, year %in% pred_vals$year, 
                  treatment == 'control') %>%
  select(-year, -treatment)

pred_vals_noy <- select(pred_vals, -year)
pred_vals_y <- select(pred_vals, year)

rodents_hel <- decostand(rodents, 'hellinger')

rodents_prda <- rda(rodents_hel ~ . + Condition(pred_vals_y$year), pred_vals_noy)

R2 <- RsquareAdj(rodents_prda)$r.squared
R2adj <- RsquareAdj(rodents_prda)$adj.r.squared

R2
```

    ## [1] 0.187431

``` r
R2adj
```

    ## [1] NA

``` r
anova(rodents_prda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)  
    ## Model     6 0.032671 1.8312  0.026 *
    ## Residual 25 0.074337                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rodents_prda, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)
    ## RDA1      1 0.014721 4.9507  0.175
    ## RDA2      1 0.010120 3.4034  0.308
    ## RDA3      1 0.004418 1.4858  0.791
    ## RDA4      1 0.002593 0.8722  0.917
    ## RDA5      1 0.000606 0.2038  1.000
    ## RDA6      1 0.000213 0.0715  0.999
    ## Residual 25 0.074337

Find the most parsimonious model...

``` r
set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_prda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                    R2.adjusted
    ## + WinterPCoAxis_1  0.173825176
    ## + SummerPCoAxis_1  0.133960638
    ## + SummerPCoAxis_2  0.132569805
    ## + WinterPCoAxis_2  0.114392435
    ## + WinterPCoAxis_3  0.018292889
    ## <none>             0.000000000
    ## + SummerPCoAxis_3 -0.009202931
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_1  1 -62.013 7.7327  0.004 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.1738252 
    ## Call: rodents_hel ~ WinterPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_2   0.3062469
    ## + SummerPCoAxis_2   0.2720064
    ## + WinterPCoAxis_3   0.2012075
    ## + SummerPCoAxis_1   0.1854104
    ## + SummerPCoAxis_3   0.1813034
    ## <none>              0.1738252
    ## 
    ##                   Df    AIC      F Pr(>F)   
    ## + WinterPCoAxis_2  1 -66.86 6.9172  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.3062469 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_2   0.3674221
    ## + WinterPCoAxis_3   0.3498341
    ## + SummerPCoAxis_1   0.3384797
    ## + SummerPCoAxis_3   0.3077858
    ## <none>              0.3062469
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + SummerPCoAxis_2  1 -69.025 3.9012  0.008 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.3674221 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + SummerPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_3   0.3776844
    ## + SummerPCoAxis_1   0.3740554
    ## + SummerPCoAxis_3   0.3714096
    ## <none>              0.3674221
    ## 
    ##                   Df     AIC      F Pr(>F)
    ## + WinterPCoAxis_3  1 -68.723 1.4782    0.2

``` r
# Most parsimonious is Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 

rod_prda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2  + Condition(pred_vals_y$year), pred_vals_noy)

pR2p <- RsquareAdj(rod_prda_pars)$r.squared
pR2adjp <- RsquareAdj(rod_prda_pars)$adj.r.squared

pR2p
```

    ## [1] 0.05830529

``` r
pR2adjp
```

    ## [1] NA

``` r
anova(rod_prda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)
    ## Model     2 0.010163 1.5217  0.163
    ## Residual 29 0.096845

``` r
anova(rod_prda_pars, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + Condition(pred_vals_y$year), data = pred_vals_noy)
    ##          Df Variance      F Pr(>F)
    ## RDA1      1 0.009390 2.8117  0.102
    ## RDA2      1 0.000773 0.2316  0.966
    ## Residual 29 0.096845

Variance partitioning
---------------------

``` r
rod_ppart <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~  WinterPCoAxis_2, ~ pred_vals_y$year, data = pred_vals_noy)
rod_ppart
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~WinterPCoAxis_1,
    ## ~WinterPCoAxis_2, ~pred_vals_y$year, data = pred_vals_noy)
    ## 
    ## Explanatory tables:
    ## X1:  ~WinterPCoAxis_1
    ## X2:  ~WinterPCoAxis_2
    ## X3:  ~pred_vals_y$year 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 5.5779 
    ##             Variance: 0.17431 
    ## No. of observations: 33 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1         1  0.19964      0.17383     TRUE
    ## [b+d+e+g] = X2         1  0.14207      0.11439     TRUE
    ## [c+e+f+g] = X3         1  0.38610      0.36630     TRUE
    ## [a+b+d+e+f+g] = X1+X2  2  0.34961      0.30625     TRUE
    ## [a+c+d+e+f+g] = X1+X3  2  0.41247      0.37330     TRUE
    ## [b+c+d+e+f+g] = X2+X3  2  0.39166      0.35110     TRUE
    ## [a+b+c+d+e+f+g] = All  3  0.44441      0.38693     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3       1               0.03583     TRUE
    ## [b] = X2 | X1+X3       1               0.01363     TRUE
    ## [c] = X3 | X1+X2       1               0.08069     TRUE
    ## [d]                    0              -0.02883    FALSE
    ## [e]                    0               0.11879    FALSE
    ## [f]                    0               0.15603    FALSE
    ## [g]                    0               0.01080    FALSE
    ## [h] = Residuals                        0.61307    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3        1               0.00700     TRUE
    ## [a+f] = X1 | X2        1               0.19185     TRUE
    ## [b+d] = X2 | X3        1              -0.01520     TRUE
    ## [b+e] = X2 | X1        1               0.13242     TRUE
    ## [c+e] = X3 | X1        1               0.19948     TRUE
    ## [c+f] = X3 | X2        1               0.23671     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_ppart, digits = 2)
```

![](narrative_files/figure-markdown_github/partial%20variance%20partitioning-1.png)

WinterPCoAxis\_1 combined with year has the largest chunk (.34); on its own, WinterPCoAxis\_1 explains an additional .1

### Years - rodents comparison

``` r
plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
abline(v = 1990, col = 'red')
```

![](narrative_files/figure-markdown_github/plot%20winter%20pcoa1%20v%20year-1.png)

``` r
colnames(rodents)
```

    ##  [1] "BA" "DM" "DO" "DS" "PB" "PE" "PF" "PH" "PI" "PL" "PM" "PP" "RF" "RM"
    ## [15] "RO"

``` r
rodents_props <- rodents/rowSums(rodents)


plot(pred_vals$year, rodents_props$DS, col = 'blue')
points(pred_vals$year, rodents_props$PE, col = 'green')
points(pred_vals$year, rodents_props$RM, col = 'purple')
points(pred_vals$year, rodents_props$DM, col = 'pink')
abline(v = 1990, col = 'red')
```

![](narrative_files/figure-markdown_github/plot%20winter%20pcoa1%20v%20year-2.png)

Regular RDA
-----------

``` r
rodents_rda <- rda(rodents_hel ~ ., pred_vals)

R2 <- RsquareAdj(rodents_rda)$r.squared
R2adj <- RsquareAdj(rodents_rda)$adj.r.squared

R2
```

    ## [1] 0.5735345

``` r
R2adj
```

    ## [1] 0.4541242

``` r
anova(rodents_rda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ year + WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     7 0.099973 4.8031  0.001 ***
    ## Residual 25 0.074337                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rodents_rda, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ year + WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3, data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.071273 23.9695  0.001 ***
    ## RDA2      1 0.011673  3.9258  0.277    
    ## RDA3      1 0.010046  3.3786  0.302    
    ## RDA4      1 0.003956  1.3305  0.845    
    ## RDA5      1 0.002476  0.8328  0.901    
    ## RDA6      1 0.000348  0.1172  1.000    
    ## RDA7      1 0.000199  0.0669  1.000    
    ## Residual 25 0.074337                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Find the most parsimonious model...

``` r
set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_rda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                    R2.adjusted
    ## + year             0.366300432
    ## + WinterPCoAxis_1  0.173825176
    ## + SummerPCoAxis_1  0.133960638
    ## + SummerPCoAxis_2  0.132569805
    ## + WinterPCoAxis_2  0.114392435
    ## + WinterPCoAxis_3  0.018292889
    ## <none>             0.000000000
    ## + SummerPCoAxis_3 -0.009202931
    ## 
    ##        Df     AIC      F Pr(>F)   
    ## + year  1 -70.766 19.497  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.3663004 
    ## Call: rodents_hel ~ year 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_3   0.3957853
    ## + SummerPCoAxis_2   0.3938742
    ## + SummerPCoAxis_1   0.3918242
    ## + SummerPCoAxis_3   0.3757239
    ## + WinterPCoAxis_1   0.3733018
    ## <none>              0.3663004
    ## + WinterPCoAxis_2   0.3511047
    ## 
    ##                   Df    AIC      F Pr(>F)  
    ## + WinterPCoAxis_3  1 -71.42 2.5128  0.052 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Most parsimonious is rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +  SummerPCoAxis_2 

rod_rda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2, pred_vals)

R2p <- RsquareAdj(rod_rda_pars)$r.squared
R2adjp <- RsquareAdj(rod_rda_pars)$adj.r.squared

R2p
```

    ## [1] 0.4554739

``` r
R2adjp
```

    ## [1] 0.3776844

``` r
anova(rod_rda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     4 0.079394 5.8552  0.001 ***
    ## Residual 28 0.094916                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rod_rda_pars, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2, data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.065300 19.2634  0.001 ***
    ## RDA2      1 0.009696  2.8604  0.197    
    ## RDA3      1 0.003906  1.1522  0.619    
    ## RDA4      1 0.000491  0.1449  0.995    
    ## Residual 28 0.094916                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
rod_part <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~ WinterPCoAxis_3, ~  WinterPCoAxis_2, ~ SummerPCoAxis_2, data = pred_vals)
rod_part
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~WinterPCoAxis_1,
    ## ~WinterPCoAxis_3, ~WinterPCoAxis_2, ~SummerPCoAxis_2, data =
    ## pred_vals)
    ## 
    ## Explanatory tables:
    ## X1:  ~WinterPCoAxis_1
    ## X2:  ~WinterPCoAxis_3
    ## X3:  ~WinterPCoAxis_2
    ## X4:  ~SummerPCoAxis_2 
    ## 
    ## No. of explanatory tables: 4 
    ## Total variation (SS): 5.5779 
    ##             Variance: 0.17431 
    ## No. of observations: 33 
    ## 
    ## Partition table:
    ##                             Df R.square Adj.R.square Testable
    ## [aeghklno] = X1              1  0.19964      0.17383     TRUE
    ## [befiklmo] = X2              1  0.04897      0.01829     TRUE
    ## [cfgjlmno] = X3              1  0.14207      0.11439     TRUE
    ## [dhijkmno] = X4              1  0.15968      0.13257     TRUE
    ## [abefghiklmno] = X1+X2       2  0.25113      0.20121     TRUE
    ## [acefghjklmno] = X1+X3       2  0.34961      0.30625     TRUE
    ## [adeghijklmno] = X1+X4       2  0.31751      0.27201     TRUE
    ## [bcefgijklmno] = X2+X3       2  0.19820      0.14475     TRUE
    ## [bdefhijklmno] = X2+X4       2  0.21116      0.15857     TRUE
    ## [cdfghijklmno] = X3+X4       2  0.25048      0.20052     TRUE
    ## [abcefghijklmno] = X1+X2+X3  3  0.41079      0.34983     TRUE
    ## [abdefghijklmno] = X1+X2+X4  3  0.34923      0.28191     TRUE
    ## [acdefghijklmno] = X1+X3+X4  3  0.42673      0.36742     TRUE
    ## [bcdefghijklmno] = X2+X3+X4  3  0.28683      0.21305     TRUE
    ## [abcdefghijklmno] = All      4  0.45547      0.37768     TRUE
    ## Individual fractions                                         
    ## [a] = X1 | X2+X3+X4          1               0.16463     TRUE
    ## [b] = X2 | X1+X3+X4          1               0.01026     TRUE
    ## [c] = X3 | X1+X2+X4          1               0.09577     TRUE
    ## [d] = X4 | X1+X2+X3          1               0.02785     TRUE
    ## [e]                          0               0.00227    FALSE
    ## [f]                          0              -0.00036    FALSE
    ## [g]                          0              -0.04129    FALSE
    ## [h]                          0               0.04045    FALSE
    ## [i]                          0               0.03332    FALSE
    ## [j]                          0               0.05285    FALSE
    ## [k]                          0              -0.01551    FALSE
    ## [l]                          0               0.01382    FALSE
    ## [m]                          0              -0.01585    FALSE
    ## [n]                          0               0.01912    FALSE
    ## [o]                          0              -0.00968    FALSE
    ## [p] = Residuals              0               0.62232    FALSE
    ## Controlling 2 tables X                                       
    ## [ae] = X1 | X3+X4            1               0.16690     TRUE
    ## [ag] = X1 | X2+X4            1               0.12334     TRUE
    ## [ah] = X1 | X2+X3            1               0.20509     TRUE
    ## [be] = X2 | X3+X4            1               0.01253     TRUE
    ## [bf] = X2 | X1+X4            1               0.00990     TRUE
    ## [bi] = X2 | X1+X3            1               0.04359     TRUE
    ## [cf] = X3 | X1+X4            1               0.09542     TRUE
    ## [cg] = X3 | X2+X4            1               0.05448     TRUE
    ## [cj] = X3 | X1+X2            1               0.14863     TRUE
    ## [dh] = X4 | X2+X3            1               0.06831     TRUE
    ## [di] = X4 | X1+X3            1               0.06118     TRUE
    ## [dj] = X4 | X1+X2            1               0.08070     TRUE
    ## Controlling 1 table X                                        
    ## [aghn] = X1 | X2             1               0.18291     TRUE
    ## [aehk] = X1 | X3             1               0.19185     TRUE
    ## [aegl] = X1 | X4             1               0.13944     TRUE
    ## [bfim] = X2 | X1             1               0.02738     TRUE
    ## [beik] = X2 | X3             1               0.03035     TRUE
    ## [befl] = X2 | X4             1               0.02600     TRUE
    ## [cfjm] = X3 | X1             1               0.13242     TRUE
    ## [cgjn] = X3 | X2             1               0.12645     TRUE
    ## [cfgl] = X3 | X4             1               0.06795     TRUE
    ## [dijm] = X4 | X1             1               0.09818     TRUE
    ## [dhjn] = X4 | X2             1               0.14028     TRUE
    ## [dhik] = X4 | X3             1               0.08612     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_part, digits = 2)
```

![](narrative_files/figure-markdown_github/regular%20variance%20partitioning-1.png)
