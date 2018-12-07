Multivar stats final project - script
================
Renata Diaz
12-05-2018

To access this code without the narrative of this document, download this repo as a .zip, set the working directory to the unzipped folder, open narrated\_script.R and run.

Download & process data
-----------------------

Download rodent and plant data from the Portal Project repo - note that functions rely on the portalr package, which is available on CRAN. I put all the raw data files in the data folder, so you don't need to run this step.

``` r
source('R/store_rodent_data.R')
source('R/store_plant_data.R')

store_rodent_data()
store_plant_data()
store_plant_species_info()
```

### Plant data

-   Extracted summer & winter plant censuses for all years
-   Restricted to winter and summer annuals
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

-   Transform raw abundance values using the Wisconsin transformation, and then created a Bray-Curtis dissimilarity matrix for each community.
-   Run PCoA on the dissimilarity matrices, separately for each community.

### Summer PCoA

``` r
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
```

    ## Warning in cmdscale(summer_dist_mat, k = nrow(summer_plants_wis) - 1, eig =
    ## T): only 26 of the first 31 eigenvalues are > 0

``` r
# Generate proportion of variance table 
eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants_wis)-1]
propVar <- eigenvalues/sum(eigenvalues)
cumVar <- cumsum(propVar)
Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)


Summer_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.0508707 0.19608563 0.1960856
    ##  [2,]   1.5954875 0.15254602 0.3486316
    ##  [3,]   1.1032454 0.10548231 0.4541140
    ##  [4,]   0.8942843 0.08550334 0.5396173
    ##  [5,]   0.7203322 0.06887162 0.6084889
    ##  [6,]   0.6103352 0.05835471 0.6668436
    ##  [7,]   0.5403830 0.05166652 0.7185101
    ##  [8,]   0.4545490 0.04345985 0.7619700
    ##  [9,]   0.3675808 0.03514474 0.7971147
    ## [10,]   0.3361137 0.03213614 0.8292509
    ## [11,]   0.2948953 0.02819521 0.8574461
    ## [12,]   0.2866688 0.02740866 0.8848547
    ## [13,]   0.2767467 0.02646001 0.9113147
    ## [14,]   0.2180342 0.02084645 0.9321612
    ## [15,]   0.1798866 0.01719913 0.9493603

``` r
# Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))
```

![](narrated_script_files/figure-markdown_github/summer%20pcoa-1.png)

``` r
# Calculate and save species scores (to look at later)
species_pc <- wascores(summer_pcoa$points[, 1:3], summer_plants_wis)
write.csv(species_pc, 'models/summer_species_scores.csv', row.names = T)
```

There seems to be an inflection point in the scree plot around axis 3. The first three axes describe a cumulative 45% of variation.

Moving forward, keeping the first 3 axes as predictor variables for the rodent community.

### Winter PCoA

``` r
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
```

    ## Warning in cmdscale(winter_dist_mat, k = nrow(winter_plants_wis) - 1, eig =
    ## T): only 24 of the first 30 eigenvalues are > 0

``` r
# Proportion of variance table 
eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants_wis)-1]
propVar <- eigenvalues/sum(eigenvalues)
cumVar <- cumsum(propVar)
winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)

winter_PCoA_Table[1:15,]
```

    ##       eigenvalues    propVar    cumVar
    ##  [1,]   2.4791232 0.22814312 0.2281431
    ##  [2,]   1.2204905 0.11231653 0.3404597
    ##  [3,]   1.1211349 0.10317326 0.4436329
    ##  [4,]   0.9427954 0.08676144 0.5303944
    ##  [5,]   0.6974427 0.06418268 0.5945770
    ##  [6,]   0.6037170 0.05555750 0.6501345
    ##  [7,]   0.5479375 0.05042435 0.7005589
    ##  [8,]   0.4965611 0.04569640 0.7462553
    ##  [9,]   0.4060608 0.03736804 0.7836233
    ## [10,]   0.3751394 0.03452248 0.8181458
    ## [11,]   0.3636433 0.03346454 0.8516103
    ## [12,]   0.3331781 0.03066095 0.8822713
    ## [13,]   0.2674315 0.02461058 0.9068819
    ## [14,]   0.2230801 0.02052911 0.9274110
    ## [15,]   0.1818489 0.01673478 0.9441458

``` r
# Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))
```

![](narrated_script_files/figure-markdown_github/winter%20pcoa-1.png)

``` r
# Calculate and save winter species scores
species_pc <- wascores(winter_pcoa$points[, 1:3], winter_plants_wis)

write.csv(species_pc, 'models/winter_species_scores.csv', row.names = T)

# Save PCoA tables for paper

write.csv(as.data.frame(Summer_PCoA_Table), 'models/summer_table.csv')
write.csv(as.data.frame(winter_PCoA_Table), 'models/winter_table.csv')
```

There seems to be an inflection point in the scree plot around axis 2 or 3. Three axes describe a cumulative 44% of variation, so moving forward with axes.

### Explore the PCoA axes

``` r
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
```

### Does seed mass predict species scores?

``` r
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
```

![](narrated_script_files/figure-markdown_github/plot%20seed%20size-1.png)

Seed size predicts the score for summer axis 1, but not the other axes.

### What is *E. ciculatum*'s score on each axis?

``` r
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
```

    ##   speciescode    saxis1     saxis2     saxis3 axis1_rank axis2_rank
    ## 1   erod.cicu 0.2403355 0.01631055 0.04430535         16         34
    ##   axis3_rank
    ## 1         34

``` r
filter(winter_scores_ranks, speciescode == 'erod.cicu')
```

    ##   speciescode    waxis1     waxis2     waxis3 axis1_rank axis2_rank
    ## 1   erod.cicu 0.2853572 0.03895151 -0.3226492         17         42
    ##   axis3_rank
    ## 1          2

*E. ciculatum* scores highly for Winter axis 3, and no other axis.

Partial RDA
-----------

Partial redundancy analysis, using combined winter and summer axes, conditioned on year, to predict the rodent community.

Restricted to years with both a winter & summer census (n = 27).

Used rodent data, summarized yearly, transformed via Hellinger transformation.

``` r
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
```

Partial RDA
===========

``` r
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
```

    ## [1] 0.3612209

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
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(year), data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     6 0.060217 4.7001  0.001 ***
    ## Residual 19 0.040571                  
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
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3 + Condition(year), data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.040712 19.0661  0.001 ***
    ## RDA2      1 0.012714  5.9542  0.006 ** 
    ## RDA3      1 0.003260  1.5268  0.797    
    ## RDA4      1 0.002051  0.9606  0.937    
    ## RDA5      1 0.001033  0.4840  0.986    
    ## RDA6      1 0.000446  0.2090  0.989    
    ## Residual 19 0.040571                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Find the most parsimonious model:
set.seed(1)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_prda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_1  0.46690627
    ## + SummerPCoAxis_1  0.22573181
    ## + WinterPCoAxis_3  0.06434779
    ## + WinterPCoAxis_2  0.06168323
    ## + SummerPCoAxis_2  0.04435527
    ## + SummerPCoAxis_3  0.01242937
    ## <none>             0.00000000
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_1  1 -63.434 23.772  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.4669063 
    ## Call: rodents_hel ~ WinterPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_3   0.5354021
    ## + WinterPCoAxis_2   0.5196228
    ## + SummerPCoAxis_2   0.5130107
    ## + SummerPCoAxis_1   0.5017403
    ## + SummerPCoAxis_3   0.4797432
    ## <none>              0.4669063
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_3  1 -66.249 4.6858  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5354021 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_2   0.5905277
    ## + SummerPCoAxis_2   0.5837467
    ## + SummerPCoAxis_1   0.5753749
    ## + SummerPCoAxis_3   0.5526490
    ## <none>              0.5354021
    ## 
    ##                   Df     AIC     F Pr(>F)   
    ## + WinterPCoAxis_2  1 -68.808 4.231  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5905277 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_2   0.6235862
    ## + SummerPCoAxis_3   0.6110886
    ## + SummerPCoAxis_1   0.5980486
    ## <none>              0.5905277
    ## 
    ##                   Df     AIC    F Pr(>F)   
    ## + SummerPCoAxis_2  1 -70.281 3.02  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6235862 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      SummerPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_3   0.6461832
    ## + SummerPCoAxis_1   0.6323416
    ## <none>              0.6235862
    ## 
    ##                   Df     AIC      F Pr(>F)  
    ## + SummerPCoAxis_3  1 -71.209 2.4051  0.034 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6461832 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      SummerPCoAxis_2 + SummerPCoAxis_3 
    ##  
    ##                   R2.adjusted
    ## <none>              0.6461832
    ## + SummerPCoAxis_1   0.6433899

``` r
# Most parsimonious is Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      SummerPCoAxis_2

# Evaluate the R2 and significance of the most parsimonious model

rod_prda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2 + Condition(year), pred_vals)

pR2p <- RsquareAdj(rod_prda_pars)$r.squared
pR2adjp <- RsquareAdj(rod_prda_pars)$adj.r.squared

pR2p
```

    ## [1] 0.3358357

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
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2 + Condition(year), data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     4 0.055986 6.5604  0.001 ***
    ## Residual 21 0.044803                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(rod_prda_pars, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2 + Condition(year), data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.040659 19.0578  0.001 ***
    ## RDA2      1 0.011742  5.5036  0.001 ***
    ## RDA3      1 0.003100  1.4530  0.416    
    ## RDA4      1 0.000485  0.2272  0.984    
    ## Residual 21 0.044803                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Only the first 2 axes - Winter 1 and 3 - are significant.

Variance partitioning
---------------------

``` r
# Run variance partitioning on significant axes of the most parsimonious RDA
rod_ppart <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~  WinterPCoAxis_3, ~ year, data = pred_vals)
rod_ppart
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~WinterPCoAxis_1,
    ## ~WinterPCoAxis_3, ~year, data = pred_vals)
    ## 
    ## Explanatory tables:
    ## X1:  ~WinterPCoAxis_1
    ## X2:  ~WinterPCoAxis_3
    ## X3:  ~year 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 4.3343 
    ##             Variance: 0.16671 
    ## No. of observations: 27 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1         1  0.48741      0.46691     TRUE
    ## [b+d+e+g] = X2         1  0.10033      0.06435     TRUE
    ## [c+e+f+g] = X3         1  0.39541      0.37123     TRUE
    ## [a+b+d+e+f+g] = X1+X2  2  0.57114      0.53540     TRUE
    ## [a+c+d+e+f+g] = X1+X3  2  0.56254      0.52608     TRUE
    ## [b+c+d+e+f+g] = X2+X3  2  0.51259      0.47197     TRUE
    ## [a+b+c+d+e+f+g] = All  3  0.62289      0.57370     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3       1               0.10174     TRUE
    ## [b] = X2 | X1+X3       1               0.04762     TRUE
    ## [c] = X3 | X1+X2       1               0.03830     TRUE
    ## [d]                    0               0.05312    FALSE
    ## [e]                    0               0.02088    FALSE
    ## [f]                    0               0.36932    FALSE
    ## [g]                    0              -0.05727    FALSE
    ## [h] = Residuals                        0.42630    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3        1               0.15486     TRUE
    ## [a+f] = X1 | X2        1               0.47105     TRUE
    ## [b+d] = X2 | X3        1               0.10074     TRUE
    ## [b+e] = X2 | X1        1               0.06850     TRUE
    ## [c+e] = X3 | X1        1               0.05918     TRUE
    ## [c+f] = X3 | X2        1               0.40762     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_ppart, digits = 2, Xnames = c('Winter1', 'Winter3', 'year'))
title(main="Variance partitioning - Partial RDA")
```

![](narrated_script_files/figure-markdown_github/partial%20variance%20partitioning-1.png)

WinterPCoAxis\_1 combined with year has the largest chunk (.369); on its own, WinterPCoAxis\_1 explains an additional .1

Regular RDA
-----------

``` r
# RDA on all variables (no condition)
rodents_rda <- rda(rodents_hel ~ ., pred_vals)

# R2 and significance test for global model
R2 <- RsquareAdj(rodents_rda)$r.squared
R2adj <- RsquareAdj(rodents_rda)$adj.r.squared
R2
```

    ## [1] 0.7566306

``` r
R2adj
```

    ## [1] 0.6669682

``` r
anova(rodents_rda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ year + WinterPCoAxis_1 + WinterPCoAxis_2 + WinterPCoAxis_3 + SummerPCoAxis_1 + SummerPCoAxis_2 + SummerPCoAxis_3, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     7 0.126134 8.4387  0.001 ***
    ## Residual 19 0.040571                  
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
    ## RDA1      1 0.090448 42.3583  0.001 ***
    ## RDA2      1 0.018913  8.8575  0.001 ***
    ## RDA3      1 0.010766  5.0418  0.006 ** 
    ## RDA4      1 0.003226  1.5108  0.798    
    ## RDA5      1 0.001655  0.7749  0.982    
    ## RDA6      1 0.000730  0.3418  0.998    
    ## RDA7      1 0.000396  0.1856  0.993    
    ## Residual 19 0.040571                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Find the most parsimonious model
set.seed(1)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_rda), 
                           R2scope = F, direction = "forward", pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: rodents_hel ~ 1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_1  0.46690627
    ## + year             0.37122615
    ## + SummerPCoAxis_1  0.22573181
    ## + WinterPCoAxis_3  0.06434779
    ## + WinterPCoAxis_2  0.06168323
    ## + SummerPCoAxis_2  0.04435527
    ## + SummerPCoAxis_3  0.01242937
    ## <none>             0.00000000
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_1  1 -63.434 23.772  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.4669063 
    ## Call: rodents_hel ~ WinterPCoAxis_1 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_3   0.5354021
    ## + year              0.5260839
    ## + WinterPCoAxis_2   0.5196228
    ## + SummerPCoAxis_2   0.5130107
    ## + SummerPCoAxis_1   0.5017403
    ## + SummerPCoAxis_3   0.4797432
    ## <none>              0.4669063
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + WinterPCoAxis_3  1 -66.249 4.6858  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5354021 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 
    ##  
    ##                   R2.adjusted
    ## + WinterPCoAxis_2   0.5905277
    ## + SummerPCoAxis_2   0.5837467
    ## + SummerPCoAxis_1   0.5753749
    ## + year              0.5737045
    ## + SummerPCoAxis_3   0.5526490
    ## <none>              0.5354021
    ## 
    ##                   Df     AIC     F Pr(>F)   
    ## + WinterPCoAxis_2  1 -68.808 4.231  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.5905277 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + year              0.6288360
    ## + SummerPCoAxis_2   0.6235862
    ## + SummerPCoAxis_3   0.6110886
    ## + SummerPCoAxis_1   0.5980486
    ## <none>              0.5905277
    ## 
    ##        Df     AIC      F Pr(>F)   
    ## + year  1 -70.661 3.3739  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.628836 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      year 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_2   0.6672563
    ## <none>              0.6288360
    ## + SummerPCoAxis_1   0.6275715
    ## + SummerPCoAxis_3   0.6268560
    ## 
    ##                   Df     AIC      F Pr(>F)   
    ## + SummerPCoAxis_2  1 -72.867 3.5402  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.6672563 
    ## Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +      year + SummerPCoAxis_2 
    ##  
    ##                   R2.adjusted
    ## + SummerPCoAxis_1   0.6696281
    ## + SummerPCoAxis_3   0.6676086
    ## <none>              0.6672563
    ## 
    ##                   Df     AIC      F Pr(>F)
    ## + SummerPCoAxis_1  1 -72.378 1.1508  0.322

``` r
# Most parsimonious is rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 +  year 

rod_rda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + year, pred_vals)

# R2 and significance for most parsimonious model
R2p <- RsquareAdj(rod_rda_pars)$r.squared
R2adjp <- RsquareAdj(rod_rda_pars)$adj.r.squared
R2p
```

    ## [1] 0.6859381

``` r
R2adjp
```

    ## [1] 0.628836

``` r
anova(rod_rda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + year, data = pred_vals)
    ##          Df Variance      F Pr(>F)    
    ## Model     4 0.114349 12.012  0.001 ***
    ## Residual 22 0.052356                  
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
    ## Model: rda(formula = rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + year, data = pred_vals)
    ##          Df Variance       F Pr(>F)    
    ## RDA1      1 0.089264 37.5089  0.001 ***
    ## RDA2      1 0.017020  7.1519  0.001 ***
    ## RDA3      1 0.007588  3.1884  0.016 *  
    ## RDA4      1 0.000478  0.2007  0.992    
    ## Residual 22 0.052356                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# The first 3 axes are significant (winter 1, winter 3, winter 2)
```

### Variance partitioning on significant axes for regular RDA

``` r
rod_part <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~ WinterPCoAxis_3, ~  WinterPCoAxis_2, data = pred_vals)
rod_part
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = rodents_hel, X = ~WinterPCoAxis_1,
    ## ~WinterPCoAxis_3, ~WinterPCoAxis_2, data = pred_vals)
    ## 
    ## Explanatory tables:
    ## X1:  ~WinterPCoAxis_1
    ## X2:  ~WinterPCoAxis_3
    ## X3:  ~WinterPCoAxis_2 
    ## 
    ## No. of explanatory tables: 3 
    ## Total variation (SS): 4.3343 
    ##             Variance: 0.16671 
    ## No. of observations: 27 
    ## 
    ## Partition table:
    ##                       Df R.square Adj.R.square Testable
    ## [a+d+f+g] = X1         1  0.48741      0.46691     TRUE
    ## [b+d+e+g] = X2         1  0.10033      0.06435     TRUE
    ## [c+e+f+g] = X3         1  0.09777      0.06168     TRUE
    ## [a+b+d+e+f+g] = X1+X2  2  0.57114      0.53540     TRUE
    ## [a+c+d+e+f+g] = X1+X3  2  0.55657      0.51962     TRUE
    ## [b+c+d+e+f+g] = X2+X3  2  0.20467      0.13840     TRUE
    ## [a+b+c+d+e+f+g] = All  3  0.63777      0.59053     TRUE
    ## Individual fractions                                   
    ## [a] = X1 | X2+X3       1               0.45213     TRUE
    ## [b] = X2 | X1+X3       1               0.07090     TRUE
    ## [c] = X3 | X1+X2       1               0.05513     TRUE
    ## [d]                    0               0.00581    FALSE
    ## [e]                    0              -0.00241    FALSE
    ## [f]                    0               0.01892    FALSE
    ## [g]                    0              -0.00996    FALSE
    ## [h] = Residuals                        0.40947    FALSE
    ## Controlling 1 table X                                  
    ## [a+d] = X1 | X3        1               0.45794     TRUE
    ## [a+f] = X1 | X2        1               0.47105     TRUE
    ## [b+d] = X2 | X3        1               0.07671     TRUE
    ## [b+e] = X2 | X1        1               0.06850     TRUE
    ## [c+e] = X3 | X1        1               0.05272     TRUE
    ## [c+f] = X3 | X2        1               0.07405     TRUE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(rod_part, digits = 2, Xnames = c('Winter1', 'Winter3','Winter2'))
title(main="Variance partitioning - unconditioned RDA")
```

![](narrated_script_files/figure-markdown_github/regular%20variance%20partitioning-1.png) Winter PCoA 1 still captures most of the variation.

### Years - rodents comparison

Plot winter PCoA1, year, and DS % abundances.

``` r
par(mfrow = c(1, 2))
rodents_props <- rodents/rowSums(rodents)

plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
title(main="Winter axis 1")
abline(v = 1990, col = 'red')

plot(pred_vals$year, rodents_props$DS, col = 'blue')
title(main = 'D. spectabilis % abundance')
abline(v = 1990, col = 'red')
```

![](narrated_script_files/figure-markdown_github/plot%20winter%20pcoa1%20v%20year-1.png)

``` r
# Important species in axis 1
# Save the ranked list of highest absolute value-scoring species for 
# winter axis 1. 
impt_species <- winter_scores_ranks %>%
  select(speciescode, axis1_rank, waxis1) %>%
  arrange(axis1_rank)

write.csv(impt_species, 'models/winter_axis_1.csv')
```
