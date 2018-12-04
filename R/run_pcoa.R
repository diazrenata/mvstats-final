run_pcoa <- function(plant_data, season) {
  
  library(vegan)
  library(ca)
  
  dist_mat <- vegdist(plant_data[,2:ncol(plant_data)], 'bray')
  
  plant_pcoa <- cmdscale(dist_mat, k = nrow(plant_data) - 1, eig = T)
  
  str(plant_pcoa)
  plant_pcoa$points
  
  # PCoordinatesA table to look at the eigenvalues
  # and the proportion of variance they capture
  
  eigenvalues <- plant_pcoa$eig[1:nrow(plant_data)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  PCoA_Table
  
  
  # Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  
  # if it's not totally intractable, I'd keep 3 axes
  # to explain a cumulative 91% of variance
  
  # Plotting the first two PCoA axes
  
  x <- plant_pcoa$points[,1]
  y <- plant_pcoa$points[,2]
  plot(x, y, xlab = 'Coord 1', ylab = 'Coord 2',
       xlim = range(x) * 1.2, 
       ylim = range(y) * 1.2, 
       type = 'n')
  text(x,y, labels = plant_data$year, cex = 0.9)
  
 pcoa_vals <- cbind(plant_data$year, plant_pcoa$points[,1:10])
  
  season_names <- rep(paste0(season, "PCoAxis_"), 10) %>%
    paste0(1:10)
  season_names <- c('year', season_names)
  colnames(pcoa_vals) <- season_names
  
  
  save(plant_pcoa,file = paste0('final-project/models/', season, '_pcoa.Rdata'))
  write.csv(pcoa_vals, paste0('final-project/models/', season, '_pcoa_vals.csv'), row.names = F)
  write.csv(PCoA_Table, paste0('final-project/models/', season, '_pcoa_table.csv'), row.names = F)
}