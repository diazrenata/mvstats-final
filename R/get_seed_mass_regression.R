get_seed_mass_regression <- function(scores_table, axis_index){
  axis <- scores_table[,axis_index]
  sizes <- scores_table$seed_mass
  
  this_model <- lm(axis ~ sizes)
  this_summ <- summary(this_model)
  
  this_coeff <- this_summ$coefficients
  this_r2 <- this_summ$r.squared
  
  outputs <- list(this_coeff, this_r2)
  
  return(outputs)
  
}