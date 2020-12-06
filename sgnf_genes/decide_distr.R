{
  library(readr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(VGAM)
  library(MASS)
  library(fitdistrplus)
  library(highcharter)
}

genemotif_count %>% dplyr::filter(cancer_type == cancer) -> df
distr = df$count

descdist(distr, boot = 100, discrete = TRUE) 
descdist(distr, boot = 100, discrete = FALSE)
