#' Plot survival time series
#' 
#' TO DO: write documentation
#' 
#' @family demography
#' @export
plotSurvivalSeries <- function(surv_data_show) {
  # surv_data_show = subset(outObs$simSurvObs,Replicate==outObs$simSurvObs$Replicate[1])

  surv_data_show$MortalitiesCertain[surv_data_show$MortalitiesCertain == 0] <- NA
  surv_data_show$Malfunctions[surv_data_show$Malfunctions == 0] <- NA
  surv_data_show$Month <- factor(surv_data_show$Month, levels = seq(1:12))
  
  if (length(unique(surv_data_show$Month)) == 1) {
    base <- ggplot(surv_data_show, 
                   aes(x = Year, y = StartTotal, group = PopulationName, 
                       colour = PopulationName)) +
      geom_line() +
      geom_point(aes(x = Year, y = MortalitiesCertain, group = PopulationName, 
                     colour = PopulationName), shape = 4) +
      geom_col(aes(x = Year, y = Malfunctions, group = PopulationName, 
                   colour = PopulationName, fill = PopulationName), 
               alpha = 0.2) +
      ylab("Number of Animals") +
      theme_bw()
  } else {
    base <- ggplot(surv_data_show, 
                   aes(x = Month, y = StartTotal, group = PopulationName, 
                       colour = PopulationName)) +
      geom_line() +
      facet_wrap(~Year) +
      geom_point(aes(x = Month, y = MortalitiesCertain, group = PopulationName, 
                     colour = PopulationName), shape = 4) +
      geom_col(aes(x = Month, y = Malfunctions, group = PopulationName, 
                   colour = PopulationName, fill = PopulationName), 
               alpha = 0.2) +
      ylab("Number of Animals") +
      theme_bw()
  }
  return(base)
}
