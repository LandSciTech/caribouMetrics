#' Plot survival time series
#' 
#' TO DO: write documentation
#' 
#' @family demography
#' @export
plotSurvivalSeries <- function(surv_data_show) {
  # surv_data_show = subset(outObs$simSurvObs,Replicate==outObs$simSurvObs$Replicate[1])

  surv_data_show$MortalitiesCertain[surv_data_show$MortalitiesCertain == 0] <- NA
  if(!hasName(surv_data_show,"Malfunctions")){surv_data_show$Malfunctions<-NA}
  surv_data_show$Malfunctions[surv_data_show$Malfunctions == 0] <- NA
  surv_data_show$Month <- factor(surv_data_show$Month, levels = seq(1:12))
  
  if (length(unique(surv_data_show$Month)) == 1) {
    base <- ggplot2::ggplot(surv_data_show, 
                            ggplot2::aes(x = Year, y = StartTotal, group = PopulationName, 
                                         colour = PopulationName)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(ggplot2::aes(x = Year, y = MortalitiesCertain, group = PopulationName, 
                                       colour = PopulationName), shape = 4) +
      ggplot2::geom_col(ggplot2::aes(x = Year, y = Malfunctions, group = PopulationName, 
                                     colour = PopulationName, fill = PopulationName), 
                        alpha = 0.2) +
      ggplot2::ylab("Number of Animals") +
      ggplot2::theme_bw()
  } else {
    base <- ggplot2::ggplot(surv_data_show, 
                            ggplot2::aes(x = Month, y = StartTotal, group = PopulationName, 
                                         colour = PopulationName)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~Year) +
      ggplot2::geom_point(ggplot2::aes(x = Month, y = MortalitiesCertain, group = PopulationName, 
                                       colour = PopulationName), shape = 4) +
      ggplot2::geom_col(ggplot2::aes(x = Month, y = Malfunctions, group = PopulationName, 
                                     colour = PopulationName, fill = PopulationName), 
                        alpha = 0.2) +
      ggplot2::ylab("Number of Animals") +
      ggplot2::theme_bw()
  }
  return(base)
}
