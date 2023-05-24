
# Copyright 2023 Daniel Eacker & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from code provided in https://doi.org/10.1002/wsb.950

#' Plot caribou Bayesian IPM results
#'
#' Plot results of the caribou Bayesian IPM with the option to include the
#' national model and observed data
#'
#' @param allRes data.frame. Summary of model results created using [getOutputTables()] 
#' @param parameter
#' @param obs
#' @param lowBound
#' @param highBound
#' @param simRange
#' @param facetVars
#' @param labFontSize
#'
#' @return
#' @export
#'
#' @examples
plotRes <- function(allRes, parameter, obs = NULL, lowBound = 0, highBound = 1,
                    simRange = NULL, facetVars = NULL, labFontSize = 14) {
  # allRes=scResults$ksDists; parameter="Recruitment";obs=scResults$obs.all;
  # lowBound=0; highBound=1;simRange=scResults$sim.all;facetVars=c("obsYears","sQuantile")
  
  pal2 = c("#EF8A62","#67A9CF")#brewer.pal(7,"RdBu")[c(2,6)]
  
  if (is.null(facetVars)) {
    titleFontSize <- 16
    # labFontSize <- 14
    breakInterval <- 1
  } else {
    titleFontSize <- 11
    labFontSize <- 10
    breakInterval <- 2
  }
  if (is.element("KSDistance", names(allRes))) {
    # plot Kolmogorov Smirnov distances
    KS <- T
    allRes$Mean <- allRes$KSDistance
  } else {
    KS <- F
  }
  
  df <- subset(allRes, allRes$Parameter == parameter)
  
  if (nrow(df) < 1) {
    stop()
  }
  
  if (!is.null(obs)) {
    pr <- parameter
    obs <- subset(obs, parameter == pr)
  }
  
  if (!KS & !is.null(simRange)) {
    pr <- parameter
    simRange <- subset(simRange, parameter == pr)
    
    df$Type <- "IPM"
    simRange$Type <- "national"
    nameSel <- c(c("Year", "Mean", "Lower 95% CRI", "Upper 95% CRI", "Type"), facetVars)
    df <- rbind(subset(df, select = nameSel), subset(simRange, select = nameSel))
    df$grp <- df$Type
    if (!is.null(facetVars)) {
      for (i in facetVars) {
        df$grp <- paste0(df$grp, df[[i]])
      }
    }
    x1 <- ggplot2::ggplot(df, ggplot2::aes(x = Year, y = Mean, fill = Type, col = Type))
  } else {
    x1 <- ggplot2::ggplot(df, ggplot2::aes(x = Year, y = Mean))
  }
  x2 <- x1 + ggplot2::theme_classic() + ggplot2::xlab("Year") +
    ggplot2::ylab(parameter) +
    ggplot2::geom_line(ggplot2::aes(x = Year, y = Mean), size = 1.75) +ggplot2::scale_color_discrete(type=pal2)+
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = labFontSize),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = labFontSize),
      axis.title.x = ggplot2::element_text(size = titleFontSize, face = "bold"),
      axis.title.y = ggplot2::element_text(size = titleFontSize, face = "bold")
    ) +
    ggplot2::scale_x_continuous(breaks = seq(
      min(df$Year, na.rm = TRUE),
      max(df$Year, na.rm = TRUE), breakInterval
    ))
  
  if (!KS) {
    x2 <- x2 + ggplot2::geom_ribbon(ggplot2::aes(ymin = `Lower 95% CRI`,
                                                 ymax = `Upper 95% CRI`),
                                    show.legend = FALSE, alpha = 0.25, colour = NA
    ) +ggplot2::scale_fill_discrete(type=pal2)+
      ggplot2::scale_y_continuous(limits = c(
        ifelse(any(df$`Lower 95% CRI` < lowBound), NA, lowBound),
        ifelse(any(df$`Upper 95% CRI` > 1), NA, highBound)
      ))
  }
  
  if (!KS & !is.null(obs)) {
    obs$Type <- "IPM"
    obs$obsError <- F
    obs$obsError[obs$type == "observed"] <- T
    x2 <- x2 + ggplot2::geom_point(data = obs,
                                   ggplot2::aes(x = Year, y = Mean,
                                                shape = obsError), col = "black",
                                   show.legend = T) +
      ggplot2::scale_shape_manual(values = c(16, 2))
  }
  
  if (!is.null(facetVars)) {
    if (length(facetVars) == 2) {
      x2 <- x2 + ggplot2::facet_grid(as.formula(paste(facetVars[1], "~", facetVars[2])),
                                     labeller = "label_both")
    } else {
      x2 <- x2 + ggplot2::facet_wrap(as.formula(paste0("~", facetVars[1])),
                                     labeller = "label_both")
    }
  }
  if (!KS & (parameter == "Population growth rate")) {
    x2 <- x2 + ggplot2::geom_hline(yintercept = 1, color = "black")
  }
  
  x2
}
