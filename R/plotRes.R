
# Copyright 2023 Daniel Eacker & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from code provided in https://doi.org/10.1002/wsb.950

#' Plot Bayesian population model results
#'
#' Plot Bayesian population model results, with (optionally) the 
#' distribution of outcomes from the national model, local observations, and true local state for comparison. 
#'
#' @param modTables list. A list of model results tables created using
#'   `[getOutputTables()]`.
#' @param parameter character. Which parameter to plot, if more than one, a list
#'   of plots is returned.
#' @param lowBound,highBound numeric. Lower and upper y axis limits
#' @param facetVars character. Optional. Vector of column names to facet by
#' @param labFontSize numeric. Optional. Label font size if there are not
#'   facets. Font size is 10 pt if facets are used.
#' @param ksDists logical. If true the `modTables$ksDists` table is used to
#'   create plots for each parameter.
#' @param legendPosition "bottom", "right", "left","top", or "none". Legend position.
#' @param breakInterval number. How many years between x tick marks?
#'
#' @return a ggplot object or list of ggplot objects if a vector of parameters
#'   was given.
#' @export
#' 
#' @family demography
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             obsAnthroSlope = 1, projAnthroSlope = 5,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#'
#' out <- caribouBayesianIPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
#'                           disturbance = simO$simDisturbance,
#'                           startYear = 2014, Nchains = 1, Niter = 100, Nburn = 10,
#'                           Nthin = 2)
#'
#' out_tbl <- getOutputTables(out, exData = simO$exData, paramTable = simO$paramTable,
#'                            simNational = getSimsNational(), getKSDists = FALSE)
#'
#' plotRes(out_tbl, parameter = "Recruitment")
plotRes <- function(modTables, parameter, lowBound = 0, highBound = 1,
                   facetVars = NULL, labFontSize = 14, ksDists = FALSE,legendPosition="right",breakInterval=1) {
  # allRes=scResults$ksDists; parameter="Recruitment";obs=scResults$obs.all;
  # lowBound=0; highBound=1;simRange=scResults$sim.all;facetVars=c("obsYears","sQuantile")
  
  if(length(parameter) > 1){
    allPlots <- lapply(parameter, function(x) {
      plotRes(modTables, x,  lowBound, highBound, facetVars, labFontSize, ksDists)
    })
    
    names(allPlots) <- parameter
    return(allPlots)
  }
  
  allRes <- modTables$rr.summary.all
  obs <- modTables$obs.all
  simRange <- modTables$sim.all
  
  pal2 = c("#EF8A62","#67A9CF")#brewer.pal(7,"RdBu")[c(2,6)]
  
  if(!is.data.frame(allRes)){
    stop("allRes must be a data.frame not a ", class(allRes), call. = FALSE)
  }
  
  exp_param_nms <- c(
    "Adult female survival", "Recruitment",
    "Adjusted recruitment", "Population growth rate", "Female population size",
    "Mean adult female survival",
    "Mean recruitment", "Mean female recruitment",
    "Median population growth rate",
    "Mean population growth rate"
  )
  
  if(!parameter %in% exp_param_nms){
    stop("parameter ", parameter, " is not one of the expected values: '",
         paste0(exp_param_nms, collapse = "', '"), call. = FALSE)
  }
  
  testTable(allRes, req_col_names = c("Year", "Parameter", "Mean", 
                                      "Lower 95% CRI", "Upper 95% CRI"),
            acc_vals = list(Parameter = exp_param_nms))
  
  if (is.null(facetVars)) {
    titleFontSize <- 16*labFontSize/14
    # labFontSize <- 14
  } else {
    titleFontSize <- 11
    labFontSize <- 10
    if(breakInterval==1){
      breakInterval <- 2
    }
  }
  if (ksDists) {
    # plot Kolmogorov Smirnov distances
    allRes <- modTables$ksDists
    allRes$Mean <- allRes$KSDistance
  } 
  
  df <- subset(allRes, allRes$Parameter == parameter)
  
  if (nrow(df) < 1) {
    stop("The parameter: ", parameter, " is not present in the data provided")
  }
  
  if (!is.null(obs)) {
    pr <- parameter
    obs <- subset(obs, parameter == pr)
  }
  
  if(!is.null(simRange)){
    pr <- parameter
    simRange <- subset(simRange, parameter == pr)
    if(nrow(simRange) == 0){
     simRange <- NULL  
    }
  }
  
  if (!ksDists & !is.null(simRange)) {
    
    df$Type <- "Bayesian"
    simRange$Type <- "national"
    nameSel <- c(c("Year", "Mean", "Lower 95% CRI", "Upper 95% CRI", "Type"), facetVars)
    df <- rbind(subset(df, select = nameSel), subset(simRange, select = nameSel))
    df$grp <- df$Type
    if (!is.null(facetVars)) {
      for (i in facetVars) {
        df$grp <- paste0(df$grp, df[[i]])
      }
    }
    x1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Year"]], y = .data[["Mean"]], 
                                           fill = .data[["Type"]], col = .data[["Type"]]))
  } else {
    x1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Year"]], y = .data[["Mean"]]))
  }
  x2 <- x1 + ggplot2::theme_classic() + ggplot2::xlab("Year") +
    ggplot2::ylab(parameter) +
    ggplot2::geom_line(ggplot2::aes(x = .data[["Year"]], y = .data[["Mean"]]), 
                       linewidth = 1.75) +
    ggplot2::scale_color_discrete(type=pal2, name = NULL)+
    ggplot2::theme(
      legend.position = legendPosition,
      axis.text.y = ggplot2::element_text(size = labFontSize),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = labFontSize),
      axis.title.x = ggplot2::element_text(size = titleFontSize, face = "bold"),
      axis.title.y = ggplot2::element_text(size = titleFontSize, face = "bold")
    ) +
    ggplot2::scale_x_continuous(breaks = seq(
      min(df$Year, na.rm = TRUE),
      max(df$Year, na.rm = TRUE), breakInterval
    ))
  
  if (!ksDists) {
    x2 <- x2 + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["Lower 95% CRI"]],
                                                 ymax = .data[["Upper 95% CRI"]]),
                                    show.legend = FALSE, alpha = 0.25, colour = NA
    ) +ggplot2::scale_fill_discrete(type=pal2, name = NULL)+
      ggplot2::scale_y_continuous(limits = c(
        ifelse(any(df$`Lower 95% CRI` < lowBound), NA, lowBound),
        ifelse(any(df$`Upper 95% CRI` > 1), NA, highBound)
      ))
  }
  
  if (!ksDists & !is.null(obs)) {
    if(nrow(obs) > 0){
      obs$Type <- "Bayesian"
      obs$obsError <- FALSE
      obs$obsError[obs$type == "observed"] <- TRUE
      obs <- filter(obs, !is.na(.data$Mean))
      x2 <- x2 + ggplot2::geom_point(data = obs,
                                     ggplot2::aes(x = .data[["Year"]], y = .data[["Mean"]],
                                                  shape = .data[["obsError"]]), 
                                     col = "black", fill = "black", inherit.aes = FALSE,
                                     show.legend = TRUE) +
        ggplot2::scale_shape_manual(values = c(16, 2), 
                                    labels = c(`TRUE` = "observed", 
                                               `FALSE` = "true\nsimulated"), 
                                    name = NULL)
    }
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
  if (!ksDists & (parameter == "Population growth rate")) {
    x2 <- x2 + ggplot2::geom_hline(yintercept = 1, color = "black")
  }
  
  x2
}
