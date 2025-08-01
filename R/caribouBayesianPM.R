# Copyright 2025 Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3

#' Bayesian population model for boreal caribou
#'
#'
#' @param survData either a path to a csv file or a survival data table in bboutools format.
#' @param recruitData either a path to a csv file or a recruitment data table in bboutools format.
#' @param disturbance either a path to a csv file or a dataframe containing the
#'   columns "Anthro","fire_excl_anthro", and "Year".
#' @param betaPriors a list of model priors. See [getPriors()]. Not used if disturbance is NA.
#' @param startYear,endYear year defining the beginning of the observation
#'   period and the end of the projection period.
#' @param niters integer. The number of iterations per chain after thinning and burn-in.
#' @param nthin integer. The number of the thinning rate.
#' @param returnSamples logical. If F returns only summaries. If T returns example trajectories.
#' @inheritParams caribouPopGrowth
#' @param inputList an optional list of inputs with names matching the above. If
#'   an argument is included in this list it will override the named argument.
#' @param niters A whole number of the number of iterations per chain after thinning and burn-in.
#'
#' @return a list with elements:
#'   * result: a list of model results:
#'     * summary: a data.frame
#'     * samples: a tibble providing the full range of MCMC trajectories from the model. 
#'     It is in a long format where "Amount" gives the value for 
#'     each metric in Anthro, fire_excl_anthro, c, survival, recruitment, X, N, 
#'     lambda, Sbar, Rbar, Xbar, Nbar, and lambda_bar, with a row for each 
#'     combination of "MetricTypeID", "Replicate", "Year" and "LambdaPercentile"
#'     * surv_data: a data.frame
#'     * recruit_data: a tibble
#'     * popInfo: a data.frame
#'   * inData: a list of data that is used as input to the jags model:
#'     * survDataIn:  survival data
#'     * disturbanceIn: disturbance data
#'     * recruitDataIn: composition data
#' @family demography
#' @export
#'
#' @examples
#' \donttest{
#'   # Note these examples take a long time to run!
#'   
#'   # Using observed survival, recruitment and disturbance data
#'   mod <- caribouBayesianPM(
#'     survData = bboudata::bbousurv_a,
#'     recruitData = bboudata::bbourecruit_a,
#'     disturbance = NULL
#'   )
#'   str(mod, max.level = 2)
#'   
#'   # Using simulated observation data
#'   scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                               obsAnthroSlope = 1, projAnthroSlope = 5,
#'                               collarCount = 20, cowMult = 5)
#'   
#'   simO <- simulateObservations(scns)
#'   
#'   out <- caribouBayesianPM(survData = simO$simSurvObs, recruitData = simO$simRecruitObs,
#'                            disturbance = simO$simDisturbance,
#'                            startYear = 2014)
#' }

caribouBayesianPM <- function(survData = bboudata::bbousurv_a,
                       recruitData = bboudata::bbourecruit_a,
                       disturbance = NULL,
                       betaPriors = "default",
                       startYear = NULL, endYear = NULL,
                       N0=1000,
                       returnSamples=F,
                       inputList = list(),
                       niters=formals(bboutools::bb_fit_survival)$niters,nthin=formals(bboutools::bb_fit_survival)$nthin,
                       ...) {

  # combine defaults in function with inputs from input list
  inputArgs <- c(
    "survData", "recruitData", "disturbance", "startYear", "endYear", "niters", "nthin"
  )
  addArgs <- inputArgs # setdiff(inputArgs,names(inp))
  inp <- list()
  for (a in addArgs) {
    if (is.element(a, names(inputList))) {
      inp[[a]] <- inputList[[a]]
    } else {
      inp[[a]] <- eval(parse(text = a))
    }
  }

  if (betaPriors[[1]] == "default") {
    betaPriors <- getPriors()
  }
  
  # Run model
  if (is.character(inp$recruitData)) {
    recruitData <- read.csv(inp$recruitData, header = T)
    recruitData$X <- NULL
  } else {
    recruitData <- inp$recruitData
  }
  if (is.character(inp$disturbance)) {
    disturbance <- read.csv(inp$disturbance)
    disturbance$X <- NULL
  } else {
    disturbance <- inp$disturbance
  }
  if (is.character(inp$survData)) {
    survData <- read.csv(inp$survData, header = T)
    survData$X <- NULL
  } else {
    survData <- inp$survData
  }

  # if decide to error when Year ranges don't match could use testTable
  if(!is.null(disturbance)){
    testTable(disturbance, c("Year", "Anthro", "fire_excl_anthro"))
  }
  #TO DO: use bboutools data test functions for survival and recruitment

  # Get start and end years from data
  if(is.null(inp$startYear)){
    inp$startYear <- min(survData$Year)
  }

  if(is.null(inp$endYear)||is.infinite(inp$endYear)){
    if(!is.null(disturbance)){
      inp$endYear <- max(disturbance$Year)
    }else{inp$endYear <- max(survData$Year);distYrs =survData$Year}
  }

  if(!is.null(disturbance)){
    disturbance <- merge(data.frame(Year = seq(inp$startYear, inp$endYear)),
                         disturbance, by = "Year", all.x = T)
    if (anyNA(select(disturbance, "Anthro", "fire_excl_anthro"))) {
      warning(
        "Years ",
        filter(disturbance, if_any(c(Anthro, fire_excl_anthro), is.na)) %>%
          pull(Year) %>% paste0(collapse = ", "),
        " have missing disturbance data. ",
        "Anthro will be filled from the next year with data and fire will be fill with 0s"
      )
      
      disturbance <- tidyr::fill(disturbance, Anthro, .direction = "downup") %>%
        mutate(
          fire_excl_anthro = tidyr::replace_na(fire_excl_anthro, 0),
          Total_dist = fire_excl_anthro + Anthro
        )
      if(anyNA(select(disturbance, "Anthro", "fire_excl_anthro"))){
        #stop("None of the disturbance data is within the requested year range",
        #     call. = FALSE)
        disturbance=NULL
      }
    }
    distYrs = disturbance$Year
  }

  ################
  # Survival data checking and fill missing yrs
  surv_data <- survData

  # check that year range is within data - model will run either way
  if (inp$endYear < max(surv_data$Year) | inp$startYear < min(surv_data$Year)) {
    warning(c("requested year range: ", inp$startYear, " - ", inp$endYear, " does not match survival data",
              c(" year range:", "  ", min(surv_data$Year), " - ", max(surv_data$Year))))
  }
  surv_data <- subset(surv_data, surv_data$Year <= inp$endYear & surv_data$Year >= inp$startYear)
  if(nrow(surv_data) == 0){
    stop("None of the survival data is within the requested year range",
         call. = FALSE)
  }
  yrs_surv_missing <- setdiff(c(inp$startYear:max(surv_data$Year)),
                             unique(surv_data$Year))
  # check that year range is within data - model will run either way
  if (length(yrs_surv_missing) > 0) {
    warning("missing years of recruitment data:", " ", 
            paste0(yrs_surv_missing, collapse = ", "), 
            call. = FALSE) 
  }

  #add missing surv yrs
  surv_data_add = expand.grid(Year=union(distYrs,surv_data$Year),Month=unique(surv_data$Month),PopulationName=unique(surv_data$PopulationName))
  surv_data=merge(surv_data,surv_data_add,all.x=T,all.y=T)
  surv_data$StartTotal[is.na(surv_data$StartTotal)]=0
  
  #dups = table(subset(surv_data,select=c(Year,Month,PopulationName)))
  

  ###################
  # Recruitment data checking and fill missing yrs
  recruit_data <- recruitData
  
  # check that year range is within data - model will run either way
  if (inp$endYear < max(recruit_data$Year) | inp$startYear < min(recruit_data$Year)) {
    warning(c("requested year range: ", inp$startYear, " - ", inp$endYear, " does not match recruitment data",
              c(" year range:", "  ", min(recruit_data$Year), " - ", max(recruit_data$Year))), 
            call. = FALSE)
  }
  
  recruit_data <- subset(recruit_data, recruit_data$Year <= inp$endYear & recruit_data$Year >= inp$startYear)
  
  if(nrow(recruit_data) == 0){
    stop("None of the recruitment data is within the requested year range",
         call. = FALSE)
  }
  
  yrs_rec_missing <- setdiff(c(inp$startYear:max(recruit_data$Year)),
                             unique(recruit_data$Year))
  # check that year range is within data - model will run either way
  if (length(yrs_rec_missing) > 0) {
    warning("missing years of recruitment data:", " ", 
            paste0(yrs_rec_missing, collapse = ", "), 
            call. = FALSE) 
  }

  #add missing recruit yrs
  recruit_data_add = expand.grid(Year=union(distYrs,recruit_data$Year),PopulationName=unique(recruit_data$PopulationName))
  recruit_data=merge(recruit_data,recruit_data_add,all.x=T,all.y=T)
  recruit_data$Month[is.na(recruit_data$Month)]=3;recruit_data$Day[is.na(recruit_data$Day)]=15

  ##################
  #fit models
  bbouResults = bbouMakeSummaryTable(surv_data, recruit_data,N0,disturbance,priors=betaPriors,
                                     return_mcmc=T,shiny_progress=F,niters=niters,nthin=nthin)
  
  #get output trajectories
  rr = getSimsInitial(bbouResults,cPars=betaPriors,skipSave=T,returnSamples=returnSamples,...)  
  
  return(list(result = rr, 
              inData = list(disturbanceIn = disturbance)))
}
