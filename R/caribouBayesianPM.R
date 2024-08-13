# Copyright 2023 Daniel Eacker & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from code provided in https://doi.org/10.1002/wsb.950

#' Bayesian population model for boreal caribou
#'
#' A Bayesian population model that integrates prior
#' information from Johnson et al.'s (2020) national analysis of
#' demographic-disturbance relationships with local demographic data to project
#' population growth.
#'
#' The model combines local observations of survival (`survData`),
#' recruitment based on calf:cow ratios (`ageRatio`) and anthropogenic
#' disturbance (`disturbance`) with prior information on
#' the relationship between disturbance and survival and recruitment from the
#' Johnson et al. (2020) national model (`getPriors()`) to reduce uncertainty and refine parameter estimates.
#'
#' For a detailed description see the
#' [vignette](https://landscitech.github.io/caribouMetrics/articles/BayesianDemographicProjection.html#integration-of-local-demographic-data-and-national-disturbance-demographic-relationships-in-a-bayesian-population-model)
#' (`vignette("BayesianDemographicProjection", package = "caribouMetrics")`).
#'
#' Note: if `survData` contains values for enter that are > 0 these rows
#' will be dropped to avoid errors when collars are added in the middle of the
#' year. This will reduce the sample size in years when new collars are added.
#'
#' @param survData either a path to a csv file or a dataframe containing the
#'   columns "Year", "event", "enter" and "exit". Enter and exit are the
#'   beginning and end of the time interval and should be a number from 1 to 12
#'   (December) or 0 (See Details). Event is 0 or 1 where 0 means the animal
#'   lived and 1 died. See [survival::Surv()] for more details.
#' @param ageRatio either a path to a csv file or a dataframe containing the
#'   columns "Year","Count", and "Class". Where class can be either "calf" or "cow"
#' @param disturbance either a path to a csv file or a dataframe containing the
#'   columns "Anthro","fire_excl_anthro", and "Year".
#' @param betaPriors a list of model priors. See [getPriors()].
#' @param startYear,endYear year defining the beginning of the observation
#'   period and the end of the projection period.
#' @param Nchains Number of chains for the MCMC algorithm.
#' @param Niter Number of iterations for the MCMC algorithm.
#' @param Nburn Length of burn-in for the MCMC algorithm.
#' @param Nthin Thinning rate for the MCMC algorithm. 
#' @param N0 Initial population size.
#' @param survAnalysisMethod Survival analysis method either "KaplanMeier" or
#'   "Exponential". The exponential method is only recommended when the number of collared animals (in survData) is small.
#' @inheritParams caribouPopGrowth
#' @param assessmentYrs Number of years over which to assess population growth rate lambda.
#' @param inputList an optional list of inputs with names matching the above. If
#'   an argument is included in this list it will override the named argument.
#' @param saveJAGStxt file path. Directory where the JAGS model txt files will
#'   be saved. Default is `tempdir()`.
#' @param quiet logical. Should jags run quietly?
#'
#' @return a list with elements:
#'   * result: an `rjags` model object see [R2jags::jags()].
#'   * inData: a list of data that is used as input to the jags model:
#'     * survDataIn:  survival data
#'     * disturbanceIn: disturbance data
#'     * ageRatioIn: composition data
#' @family demography
#' @export
#'
#' @examples
#' # Using observed survival, recruitment and disturbance data
#' mod <- caribouBayesianPM(
#'   survData = system.file("extdata/simSurvData.csv",
#'                          package = "caribouMetrics"),
#'   ageRatio = system.file("extdata/simAgeRatio.csv",
#'                          package = "caribouMetrics"),
#'   disturbance = system.file("extdata/simDisturbance.csv",
#'                             package = "caribouMetrics"),
#'   Nchains = 1, Niter = 100, Nburn = 10, Nthin = 2
#' )
#' str(mod, max.level = 2)
#'
#' # Using simulated observation data
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             obsAnthroSlope = 1, projAnthroSlope = 5,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#'
#' out <- caribouBayesianPM(survData = simO$simSurvObs, ageRatio = simO$ageRatioOut,
#'                           disturbance = simO$simDisturbance,
#'                           startYear = 2014, Nchains = 1, Niter = 100, Nburn = 10,
#'                           Nthin = 2)

caribouBayesianPM <- function(survData = system.file("extdata/simSurvData.csv",
                                              package = "caribouMetrics"),
                       ageRatio = system.file("extdata/simAgeRatio.csv",
                                                   package = "caribouMetrics"),
                       disturbance = system.file("extdata/simDisturbance.csv",
                                                 package = "caribouMetrics"),
                       betaPriors = "default",
                       startYear = NULL, endYear = NULL, Nchains = 4,
                       Niter = 15000, Nburn = 10000, Nthin = 2, N0 = 1000,
                       survAnalysisMethod = "KaplanMeier", adjustR = FALSE,
                       assessmentYrs = 1,
                       inputList = list(), saveJAGStxt = tempdir(),
                       quiet = TRUE) {
  # survData=oo$simSurvObs;ageRatio=oo$ageRatioOut;disturbance=oo$simDisturbance;
  # betaPriors=betaPriors;startYear = minYr;endYear=maxYr;N0=cs$N0;survAnalysisMethod = "KaplanMeier"; adjustR=F
  # Nchains = 2;Niter = 20000;Nburn = 10000;Nthin = 1;assessmentYrs = 3;inputList=list();saveJAGStxt=tempdir();quiet=F

  # combine defaults in function with inputs from input list
  inputArgs <- c(
    "survData", "ageRatio", "disturbance", "startYear", "endYear",
    "Nchains", "Niter", "Nburn", "Nthin", "N0", "survAnalysisMethod", "adjustR",
    "assessmentYrs"
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
  if (is.character(inp$ageRatio)) {
    ageRatio <- read.csv(inp$ageRatio, header = T)
    ageRatio$X <- NULL
  } else {
    ageRatio <- inp$ageRatio
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
  testTable(disturbance, c("Year", "Anthro", "fire_excl_anthro"))
  testTable(ageRatio, c("Year", "Count", "Class"))
  testTable(survData, c("id", "Year", "event", "enter", "exit"))
  
  # Get start and end years from data
  if(is.null(inp$startYear)){
    inp$startYear <- min(survData$Year)
  }
  
  if(is.null(inp$endYear)){
    inp$endYear <- max(disturbance$Year)
  }

  disturbance <- merge(data.frame(Year = seq(inp$startYear, inp$endYear)),
                       disturbance, by = "Year", all.x = T)
  if (anyNA(select(disturbance, "Anthro", "fire_excl_anthro"))) {
    warning(
      "Years ",
      filter(disturbance, if_any(c(.data$Anthro, .data$fire_excl_anthro), is.na)) %>%
        pull(.data$Year) %>% paste0(collapse = ", "),
      " have missing disturbance data. ",
      "Anthro will be filled from the next year with data and fire will be fill with 0s"
    )

    disturbance <- tidyr::fill(disturbance, .data$Anthro, .direction = "downup") %>%
      mutate(
        fire_excl_anthro = tidyr::replace_na(.data$fire_excl_anthro, 0),
        Total_dist = .data$fire_excl_anthro + .data$Anthro
      )
    if(anyNA(select(disturbance, "Anthro", "fire_excl_anthro"))){
      stop("None of the disturbance data is within the requested year range",
           call. = FALSE)
    }
  }

  # TODO: Give a better name. survData is overwritten later so it is confusing
  data <- survData

  # check that year range is within data - model will run either way
  data$Year <- as.numeric(data$Year)

  # check that year range is within data - model will run either way
  if (inp$endYear < max(data$Year) | inp$startYear < min(data$Year)) {
    warning(c("requested year range does not match survival data",
              c(" year range:", "  ", min(data$Year), " - ", max(data$Year))))
  }

  data <- subset(data, data$Year <= inp$endYear & data$Year >= inp$startYear)

  if(nrow(data) == 0){
    stop("None of the survival data is within the requested year range",
         call. = FALSE)
  }

  data$id <- factor(data$id)

  test1 <- length(c(inp$startYear:max(data$Year)))
  test2 <- length(unique(data$Year))

  # check that year range is within data - model will run either way
  if (test1 > test2) {
    warning(c("missing years of survival data.",
              "Start model at beginning of consecutive years.",
              " ", c("Years of survival data:", "  ",
                     list(sort(unique(data$Year))))))
  }

  data <- data[order(data$exit), ]
  list_data1 <- split(data, data$Year)
  nYears <- length(levels(as.factor(data$Year)))
  n.ind <- numeric(nYears)

  for (i in 1:nYears) {
    n.ind[i] <- length(list_data1[[i]]$exit)
  }

  # check for low sample size - model will run either way
  if (any(n.ind < 20)) {
    warning(c("warning, low sample size of adult females in at least one year"))
  }

  # get KM estimates to use for adult female survival

  # Note: biased results from years with <12 months of observations.
  # And problems with adding animals part way through the year, so omitting those
  #data = data.frame(id=1,Year=2023,event=NA,enter=NA,exit=NA)
  dSubset <- subset(data, data$enter == 0) # ;dSubset=subset(dSubset,!((exit<12)&(event==0)))
  if (nrow(dSubset) == 0) {
    #stop("Collars not present at the start of a year are omitted from survival ",
    #     "analysis in that year. Please ensure there is at least one year with",
    #     " a collar in the first month.")
    warning("Collars not present at the start of a year are omitted from survival ",
         "analysis in that year. There are no years with collars in the first month.")
    inp$survAnalysisMethod <- "Exponential"
    survData = merge(data.frame(X1=0,X2=0,X3=0,X4=0,X5=0,X6=0,X7=0,X8=0,X9=0,X10=0,X11=0,X12=0,X13=0),data)
  }else{
    if (nrow(dSubset) == 1) {
      warning("Switching to exponential survival analysis method because there is",
              " only one collared animal.")
      inp$survAnalysisMethod <- "Exponential"
    }
    if (sum(dSubset$event, na.rm = T) == 0) {
      warning("Switching to exponential survival analysis method because there ",
              "are no recorded deaths.")
      inp$survAnalysisMethod <- "Exponential"
    }
    
    if (inp$survAnalysisMethod == "KaplanMeier") {
      survData <- getKMSurvivalEstimates(dSubset)
      # omitting years with less than 12 months of observations of collared animals
      sumDat <- dSubset %>%
        group_by(.data$Year) %>%
        summarise(minEnter = min(.data$enter), maxExit = max(.data$exit))
      includeYrs <- subset(sumDat, sumDat$minEnter == 0 & sumDat$maxExit == 12)
      survData$Year <- as.numeric(gsub("as.factor(Year)=", "",
                                       as.character(survData$Var1), fixed = T))
      survData <- merge(survData, includeYrs)
      if (nrow(survData) == 0) {
        warning("Years with less than 12 months of collar data are omitted from",
                " survival analysis. Please ensure there is 12 months of collar ",
                "data in at least one year.")
        warning("Switching to exponential survival analysis method because there ",
                "are no years with 12 months of collar data.")
        inp$survAnalysisMethod <- "Exponential"
      }
    }
    if (inp$survAnalysisMethod == "KaplanMeier") {
      message("using Kaplan-Meier survival model")
      if (any(survData$surv == 1)) {
        # which years does survival equal 1
        survOne <- which(unlist(lapply(split(survData, survData$Year),
                                       function(x) sum(x$event))) == 0)
        yearsOne <- as.numeric(names(survOne))
        data.sub <- data[data$Year %in% yearsOne, ]
        nriskYears <- data.frame(table(data.sub$Year))
        
        binLikFile <- file.path(saveJAGStxt, "binLik.txt")
        
        # Specify model
        sink(binLikFile)
        cat("
	   model{
	     for(i in 1:N){
	      lived[i] ~ dbin(s[i], atrisk[i])
	      s[i] ~ dbeta(1,1) # vague prior
	      }
	     }
	      ", fill = TRUE)
        sink()
        
        data1 <- list(N = nrow(nriskYears), lived = nriskYears$Freq,
                      atrisk = nriskYears$Freq)
        params <- c("s")
        inits <- function() {
          list(s = runif(nrow(nriskYears), 0.80, 0.99))
        }
        
        # run model in JAGS
        out1 <- R2jags::jags(
          data = data1, inits = inits, parameters.to.save = params,
          model.file = binLikFile, n.chains = 2, n.iter = 5000,
          n.burnin = 1000, n.thin = 1, quiet = quiet
        )
        # create standard deviation variable from survData$tau above
        survData$tau <- 1 / (survData$se^2)
        survData$tau[survOne] <- 1 / (out1$BUGSoutput$sd$s^2)
      } else {
        survData$tau <- 1 / (survData$se^2)
      }
      
      surv_id <- which(survData$surv != 1)
      nSurv <- length(surv_id)
      survData$Var1 <- as.character(survData$Var1)
    } else {
      message("expanding survival record")
      dExpand <- apply(subset(dSubset, select = c("id", "Year", "event", "enter", "exit")),
                       1, expandSurvivalRecord)
      survData <- do.call(rbind, dExpand)
    }
  }


  # split data into calf and cow recruitment data

  calf.cnt <- subset(ageRatio, ageRatio$Class == "calf")
  calf.cnt$Class <- factor(calf.cnt$Class)
  cow.cnt <- subset(ageRatio, ageRatio$Class == "cow")
  cow.cnt$Class <- factor(cow.cnt$Class)

  # deal with missing years of data between year ranges
  Years2 <- data.frame(sort(unique(data$Year)))
  names(Years2) <- "Year"
  y1 <- merge(Years2, calf.cnt, by = "Year", all = TRUE)
  data3 <- y1[, 1:3]

  y2 <- merge(Years2, cow.cnt, by = "Year", all = TRUE)
  data4 <- y2[, 1:3]

  data3$Count <- ifelse(data3$Count > data4$Count, NA, data3$Count)
  data4$Count <- ifelse(data3$Count > data4$Count, NA, data4$Count)

  data3$Count[(data3$Count==1)&(data4$Count==1)]=NA #JAGS fails in this case
  data4$Count[(data3$Count==1)&(data4$Count==1)]=NA #JAGS fails in this case
  
  xCalf <- which(is.na(data3$Count) == T)
  xCow <- which(is.na(data4$Count) == T)
  
  Years4 <- levels(as.factor(data4$Year))[xCalf]

  if (any(is.na(data3$Count) == T)) {
    warning("missing composition data; missing years:", " ", list(Years4))
  }

  t.pred <- max(inp$endYear - max(data3$Year), 0)

  # also add missing history and NA for projection period
  missingSurvYrs <- setdiff(seq(inp$startYear, inp$endYear), survData$Year)
  if (length(missingSurvYrs) > 0) {
    survAddBit <- survData[1, ]
    if (inp$survAnalysisMethod == "KaplanMeier") {
      survAddBit[1, ] <- NA
      survAddBit$Var1 <- NULL
      survAddBit$Year <- NULL
      survAddBit <- merge(survAddBit, data.frame(Var1 = missingSurvYrs,
                                                 Year = missingSurvYrs))
    } else {
      survAddBit[1:ncol(survAddBit)] <- NA
      survAddBit$Year <- NULL
      survAddBit <- merge(survAddBit, data.frame(Year = missingSurvYrs))
    }
    survDatat <- rbind(survData, survAddBit)
    survDatat <- survDatat[order(survDatat$Year), ]
  } else {
    survDatat <- survData
  }

  missingRecYrs <- setdiff(seq(inp$startYear, inp$endYear), data3$Year)

  if (length(missingRecYrs) > 0) {
    dat3Bit <- data3[1, ]
    dat3Bit[, 2:3] <- NA
    dat3Bit$Year <- NULL
    dat3Bit <- merge(dat3Bit, data.frame(Year = missingRecYrs))
    data3t <- rbind(data3, dat3Bit)
    data3t <- data3t[order(data3t$Year), ]
    data4t <- rbind(data4, dat3Bit)
    data4t <- data4t[order(data4t$Year), ]
  } else {
    data3t <- data3
    data4t <- data4
  }

  if (inp$adjustR) {
    adjustString <- "Rfemale[k] <- (composition.bias*R[k]/2)/(1+(composition.bias*R[k]/2))"
  } else {
    adjustString <- "Rfemale[k] <- composition.bias*R[k]/2"
  }

  if (inp$survAnalysisMethod == "KaplanMeier") {
    survString <- "Surv[surv_id[k]] ~ dnorm(S.annual.KM[surv_id[k]], tau[surv_id[k]])T(0,1)"
  } else {
    survString <- paste(c("for(t in 1:12){", "surv[surv_id[k],t+1] ~ dbern(S.annual.KM[survYr[surv_id[k]]]^(1/12)*surv[surv_id[k],t])", "}"), collapse = "\n")
  }
  
  if(is.na(betaPriors$bias.Prior1)){
    biasString <- "composition.bias <- bias.Prior1+0*bias.Prior2"
  }else{
    if(betaPriors$bias.Prior2==0){
      biasString <- "composition.bias <- exp(bias.Prior1+0*bias.Prior2)"
    }else{
      biasString <- paste0(c("lbias~dnorm(bias.Prior1,pow(bias.Prior2,-2))","composition.bias <- exp(lbias)"),collapse="\n")
    }
  }

  jagsTemplate <- paste(readLines(system.file("templates/JAGS_template.txt",
                                              package = "caribouMetrics"
  )), collapse = "\n")
  jagsTemplate <- gsub("_survString_", survString, jagsTemplate, fixed = T)
  jagsTemplate <- gsub("_adjustString_", adjustString, jagsTemplate, fixed = T)
  jagsTemplate <- gsub("_biasString_", biasString, jagsTemplate, fixed = T)
  
  jagsFile <- file.path(saveJAGStxt, "JAGS_run.txt")

  sink(jagsFile)
  cat(jagsTemplate, fill = TRUE)
  sink()

  sp.data <- list(
    anthro = disturbance$Anthro,
    fire = disturbance$fire_excl_anthro,
    beta.Saf.Prior1 = betaPriors$beta.Saf.Prior1,
    beta.Saf.Prior2 = betaPriors$beta.Saf.Prior2,
    beta.Rec.anthro.Prior1 = betaPriors$beta.Rec.anthro.Prior1,
    beta.Rec.anthro.Prior2 = betaPriors$beta.Rec.anthro.Prior2,
    beta.Rec.fire.Prior1 = betaPriors$beta.Rec.fire.Prior1,
    beta.Rec.fire.Prior2 = betaPriors$beta.Rec.fire.Prior2,
    l.Saf.Prior1 = betaPriors$l.Saf.Prior1,
    l.Saf.Prior2 = betaPriors$l.Saf.Prior2,
    l.R.Prior1 = betaPriors$l.R.Prior1,
    l.R.Prior2 = betaPriors$l.R.Prior2,
    sig.Saf.Prior1 = betaPriors$sig.Saf.Prior1,
    sig.Saf.Prior2 = betaPriors$sig.Saf.Prior2,
    sig.R.Prior1 = betaPriors$sig.R.Prior1,
    sig.R.Prior2 = betaPriors$sig.R.Prior2,
    bias.Prior1 = betaPriors$bias.Prior1,
    bias.Prior2 = betaPriors$bias.Prior2,
    Ninit = inp$N0,
    assessmentYrs = inp$assessmentYrs,
    nCounts = length(which(is.na(data3t$Count) == FALSE)),
    count_id = which(is.na(data3t$Count) == FALSE),
    nYears = inp$endYear - inp$startYear + 1,
    calves = round(data3t$Count),
    CountAntlerless = round(data4t$Count)
  )

  if (inp$survAnalysisMethod == "KaplanMeier") {
    sp.data <- c(sp.data, list(
      Surv = survDatat$surv, tau = survDatat$tau,
      nSurvs = length(which(is.na(survDatat$surv) == FALSE)),
      surv_id = which(is.na(survDatat$surv) == FALSE)
    ))
  } else {
    sp.data <- c(sp.data, list(
      surv = survDatat[, 1:13], survYr = survDatat$Year - inp$startYear + 1,
      nSurvs = length(which(is.na(survDatat[, 1]) == FALSE)),
      surv_id = which(is.na(survDatat$Year) == FALSE)
    ))
  }

  sp.params <- c("S.annual.KM", "R", "Rfemale", "pop.growth",
                 "fpop.size", "l.R", "l.Saf",
                 "beta.Rec.anthro", "beta.Rec.fire", "beta.Saf","composition.bias")

  
  
  rr.surv <- try(R2jags::jags(
    data = sp.data, parameters.to.save = sp.params,
    model.file = jagsFile,
    n.chains = inp$Nchains, n.iter = inp$Niter, n.burnin = inp$Nburn,
    n.thin = inp$Nthin, quiet = quiet
  ))
  
  ageRatioIn <- rbind(
    mutate(data3t, Class = "calf"), 
    mutate(data4t, Class = "cow")
  ) %>% arrange(.data$Year)

  return(list(result = rr.surv, 
              inData = list(survDataIn = survDatat, 
                            disturbanceIn = disturbance, 
                            ageRatioIn = ageRatioIn)))
}
