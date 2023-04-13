#' Run recruitment and mortality model
#'
#'
#'
#' @param survData either a path to a csv file or a dataframe containing the
#'   columns "Year", "event", "enter" and "exit". Enter and exit are the
#'   beginning and end of the time interval and should be a number from 1 to 12
#'   (December) or 0. Event is 0 or 1 where 0 means the animal lived and 1 died.
#'   See [survival::Surv()] for more details.
#' @param ageRatio.herd either a path to a csv file or a dataframe containing
#'   the columns "Year","Count", and "Class". Where class can be either calf or
#'   cow
#' @param disturbance either a path to a csv file or a dataframe containing the
#'   columns "Anthro","fire_excl_anthro", and "Year"
#' @param betaPriors a list of model priors see [getPriors()]
#' @param startYear,endYear year defining the beginning of the observation
#'   period and the end of the projection period.
#' @param Nchains Number of chains for the Bayesian model
#' @param Niter Number of iterations for the Bayesian model
#' @param Nburn Length of burn-in for the Bayesian model
#' @param Nthin Thinning rate for the Bayesian model
#' @param N0 Initial population size
#' @param survAnalysisMethod Survival analysis method either "KaplanMeier" or
#'   "Exponential". Exponential is only recommended for small sample sizes
#' @inheritParams popGrowthJohnson
#' @param assessmentYrs Number of years over which to assess lambda (growth rate)
#' @param inpFixed an optional list of inputs with names matching the above, if
#'   an argument is included in this list it will override the named argument
#' @param saveJAGStxt file path. Directory where the JAGS model txt files will
#'   be saved. Default is `tempdir()`
#' @param quiet logical. Should jags be run quietly?
#'
#' @return a list with elements:
#'   * result: an `rjags` model object see [R2jags::jags()].
#'   * survInput: a data.frame with the input data used in the model.
#'
#' @export
#'
#' @examples
#' # this uses example data shipped with the package
#' runRMModel(startYear = 2009, Nchains = 1, Niter = 100, Nburn = 10, Nthin = 2)
#' 
runRMModel <- function(survData = system.file("extdata/simSurvData.csv",
                                              package = "caribouMetrics"),
                       ageRatio.herd = system.file("extdata/simAgeRatio.csv",
                                                   package = "caribouMetrics"),
                       disturbance = system.file("extdata/simDisturbance.csv",
                                                 package = "caribouMetrics"),
                       betaPriors = "default",
                       startYear = 1998, endYear = 2023, Nchains = 4,
                       Niter = 15000, Nburn = 10000, Nthin = 2, N0 = 1000,
                       survAnalysisMethod = "KaplanMeier", adjustR = F,
                       assessmentYrs = 1,
                       inpFixed = list(), saveJAGStxt = tempdir(),
                       quiet = TRUE) {
  # survData=oo$simSurvObs;ageRatio.herd=oo$ageRatioOut;disturbance=oo$simDisturbance;
  # betaPriors=betaPriors;startYear = minYr;endYear=maxYr;N0=cs$N0;survAnalysisMethod = "KaplanMeier"
  # Nchains = 2;Niter = 20000;Nburn = 10000;Nthin = 1;assessmentYrs = 3;inpFixed=list()

  # combine defaults in function with inputs from input list
  inputArgs <- c(
    "survData", "ageRatio.herd", "disturbance", "startYear", "endYear",
    "Nchains", "Niter", "Nburn", "Nthin", "N0", "survAnalysisMethod", "adjustR",
    "assessmentYrs"
  )
  addArgs <- inputArgs # setdiff(inputArgs,names(inp))
  inp <- list()
  for (a in addArgs) {
    if (is.element(a, names(inpFixed))) {
      inp[[a]] <- inpFixed[[a]]
    } else {
      inp[[a]] <- eval(parse(text = a))
    }
  }

  if (betaPriors[[1]] == "default") {
    betaPriors <- getPriors()
  }

  # Run model
  if (is.character(inp$ageRatio.herd)) {
    ageRatio.herd <- read.csv(inp$ageRatio.herd, header = T)
    ageRatio.herd$X <- NULL
  } else {
    ageRatio.herd <- inp$ageRatio.herd
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
  testTable(ageRatio.herd, c("Year", "Count", "Class"))
  testTable(survData, c("Year", "event", "enter", "exit"))

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
  dSubset <- subset(data, enter == 0) # ;dSubset=subset(dSubset,!((exit<12)&(event==0)))
  if (nrow(dSubset) == 0) {
    stop("Collars not present at the start of a year are omitted from survival ",
         "analysis in that year. Please ensure there is at least one year with",
         " a collar in the first month.")
  }

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
      group_by(Year) %>%
      summarise(minEnter = min(enter), maxExit = max(exit))
    includeYrs <- subset(sumDat, minEnter == 0 & maxExit == 12)
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
      nriskYears <- data.frame(with(data.sub, table(Year)))

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
    dExpand <- apply(subset(dSubset, select = c(id, Year, event, enter, exit)),
                     1, expandSurvivalRecord)
    survData <- do.call(rbind, dExpand)
  }

  # split data into calf and cow recruitment data
  calf.cnt <- subset(ageRatio.herd, ageRatio.herd$Class == "calf")
  calf.cnt$Class <- factor(calf.cnt$Class)
  cow.cnt <- subset(ageRatio.herd, ageRatio.herd$Class == "cow")
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
    adjustString <- "Rfemale[k] <- (RT[k]/2)/(1+(RT[k]/2))"
  } else {
    adjustString <- "Rfemale[k] <- RT[k]/2"
  }

  if (inp$survAnalysisMethod == "KaplanMeier") {
    survString <- "Surv[surv_id[k]] ~ dnorm(S.annual.KM[surv_id[k]], tau[surv_id[k]])"
  } else {
    survString <- paste(c("for(t in 1:12){", "surv[surv_id[k],t+1] ~ dbern(S.annual.KM[survYr[surv_id[k]]]^(1/12)*surv[surv_id[k],t])", "}"), collapse = "\n")
  }

  jagsTemplate <- paste(readLines(system.file("templates/JAGS_template.txt",
                                              package = "caribouMetrics"
  )), collapse = "\n")
  jagsTemplate <- gsub("_survString_", survString, jagsTemplate, fixed = T)
  jagsTemplate <- gsub("_adjustString_", adjustString, jagsTemplate, fixed = T)

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
                 "fpop.size", "var.R.real", "l.R", "l.Saf",
                 "beta.Rec.anthro", "beta.Rec.fire", "beta.Saf")
  rr.surv <- try(R2jags::jags(
    data = sp.data, parameters.to.save = sp.params,
    model.file = jagsFile,
    n.chains = inp$Nchains, n.iter = inp$Niter, n.burnin = inp$Nburn,
    n.thin = inp$Nthin, quiet = quiet
  ))

  return(list(result = rr.surv, survInput = survDatat))
}
