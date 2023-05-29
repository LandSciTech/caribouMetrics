# internal functions related to caribou demographics


# Helpers for caribouBayesianIPM --------------------------------------------------
# Copyright 2023 Daniel Eacker & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from code provided in https://doi.org/10.1002/wsb.950

getKMSurvivalEstimates <- function(dSubset) {
  sModel <- survival::survfit(survival::Surv(enter, exit, event) ~ as.factor(Year),
                              conf.type = "log-log", data = dSubset)
  reg.out <- summary(sModel)
  # TODO: Josie can this be deleted?
  if (0) {
    # check <- data.frame(strata = reg.out$strata, survival = reg.out$surv,
    #                     time = reg.out$time)
    # check$type <- "est"
    # check$Year <- as.numeric(gsub("as.factor(Year)=", "", as.character(check$strata),
    #                               fixed = T))
    # check$strata <- NULL
    # tt <- subset(oo$exData, select = c(Year, survival))
    # tt$time <- 12
    # tt$type <- "truth"
    # check <- rbind(check, tt)
    # base <- ggplot2::ggplot(check,
    #                         ggplot2::aes(x = time, y = survival, colour = type,
    #                                      shape = type)) +
    #   ggplot2::geom_point() +
    #   ggplot2::facet_wrap(~Year)
    # print(base)
  }

  if (is.null(reg.out$strata)) {
    reg.out$strata <- as.character(dSubset$Year[1])
  }
  data5 <- data.frame(reg.out$strata, reg.out$surv)
  data.se <- data.frame(reg.out$strata, reg.out$std.err)
  data.l <- data.frame(reg.out$strata, reg.out$lower)
  data.u <- data.frame(reg.out$strata, reg.out$upper)
  data6 <- data.frame(table(data5$reg.out.strata))
  num <- data6$Freq
  if (inherits(data5$reg.out.strata, "character")) {
    data5$reg.out.strata <- as.factor(data5$reg.out.strata)
  }
  levs <- data.frame(levels(data5$reg.out.strata))
  names(levs) <- c("reg.out.strata")
  data7 <- subset(data6, data6$Freq > 0)
  survs <- numeric(length(data7$Freq))
  se <- numeric(length(data7$Freq))
  lower <- numeric(length(data7$Freq))
  upper <- numeric(length(data7$Freq))
  index <- cumsum(data7$Freq)

  for (i in 1:length(survs)) {
    survs[i] <- data5$reg.out.surv[index[i]]
    se[i] <- data.se$reg.out.std.err[index[i]]
    lower[i] <- data.l$reg.out.lower[index[i]]
    upper[i] <- data.u$reg.out.upper[index[i]]
  }

  indexS <- which(data6$Freq != 0)
  data6$surv <- numeric(length(data6$Freq))
  data6$se <- numeric(length(data6$Freq))
  data6$lower <- numeric(length(data6$Freq))
  data6$upper <- numeric(length(data6$Freq))
  data6$surv <- ifelse(data6$Freq < 1, 1, NA)
  data6$lower <- ifelse(data6$Freq < 1, NA, 0)
  data6$upper <- ifelse(data6$Freq < 1, NA, 0)

  for (i in 1:length(data6$Freq)) {
    data6$surv[indexS[i]] <- survs[i]
    data6$se[indexS[i]] <- se[i]
    data6$lower[indexS[i]] <- lower[i]
    data6$upper[indexS[i]] <- upper[i]
  }

  survData <- data6

  return(survData)
}

# Takes enter exit table and converts it to a larger table of 1s and 0s for all
# months. This is large but can run when there is very little info and not
# enough for KM. Makes more assumptions about survival in every month being the
# same
expandSurvivalRecord <- function(crow) {
  # crow=subset(dSubset,exit==12)[1,]
  crow <- as.numeric(crow)
  mnths <- c(rep(1, crow[5] - crow[4]), rep(!crow[3], 12 - crow[5] + 1))
  mnths <- data.frame(t(mnths))
  mnths$Year <- crow[2]
  mnths$enter <- crow[4]
  mnths$exit <- crow[5]
  mnths$event <- crow[3]
  return(mnths)
}

# Helpers for simulateObservations -----------------------------------------

simTrajectory <- function(numYears, covariates, survivalModelNumber = "M1",
                          recruitmentModelNumber = "M4",
                          popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                          recSlopeMultiplier = 1, sefSlopeMultiplier = 1,
                          recQuantile = 0.5, sefQuantile = 0.5,
                          stepLength = 1, N0 = 1000, adjustR = T) {
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4";
  # recSlopeMultiplier=1;sefSlopeMultiplier=1;recQuantile=0.5;sefQuantile=0.5
  # stepLength=1;N0=1000

  # alter coefficients
  growthTab <- popGrowthTable

  growthTab$Value[(growthTab$Coefficient == "Anthro") &
                    (growthTab$responseVariable == "recruitment")] <-
    recSlopeMultiplier * growthTab$Value[(growthTab$Coefficient == "Anthro") &
                                           (growthTab$responseVariable == "recruitment")]

  growthTab$Value[(growthTab$Coefficient == "Anthro") &
                    (growthTab$responseVariable == "femaleSurvival")] <-
    sefSlopeMultiplier * growthTab$Value[(growthTab$Coefficient == "Anthro") &
                                           (growthTab$responseVariable == "femaleSurvival")]

  popGrowthParsSmall <- demographicCoefficients(
    2,
    modelVersion = "Johnson",
    survivalModelNumber = survivalModelNumber,
    recruitmentModelNumber = recruitmentModelNumber,
    populationGrowthTable = growthTab,
    useQuantiles = c(recQuantile, recQuantile)
  )
  # manually set quantiles for example population
  popGrowthParsSmall$coefSamples_Survival$quantiles <- sefQuantile

  # Only use precision if included in the table for this model number for both
  # rec and surv
  usePrec <- "Precision" %in% names(popGrowthParsSmall$coefSamples_Survival$coefValues) &
    "Precision" %in% names(popGrowthParsSmall$coefSamples_Recruitment$coefValues)
  # at each time,  sample demographic rates and project, save results
  # TODO: SE thinks this can be done all at once with a table of demographic rates 
  pars <- data.frame(N0 = N0)
  for (t in 1:numYears) {
    # t=1
    covs <- subset(covariates, time == t)

    rateSamples <- demographicRates(
      covTable = covs,
      popGrowthPars = popGrowthParsSmall,
      ignorePrecision = !usePrec,
      returnSample = TRUE
    )
    if (is.element("N", names(pars))) {
      pars <- subset(pars, select = c("replicate", "N"))
      names(pars)[names(pars) == "N"] <- "N0"
    }
    pars <- merge(pars, rateSamples)
    pars <- cbind(
      pars,
      caribouPopGrowth(pars$N0,
                       R_bar = pars$R_bar, S_bar = pars$S_bar,
                       numSteps = stepLength, K = FALSE, l_R = 1e-06, adjustR = adjustR, 
                       progress = FALSE
      )
    )

    # add results to output set
    fds <- subset(pars, select = c("replicate", "Anthro", "fire_excl_anthro",
                                   "S_t", "R_t", "N",
                                   "lambda", "n_recruits", "surviving_adFemales"))
    fds$replicate <- as.numeric(gsub("V", "", fds$replicate))
    names(fds) <- c("Replicate", "Anthro", "fire_excl_anthro", "survival",
                    "recruitment", "N", "lambda", "n_recruits", "n_cows")
    # apparent number of calves per cow, not actual, from unadjusted R_t
    fds$n_recruits <- fds$recruitment * fds$n_cows
    fds <- tidyr::pivot_longer(fds, !.data$Replicate, names_to = "MetricTypeID",
                               values_to = "Amount")
    fds$Timestep <- t * stepLength
    if (t == 1) {
      popMetrics <- fds
    } else {
      popMetrics <- rbind(popMetrics, fds)
    }
  }

  popMetrics$MetricTypeID <- as.character(popMetrics$MetricTypeID)
  popMetrics$Replicate <- paste0("x", popMetrics$Replicate)
  return(subset(popMetrics, popMetrics$Replicate == "x1"))
}

simCalfCowRatios <- function(cowCounts, minYr, exData) {
  # assume info from only one herd
  ageRatioOut <- subset(cowCounts, (cowCounts$Year >= minYr),
                        select = c("Year", "Class", "Count"))
  ageRatioOut <- tidyr::pivot_wider(ageRatioOut, id_cols = c("Year"),
                                    names_from = "Class", values_from = "Count")
  ageRatioOut <- merge(ageRatioOut,
                       subset(exData, select = c("Year", "n_recruits", "n_cows")))
  # rbinom needs n_recruits to be <= n_cows and n_cows not 0
  n_recs <- pmin(ageRatioOut$n_cows, ageRatioOut$n_recruits)
  ageRatioOut$calf <- ifelse(ageRatioOut$n_cows == 0, 0,
                             rbinom(
                               n = nrow(ageRatioOut), size = ageRatioOut$cow,
                               prob = n_recs / ageRatioOut$n_cows
                             )
  )
  ageRatioOut <- subset(ageRatioOut, select = c("Year", "calf", "cow"))
  ageRatioOut <- tidyr::pivot_longer(ageRatioOut, cols = c(.data$calf, .data$cow),
                                     names_to = "Class", values_to = "Count")
  return(ageRatioOut)
}

simSurvivalData <- function(freqStartsByYear, exData, collarNumYears, collarOffTime,
                            collarOnTime, topUp = FALSE) {
  # topUp=T
  # for simplicity, ignore variation in survival probability among months
  initYear <- min(exData$Year)
  freqStartsByYear <- subset(freqStartsByYear, 
                             (freqStartsByYear$Year >= initYear) & 
                               (freqStartsByYear$numStarts > 0))
  
  freqStartsByYear <- freqStartsByYear[order(freqStartsByYear$Year), ]
  survivalSeries <- subset(exData, select = c("survival", "Year"))
  
  # freqStartsByYear$numStarts=10000;collarNumYears=2
  
  # if(initYear<min(survivalSeries$Year)){
  #  missingYrs = seq(initYear,min(survivalSeries$Year)-1)
  #  addBit = subset(survivalSeries,Year==min(survivalSeries$Year))
  #  addBit$Year=NULL;addBit = merge(addBit, data.frame(Year=missingYrs))
  #  survivalSeries=rbind(survivalSeries,addBit)
  # }
  
  animalID <- 1
  simSurvObs <- data.frame(id = NA, Year = NA, event = NA, enter = NA, exit = NA)
  simSurvObs <- subset(simSurvObs, !is.na(id))
  # collarNumYears=4
  for (k in 1:nrow(freqStartsByYear)) {
    # k =1
    if (is.na(freqStartsByYear$numStarts[k]) | (freqStartsByYear$numStarts[k] <= 0)) {
      next
    }
    startYear <- freqStartsByYear$Year[k]
    if (topUp) {
      collarsExisting <- nrow(subset(simSurvObs, (simSurvObs$enter == 0) & 
                                       (simSurvObs$Year == startYear)))
      nstarts <- max(0, freqStartsByYear$numStarts[k] - collarsExisting)
    } else {
      nstarts <- freqStartsByYear$numStarts[k]
    }
    if (nstarts == 0) {
      next
    }
    for (n in 1:nstarts) {
      # n=1
      addS <- simSurvivalObs(animalID, startYear = startYear,
                             collarNumYears = collarNumYears,
                             collarOffTime = collarOffTime,
                             collarOnTime = collarOnTime,
                             survivalSeries = survivalSeries)
      animalID <- animalID + 1
      simSurvObs <- rbind(simSurvObs, addS)
    }
  }
  
  # 1-sum(simSurvObs$event)/nrow(simSurvObs)
  # exData
  
  simSurvObs <- subset(simSurvObs, is.element(simSurvObs$Year, exData$Year))
  simSurvObs <- simSurvObs[order(simSurvObs$Year), ]
  
  addBit <- unique(subset(freqStartsByYear,
                          select = setdiff(names(freqStartsByYear),
                                           c("numStarts", names(simSurvObs)))))
  if (nrow(addBit) > 1) {
    stop()
  } else if (nrow(addBit) > 0) {
    simSurvObs <- merge(simSurvObs, addBit)
  }
  
  return(simSurvObs)
}

simSurvivalObs <- function(animalID, startYear, collarNumYears, collarOffTime,
                           collarOnTime, survivalSeries) {
  # animalID =  1; startYear = 2016
  for (i in startYear:min((startYear + collarNumYears - 1), max(survivalSeries$Year))) {
    # i = startYear
    surv <- survivalSeries$survival[survivalSeries$Year == i]^(1 / 12) # monthly survival
    
    if (i == startYear) {
      enter <- collarOnTime - 1
    } else {
      enter <- 0
    }
    
    if (i == (startYear + collarNumYears - 1)) {
      exit <- collarOffTime
    } else {
      exit <- 12
    }
    
    event <- 0
    for (j in (enter + 1):exit) {
      die <- rbinom(1, 1, 1 - surv)
      if (die) {
        event <- 1
        exit <- j
        break
      }
    }
    
    addBit <- data.frame(id = animalID, Year = i, event = die, enter = enter,
                         exit = exit)
    
    if (i == startYear) {
      survObs <- addBit
    } else {
      survObs <- rbind(survObs, addBit)
    }
    if (die) {
      break
    }
  }
  
  return(survObs)
}

# Helpers for runScnSet and App -------------------------------------------
# Copyright 2023 Daniel Eacker & Her Majesty the Queen in Right of Canada as represented by the Minister of the Environment
# License GPL-3
#NOTICE: This function has been modified from code provided in https://doi.org/10.1002/wsb.950

getSumStats <- function(param, rrSurvMod, startYear, endYear, doSummary = T) {
  # param = "pop.growth";doSummary=T
  paramNames <- data.frame(
    parameter = c(
      "S.annual.KM", "R", "Rfemale", "pop.growth", "fpop.size",
      "meanAFsurv", "meanR", "meanRfemale",
      "medianLambda", "meanLambda"
    ),
    name = c(
      "Adult female survival", "Recruitment",
      "Female-only recruitment", "Population growth rate", "Female population size",
      "Mean adult female survival",
      "Mean recruitment", "Mean female recruitment",
      "Median population growth rate",
      "Mean population growth rate"
    )
  )

  paramNames <- subset(paramNames, is.element(paramNames$parameter, rrSurvMod$parameters.to.save))

  if (grepl("mean|median", param)) {
    yr <- NA_integer_
  } else {
    yr <- startYear:endYear
  }

  if (!param %in% paramNames$parameter) {
    stop(
      "param ", param, "is not recognized\n",
      "accepted params are: ", paramNames$parameter
    )
  }

  if (doSummary) {
    lower.cri <- apply(
      rrSurvMod$BUGSoutput$sims.list[[param]], 2,
      function(x) {
        quantile(x, 0.025)
      }
    )
    upper.cri <- apply(
      rrSurvMod$BUGSoutput$sims.list[[param]], 2,
      function(x) {
        quantile(x, 0.975)
      }
    )
    probViable <- apply(
      rrSurvMod$BUGSoutput$sims.list[[param]], 2,
      function(x) {
        mean(x > 0.99)
      }
    )

    results <- data.frame(
      Year = yr,
      Parameter = subset(paramNames, paramNames$parameter == param,
                         select = "name", drop = T
      ),
      Mean = round(rrSurvMod$BUGSoutput$mean[[param]], digits = 3),
      SD = round(rrSurvMod$BUGSoutput$sd[[param]], digits = 3),
      `Lower 95% CRI` = round(lower.cri, digits = 3),
      `Upper 95% CRI` = round(upper.cri, digits = 3),
      probViable = round(probViable, digits = 3),
      check.names = FALSE
    )
    if(!grepl("growth rate", results$Parameter[1])){
      results$probViable <- NA
    }
  } else {
    # rrSurvMod= result
    wideRes <- data.frame(rrSurvMod$BUGSoutput$sims.list[[param]])
    names(wideRes) <- yr

    results <- wideRes %>%
      tidyr::pivot_longer(cols = names(wideRes), names_to = "Year",
                          values_to = "Value")
    results$Year <- as.numeric(results$Year)
    results$Parameter <- subset(paramNames, paramNames$parameter == param,
                                select = "name", drop = T
    )
    results <- as.data.frame(results)
  }
  return(results)
}



movingAveGrowthRate <- function(obs, assessmentYrs) {
  # obs=obsLam
  if (assessmentYrs == 1) {
    return(obs)
  }
  obsOut <- obs
  for (k in assessmentYrs:nrow(obsOut)) {
    # k=3
    obsOut$Mean[k] <- mean(obs$Mean[(k - assessmentYrs + 1):k])
  }
  obsOut
}

getKSDist <- function(Value, type) {
  # sampleBit=subset(allSamples,(Year==2017)&(Parameter==allSamples$Parameter[1]))
  # Value=sampleBit$Value;type=sampleBit$type

  if (length(Value[type == "national"]) == 0) {
    out <- data.frame(KSDistance = NA, KSpvalue = NA)
    return(out)
  }
  res <- ks.test(Value[type == "local"], Value[type == "national"])

  out <- data.frame(KSDistance = res$statistic, KSpvalue = res$p.value)
  return(out)
}


# General helpers ---------------------------------------------------------

# test a popGrowthTable has the right format
testPopGrowthTable <- function(df) {
  # required columns
  missed_nms <- setdiff(
    c("responseVariable", "Coefficient", "Value"),
    colnames(df)
  )

  # need stdErr or CI
  if (!"StdErr" %in% colnames(df)) {
    missed_nms <- c(
      missed_nms,
      setdiff(c("lowerCI", "upperCI"), colnames(df))
    )
    df <- mutate(df, StdErr = NA_real_)
  } else if (!all(c("lowerCI", "upperCI") %in% colnames(df))) {
    df <- mutate(df, lowerCI = NA_real_, upperCI = NA_real_)
  }

  if (length(missed_nms) > 0) {
    stop("The model coefficient file loaded is missing the columns ",
         paste(missed_nms, collapse = ", "),
         call. = FALSE
    )
  }

  # should only give one model number in table because no method to select a mod num
  if (!is.null(df$ModelNumber)) {
    nmod <- df %>%
      group_by(.data$responseVariable) %>%
      summarise(nmod = n_distinct(.data$ModelNumber)) %>%
      pull(.data$nmod)

    if (any(nmod > 1)) {
      stop("The model coefficient file loaded contains more than one model per response variable",
           call. = FALSE
      )
    }
  }

  # add expected columns that should never change
  df <- mutate(df,
               modelVersion = "Johnson",
               ModelNumber = ifelse(.data$responseVariable == "recruitment", "M4", "M1"),
               Type = "National"
  )

  # expected values
  diff_res <- setdiff(unique(df$responseVariable),
                      unique(caribouMetrics::popGrowthTableJohnsonECCC$responseVariable))

  if (!setequal(unique(caribouMetrics::popGrowthTableJohnsonECCC$responseVariable),
                unique(df$responseVariable))) {
    stop("The model coefficient file loaded contains unrecognized responseVariable: ",
         paste(diff_res, collapse = ", "),
         call. = FALSE
    )
  }

  diff_coef <- setdiff(unique(df$Coefficient), c(
    "Intercept", "Precision",
    "Anthro", "fire_excl_anthro"
  ))

  if (length(diff_coef) > 0) {
    stop("The model coefficient file loaded contains unrecognized Coefficient: ",
         paste(diff_coef, collapse = ", "),
         call. = FALSE
    )
  }

  testStdCI <- df %>% mutate(stdOrCI = !is.na(.data$StdErr) | (!is.na(.data$lowerCI) & !is.na(.data$upperCI)))

  if (!all(testStdCI$stdOrCI)) {
    stop("The model coefficient file loaded is missing StdErr or lowerCI and upperCI for:\n",
         "femaleSurvival: ",
         testStdCI %>% filter(!.data$stdOrCI, .data$responseVariable == "femaleSurvival") %>%
           pull(.data$Coefficient) %>% paste0(collapse = ", "), "\n",
         "recruitment: ",
         testStdCI %>% filter(!.data$stdOrCI, .data$responseVariable == "recruitment") %>%
           pull(.data$Coefficient) %>% paste0(collapse = ", "),
         call. = FALSE
    )
  }

  return(df)
}


#' Test table
#'
#' Test has expected column names and optionally that certain columns have
#' expected values
#'
#' @param df data.frame. The table to test
#' @param req_col_names character. Required column names. A vector of column
#'   names that must be present in `df`. Other columns are allowed
#' @param req_vals list.  A named list where the name is a column name and the
#'   value is a vector of required values. Values in the list and not in the
#'   column will throw an error
#' @param acc_vals list. A named list where the name is a column name and the
#'   value is vector of accepted values. Any values not in the list will throw
#'   an error.
#'
#' @return throws an error if failed otherwise invisible NULL
#'
#' @noRd
testTable <- function(df, req_col_names, req_vals = NULL, acc_vals = NULL){
  df_name <- deparse(substitute(df))
  missing_cols <- setdiff(req_col_names, colnames(df))
  if(length(missing_cols) > 0){
    stop(df_name, " is missing expected columns: ",
         paste0(missing_cols, collapse = ", "), call. = FALSE)
  }

  if(!is.null(req_vals)){
    Map(function(x, nm){
      missing_vals <- setdiff(x, df[[nm]])
      if(length(missing_vals) > 0){
        stop(df_name, "$", nm, " is missing expected values: ",
             paste0(missing_vals, collapse = ", "), call. = FALSE)
      }
    }, req_vals, names(req_vals))
  }

  if(!is.null(acc_vals)){
    Map(function(x, nm){
      wrong_vals <- setdiff(df[[nm]], x)
      if(length(wrong_vals) > 0){
        stop(df_name, "$", nm, " contains unexpected values: ",
             paste0(wrong_vals, collapse = ", "), ".\n Expected values: ",
             paste0(x, collapse = ", "), call. = FALSE)
      }
    }, acc_vals, names(acc_vals))
  }
  return(invisible(NULL))
}

# saves the cached version of the national sims to a local file so it doesn't
# need to be re-run every time you restart
savePersistentCache <- function(env = cacheEnv){
  obj_nms <- ls(envir = env)
  lapply(obj_nms, function(x){
    obj <- get(x, envir=env)
    saveRDS(obj, paste0("inst/extdata/", x, ".rds"))
  })
  return(invisible())
}
