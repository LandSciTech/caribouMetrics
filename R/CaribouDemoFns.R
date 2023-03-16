getKMSurvivalEstimates <- function(dSubset) {
  sModel <- survival::survfit(survival::Surv(enter, exit, event) ~ as.factor(Year), conf.type = "log-log", data = dSubset)
  reg.out <- summary(sModel)
  if (0) {
    check <- data.frame(strata = reg.out$strata, survival = reg.out$surv, time = reg.out$time)
    check$type <- "est"
    check$Year <- as.numeric(gsub("as.factor(Year)=", "", as.character(check$strata), fixed = T))
    check$strata <- NULL
    tt <- subset(oo$exData, select = c(Year, survival))
    tt$time <- 12
    tt$type <- "truth"
    check <- rbind(check, tt)
    base <- ggplot2::ggplot(check, ggplot2::aes(x = time, y = survival, colour = type, shape = type)) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~Year)
    print(base)
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
  if (class(data5$reg.out.strata) == "character") {
    data5$reg.out.strata <- as.factor(data5$reg.out.strata)
  }
  levs <- data.frame(levels(data5$reg.out.strata))
  names(levs) <- c("reg.out.strata")
  data7 <- subset(data6, Freq > 0)
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

# install and load required packages
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

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


simTrajectory <- function(numYears, covariates, survivalModelNumber = "M1", recruitmentModelNumber = "M4",
                          popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                          recSlopeMultiplier = 1, sefSlopeMultiplier = 1, recQuantile = 0.5, sefQuantile = 0.5,
                          stepLength = 1, N0 = 1000, adjustR = T) {
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4";
  # recSlopeMultiplier=1;sefSlopeMultiplier=1;recQuantile=0.5;sefQuantile=0.5
  # stepLength=1;N0=1000

  # alter coefficients
  growthTab <- popGrowthTable
  growthTab$Value[(growthTab$Coefficient == "Anthro") & (growthTab$responseVariable == "recruitment")] <- recSlopeMultiplier * growthTab$Value[(growthTab$Coefficient == "Anthro") & (growthTab$responseVariable == "recruitment")]
  growthTab$Value[(growthTab$Coefficient == "Anthro") & (growthTab$responseVariable == "femaleSurvival")] <- sefSlopeMultiplier * growthTab$Value[(growthTab$Coefficient == "Anthro") & (growthTab$responseVariable == "femaleSurvival")]

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
  # TO DO: across full covariate range, note proportion of this scenario that is outside distribution of observations across the country

  # Only use precision if included in the table for this model number for both rec and surv
  usePrec <- "Precision" %in% names(popGrowthParsSmall$coefSamples_Survival$coefValues) &
    "Precision" %in% names(popGrowthParsSmall$coefSamples_Recruitment$coefValues)
  # at each time,  sample demographic rates and project, save results
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
      pars <- subset(pars, select = c(replicate, N))
      names(pars)[names(pars) == "N"] <- "N0"
    }
    pars <- merge(pars, rateSamples)
    pars <- cbind(
      pars,
      popGrowthJohnson(pars$N0,
        R_bar = pars$R_bar, S_bar = pars$S_bar,
        numSteps = stepLength, K = F, l_R = 1e-06, adjustR = adjustR
      )
    )

    # add results to output set
    fds <- subset(pars, select = c(replicate, Anthro, fire_excl_anthro, S_t, R_t, N, lambda, n_recruits, surviving_adFemales))
    fds$replicate <- as.numeric(gsub("V", "", fds$replicate))
    names(fds) <- c("Replicate", "Anthro", "fire_excl_anthro", "survival", "recruitment", "N", "lambda", "n_recruits", "n_cows")
    fds$n_recruits <- fds$recruitment * fds$n_cows # apparent number of calves per cow, not actual, from unadjusted R_t
    fds <- tidyr::pivot_longer(fds, !Replicate, names_to = "MetricTypeID", values_to = "Amount")
    fds$Timestep <- t * stepLength
    if (t == 1) {
      popMetrics <- fds
    } else {
      popMetrics <- rbind(popMetrics, fds)
    }
  }

  popMetrics$MetricTypeID <- as.character(popMetrics$MetricTypeID)
  popMetrics$Replicate <- paste0("x", popMetrics$Replicate)
  return(subset(popMetrics, Replicate == "x1"))
}

simCalfCowRatios <- function(cowCounts, minYr, exData) {
  ageRatioOut <- subset(cowCounts, (Year >= minYr), select = c(Year, Class, Count)) # assume info from only one herd
  ageRatioOut <- tidyr::pivot_wider(ageRatioOut, id_cols = c("Year"), names_from = "Class", values_from = "Count")
  ageRatioOut <- merge(ageRatioOut, subset(exData, select = c("Year", "n_recruits", "n_cows")))
  # rbinom needs n_recruits to be <= n_cows and n_cows not 0
  n_recs <- pmin(ageRatioOut$n_cows, ageRatioOut$n_recruits)
  ageRatioOut$calf <- ifelse(ageRatioOut$n_cows == 0, 0,
    rbinom(
      n = nrow(ageRatioOut), size = ageRatioOut$cow,
      prob = n_recs / ageRatioOut$n_cows
    )
  )
  ageRatioOut <- subset(ageRatioOut, select = c(Year, calf, cow))
  ageRatioOut <- tidyr::pivot_longer(ageRatioOut, cols = c(calf, cow), names_to = "Class", values_to = "Count")
  return(ageRatioOut)
}

simSurvivalData <- function(freqStartsByYear, exData, collarNumYears, collarOffTime, collarOnTime, topUp = F) {
  # topUp=T
  # for simplicity, ignore variation in survival probability among months
  initYear <- min(exData$Year)
  freqStartsByYear <- subset(freqStartsByYear, (Year >= initYear) & (numStarts > 0))

  freqStartsByYear <- freqStartsByYear[order(freqStartsByYear$Year), ]
  survivalSeries <- subset(exData, select = c(survival, Year))

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
      collarsExisting <- nrow(subset(simSurvObs, (enter == 0) & (Year == startYear)))
      nstarts <- max(0, freqStartsByYear$numStarts[k] - collarsExisting)
    } else {
      nstarts <- freqStartsByYear$numStarts[k]
    }
    if (nstarts == 0) {
      next
    }
    for (n in 1:nstarts) {
      # n=1
      addS <- simSurvivalObs(animalID, startYear = startYear, collarNumYears = collarNumYears, collarOffTime = collarOffTime, collarOnTime = collarOnTime, survivalSeries = survivalSeries)
      animalID <- animalID + 1
      simSurvObs <- rbind(simSurvObs, addS)
    }
  }

  # 1-sum(simSurvObs$event)/nrow(simSurvObs)
  # exData

  simSurvObs <- subset(simSurvObs, is.element(Year, exData$Year))
  simSurvObs <- simSurvObs[order(simSurvObs$Year), ]

  addBit <- unique(subset(freqStartsByYear, select = setdiff(names(freqStartsByYear), c("numStarts", names(simSurvObs)))))
  if (nrow(addBit) > 1) {
    stop()
  } else if (nrow(addBit) > 0) {
    simSurvObs <- merge(simSurvObs, addBit)
  }

  return(simSurvObs)
}

simSurvivalObs <- function(animalID, startYear, collarNumYears, collarOffTime, collarOnTime, survivalSeries) {
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

    addBit <- data.frame(id = animalID, Year = i, event = die, enter = enter, exit = exit)

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


# Tables ------------------------------------------------------------------

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

  paramNames <- subset(paramNames, is.element(parameter, rrSurvMod$parameters.to.save))

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
        select = name, drop = T
      ),
      Mean = round(rrSurvMod$BUGSoutput$mean[[param]], digits = 3),
      SD = round(rrSurvMod$BUGSoutput$sd[[param]], digits = 3),
      `Lower 95% CRI` = round(lower.cri, digits = 3),
      `Upper 95% CRI` = round(upper.cri, digits = 3),
      probViable = round(probViable, digits = 3),
      check.names = FALSE
    )
  } else {
    # rrSurvMod= result
    wideRes <- data.frame(rrSurvMod$BUGSoutput$sims.list[[param]])
    names(wideRes) <- yr

    results <- wideRes %>% tidyr::pivot_longer(cols = names(wideRes), names_to = "Year", values_to = "Value")
    results$Year <- as.numeric(results$Year)
    results$Parameter <- subset(paramNames, paramNames$parameter == param,
      select = name, drop = T
    )
  }
  return(results)
}

tabAllRes <- function(rrSurvMod, startYear, endYear, doSummary = T) {
  # rrSurvMod=rr.surv;startYear= minYr;endYear= maxYr
  # rrSurvMod=result;doSummary=T

  allParams <- c(
    "S.annual.KM", "R", "Rfemale", "pop.growth", "fpop.size",
    "meanAFsurv", "meanR", "meanRfemale",
    "medianLambda", "meanLambda"
  )
  allParams <- allParams[is.element(allParams, rrSurvMod$parameters.to.save)]

  allResults <- lapply(allParams, getSumStats, rrSurvMod, startYear, endYear, doSummary = doSummary)

  allResults <- do.call(rbind, allResults)

  allResults <- allResults[order(allResults$Year), ]
  allResults <- allResults[order(allResults$Parameter), ]
  row.names(allResults) <- 1:length(allResults$Year)
  allResults
}

getParamsFromEacker <- function(path) {
  ###############
  # Use Eacker example data for sample sizes in each year.
  survData <- paste0(path, "/tte_caribouFEMALES.csv")
  ageRatio.herd <- paste0(path, "/ageRatio.herd.csv")
  ageRatio.herd2 <- read.csv(ageRatio.herd, header = T)
  tte_caribou2 <- read.csv(survData, header = T)

  # need table of observed number of cows each year as input for simulating calf:cow ratios. Use Eaker as example.
  cowCounts <- subset(ageRatio.herd2, Class == "cow")
  write.csv(cowCounts, "tabs/cowCounts.csv")
  yrRange <- max(tte_caribou2$Year) - min(tte_caribou2$Year)

  # get survival sampling parameters from example data
  animalStarts <- subset(tte_caribou2, select = c(id, Year)) %>%
    group_by(id) %>%
    summarise(startYear = min(Year))
  freqStartsByYear <- as.data.frame(table(animalStarts$startYear))
  names(freqStartsByYear) <- c("Year", "numStarts")
  freqStartsByYear$Year <- as.numeric(as.character(freqStartsByYear$Year))
  sm <- tte_caribou2 %>%
    group_by(id) %>%
    summarize(startYear = min(Year), endYear = max(Year), numYears = max(Year) - min(Year), died = sum(event))
  collarNumYears <- median(subset(sm, !died)$numYears) # for simplicity, pick a single number of years that collars remain on

  addInfo <- unique(subset(tte_caribou2, select = c(HerdDescription, HerdCode, Range_ID, RangeDescription, RangeCode)))
  # TO DO: remove requirement for additional herd ID and range ID columns in UI code
  freqStartsByYear <- merge(freqStartsByYear, addInfo)
  write.csv(freqStartsByYear, "tabs/freqStartsByYear.csv")

  offSet <- subset(sm, !died, select = c(id, endYear))
  names(offSet) <- c("id", "Year")
  offSet <- merge(offSet, tte_caribou2, all.x = T)
  collarOffTime <- median(offSet$exit) # for simplicity, pick a single month that collars fall off

  onSet <- subset(sm, select = c(id, startYear))
  names(onSet) <- c("id", "Year")
  onSet <- merge(onSet, tte_caribou2, all.x = T)
  collarOnTime <- median(onSet$enter) # for simplicity, pick a single month that collars are put on
  # TO DO: allow users to set freqStartsByYear, collarNumYears, collarOffTime, and collarOnTime as parameters OR
  # provide a file formatted as in Eacker from which these parameters can be derived.

  # freqStartsByYear$numStarts=30
  return(list(cowCounts = cowCounts, freqStartsByYear = freqStartsByYear, collarOnTime = collarOnTime, collarOffTime = collarOffTime, collarNumYears = collarNumYears))
}



getOutputTables <- function(result, startYear, endYear, survInput, oo, simBig, getKSDists) {
  # result=out$result;startYear=minYr;endYear=maxYr;survInput=out$survInput;oo=oo;simBig=simBig

  # get summary info for plots
  rr.summary <- tabAllRes(result, startYear, endYear)

  if (!is.element("surv", names(survInput))) {
    if (sum(survInput$event, na.rm = T) > 0) {
      obsSurv <- getKMSurvivalEstimates(survInput)
    } else {
      obsSurv <- unique(subset(survInput, !is.na(enter), select = c(Year)))
      obsSurv$surv <- NA
      obsSurv$Var1 <- obsSurv$Year
    }
  } else {
    obsSurv <- survInput
  }

  obsSurv$Mean <- obsSurv$surv
  obsSurv$Year <- as.numeric(gsub("as.factor(Year)=", "", obsSurv$Var1, fixed = T))
  obsSurv <- subset(obsSurv, Year > 1000)

  obsSurv$parameter <- "Adult female survival"
  obsSurv$type <- "observed"

  trueSurv <- subset(oo$exData, select = c(Year, survival))
  names(trueSurv) <- c("Year", "Mean")
  trueSurv$parameter <- "Adult female survival"
  trueSurv$type <- "true"

  obsRec <- subset(oo$ageRatioOut, select = c(Year, Count, Class))
  obsRec <- tidyr::pivot_wider(obsRec, id_cols = c("Year"), names_from = "Class", values_from = "Count")
  obsRec$Mean <- obsRec$calf / obsRec$cow
  obsRec$parameter <- "Recruitment"
  obsRec$type <- "observed"

  trueRec <- subset(oo$exData, select = c(Year, recruitment))
  names(trueRec) <- c("Year", "Mean")
  trueRec$parameter <- "Recruitment"
  trueRec$type <- "true"

  obsLam <- subset(oo$exData, select = c(Year, lambda))
  names(obsLam) <- c("Year", "Mean")
  obsLam <- movingAveGrowthRate(obsLam, oo$cs$assessmentYrs)
  obsLam$parameter <- "Population growth rate"
  obsLam$type <- "true"

  obsSize <- subset(oo$exData, select = c(Year, N))
  names(obsSize) <- c("Year", "Mean")
  obsSize$Year <- obsSize$Year + 1 # pop size returned from Bayesian model is at the start of the year, not the end.
  obsSize$parameter <- "Female population size"
  obsSize$type <- "true"

  obsAll <- rbind(obsLam, obsSize, subset(obsRec, select = names(obsLam)), trueRec, subset(obsSurv, select = names(obsLam)), trueSurv)

  simBigO <- subset(simBig$summary, select = c(Anthro, Mean, lower, upper, parameter))
  names(simBigO) <- c("Anthro", "Mean", "Lower 95% CRI", "Upper 95% CRI", "parameter")

  # combine cs and simDisturbance and add to all output tables, nest params in a list
  dist_params <- merge(oo$simDisturbance, oo$cs)

  rr.summary <- merge(rr.summary, dist_params)
  simBigO <- merge(simBigO, dist_params)
  obsAll <- merge(obsAll, dist_params)
  rr.summary.all <- rr.summary
  sim.all <- simBigO
  obs.all <- obsAll

  if (getKSDists) {
    # get Kolmogorov smirnov distance between samples at each point

    variables <- unique(simBig$summary$parameter)
    anthroPts <- unique(subset(rr.summary, select = c(Year, Anthro)))
    # TO DO: make this step faster
    bmSamples <- tabAllRes(result, startYear, endYear, doSummary = F)
    bmSamples$type <- "local"

    simSamples <- merge(anthroPts, simBig$samples)
    simSamples$Anthro <- NULL
    simSamples$type <- "national"

    allSamples <- rbind(subset(bmSamples, is.element(Parameter, unique(simSamples$Parameter))), simSamples)

    ksDists <- allSamples %>%
      group_by(Year, Parameter) %>%
      group_modify(~ {
        getKSDist(.x$Value, .x$type)
      })
  } else {
    ksDists <- unique(subset(rr.summary, select = c(Year, Parameter)))
    ksDists$KSDistance <- NA
    ksDists$KSpvalue <- NA
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all, obs.all = obs.all, ksDists = ksDists))
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

makeInterceptPlots <- function(scResults, addBit = "", facetVars = c("P", "sQ"), loopVars = NULL,
                               whichPlots = c("Adult female survival", "Population growth rate", "Recruitment", "Female population size"),
                               survLow = 0.6, type = "png", useNational = T) {
  # facetVars=c("lre","sre");loopVars="srv";scResults=scResultsHigh

  if (!is.null(loopVars)) {
    loopSet <- unique(subset(scResults$rr.summary.all, select = loopVars))
    loopSet$dummy <- 1
  } else {
    loopSet <- data.frame(dummy = 1)
  }


  for (l in 1:nrow(loopSet)) {
    # l = 1
    crow <- loopSet[l, ]

    aa <- ""
    for (n in names(crow)) {
      if (n == "dummy") {
        next
      }
      aa <- paste0(aa, n, crow[[n]])
    }

    addBitO <- paste0(addBit, aa)

    if (useNational) {
      simRange <- merge(scResults$sim.all, crow)
    } else {
      simRange <- NULL
    }

    if (is.element("Adult female survival", whichPlots)) {
      if (type == "png") {
        png(here::here(paste0("figs/Surv", addBitO, ".png")),
          height = 6, width = 7.48, units = "in", res = 600
        )
      } else {
        pdf(paste0("figs/Surv", addBitO, ".pdf"), width = 10, height = 7)
      }
      print(plotRes(merge(scResults$rr.summary.all, crow), "Adult female survival",
        obs = merge(scResults$obs.all, crow),
        lowBound = survLow, simRange = simRange, facetVars = facetVars
      ))
      dev.off()
    }

    if (is.element("Population growth rate", whichPlots)) {
      if (type == "png") {
        png(here::here(paste0("figs/Lambda", addBitO, ".png")),
          height = 6, width = 7.48, units = "in", res = 600
        )
      } else {
        pdf(paste0("figs/Lambda", addBitO, ".pdf"), width = 10, height = 7)
      }
      print(plotRes(merge(scResults$rr.summary.all, crow), "Population growth rate",
        obs = merge(scResults$obs.all, crow),
        lowBound = 0, simRange = simRange, facetVars = facetVars
      ))
      dev.off()
    }

    if (is.element("Recruitment", whichPlots)) {
      if (type == "png") {
        png(here::here(paste0("figs/Rec", addBitO, ".png")),
          height = 6, width = 7.48, units = "in", res = 600
        )
      } else {
        pdf(paste0("figs/Rec", addBitO, ".pdf"), width = 10, height = 7)
      }
      print(plotRes(merge(scResults$rr.summary.all, crow), "Recruitment",
        obs = merge(scResults$obs.all, crow),
        lowBound = 0, simRange = simRange, facetVars = facetVars
      ))
      dev.off()
    }

    if (is.element("Female population size", whichPlots)) {
      if (type == "png") {
        png(here::here(paste0("figs/FPOP", addBitO, ".png")),
          height = 6, width = 7.48, units = "in", res = 600
        )
      } else {
        pdf(paste0("figs/FPOP", addBitO, ".pdf"), width = 10, height = 7)
      }
      print(plotRes(merge(scResults$rr.summary.all, crow), "Female population size",
        obs = merge(scResults$obs.all, crow),
        lowBound = 0, facetVars = facetVars
      ))
      dev.off()
    }
  }
}



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
      group_by(responseVariable) %>%
      summarise(nmod = n_distinct(ModelNumber)) %>%
      pull(nmod)

    if (any(nmod > 1)) {
      stop("The model coefficient file loaded contains more than one model per response variable",
        call. = FALSE
      )
    }
  }

  # add expected columns that should never change
  df <- mutate(df,
    modelVersion = "Johnson",
    ModelNumber = ifelse(responseVariable == "recruitment", "M4", "M1"),
    Type = "National"
  )

  # expected values
  diff_res <- setdiff(unique(df$responseVariable), unique(caribouMetrics::popGrowthTableJohnsonECCC$responseVariable))

  if (!setequal(unique(caribouMetrics::popGrowthTableJohnsonECCC$responseVariable), unique(df$responseVariable))) {
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

  testStdCI <- df %>% mutate(stdOrCI = !is.na(StdErr) | (!is.na(lowerCI) & !is.na(upperCI)))

  if (!all(testStdCI$stdOrCI)) {
    stop("The model coefficient file loaded is missing StdErr or lowerCI and upperCI for:\n",
      "femaleSurvival: ",
      testStdCI %>% filter(!stdOrCI, responseVariable == "femaleSurvival") %>%
        pull(Coefficient) %>% paste0(collapse = ", "), "\n",
      "recruitment: ",
      testStdCI %>% filter(!stdOrCI, responseVariable == "recruitment") %>%
        pull(Coefficient) %>% paste0(collapse = ", "),
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
         paste0(missing_cols, collapse = ", "))
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

# Plots -------------------------------------------------------------------

plotRes <- function(allRes, parameter, obs = NULL, lowBound = 0, highBound = 1, simRange = NULL, facetVars = NULL) {
  # allRes=scResults$ksDists; parameter="Recruitment";obs=scResults$obs.all;lowBound=0; highBound=1;simRange=scResults$sim.all;facetVars=c("P","sQ")

  if (is.null(facetVars)) {
    titleFontSize <- 16
    labFontSize <- 14
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

    df$Type <- "local"
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
  x2 <- x1 + ggplot2::theme_classic() + ggplot2::xlab("Year") + ggplot2::ylab(parameter) +
    ggplot2::geom_line(ggplot2::aes(x = Year, y = Mean), size = 1.75) +
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
    x2 <- x2 + ggplot2::geom_ribbon(ggplot2::aes(ymin = `Lower 95% CRI`, ymax = `Upper 95% CRI`),
      show.legend = FALSE, alpha = 0.25, colour = NA
    ) +
      ggplot2::scale_y_continuous(limits = c(
        ifelse(any(df$`Lower 95% CRI` < lowBound), NA, lowBound),
        ifelse(any(df$`Upper 95% CRI` > 1), NA, highBound)
      ))
  }

  if (!KS & !is.null(obs)) {
    obs$Type <- "local"
    obs$obsError <- F
    obs$obsError[obs$type == "observed"] <- T
    x2 <- x2 + ggplot2::geom_point(data = obs, ggplot2::aes(x = Year, y = Mean, shape = obsError), col = "black", show.legend = T) +
      ggplot2::scale_shape_manual(values = c(16, 2))
  }

  if (!is.null(facetVars)) {
    if (length(facetVars) == 2) {
      x2 <- x2 + ggplot2::facet_grid(as.formula(paste(facetVars[1], "~", facetVars[2])), labeller = "label_both")
    } else {
      x2 <- x2 + ggplot2::facet_wrap(as.formula(paste0("~", facetVars[1])), labeller = "label_both")
    }
  }
  if (!KS & (parameter == "Population growth rate")) {
    x2 <- x2 + ggplot2::geom_hline(yintercept = 1, color = "black")
  }

  x2
}



runScnSet <- function(scns, ePars, simBig, survAnalysisMethod = "KaplanMeier", getKSDists = T, printProgress = F) {
  # ePars=eParsIn;survAnalysisMethod="KaplanMeier";getKSDists=T;printProgress=F
  scns <- fillDefaults(scns)
  errorLog <- list()
  for (p in 1:nrow(scns)) {
    # p=1
    cs <- scns[p, ]
    if (printProgress) {
      print(paste0(c(p, scns[p, ]), collapse = " "))
    }

    if (is.element("cw", names(cs))) {
      ePars$cowCounts$Count <- cs$cw
    }
    oo <- simulateObservations(cs, cowCounts = ePars$cowCounts, freqStartsByYear = ePars$freqStartsByYear, collarNumYears = ePars$collarNumYears, collarOffTime = ePars$collarOffTime, collarOnTime = ePars$collarOnTime)
    betaPriors <- getPriors(cs)
    minYr <- min(oo$exData$Year)
    maxYr <- max(oo$simDisturbance$Year)
    out <- try(runRMModel(
      survData = oo$simSurvObs, ageRatio.herd = oo$ageRatioOut, disturbance = oo$simDisturbance,
      betaPriors = betaPriors, startYear = minYr, endYear = maxYr, N0 = cs$N0, survAnalysisMethod = survAnalysisMethod, adjustR = cs$adjustR, assessmentYrs = cs$assessmentYrs
    ))
    if (inherits(out, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all, obs.all = obs.all, ksDists = ksDists, errorLog = errorLog), "temp.Rds")
      next
    }

    if (inherits(out$result, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out$result)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all, obs.all = obs.all, ksDists = ksDists, errorLog = errorLog), "temp.Rds")
      next
    }

    outTabs <- getOutputTables(result = out$result, startYear = minYr, endYear = maxYr, survInput = out$survInput, oo = oo, simBig = simBig, getKSDists = getKSDists)

    if (p == 1) {
      rr.summary.all <- outTabs$rr.summary.all
      sim.all <- outTabs$sim.all
      obs.all <- outTabs$obs.all
      ksDists <- merge(outTabs$ksDists, cs)
    } else {
      rr.summary.all <- rbind(rr.summary.all, outTabs$rr.summary.all)
      sim.all <- rbind(sim.all, outTabs$sim.all)
      obs.all <- rbind(obs.all, outTabs$obs.all)
      ksDists <- rbind(ksDists, merge(outTabs$ksDists, cs))
    }
  }
  if (length(errorLog) > 0) {
    print(errorLog)
  }
  return(list(rr.summary.all = rr.summary.all, sim.all = sim.all, obs.all = obs.all, ksDists = ksDists, errorLog = errorLog))
}
