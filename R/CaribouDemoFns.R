




# TODO: These are not used any where else. Remove?
simSurvivalData <- function(freqStartsByYear, exData, collarNumYears, collarOffTime,
                            collarOnTime, topUp = F) {
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

  simSurvObs <- subset(simSurvObs, is.element(Year, exData$Year))
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

getParamsFromEacker <- function(path) {
  ###############
  # Use Eacker example data for sample sizes in each year.
  survData <- paste0(path, "/tte_caribouFEMALES.csv")
  ageRatio.herd <- paste0(path, "/ageRatio.herd.csv")
  ageRatio.herd2 <- read.csv(ageRatio.herd, header = T)
  tte_caribou2 <- read.csv(survData, header = T)

  # need table of observed number of cows each year as input for simulating
  # calf:cow ratios. Use Eaker as example.
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
    summarize(startYear = min(Year), endYear = max(Year),
              numYears = max(Year) - min(Year), died = sum(event))
  # for simplicity, pick a single number of years that collars remain on
  collarNumYears <- median(subset(sm, !died)$numYears)

  addInfo <- unique(subset(tte_caribou2,
                           select = c(HerdDescription, HerdCode, Range_ID,
                                      RangeDescription, RangeCode)))
  # TO DO: remove requirement for additional herd ID and range ID columns in UI code
  freqStartsByYear <- merge(freqStartsByYear, addInfo)
  write.csv(freqStartsByYear, "tabs/freqStartsByYear.csv")

  offSet <- subset(sm, !died, select = c(id, endYear))
  names(offSet) <- c("id", "Year")
  offSet <- merge(offSet, tte_caribou2, all.x = T)
  # for simplicity, pick a single month that collars fall off
  collarOffTime <- median(offSet$exit)

  onSet <- subset(sm, select = c(id, startYear))
  names(onSet) <- c("id", "Year")
  onSet <- merge(onSet, tte_caribou2, all.x = T)
  # for simplicity, pick a single month that collars are put on
  collarOnTime <- median(onSet$enter)
  # TO DO: allow users to set freqStartsByYear, collarNumYears, collarOffTime,
  # and collarOnTime as parameters OR provide a file formatted as in Eacker from
  # which these parameters can be derived.

  # freqStartsByYear$numStarts=30
  return(list(cowCounts = cowCounts, freqStartsByYear = freqStartsByYear,
              collarOnTime = collarOnTime, collarOffTime = collarOffTime,
              collarNumYears = collarNumYears))
}

# TODO: used in paper only. move there?
makeInterceptPlots <- function(scResults, addBit = "", facetVars = c("P", "sQ"),
                               loopVars = NULL,
                               whichPlots = c("Adult female survival",
                                              "Population growth rate",
                                              "Recruitment",
                                              "Female population size"),
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









