# internal functions related to caribou demographics

# Helpers for table format conversion --------------------------------------------------
#' Format trajectory tables
#'
#' @param pars 
#'
#' @returns convertTrajectories: formatted tables
#' @export
#'
#' @rdname simulateTrajectoriesFromPosterior
#' 

convertTrajectories<-function(pars){
  #converts output from simPopsOverTime to alternate form
  #pars = trajectories
  if(!is.element("lamPercentile",names(pars))){
    pars$lamPercentile=NA
  }
  if(!is.element("c",names(pars))){pars$c=NA}
  
  nameChange <- data.frame(inName=c("id","lamPercentile", "Year","PopulationName","Anthro", "fire_excl_anthro","c",
                                    "S_t", "R_t", "X_t", "N",
                                    "lambda","S_bar","R_bar","X_bar","N_bar","lambdaE_bar"),
                           outName=c("Replicate","LambdaPercentile","Year", "PopulationName","Anthro", "fire_excl_anthro","c", 
                                     "survival","recruitment","X", "N", "lambda","Sbar","Rbar","Xbar","Nbar","lambda_bar"))
  
  if(!is.element("lambdaE_bar",names(pars))){
    pars$lambdaE_bar = pars$lambdaE
  }
  nameChange <-subset(nameChange,is.element(inName,names(pars)))
  fds <- subset(pars, select = nameChange$inName)
  names(fds) <- nameChange$outName
  
  if(is.element("Anthro", colnames(fds))){
    fds$AnthroID = round(fds$Anthro);fds$fire_excl_anthroID=round(fds$fire_excl_anthro)
  }
  
  fds$Timestep = as.numeric(fds$Year)
  fds$Year=as.numeric(as.character(fds$Year))
  fds <- tidyr::pivot_longer(fds, !any_of(c("Replicate","LambdaPercentile","Year","Timestep","PopulationName","AnthroID","fire_excl_anthroID")), names_to = "MetricTypeID",
                             values_to = "Amount")
  fds$MetricTypeID <- as.character(fds$MetricTypeID)
  fds$Replicate <- paste0("x", fds$Replicate)
  
  if(!is.element("lamPercentile",names(pars))){
    fds$LambdaPercentile=NULL
  }
  return(fds)
}

#' Get 95% prediction intervals from trajectories
#'
#' @param pars 
#' @param returnSamples 
#'
#' @returns summarizeTrajectories:
#' @export
#' @family demography
#'
#' @rdname simulateTrajectoriesFromPosterior
summarizeTrajectories <- function(pars,returnSamples=T){

  if(is.element("AnthroID",names(pars))){  
    simSum <- pars  %>%
      group_by(Year,PopulationName,MetricTypeID,AnthroID,fire_excl_anthroID) %>%
      summarize(Mean = mean(Amount,na.rm=T), lower = quantile(Amount, 0.025,na.rm=T),
                upper = quantile(Amount, 0.975,na.rm=T),probViable=mean(Amount > 0.99,na.rm=T))
  }else{
    simSum <- pars  %>%
      group_by(Year,PopulationName,MetricTypeID) %>%
      summarize(Mean = mean(Amount,na.rm=T), lower = quantile(Amount, 0.025,na.rm=T),
                upper = quantile(Amount, 0.975,na.rm=T),probViable=mean(Amount > 0.99,na.rm=T))
  }  
  names = data.frame(MetricTypeID = c("survival","recruitment","X", "lambda","N","c",
                                      "Sbar","Rbar","Xbar","lambda_bar"),
                     Parameter = c("Adult female survival","Recruitment","Adjusted recruitment",
                                   "Population growth rate","Female population size","c",
                                   "Expected survival","Expected recruitment","Expected adjusted recruitment","Expected growth rate"
                                   ))
  simSum=merge(simSum,names)
  if (returnSamples){
    simBig <- list(summary = simSum, samples = pars)
  } else {
    simBig <- list(summary = simSum)
  }

  return(simBig)
}

# Helpers for simulateObservations -----------------------------------------

simTrajectory <- function(numYears, covariates, survivalModelNumber = "M1",
                          recruitmentModelNumber = "M4",
                          popGrowthTable = caribouMetrics::popGrowthTableJohnsonECCC,
                          recSlopeMultiplier = 1, sefSlopeMultiplier = 1,
                          rQuantile = NULL, sQuantile = NULL,
                          stepLength = 1, N0 = 1000,cowMult=1,
                          qMin=0,qMax=0,uMin=0,uMax=0,zMin=0,zMax=0,interannualVar = eval(formals(caribouPopGrowth)$interannualVar)) {
  # survivalModelNumber = "M1";recruitmentModelNumber = "M4";
  # recSlopeMultiplier=1;sefSlopeMultiplier=1;recQuantile=0.5;sefQuantile=0.5
  # stepLength=1;N0=1000
  
  if(is.null(rQuantile)||is.na(rQuantile)){rQuantile<-runif(1)}
  if(is.null(sQuantile)||is.na(sQuantile)){sQuantile<-runif(1)} 
  
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
  
  popGrowthParsSmall <- getNationalCoefficients(
    2,
    modelVersion = "Johnson",
    survivalModelNumber = survivalModelNumber,
    recruitmentModelNumber = recruitmentModelNumber,
    populationGrowthTable = growthTab,
    useQuantiles = c(rQuantile, rQuantile)
  )
  # set quantiles for example population
  popGrowthParsSmall$coefSamples_Survival$quantiles <- sQuantile
  
  # Only use precision if included in the table for this model number for both
  # rec and surv
  usePrec <- "Precision" %in% names(popGrowthParsSmall$coefSamples_Survival$coefValues) &
    "Precision" %in% names(popGrowthParsSmall$coefSamples_Recruitment$coefValues)
  # at each time,  sample demographic rates and project, save results
  
  
  pars <- data.frame(N0 = N0)
  
  # sample rates with covariates from each timestep
  rateSamples <- estimateNationalRates(
    covTable = covariates,
    popGrowthPars = popGrowthParsSmall,
    ignorePrecision = !usePrec,
    returnSample = TRUE
  )[1:nrow(covariates),] # only using the first replicate but doesn't work with just 1
  
  #set bias correction term for each example population - constant over time.
  bc = unique(subset(rateSamples,select=replicate))
  nr=nrow(bc)
  c = compositionBiasCorrection(q=runif(nr,qMin,qMax),w=cowMult,
                                   u=runif(nr,uMin,uMax),z=runif(nr,zMin,zMax))

  popMetrics <- simPopsOverTime(N0, numSteps = numYears, R_samp = rateSamples$R_bar,
                              S_samp = rateSamples$S_bar, 
                              interannualVar = interannualVar,
                              onePop = TRUE,
                              stepLength = stepLength,
                              K = FALSE,
                              l_R = 1e-06,
                              c = c, 
                              progress = FALSE)
  popMetrics <- merge(popMetrics, rateSamples)
  
  popMetrics <- convertTrajectories(popMetrics)
  return(popMetrics)
}

simSurvivalData <- function(freqStartsByYear, exData, collarNumYears, collarOffTime,
                            collarOnTime, caribouYearStart,topUp = FALSE,forceMonths=FALSE) {
  #Note: If collarOffTime and collarOnTime both equal caribouYearStart simulation will be faster because we can ignore variation in number of collars at the start of each month.
  # topUp=T;caribouYearStart=4
  # for simplicity, ignore variation in survival probability among months
  
  options(dplyr.summarise.inform = FALSE)
  
  if(!forceMonths&&(collarOnTime==caribouYearStart)&&(collarOffTime==caribouYearStart)){
    nMonths = 1
  }else{
    nMonths = 12
  }
  
  survivalSeries <- subset(exData, select = c("survival", "N","Year","PopulationName","Replicate"))
  if(nMonths>1){
    survivalSeries <- merge(survivalSeries,data.frame(Month=seq(1:nMonths)))
  }else{
    survivalSeries$Month=caribouYearStart
  }
  survivalSeries$Year[survivalSeries$Month<caribouYearStart]= survivalSeries$Year[survivalSeries$Month<caribouYearStart]+1
  
  initYear <- min(survivalSeries$Year)
  freqStartsByYear <- subset(freqStartsByYear,
                             (freqStartsByYear$Year >= initYear))
  freqStartsByYear$Month = collarOnTime
  
  survivalSeries = merge(survivalSeries,freqStartsByYear,all.x=T)
  survivalSeries$numStarts[is.na(survivalSeries$numStarts)]=0
  
  startYrs = sort(unique(freqStartsByYear$Year))
  
  firstStep=T
  
  for (sy in startYrs) {
    #sy = startYrs[2]
    y = sy
    cMonth = caribouYearStart
    if(nMonths==1){prevMonth=cMonth}else{prevMonth = cMonth-1}
    prevYear = y
    for (yId in seq(sy,min(sy+collarNumYears,max(survivalSeries$Year)))){
      if ((y == sy+collarNumYears)&(cMonth==collarOffTime)){break}
      for(mId in 1:nMonths){
        #print(paste(sy, y,cMonth))
        cInfo = subset(survivalSeries,(Month==cMonth)&(Year==y))
        if(nrow(cInfo)==0){break}
        cInfo$startYr = sy
        
        if(!firstStep){
          cInfo$Prevs = NULL;cInfo$PrevsAll=NULL
          prevs = subset(cAll,(Month==prevMonth)&(Year==prevYear)&(startYr==sy)) %>%
            group_by(PopulationName, Replicate) %>%
            summarise(Prevs = sum(StartTotal,na.rm=T)-sum(MortalitiesCertain,na.rm=T))
          if((cMonth==prevMonth)&(y==prevYear)){
            prevsAll = subset(cAll,(Month==prevMonth)&(Year==prevYear)) %>%
              group_by(PopulationName, Replicate) %>%
              summarise(PrevsAll = sum(StartTotal,na.rm=T))
          }else{
            prevsAll = subset(cAll,(Month==prevMonth)&(Year==prevYear)) %>%
              group_by(PopulationName, Replicate) %>%
              summarise(PrevsAll = sum(StartTotal,na.rm=T)-sum(MortalitiesCertain,na.rm=T))
          }
          
          cInfo = merge(cInfo,prevs,all.x=T)
          cInfo = merge(cInfo,prevsAll,all.x=T)
          cInfo$Prevs[is.na(cInfo$Prevs)]=0;cInfo$PrevsAll[is.na(cInfo$PrevsAll)]=0
        }else{
          cInfo$Prevs=0;cInfo$PrevsAll=0
        }
        
        if ((y == sy+collarNumYears)&(cMonth==collarOffTime)){break}
        
        if(y==sy){
          if (topUp) {
            cInfo$StartTotal=pmax((cInfo$numStarts-cInfo$PrevsAll),cInfo$Prevs)
          } else {
            cInfo$StartTotal=cInfo$Prevs+cInfo$numStarts
          }
        }else{
          cInfo$StartTotal=cInfo$Prevs
          #if(sum(cInfo$StartTotal)==0){break}
        }
        
        if(any(cInfo$StartTotal[!is.na(cInfo$N)]>cInfo$N[!is.na(cInfo$N)])){
          warning("Target number of collars exceeds population size. Adjusting number of collars for consistency.")
          cInfo$StartTotal[cInfo$StartTotal>cInfo$N] = cInfo$N[cInfo$StartTotal>cInfo$N]
        }
        
        cInfo$MortalitiesCertain = rbinom(nrow(cInfo),cInfo$StartTotal,prob=(1-cInfo$survival^(1/nMonths)))
        
        if(!firstStep){
          cAll = rbind(cAll,cInfo)
        }else{
          cAll = cInfo
        }
        
        prevMonth = cMonth;prevYear = y;firstStep=F
        
        if(nMonths==1){
          y=y+1
        }else{
          if(cMonth==12){cMonth=1;y=y+1}else{cMonth=cMonth+1}
        }
        #print(cInfo)
      }
    }
  }

  
  simSurvs <- cAll %>%
    group_by(PopulationName, Replicate, Year, Month) %>%
    summarise(StartTotal=sum(StartTotal,na.rm=T),MortalitiesCertain=sum(MortalitiesCertain,na.rm=T), survival=mean(survival,na.rm=T))
  simSurvs$Malfunctions = 0
  simSurvs$MortalitiesUncertain = 0
  
  #plot(plotSurvivalSeries(subset(simSurvs,Replicate==simSurvObs$Replicate[1])))
  if(topUp){
    if(any(simSurvs$StartTotal>max(freqStartsByYear$numStarts))){stop("Error in simSurvivalData: too many collars")}
  }
  if(nrow(simSurvs)==0){
    stop("TO DO: deal with no sampling case")
  }
  
  simSurvs$MortalitiesCertain[simSurvs$StartTotal==0]=NA
  
  return(simSurvs)
}

simCalfCowRatios <- function(cowCounts, exData) {
  # assume info from only one herd
  
  recruitmentSeries <- subset(exData, select = c("recruitment", "N","Year","PopulationName","Replicate"))

  simRecruitObs <- merge(cowCounts,recruitmentSeries,all.x=T)
  
  if(any(simRecruitObs$Cows>simRecruitObs$N,na.rm=T)){
    warning("The expected number of cows in composition survey exceeds population size. Adjusting cows in survey for consistency.")
    simRecruitObs$Cows = pmin(simRecruitObs$Cows,simRecruitObs$N)
  }
  apparentCows <- simRecruitObs$Cows
  if(is.element("UnknownAdults",names(simRecruitObs))){
    apparentCows = apparentCows + simRecruitObs$UnknownAdults*0.65
  }
  if(is.element("Yearlings",names(simRecruitObs))){
    apparentCows = apparentCows + simRecruitObs$Yearlings*0.5
  }
  
  #apparent number of calves (M+F) from apparent number of cows using apparent recruitment rate
  # removing NAs and then putting them back to avoid warning in rbinom
  na_cows <- which(is.na(apparentCows))
  apparentCows[na_cows] <- 0

  simRecruitObs$Calves <-  rbinom(
    n = nrow(simRecruitObs), size = round(apparentCows),
    prob = simRecruitObs$recruitment
  )
  
  simRecruitObs$Calves[na_cows] <- NA_integer_
       
  simRecruitObs$recruitment = NULL;simRecruitObs$N=NULL;simRecruitObs$StartTotal=NULL
  
  simRecruitObs$Calves[simRecruitObs$Cows==0]=NA
  return(simRecruitObs)
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

# saves the cached version of the initial sims to a local file so it doesn't
# need to be re-run every time you restart
savePersistentCache <- function(env = cacheEnv){
  obj_nms <- ls(envir = env)
  lapply(obj_nms, function(x){
    obj <- get(x, envir=env)
    saveRDS(obj, paste0("results/", x, ".rds"))
  })
  return(invisible())
}
