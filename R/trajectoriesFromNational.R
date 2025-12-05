# this creates an environment where we can store objects that will be available
# to multiple functions/multiple function calls. Does not persist across
# sessions but it only take ~ 20s so once per session is probably ok.
# See explanation here: https://r-pkgs.org/data.html#sec-data-state
cacheEnv <- new.env(parent = emptyenv())

# Not supposed to save files to user computer on CRAN so for users the cache is
# only preserved within a session but for dev I have added this "persistent
# cache" use savePersistentCache function to update/create it after having run
# trajectoriesFromNational
if(file.exists("results/simsInitial.rds")){
  simsInitial <- readRDS( "results/simsInitial.rds")
  bayesianResults <- readRDS( "results/bayesianResults.rds")
  
  assign("simsInitial", simsInitial, envir = cacheEnv)
  assign("bayesianResults", bayesianResults, envir = cacheEnv)
}


#' Get a set of simulation results from the national demographic model
#' 
#' Simulate demograhic rates based on the National demographic - disturbance model
#' If a disturbance scenario containing Years is supplied trajectories will show
#' growth of a population over time based on the National demographic - disturbance model
#'
#' @param disturbance data frame with Anthro, fire_excl_anthro and Year numeric
#'   columns. Anthro and fire_excl_anthro are vectors of numbers between 0 and 100
#'   representing the percentage of the landscape covered by anthropogenic
#'   disturbance buffered by 500 m, and the percentage covered by fire that does
#'   not overlap anthropogenic disturbance.
#' @inheritParams getNationalCoefficients
#' @inheritParams caribouPopGrowth
#' @param N0 initial population size
#' @param cPars optional. Parameters for calculating composition survey bias term.
#' @param doSummary logical. Default TRUE. If FALSE returns unprocessed outcomes from caribouPopGrowth. 
#'  If TRUE returns summaries and (if returnSamples = T) sample trajectories from prepareTrajectories.
#' @param returnSamples logical. If FALSE returns only summaries. If TRUE
#'   returns example trajectories as well. By default summaries are not returned
#'   unless the disturbance data provided contains a column named "Year".
#' @param numSteps numeric. Number of steps to run `caribouPopGrowth()` at each 
#'   disturbance level. 
#'
#' @return Output from caribouPopGrowth function.
#' 
#' @family demography
#' @export
#'
#' @examples
#' trajectoriesFromNational()
trajectoriesFromNational <- function(replicates = 1000, N0 = 1000,
                            useQuantiles  = NULL,
                            populationGrowthTable  = NULL,
                            cPars = subset(getScenarioDefaults(),select=-iAnthro),
                            interannualVar = eval(formals(caribouPopGrowth)$interannualVar),
                            disturbance = NULL,
                            skipSave = FALSE,
                            forceUpdate = FALSE,
                            doSummary = TRUE,
                            returnSamples = "default",
                            numSteps = 1) {
  # replicates=1000;N0=1000;Anthro=seq(0,100,by=1);fire_excl_anthro=0;
  # useQuantiles =NULL

  # from trajectoriesFromAny #===========================================
  doSave <- FALSE

  hasAnthro <- is.element("iAnthro",names(cPars))
  cPars <- getScenarioDefaults(cPars)

  if(!skipSave){
    check <- as.list(match.call())

    saveName <- "simsInitial"

    if (length(check) == 1) {
      if (exists(saveName, envir=cacheEnv)) {
        message("Using saved object")
        return(get(saveName, envir=cacheEnv))
      } else {
        doSave <- TRUE
      }
    }
    check$forceUpdate <- NULL

    if (forceUpdate & (length(check) == 1)) {
      doSave <- TRUE
    }
  }
  hasYear <- T
  if(is.null(disturbance)){
    if(hasAnthro){
      distPars = unique(subset(cPars,select=c(iAnthro,iFire,preYears,obsYears,projYears,obsAnthroSlope,projAnthroSlope,preYears,startYear)))
      first<-T
      for(r in 1:nrow(distPars)){
        #r=60
        cr <- distPars[r,]
        covariates <- simCovariates(cr$iAnthro, cr$iFire,
                                    cr$preYears+cr$obsYears + cr$projYears,
                                    cr$obsAnthroSlope, cr$projAnthroSlope,
                                    cr$obsYears + cr$preYears + 1)
        covariates$Year <- cr$startYear + covariates$time - 1
        covariates$fire_excl_anthro=round(covariates$fire_excl_anthro)
        if(first){
          covTableObs <- covariates
          first=F
        }else{
          covTableObs <- unique(rbind(covTableObs, covariates))
        }
      }
    }else{
      if(!is.element("Anthro",names(cPars))){
        covTableObs <- expand.grid(Anthro=seq(0,100,by=1),fire_excl_anthro=0,Year=NA)
        covTableObs$Year <- covTableObs$Anthro
        hasYear <- F
      }else{
        covTableObs <- unique(subset(cPars, select = c("Year","Anthro","fire_excl_anthro")))
      }
    }
  }else {
    if(!is.element("Year",names(disturbance))){
      hasYear <- F
      disturbance$Year <- disturbance$Anthro
    }
    covTableObs <- disturbance %>% select(Year, Anthro, fire_excl_anthro)
  }
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  if(length(N0)>1){
    stop("Specify a single initial population size for trajectories from national model.")
  }

  # original trajectoriesFromNational #============================================

  covTableObs$Total_dist <- covTableObs$Anthro + covTableObs$fire_excl_anthro

  if (is.null(populationGrowthTable )) {
    populationGrowthTable  <- caribouMetrics::popGrowthTableJohnsonECCC
  }
  if (is.null(useQuantiles )) {
    popGrowthPars <- getNationalCoefficients(
      replicates,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- estimateNationalRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = TRUE, useQuantiles = FALSE)
  } else {
    popGrowthPars <- getNationalCoefficients(
      replicates, useQuantiles = useQuantiles,
      populationGrowthTable = populationGrowthTable
    )
    rateSamplesAll <- estimateNationalRates(covTable = covTableObs,
                                       popGrowthPars = popGrowthPars,
                                       returnSample = T)
  }

  bc = unique(subset(rateSamplesAll,select=replicate));nr=nrow(bc)
  bc$c = compositionBiasCorrection(q=runif(nr,cPars$qMin,cPars$qMax),w=cPars$cowMult,u=runif(nr,cPars$uMin,cPars$uMax),
                                   z=runif(nr,cPars$zMin,cPars$zMax))
  rateSamplesAll$c = NULL; rateSamplesAll= merge(rateSamplesAll, bc)

  #print(paste("trajectoriesFromNational",mean(bc$c)))
  pars <- merge(data.frame(N0 = N0, PopulationName = "National"), rateSamplesAll)

  if(hasYear){
    R_dat <- pars %>% select(Year, replicate, R_bar) %>% 
      pivot_wider(names_from = Year, values_from = R_bar) %>% 
      tibble::column_to_rownames("replicate")
    S_dat <- pars %>% select(Year, replicate, S_bar) %>% 
      pivot_wider(names_from = Year, values_from = S_bar) %>% 
      tibble::column_to_rownames("replicate")
    
    out <- simPopsOverTime(
      N0, numSteps = n_distinct(pars$Year), R_samp = R_dat,
      S_samp = S_dat, dynamicRates = TRUE, stepLength = numSteps,
      c = unique(pars$c), interannualVar = interannualVar, 
      progress = FALSE, K = FALSE
    ) %>% mutate(replicate = factor(id, labels = rownames(R_dat)), 
                   Year = factor(time, labels = colnames(R_dat)) %>%
                     as.character() %>% as.numeric(), .keep = "unused") 
    
    pars <- full_join(pars %>% select(-N0), out, by = c("replicate", "Year"))
  } else {
    pars <- cbind(subset(pars,select=-N0), 
                  caribouPopGrowth(pars$N0, R_bar = pars$R_bar,
                                   S_bar = pars$S_bar, numSteps = numSteps,
                                   K = FALSE, c = pars$c,
                                   interannualVar=interannualVar, progress = FALSE))
  }
  
  
  names(pars)[names(pars)=="replicate"]= "id"
  
  if(returnSamples == "default"){
      returnSamples <- hasYear
  } else if (!hasYear & isTRUE(returnSamples)){
    warning("returnSamples is set to FALSE when Year is not included in the disturbance scenario")
  }

  if(doSummary){
    if(!hasYear){
      simBig <- prepareTrajectories(pars, returnSamples = FALSE)
      simBig$summary$Year = NULL
      simBig$summary <- subset(simBig$summary,MetricTypeID!="N")
    }else{
      simBig <- prepareTrajectories(pars, returnSamples = returnSamples)
    }
  }else {
    simBig <- pars
  }

  # Note this must be the last thing before return or the first and cached
  # results won't match
  if (doSave) {
    message("Updating cached initial simulations.")
    assign(saveName, simBig, envir = cacheEnv)
  }
  return(simBig)
}
