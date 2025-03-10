# this creates an environment where we can store objects that will be available
# to multiple functions/multiple function calls. Does not persist across
# sessions but it only take ~ 20s so once per session is probably ok.
# See explanation here: https://r-pkgs.org/data.html#sec-data-state
cacheEnv <- new.env(parent = emptyenv())

# Not supposed to save files to user computer on CRAN so for users the cache is
# only preserved within a session but for dev I have added this "persistent
# cache" use savePersistentCache function to update/create it after having run
# getSimsInitial
if(file.exists("inst/extdata/simsInitiald.rds")){
  simsInitial <- readRDS( "inst/extdata/simsInitial.rds")
  bbouResults <- readRDS( "inst/extdata/bbouResults.rds")

  assign("simsInitial", simsInitial, envir = cacheEnv)
  assign("bbouResults", bbouResults, envir = cacheEnv)
}

#' Get a set of simulation results from fitted demographic models
#'
#' @param bbouResults Fitted bboutools model and summary table created by bbouMakeSummaryTable(), 
#'   or a path to those results.
#' @param N0 initial population size. If NULL, will use information in bbouResults.
#' @param cPars ptional. Parameters for calculating composition survey bias term.
#' @param forceUpdate logical. If the default inputs are used the result is
#'   cached. Set `forceUpdate` to TRUE to ensure the simulations are re-run.
#' @param skipSave logical. If F default consider saved results. Set T to ignore the saved file.
#' @param returnSamples logical. If T default return example trajectories. Set F to return only summaries.
#' @inheritParams caribouPopGrowth
#'
#' @return a list with two elements:
#'  * summary: a tibble with a summary of parameter values for each scenario.
#'    Column names are year, PopulationName, Mean, lower, upper, Parameter.
#'  * samples: a tibble with parameter values for each scenario and replicate
#'    4 rows per replicate \* scenario. Column names are year, PopulationName,  Parameter and Value
#'    
#' 
#' @family demography
#' @export
#'
#' @examples
#' getSimsInitial()
getSimsInitial <- function(bbouResults, N0=NULL,
                            cPars=getScenarioDefaults(), forceUpdate = F,skipSave=F,returnSamples=T,...) {
  doSave <- FALSE

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
  
  #bbouResults = bbouResultFile
  if(is.character(bbouResults) && (length(bbouResults) == 1) && file.exists(bbouResults)){
    bbouResults <- readRDS(bbouResults)
  }
  
  if(is.null(N0)){
    N0 <- subset(bbouResults$parTab,select=c(pop_name,N0))
  }
  if(length(N0)==1){
    N0 = data.frame(N0=N0)
    N0 = merge(N0,subset(bbouResults$parTab,select=c(pop_name)))
  }
  
  surv_pred <- bb_predict_survival (bbouResults$surv_fit,year=T,month=F,conf_level=F)
  rec_pred <- bb_predict_calf_cow_ratio(bbouResults$recruit_fit,year=T,conf_level=F)
  
  nr <- dim(surv_pred$samples)[1]*dim(surv_pred$samples)[2]
  
  popInfo <- merge(data.frame(id=seq(1:nr)),N0)
  popInfo$c <- compositionBiasCorrection(q=runif(nrow(popInfo),cPars$qMin,cPars$qMax),w=cPars$cowMult,u=runif(nr,cPars$uMin,cPars$uMax),
                                       z=runif(nr,cPars$zMin,cPars$zMax))
  pars <- caribouPopSimMCMC(popInfo,rec_pred,surv_pred,progress=F,...)
  
  pars$Anthro=NA;pars$fire_excl_anthro=NA
  
  popInfo$PopulationName <- popInfo$pop_name
  
  pars <- merge(pars,subset(popInfo,select=c(-N0,-pop_name)))
  
  #get the lambda percentile for each id - to allow users to select extreme examples
  simSum <- pars  %>%
    group_by(id) %>%
    summarize(MeanLam = mean(lambdaTrue,na.rm=T))
  simSum <-simSum[order(simSum$MeanLam),]
  simSum$lamPercentile <- round(100*seq(1:nrow(simSum))/nrow(simSum))
  simSum$MeanLam=NULL
  pars <- merge(pars,simSum)

  pars <- convertTrajectories(pars)
  simBig <- summarizeCaribouPopSim(pars,returnSamples=returnSamples)

  if (doSave) {
    message("Updating cached initial simulations.")
    assign(saveName, simBig, envir = cacheEnv)
  }
  simBig$surv_data = bbouResults$surv_fit$data
  simBig$recruit_data = bbouResults$recruit_fit$data
  return(simBig)
}
