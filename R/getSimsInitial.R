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

#' Get a set of simulation results from fitted demographic models (if bbouResults argument is provided) or from national demographic disturbance model.
#'
#' @param bbouResults Optional. Fitted bboutools model and summary table created by bbouMakeSummaryTable(), 
#'   or a path to those results. If not specified trajectories will be from national demographic disturbance model.
#'   To specify a disturbance table for the national model set bbouResults=disturbanceTable. 
#'   disturbanceTable should be a data.frame with columns Anthro, fire_excl_anthro, and Year.
#' @param N0 initial population size. If NULL, will use information in bbouResults.
#' @param replicates Number of replicate trajectories. Default "all"   
#' @param cPars optional. Parameters for calculating composition survey bias term and/or disturbance scenario. Note that cPars can specify multiple disturbance scenarios, but only one unique composition bias scenario.
#' @param forceUpdate logical. If the default inputs are used the result is
#'   cached. Set `forceUpdate` to TRUE to ensure the simulations are re-run.
#' @param skipSave logical. If F default consider saved results. Set T to ignore the saved file.
#' @param returnSamples logical. If T default return example trajectories. Set F to return only summaries.
#' @inheritParams caribouPopGrowth
#'
#' @return a list with two elements:
#'  * summary: a tibble with a summary of parameter values for each scenario.
#'    Column names are Year or Anthro, MetricTypeID, PopulationName, Mean, lower, upper, probViable and Parameter
#'  * samples: a tibble with parameter values for each scenario and replicate
#'    1 row per MetricTypeID per replicate \* scenario. Column names are Year or Anthro, Replicate, PopulationName, MetricTypeID and Amount
#'    
#' 
#' @family demography
#' @export
#'
#' @examples
#' getSimsInitial()
getSimsInitial <- function(bbouResults=NULL, N0=NULL, replicates = "all",
                            cPars=subset(getScenarioDefaults(),select=-iAnthro), forceUpdate = F,skipSave=F,returnSamples=T,...) {
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
  rmSamples<-F

  if(is.null(bbouResults)){
    if(hasAnthro){
      distPars = unique(subset(cPars,select=c(iAnthro,iFire,preYears,obsYears,projYears,obsAnthroSlope,projAnthroSlope,preYears)))
      
      first<-T
      for(r in 1:nrow(distPars)){
        #r=1
        cr <- cPars[r,]
        covariates <- simCovariates(cr$iAnthro, cr$iFire, 
                                    cr$preYears+cr$obsYears + cr$projYears, 
                                    cr$obsAnthroSlope, cr$projAnthroSlope, 
                                    cr$obsYears + cr$preYears + 1)
        covariates$Year <- cr$startYear + covariates$time - 1
        covariates$fire_excl_anthro=round(covariates$fire_excl_anthro)
        if(first){
          simDisturbance <- covariates
          first=F
        }else{
          simDisturbance <- unique(rbind(simDisturbance,covariates))
        }
      }
      bbouResults <- list(parTab=subset(simDisturbance,select=c(Anthro,fire_excl_anthro,Year)))
    }else{
      warning("To create sample trajectories from the national model a disturbance scenario must be specified in bbouResults$parTab or cPars arguments.")
      rmSamples <- T
      bbouResults <- list(parTab = eval(formals(getSimsNational)$covTableObs))
      bbouResults$parTab$Year=bbouResults$parTab$Anthro
    }
  }
  if(is.element("Anthro",names(bbouResults))){
    bbouResults=list(parTab=bbouResults)
  }
    
  #bbouResults = bbouResultFile
  if(is.character(bbouResults) && (length(bbouResults) == 1) ){
    if(file.exists(bbouResults)){
      bbouResults <- readRDS(bbouResults)
    }else{
      stop(paste("bbouResults file not found,",bbouResults))
    }
  }

  if(!is.element("pop_name",names(bbouResults$parTab))){bbouResults$parTab$pop_name=NA}
  
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  
  if(is.null(N0)){
    if(is.element("N0",names(bbouResults$parTab))){
      N0 <- subset(bbouResults$parTab,select=c(pop_name,N0))
    }else{
      N0 <- eval(formals(getSimsNational)$N0)
    }
  }

  if(length(N0)==1){
    if(class(bbouResults$parTab) == "list"){
      bbouResults$parTab <- as.data.frame(bbouResults$parTab)
    }
    N0 = merge(data.frame(N0=N0),subset(bbouResults$parTab,select=c(pop_name)))
  }
  N0$PopulationName = N0$pop_name

  if(!is.null(bbouResults$surv_fit)){
    if(is.element("bboufit",class(bbouResults$surv_fit))){
      nr <- dim(bbouResults$surv_fit$samples$b0)[1]*dim(bbouResults$surv_fit$samples$b0)[2]
    }else{
      if(sum(grepl("Sbar",colnames(bbouResults$surv_fit$samples[[1]]),fixed=T))>0){
        divBy=2
      }else{
        divBy=1
      }
      nr <- dim(bbouResults$surv_fit$samples[[1]])[1]*dim(bbouResults$surv_fit$samples[[1]])[2]*length(bbouResults$surv_fit$samples)/divBy
    }
    
    popInfo <- merge(data.frame(id=seq(1:nr)),N0)
    popInfo$c <- compositionBiasCorrection(q=runif(nrow(popInfo),ccPars$qMin,ccPars$qMax),w=ccPars$cowMult,u=runif(nr,ccPars$uMin,ccPars$uMax),
                                           z=runif(nr,ccPars$zMin,ccPars$zMax))
    #print(paste("getSimsInitial",mean(popInfo$c)))
    parsBar <- caribouPopSimMCMC(popInfo,bbouResults$recruit_fit,bbouResults$surv_fit,progress=F,
                                 correlateRates=ccPars$correlateRates,returnExpected=T,...)
    pars <- caribouPopSimMCMC(popInfo,bbouResults$recruit_fit,bbouResults$surv_fit,progress=F,
                              correlateRates=ccPars$correlateRates,...)
    nrow(pars);nrow(parsBar)
    pars <- merge(pars,parsBar)
    nrow(pars)

    popInfo$PopulationName <- popInfo$pop_name
    
    pars <- merge(pars,subset(popInfo,select=c(-N0,-pop_name)))
    
  }else{
    if(replicates == "all"){replicates = formals(getSimsNational)$replicates}
    covTableObs <- unique(subset(bbouResults$parTab, select=c(Anthro,fire_excl_anthro,Year)))
    N0s <- unique(N0$N0)
    if(length(N0s)>1){
      stop("Specify a single initial population size for trajectories from national model.")
    }
    pars<- getSimsNational(replicates = max(replicates,2),N0 = N0s,covTableObs = covTableObs,cPars=ccPars,...)
    
    if(replicates==1){
      pars=subset(pars,id==pars$id[1])
    }
    #pars$Year <- pars$Anthro
    pars$PopulationName <- "A"
  }  
  
  if(!is.element("Anthro",names(pars))){pars$Anthro=NA}
  if(!is.element("fire_excl_anthro",names(pars))){pars$fire_excl_anthro=NA}
  
  #get the lambda percentile for each id - to allow users to select extreme examples
  simSum <- pars  %>%
    group_by(id) %>%
    summarize(MeanLam = mean(lambda,na.rm=T))
  simSum <-simSum[order(simSum$MeanLam),]
  simSum$lamPercentile <- round(100*seq(1:nrow(simSum))/nrow(simSum))
  simSum$MeanLam=NULL
  pars <- merge(pars,simSum)

  pars <- convertTrajectories(pars)
  
  simBig <- summarizeCaribouPopSim(pars,returnSamples=returnSamples)
  
  if(is.element("surv_fit",names(bbouResults))){
    simBig$surv_data = bbouResults$surv_fit$data
    simBig$recruit_data = bbouResults$recruit_fit$data
    simBig$popInfo = popInfo
  }
  
  if(max(simBig$summary$Year)<=100){simBig$summary$Year=NULL}
  if(max(simBig$samples$Year)<=100){simBig$summary$Year=NULL}

  if (doSave) {
    message("Updating cached initial simulations.")
    assign(saveName, simBig, envir = cacheEnv)
  }
  
  if(rmSamples){
    simBig$samples<-NULL
  }
  return(simBig)
}
