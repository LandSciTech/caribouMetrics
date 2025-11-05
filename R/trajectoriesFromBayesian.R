trajectoriesFromBayesian <- function(bayesianResults, N0 = NULL,
                                     cPars=subset(getScenarioDefaults(),select=-iAnthro),
                                     returnSamples = TRUE, doSummary = TRUE, ...){
  if(is.element("Anthro",names(bayesianResults))){
    message("Anthro from bayesianResults")
    bayesianResults=list(parTab=bayesianResults)
  }
  
  #bayesianResults = bbouResultFile
  if(is.character(bayesianResults) && (length(bayesianResults) == 1) ){
    if(file.exists(bayesianResults)){
      bayesianResults <- readRDS(bayesianResults)
    }else{
      stop(paste("bayesianResults file not found,",bayesianResults))
    }
  }
  
  if(!is.element("pop_name",names(bayesianResults$parTab))){bayesianResults$parTab$pop_name=NA}
  
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  
  if(is.null(N0)){
    if(is.element("N0",names(bayesianResults$parTab))){
      N0 <- unique(subset(bayesianResults$parTab,select=c(pop_name,N0)))
    }else{
      N0 <- eval(formals(getSimsNational)$N0)
    }
  }
  
  if(length(N0)==1){
    if(class(bayesianResults$parTab) == "list"){
      bayesianResults$parTab <- as.data.frame(bayesianResults$parTab)
    }
    N0 = merge(data.frame(N0=N0),
               unique(subset(bayesianResults$parTab, select=c(pop_name))))
  }
  N0$PopulationName = N0$pop_name
  
  if(is.element("bboufit",class(bayesianResults$surv_fit))){
    nr <- dim(bayesianResults$surv_fit$samples$b0)[1]*dim(bayesianResults$surv_fit$samples$b0)[2]
  }else{
    if(sum(grepl("Sbar",colnames(bayesianResults$surv_fit$samples[[1]]),fixed=T))>0){
      divBy=2
    }else{
      divBy=1
    }
    nr <- dim(bayesianResults$surv_fit$samples[[1]])[1]*dim(bayesianResults$surv_fit$samples[[1]])[2]*length(bayesianResults$surv_fit$samples)/divBy
  }
  
  popInfo <- merge(data.frame(id=seq(1:nr)),N0)
  popInfo$c <- compositionBiasCorrection(q=runif(nrow(popInfo),ccPars$qMin,ccPars$qMax),w=ccPars$cowMult,u=runif(nr,ccPars$uMin,ccPars$uMax),
                                         z=runif(nr,ccPars$zMin,ccPars$zMax))
  #print(paste("getSimsInitial",mean(popInfo$c)))
  parsBar <- caribouPopSimMCMC(popInfo,
                               bayesianResults$recruit_fit,
                               bayesianResults$surv_fit,
                               progress=F,
                               correlateRates=ccPars$correlateRates,
                               returnExpected=T,
                               ...)
  pars <- caribouPopSimMCMC(popInfo,
                            bayesianResults$recruit_fit,
                            bayesianResults$surv_fit,
                            progress=F,
                            correlateRates=ccPars$correlateRates,
                            ...)

  pars <- merge(pars,parsBar)
  nrow(pars)
  
  pi <- bayesianResults$parTab
  pi$PopulationName <- pi$pop_name
  pi$R_bar=NULL;pi$S_bar=NULL;pi$N0=NULL;pi$pop_name=NULL
  pars <- merge(pars,pi)
  
  if(doSummary){
    simBig <- prepareTrajectories(pars, returnSamples = returnSamples)
    
    if(max(simBig$summary$Year)<=100){simBig$summary$Year=NULL}
    if(max(simBig$samples$Year)<=100){simBig$samples$Year=NULL}
    
  }else {
    simBig <- pars
  }
  
  if(is.element("surv_fit",names(bayesianResults))){
    simBig$surv_data = bayesianResults$surv_fit$data
    simBig$recruit_data = bayesianResults$recruit_fit$data
    simBig$popInfo = popInfo
  }
  return(simBig)
}