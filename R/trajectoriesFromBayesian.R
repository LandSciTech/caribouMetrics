#' Get trajectories from a Bayesian model result
#' 
#' 
#'
#' @param bayesianResults 
#' @param N0 
#' @param cPars 
#' @param returnSamples 
#' @param doSummary 
#' @param ... 
#'
#' @returns
#' @export
#' @family demography
#' @examples
trajectoriesFromBayesian <- function(bayesianResults, N0 = NULL,
                                     cPars=subset(getScenarioDefaults(),select=-iAnthro),
                                     returnSamples = TRUE, doSummary = TRUE, ...){
  cPars <- getScenarioDefaults(cPars)
  
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
  
  if(!is.element("PopulationName",names(bayesianResults$parTab))){bayesianResults$parTab$PopulationName=NA}
  
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  
  if(is.null(N0)){
    if(is.element("N0",names(bayesianResults$parTab))){
      if(class(bayesianResults$parTab)=="list"){
        N0 <- unique(subset(bayesianResults$parTab$N0,select=c(PopulationName,N0)))
      }else{
        N0 <- unique(subset(bayesianResults$parTab,select=c(PopulationName,N0)))
      }
    }else{
      N0 <- eval(formals(trajectoriesFromNational)$N0)
    }
  }
  
  if(length(N0)==1){
    if(class(bayesianResults$parTab) == "list"){
      bayesianResults$parTab <- as.data.frame(bayesianResults$parTab)
    }
    N0 = merge(data.frame(N0=N0),
               unique(subset(bayesianResults$parTab, select=c(PopulationName))))
  }
  
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
  
  Nnames <- intersect(c("N0","N.sd","N.lower","N.upper"),names(N0))
  Nuse <- unique(subset(N0,select=Nnames))
  if(nrow(Nuse)>1){
    stop("Handle this case - variation in initial N among populations")
  }else{
    popInfo <- merge(data.frame(id=seq(1:nr)),Nuse)
  }
  
  if(length(intersect(c("N.sd","N.lower"),names(Nuse)))>0){
    #model variation in N0
    popInfo$mean = popInfo$N0
    if(is.element("N.sd",names(popInfo))){popInfo$N0 = rnorm(nrow(popInfo),mean=popInfo$mean,sd=popInfo$N.sd)}
    if(is.element("N.lower",names(popInfo))){popInfo$N0=pmax(popInfo$N.lower,popInfo$N0)}
    if(is.element("N.upper",names(popInfo))){popInfo$N0=pmin(popInfo$N.upper,popInfo$N0)}
    popInfo$N0 <- round(popInfo$N0)
    popInfo <- subset(popInfo,select=setdiff(names(popInfo),c("mean","N.sd","N.lower","N.upper")))
  }
  
  popInfo$c <- compositionBiasCorrection(q=runif(nrow(popInfo),ccPars$qMin,ccPars$qMax),w=ccPars$cowMult,u=runif(nr,ccPars$uMin,ccPars$uMax),
                                         z=runif(nr,ccPars$zMin,ccPars$zMax))

  parsBar <- simulateTrajectoriesFromPosterior(popInfo=popInfo,
                               rec_pred=bayesianResults$recruit_fit,
                               surv_pred=bayesianResults$surv_fit,
                               progress=F,
                               correlateRates=ccPars$correlateRates,
                               returnExpected=T,
                               ...)
  pars <- simulateTrajectoriesFromPosterior(popInfo=popInfo,
                            rec_pred=bayesianResults$recruit_fit,
                            surv_pred=bayesianResults$surv_fit,
                            progress=F,
                            correlateRates=ccPars$correlateRates,
                            ...)

  pars <- merge(pars,parsBar)
  nrow(pars)
  
  if(class(bayesianResults$parTab)=="list"){
    pi <- bayesianResults$parTab$N0
  }else{
    pi <- bayesianResults$parTab
  }
  pi$R_bar=NULL;pi$S_bar=NULL;pi$N0=NULL
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
