#' Get trajectories from a Bayesian model result
#' 
#' 
#'
#' @param bayesianResults A result from `estimateBayesianRates`
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
#' surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
#' recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
#' bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
#'                                          return_mcmc = TRUE)
#' 
#' trajB <- trajectoriesFromBayesian(bbouInformative)
#' str(trajB, max.level = 1)
#' plotTrajectories(trajB)
#' 
trajectoriesFromBayesian <- function(bayesianResults, N0 = NULL,
                                     cPars=demographyDefaults(),
                                     returnSamples = TRUE, doSummary = TRUE, ...){
  cPars <- demographyDefaults(cPars)
  
  if(hasName(bayesianResults,"Anthro")){
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
  
  if(!hasName(bayesianResults$parTab,"PopulationName")){bayesianResults$parTab$PopulationName="A"}
  
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  
  if(is.null(N0)){
    if(hasName(bayesianResults$parTab,"N0")){
      if(class(bayesianResults$parTab)=="list"){
        N0 <- getN0Pars(bayesianResults$parTab$N0)
      }else{
        N0 <- getN0Pars(bayesianResults$parTab)
      }
    }else{
      N0 <- eval(formals(trajectoriesFromNational)$N0)
    }
  }
  
  Nuse <- getN0Pars(unique(N0),popNames = unique(bayesianResults$parTab$PopulationName))
  
  if(nrow(Nuse)>length(unique(Nuse$PopulationName))){
    stop("Expecting a single value or distribution of N0 for each population.")
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
  
  popInfo <- merge(data.frame(id=seq(1:nr/length(unique(Nuse$PopulationName)))),Nuse)

  popInfo <- addN0Variation(popInfo)

  popInfo$c <- compositionBiasCorrection(q=runif(nrow(popInfo),ccPars$qMin,ccPars$qMax),w=ccPars$cowMult,
                                         u=runif(nrow(popInfo),ccPars$uMin,ccPars$uMax),
                                         z=runif(nrow(popInfo),ccPars$zMin,ccPars$zMax))
  
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

  if(class(bayesianResults$parTab)=="list"){
    pi <- bayesianResults$parTab$N0
  }else{
    pi <- bayesianResults$parTab
  }
  pi$R_bar=NULL;pi$S_bar=NULL;pi$N0=NULL; pi$N.sd=NULL; pi$N.lower=NULL; pi$N.upper=NULL
  pars <- merge(pars,pi)
  
  if(max(table(subset(pars,select=c("PopulationName","Year","id"))))>1){stop("Error in trajectoriesFromBayesian: trajectories are not uniquely id'd")}
  
  if(doSummary){
    simBig <- prepareTrajectories(pars, returnSamples = returnSamples)
    
    if(max(simBig$summary$Year)<=100){simBig$summary$Year=NULL}
    if(returnSamples){
      if(max(simBig$samples$Year)<=100){simBig$samples$Year=NULL}
    }
    
  }else {
    simBig <- pars
  }
  
  if(hasName(bayesianResults,"surv_fit")){
    simBig$surv_data = convertBbouData(bayesianResults$surv_fit$data)
    simBig$recruit_data = convertBbouData(bayesianResults$recruit_fit$data)
    simBig$popInfo = popInfo
  }
  return(simBig)
}
