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
  
  if(!hasName(bayesianResults$parTab,"PopulationName")){bayesianResults$parTab$PopulationName=NA}
  
  ccPars = unique(subset(cPars,select=c(qMin,qMax,uMin,uMax,zMin,zMax,cowMult,correlateRates)))
  if(nrow(ccPars)>1){
    stop("Do not include more than one composition bias scenario in cPars")
  }
  
  if(is.null(N0)){
    if(hasName(bayesianResults$parTab,"N0")){
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
  
  Nnames <- intersect(c("N0","N.sd","N.lower","N.upper","PopulationName"),names(N0))
  Nuse <- unique(subset(N0,select=Nnames))
  if(!hasName(Nuse,"PopulationName")){Nuse$PopulationName="A"}
  
  popInfo <- merge(data.frame(id=seq(1:nr/length(unique(Nuse$PopulationName)))),Nuse)

  if(length(intersect(c("N.sd","N.lower"),names(Nuse)))>0){
    #model variation in N0
    popInfo$mean = popInfo$N0
    if(hasName(popInfo,"N.sd")){popInfo$N0 = rnorm(nrow(popInfo),mean=popInfo$mean,sd=popInfo$N.sd)}
    if(hasName(popInfo,"N.lower")){popInfo$N0=pmax(popInfo$N.lower,popInfo$N0)}
    if(hasName(popInfo,"N.upper")){popInfo$N0=pmin(popInfo$N.upper,popInfo$N0)}
    popInfo$N0 <- round(popInfo$N0)
    popInfo <- subset(popInfo,select=setdiff(names(popInfo),c("mean","N.sd","N.lower","N.upper")))
  }
  
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
  pi$R_bar=NULL;pi$S_bar=NULL;pi$N0=NULL
  pars <- merge(pars,pi)
  
  if(max(table(subset(pars,select=c("PopulationName","Year","id"))))>1){stop("Error in trajectoriesFromBayesian: trajectories are not uniquely id'd")}
  
  if(doSummary){
    simBig <- prepareTrajectories(pars, returnSamples = returnSamples)
    
    if(max(simBig$summary$Year)<=100){simBig$summary$Year=NULL}
    if(max(simBig$samples$Year)<=100){simBig$samples$Year=NULL}
    
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
