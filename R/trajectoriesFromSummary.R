#' Projections of population growth from demographic model summaries.
#'
#' @param replicates 
#' @param Rbar,Sbar Mean and standard deviation of R_bar and S_bar over time. See `estimateBayesianRates()$parList` for expected form.
#' @param Riv,Siv Parameters defining the distribution of interannual variation. See `estimateBayesianRates()$parList` for expected form.
#' @param addl_params TO DO explain population mgmt parameters
#' @param type The distribution of interannual variation varies between "beta" or "bbou" model types. 
#' @param doSummary logical. Default TRUE. If FALSE returns unprocessed outcomes from caribouPopGrowth. 
#'  If TRUE returns summaries and (if returnSamples = T) sample trajectories from prepareTrajectories.
#' @param returnSamples logical. If FALSE returns only summaries. If TRUE
#'   returns example trajectories as well. 
#' @param ... Additional arguments passed to `caribouPopGrowth`
#' @return a data.frame
#' @family demography
#'
#' @export
#' 
#' @examples
#'
#'
trajectoriesFromSummary <- function(replicates, N0, Rbar, Sbar, Riv, Siv,  
                  type = "beta", addl_params = list(), doSummary = T, returnSamples = T,
                  nthin=formals(bboutools::bb_fit_survival)$nthin,...){

  # TO DO: check for expected elements of input args
  #Finish documenting trajectoriesFromSummary
  #Tests for trajectoriesFromSummary
  #Update vignettes to use trajectoriesFromSummary, and add variation over time + population mgmt example to vignettes.
  #Think about reducing compute requirements for disturbance projection updating.
  #Eventually - update app to include disturbance, and remove trajectoriesFromSummaryForApp and associated parTab.
  
  if(type=="beta"){
    results <- ratesFromBetaSummary(N0, Rbar, Sbar, Riv, Siv, addl_params, replicates, nthin)
  }else{
    stop("Handle bbou case")
  }
  
  rr <- trajectoriesFromBayesian(results, cPars = addl_params,returnSamples = returnSamples, 
                                 doSummary = doSummary, ...)
  
  return(rr)
  
}

ratesFromBetaSummary <- function(N0, Rbar, Sbar, Riv, Siv, addl_params, replicates, nthin){
  #Assumes gaussian distributed variation in means, beta distributed interannual variation
  #Adapted from Shimoda QC workflow
  
  nc <- 2      # number of chains
  niters <- round(replicates/nc)
  ni <- niters * nthin   # number of samples for each chain
  nb <- ni / 2    # number of samples to discard as burnin

  r_priors <- c(cv_min= Siv$S_cv_min,cv_max = Siv$S_cv_max)
  pm.r = addl_params$pm.s
  params = c("Sbar","Survival") 
  surv_fit <- rateFromBetaSummary(Sbar,Siv,addl_params,r_priors,pm.r,params,nc,nthin,ni,nb)
  r_priors <- c(cv_min= Riv$R_cv_min,cv_max = Riv$R_cv_max)
  pm.r = addl_params$pm.r
  params = c("Rbar","Recruitment") 
  recruit_fit <- rateFromBetaSummary(Rbar,Riv,addl_params,r_priors,pm.r,params,nc,nthin,ni,nb)

  results <- list(parTab=N0,surv_fit=surv_fit,recruit_fit=recruit_fit)
  
  return(results)
}

rateFromBetaSummary<- function(Rbar,Riv,addl_params,r_priors,pm.r,params,nc,nt,ni,nb,fname="PopDynMod_"){
  #####
  # Model - assign model file name (needed to match with the file name in the model object in jags.model())
  
  ffname = tempfile(pattern = fname, fileext = ".txt")
  sink(ffname)
  cat("

model {

for (t in 1:nAnnual) {
  pm.period[t]<-ifelse(t > pm.start && t <= pm.end, 1,0) # pm treated period (to set return 1, to turn off return 0)
  for (k in 1:nPops) {
	  mu.R[t,k]<- min(R[t,k] + pm.r[t]*pm.period[t],0.99) # Anything above 1 will be replaced by 1 (min-max truncation)
    sig.R[t,k] <- min(cv.R[k]*mu.R[t,k],0.99*pow(mu.R[t,k]*(1-mu.R[t,k]),0.5)) 
    alpha[t,k] <- ((1-mu.R[t,k])/pow(sig.R[t,k],2) - 1/mu.R[t,k]) * pow(mu.R[t,k],2)
    beta[t,k] <- alpha[t,k] * (1/mu.R[t,k] - 1)
    Rbar[t,k] <- alpha[t,k]/(alpha[t,k]+beta[t,k])
    Recruitment[t,k] ~ dbeta(alpha[t,k],beta[t,k])T(0.01,0.99)
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  R[Annual[i],PopulationID[i]]<-max(r[i],0.01)
  r[i] <- rmu[i] + rsd[i]*qq[PopulationID[i]]  
}

for (t in 1:nAnnual) {
	pm.r[t]~dnorm(pmr.mu[t], pmr.tau[t])
	pmr.tau[t]<-1/pow(pmr.sd[t],2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
  cv.R[k]~dunif(cv_min,cv_max)
}

}

", fill = TRUE)
  sink()
  
  ### Define data, parameters, initials and settings
  # Setting the data that you want to pass the model objects
  data <- Rbar
  data <- data[order(data$Annual,data$PopulationName),]
  data$PopulationID <- as.factor(data$PopulationName)
  nAnnual <- length(unique(data$Annual))
  
  # Setting the population management timing and duration
  pm.start <- addl_params$pm.startYear-min(data$Year) # vector location
  pm.end <- addl_params$pm.endYear-min(data$Year)  # vector location
  
  datal = list(
    nObs = nrow(data),
    nAnnual=nAnnual,
    Annual = as.integer(data$Annual),
    nPops=length(unique(data$PopulationName)),
    PopulationID = as.integer(data$PopulationID),
    rmu=data$mean,
    rsd=data$sd,
    pm.start=pm.start,
    pm.end=pm.end,
    pmr.mu = rep(pm.r[,2],nAnnual),
    pmr.sd = rep(pm.r[,3],nAnnual))
  datal <- c(datal,r_priors)
  
  inits = parallel.seeds("base::BaseRNG", nc) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain
  
  return(jagsRunAndSummarize(data,datal,params,ffname,inits,nc,ni,nb,nt))
}



