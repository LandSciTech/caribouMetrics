#' Projections of population growth from demographic model summaries.
#'
#' @param replicates 
#' @param Rbar,Sbar Mean and standard deviation of R_bar and S_bar over time. See `estimateBayesianRates()$parList` for expected form.
#' @param Riv,Siv Parameters defining the distribution of interannual variation. See `estimateBayesianRates()$parList` for expected form.
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
                  type = "beta", cPars = subset(getScenarioDefaults(),select=-iAnthro), 
                  doSummary = T, returnSamples = T,
                  nthin=formals(bboutools::bb_fit_survival)$nthin,varPersists=T,...){

  # TO DO:
  # See https://github.com/LandSciTech/caribouMetrics/issues/146
  
  if(type=="beta"){
    results <- ratesFromBetaSummary(Rbar, Sbar, Riv, Siv, replicates, nthin,varPersists)
  }else{
    results <- ratesFromLogisticSummary(Rbar, Sbar, Riv, Siv, replicates, nthin,varPersists)
  }
  rr <- trajectoriesFromBayesian(results, N0 = N0,returnSamples = returnSamples, 
                                 doSummary = doSummary,cPars=cPars, ...)
  return(rr)
}

ratesFromLogisticSummary <- function(Rbar, Sbar, Riv, Siv, replicates, nthin, varPersists){
  #Assumes gaussian distributed variation in means, gaussian random effect of year and log link.
  #Adapted from Shimoda QC workflow
  
  nc <- 2      # number of chains
  niters <- round(replicates/nc)
  ni <- niters * nthin   # number of samples for each chain
  nb <- ni / 2    # number of samples to discard as burnin
  
  s_priors <- list(iv_mean= Siv$S_iv_mean,iv_shape = Siv$S_iv_shape)
  sparams = c("Sbar","Survival") 
  r_priors <- list(iv_mean= Riv$R_iv_mean,iv_shape = Riv$R_iv_shape)
  params = c("Rbar","Recruitment") 
  if(varPersists){
    surv_fit <- rateFromLogisticSummary(Sbar,s_priors,sparams,nc,nthin,ni,nb)
    recruit_fit <- rateFromLogisticSummary(Rbar,r_priors,params,nc,nthin,ni,nb)
  }else{
    stop("varPersists = F is not an option for bbou logistic models")
  }
  
  parTab <- subset(Rbar,select=intersect(c("Year","PopulationName","Anthro","Fire_excl_anthro","PopulationID"),names(Rbar)))
  results <- list(parTab=parTab,surv_fit=surv_fit,recruit_fit=recruit_fit)
  
  return(results)
}

rateFromLogisticSummary<- function(Rbar,r_priors,params,nc,nt,ni,nb,fname="RatesLogistic_"){
  
  if(!is.element("adjust.mu",names(Rbar))){Rbar$adjust.mu = 0}
  if(!is.element("adjust.sd",names(Rbar))){Rbar$adjust.sd = 0}
      
  #####
  # Model - assign model file name (needed to match with the file name in the model object in jags.model())
  
  ffname = tempfile(pattern = fname, fileext = ".txt")
  sink(ffname)
  if(r_priors$iv_mean>0){
    cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  logit(Rbar[Annual[i],PopulationID[i]]) <- logit(rmu[i]+pm.r[i]) + rsd[i]*qq[PopulationID[i]]
  logit(Recruitment[Annual[i],PopulationID[i]]) <- logit(rmu[i]+pm.r[i]) + rsd[i]*qq[PopulationID[i]] + iv[i]  
  iv[i] ~ dnorm(0,1/pow(iv.sd[PopulationID[i]],2))
  pm.r[i]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
  iv.sd[k]~dgamma(iv_shape,iv_shape/iv_mean)
}

}

", fill = TRUE)
  }else{
    cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  logit(Rbar[Annual[i],PopulationID[i]]) <- logit(rmu[i]+pm.r[i]) + rsd[i]*qq[PopulationID[i]]
  logit(Recruitment[Annual[i],PopulationID[i]]) <- logit(rmu[i]+pm.r[i]) + rsd[i]*qq[PopulationID[i]]  
  pm.r[i]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
}

}

", fill = TRUE)
    
  }
sink()

### Define data, parameters, initials and settings
# Setting the data that you want to pass the model objects
data <- Rbar
data <- data[order(data$Annual,data$PopulationName),]
data$PopulationID <- as.factor(data$PopulationName)
nAnnual <- length(unique(data$Annual))
if(min(as.integer(data$Annual))>1){
  data$Annual = as.factor(data$Year)
}

datal = list(
  nObs = nrow(data),
  nAnnual=nAnnual,
  Annual = as.integer(data$Annual),
  nPops=length(unique(data$PopulationName)),
  PopulationID = as.integer(data$PopulationID),
  rmu=data$mean,
  rsd=data$sd,
  pmr.mu = data$adjust.mu,
  pmr.sd = data$adjust.sd)
if(r_priors$iv_mean>0){
  datal <- c(datal,r_priors)
}
inits = rjags::parallel.seeds("base::BaseRNG", nc) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain

return(jagsRunAndSummarize(data,datal,params,ffname,inits,nc,ni,nb,nt))
}


ratesFromBetaSummary <- function(Rbar, Sbar, Riv, Siv, replicates, nthin, varPersists){
  #Assumes gaussian distributed variation in means, beta distributed interannual variation
  #Adapted from Shimoda QC workflow
  
  nc <- 2      # number of chains
  niters <- round(replicates/nc)
  ni <- niters * nthin   # number of samples for each chain
  nb <- ni / 2    # number of samples to discard as burnin

  s_priors <- list(cv_min= Siv$S_cv_min,cv_max = Siv$S_cv_max)
  sparams = c("Sbar","Survival") 
  r_priors <- list(cv_min= Riv$R_cv_min,cv_max = Riv$R_cv_max)
  params = c("Rbar","Recruitment") 
  if(varPersists){
    surv_fit <- rateFromBetaSummary(Sbar,s_priors,sparams,nc,nthin,ni,nb)
    recruit_fit <- rateFromBetaSummary(Rbar,r_priors,params,nc,nthin,ni,nb)
  }else{
    #treats all variation as interannual variation
    surv_fit <- rateFromBetaSummaryYS(Sbar,s_priors,sparams,nc,nthin,ni,nb)
    recruit_fit <- rateFromBetaSummaryYS(Rbar,r_priors,params,nc,nthin,ni,nb)
  }
  
  parTab <- subset(Rbar,select=intersect(c("Year","PopulationName","Anthro","Fire_excl_anthro","PopulationID"),names(Rbar)))
  results <- list(parTab=parTab,surv_fit=surv_fit,recruit_fit=recruit_fit)
  
  return(results)
}

rateFromBetaSummary<- function(Rbar,r_priors,params,nc,nt,ni,nb,fname="RatesBeta_"){
  
  if(!is.element("adjust.mu",names(Rbar))){
    Rbar$adjust.mu = 0
  }
  
  if(!is.element("adjust.sd",names(Rbar))){
    Rbar$adjust.sd = 0
  }
  
  #####
  # Model - assign model file name (needed to match with the file name in the model object in jags.model())
  
  ffname = tempfile(pattern = fname, fileext = ".txt")
  sink(ffname)
  if(r_priors$cv_max>0){
    cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
	  mu.R[t,k]<- min(R[t,k] + pm.r[t,k],0.999) # Anything above 1 will be replaced by 1 (min-max truncation)
    sig.R[t,k] <- min(cv.R[k]*mu.R[t,k],0.999*pow(mu.R[t,k]*(1-mu.R[t,k]),0.5)) 
    alpha[t,k] <- ((1-mu.R[t,k])/pow(sig.R[t,k],2) - 1/mu.R[t,k]) * pow(mu.R[t,k],2)
    beta[t,k] <- alpha[t,k] * (1/mu.R[t,k] - 1)
    Rbar[t,k] <- alpha[t,k]/(alpha[t,k]+beta[t,k])
    Recruitment[t,k] ~ dbeta(alpha[t,k],beta[t,k])T(0.001,0.999)
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  R[Annual[i],PopulationID[i]]<-max(r[i],0.01)
  r[i] <- rmu[i] + rsd[i]*qq[PopulationID[i]]  
  pm.r[Annual[i],PopulationID[i]]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
  cv.R[k]~dunif(cv_min,cv_max)
}

}

", fill = TRUE)
  }else{
    cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
	  Rbar[t,k]<- min(R[t,k] + pm.r[t,k],1) # Anything above 1 will be replaced by 1 (min-max truncation)
    Recruitment[t,k] <- Rbar[t,k]
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  R[Annual[i],PopulationID[i]]<-max(r[i],0.01)
  r[i] <- rmu[i] + rsd[i]*qq[PopulationID[i]]  
  pm.r[Annual[i],PopulationID[i]]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
}

}

", fill = TRUE)
    
  }
sink()

### Define data, parameters, initials and settings
# Setting the data that you want to pass the model objects
data <- Rbar
data <- data[order(data$Annual,data$PopulationName),]
data$PopulationID <- as.factor(data$PopulationName)
nAnnual <- length(unique(data$Annual))
if(min(as.integer(data$Annual))>1){
  data$Annual = as.factor(data$Year)
}

datal = list(
  nObs = nrow(data),
  nAnnual=nAnnual,
  Annual = as.integer(data$Annual),
  nPops=length(unique(data$PopulationName)),
  PopulationID = as.integer(data$PopulationID),
  rmu=data$mean,
  rsd=data$sd,
  pmr.mu = data$adjust.mu,
  pmr.sd = data$adjust.sd)
if(r_priors$cv_max>0){
  datal <- c(datal,r_priors)
}
inits = parallel.seeds("base::BaseRNG", nc) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain

return(jagsRunAndSummarize(data,datal,params,ffname,inits,nc,ni,nb,nt))
}


rateFromBetaSummaryYS<- function(Rbar,r_priors,params,nc,nt,ni,nb,fname="PopDynMod_"){
  #####
  # Model - assign model file name (needed to match with the file name in the model object in jags.model())
  
  ffname = tempfile(pattern = fname, fileext = ".txt")
  sink(ffname)
  if(r_priors$cv_max>0){
  cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
	  mu.R[t,k]<- min(R[t,k] + pm.r[t,k],0.999) # Anything above 1 will be replaced by 1 (min-max truncation)
    sig.R[t,k] <- min(cv.R[k]*mu.R[t,k],0.999*pow(mu.R[t,k]*(1-mu.R[t,k]),0.5)) 
    alpha[t,k] <- ((1-mu.R[t,k])/pow(sig.R[t,k],2) - 1/mu.R[t,k]) * pow(mu.R[t,k],2)
    beta[t,k] <- alpha[t,k] * (1/mu.R[t,k] - 1)
    Rbar[t,k] <- alpha[t,k]/(alpha[t,k]+beta[t,k])
    Recruitment[t,k] ~ dbeta(alpha[t,k],beta[t,k])T(0.001,0.999)
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  R[Annual[i],PopulationID[i]]<-max(r[i],0.01)
  #r[i] <- rmu[i] + rsd[i]*qq[PopulationID[i]]  
  r[i]~dnorm(rmu[i],rtau[i])
  rtau[i]<-1/pow(rsd[i],2)
  pm.r[Annual[i],PopulationID[i]]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
  cv.R[k]~dunif(cv_min,cv_max)
}

}

", fill = TRUE)
  }else{
    cat("

model {

for (t in 1:nAnnual) {
  for (k in 1:nPops) {
	  Rbar[t,k]<- min(R[t,k] + pm.r[t,k],1) # Anything above 1 will be replaced by 1 (min-max truncation)
    Recruitment[t,k] <- Rbar[t,k]
    Sbar[t,k] <- Rbar[t,k]
    Survival[t,k] <- Recruitment[t,k]
}}

for(i in 1:nObs) {
  R[Annual[i],PopulationID[i]]<-max(r[i],0.01)
  #r[i] <- rmu[i] + rsd[i]*qq[PopulationID[i]]  
  r[i]~dnorm(rmu[i],rtau[i])
  rtau[i]<-1/pow(rsd[i],2)
  pm.r[Annual[i],PopulationID[i]]<-pm.rr[i]*ifelse(max(pmr.mu[i],pmr.sd[i]) > 0, 1,0)
  pm.rr[i]~dnorm(pmr.mu[i], pmr.tau[i])
	pmr.tau[i]<-1/pow(max(pmr.sd[i],0.0000000001),2)
}

for (k in 1:nPops) {
  qq[k]~dnorm(0,1)
}

}

", fill = TRUE)
    
  }
  sink()
  
  ### Define data, parameters, initials and settings
  # Setting the data that you want to pass the model objects
  data <- Rbar
  data <- data[order(data$Annual,data$PopulationName),]
  data$PopulationID <- as.factor(data$PopulationName)
  nAnnual <- length(unique(data$Annual))
  if(min(as.integer(data$Annual))>1){
    data$Annual = as.factor(data$Year)
  }
  
  datal = list(
    nObs = nrow(data),
    nAnnual=nAnnual,
    Annual = as.integer(data$Annual),
    nPops=length(unique(data$PopulationName)),
    PopulationID = as.integer(data$PopulationID),
    rmu=data$mean,
    rsd=data$sd,
    pmr.mu = data$adjust.mu,
    pmr.sd = data$adjust.sd)
  if(r_priors$cv_max>0){
    datal <- c(datal,r_priors)
  }
  inits = rjags::parallel.seeds("base::BaseRNG", nc) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain
  
  return(jagsRunAndSummarize(data,datal,params,ffname,inits,nc,ni,nb,nt))
}



