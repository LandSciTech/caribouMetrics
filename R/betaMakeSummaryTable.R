betaMakeSummaryTable <- function(surv_data, recruit_data, disturbance, ...){
  #Note: using bboutools to check and structure the data without fitting the models...
  recruit_fit_in <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)
  recruit_fit <- betaRecruitment(recruit_fit_in,disturbance)
  
  surv_fit_in <- bboutools::bb_fit_survival(surv_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)
  surv_fit <- betaSurvival(surv_fit_in,disturbance)
  
  return(list(parTab=NULL,surv_fit=surv_fit,recruit_fit=recruit_fit))
}

betaSurvival <-function(surv_fit,disturbance){
  data <- surv_fit$data
  data <- as.data.frame(data)
  data <- subset(data,is.element(Annual,disturbance$Year))
  data$Annual <- as.factor(as.character(data$Annual))
  disturbance <- merge(disturbance,unique(subset(data,select=c(Annual,Year,PopulationID))))
  anthro <- spread(subset(disturbance,select=c(Annual,PopulationID,Anthro)), PopulationID, Anthro)
  if(ncol(anthro)==2){
    anthro = matrix(anthro[,2:ncol(anthro)],ncol= ncol(anthro)-1)
  }else{
    anthro = anthro[,2:ncol(anthro)]
  }
  
  datal <- list(nObs = nrow(data),
                StartTotal = data$StartTotal,
                Mortalities = data$Mortalities,
                nAnnual = length(unique(data$Annual)),
                Annual = as.integer(data$Annual),
                nPops = length(unique(data$PopulationName)),
                PopulationID = as.integer(data$PopulationID),
                anthro = anthro
  )

  priors <- getPriors()
  surv_priors <- priors[c("l.Saf.Prior1","l.Saf.Prior2","beta.Saf.Prior1","beta.Saf.Prior2","sig.Saf.Prior1","sig.Saf.Prior2")] # TO DO add composition bias: "bias.Prior1","bias.Prior2")]
  names(surv_priors)<- c("b0_mu","b0_sd","b1_mu","b1_sd","sig.S.Prior1","sig.S.Prior2")

  datal <- c(datal,surv_priors)
  
  ###### Survival Model ######
  sink("survival_model.txt")  # assign model file name #
  cat("

model {
  for(i in 1:nObs) {
    eSurvival[i] <- AnnualSurvival[Annual[i],PopulationID[i]]
    Mortalities[i] ~ dbin(1 - eSurvival[i], StartTotal[i])
  }

  for (i in 1:nAnnual) {
    for (k in 1:nPops) {
      mu.S[i,k] <- max(0.01,min(0.99,(46*exp(b0[k] + b1[k]*anthro[i,k])-0.5)/45))
      sig.S[i,k] <- min(cv.S[k]*mu.S[i,k],0.99*pow(mu.S[i,k]*(1-mu.S[i,k]),0.5))
      alpha[i,k] <- ((1-mu.S[i,k])/pow(sig.S[i,k],2) - 1/mu.S[i,k]) * pow(mu.S[i,k],2)
      beta[i,k] <- alpha[i,k] * (1/mu.S[i,k] - 1)
      Survival[i,k]~ dbeta(alpha[i,k],beta[i,k])T(0.01,0.99)
      AnnualSurvival[i,k] <- pow(Survival[i,k],1/12)
    }
  }

  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
    cv.S[k]~dunif(sig.S.Prior1,sig.S.Prior2)
  }

}

", fill = TRUE)
  sink()
  
  # Setting parameters - setting parameters that you want to monitor
  params = c("Survival")
  
  # Setting initial values - not assigned
  inits1 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits2 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits3 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits <- list(inits1,inits2,inits3)
  
  # option2: seed setting
  # inits = parallel.seeds("base::BaseRNG", 2) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain
  
  # MCMC settings - (bboutools default: 1000 MCMC samples from 3 chains, number of )
  nc <- 3      # number of chains
  ni <- 1000   # number of samples for each chain
  nb <- 500    # number of samples to discard as burnin
  nt <- 10     # thinning rate
  
  #################################################################################################################
  ### Running JAGS
  # Create a model object - this compiles and initialize the model (if adaptation is required then a prgress bar made of '+' signs will be printed)
  model.fit <- rjags::jags.model(file="survival_model.txt", data=datal, n.adapt=nb, n.chains = nc)
  
  # To get samples from the posterior distribution of the parameters
  update(model.fit, n.iter=ni)
  model.samples <- rjags::coda.samples(model.fit, params, n.iter=ni, thin = nt)
  surv_data <- subset(data,is.element(Year,disturbance$Year))
  return(list(data=data,samples=model.samples)) 
}

betaRecruitment <- function(rec_fit, disturbance){
  rec_data <- rec_fit$data
  data <- as.data.frame(rec_data)
  data <- merge(data,disturbance)
  data$Annual <- as.factor(as.character(data$Annual))
  
  data$Cows[is.na(data$Cows)]=0
  data$CowsBulls[is.na(data$CowsBull)]=0
  data$UnknownAdults[is.na(data$UnknownAdults)]=0
  data$Yearlings[is.na(data$Yearlings)]=0
  
  anthro <- spread(subset(data,select=c(Annual,PopulationID,Anthro)), PopulationID, Anthro)
  fire <- spread(subset(data,select=c(Annual,PopulationID,fire_excl_anthro)), PopulationID, fire_excl_anthro)
  if(ncol(anthro)==2){
    anthro = matrix(anthro[,2:ncol(anthro)],ncol= ncol(anthro)-1)
    fire = matrix(fire[,2:ncol(fire)],ncol= ncol(fire)-1)
  }else{
    anthro = anthro[,2:ncol(anthro)]
    fire = fire[,2:ncol(fire)]
  }
  
  datal <- list(
    nAnnual = length(unique(data$Annual)),
    nPops = length(unique(data$PopulationName)),
    nObs = nrow(data),
    Cows = data$Cows,
    CowsBulls = data$CowsBulls,
    UnknownAdults = data$UnknownAdults,
    Yearlings = data$Yearlings,
    Calves = data$Calves,
    Annual = as.integer(data$Annual),
    PopulationID = as.integer(data$PopulationID),
    anthro = anthro,
    fire = fire,
    adult_female_proportion =0.65,
    sex_ratio=0.5
  )

  # Setting priors from caribouMetrics
  priors <- getPriors()
  
  rec_priors = priors[c("l.R.Prior1","l.R.Prior2","beta.Rec.anthro.Prior1","beta.Rec.anthro.Prior2","beta.Rec.fire.Prior1","beta.Rec.fire.Prior2","sig.R.Prior1","sig.R.Prior2")]
  names(rec_priors)<- c("b0_mu","b0_sd","b1_mu","b1_sd","b2_mu","b2_sd","sig.R.Prior1","sig.R.Prior2")
  
  datal <- c(datal,rec_priors)
  
  # Currently set "adult_female_proportion" to be constant - potential to consider to be stochastic (beta distribution)
  # adult_female_proportion_alpha=rec_priors[["adult_female_proportion_alpha"]]
  # adult_female_proportion_beta=rec_priors[["adult_female_proportion_beta"]]
  
  ###### Recruitment Model ######
  sink("recruit_model.txt")  # observation model can be connected by CaribouYear[i] 
  cat("

model {
  for (i in 1:nObs){
    FemaleYearlings[i] ~ dbin(sex_ratio, Yearlings[i])
    Cows[i] ~ dbin(adult_female_proportion, CowsBulls[i])
    OtherAdultsFemales[i] ~ dbin(adult_female_proportion, UnknownAdults[i])
    AdultsFemales[i] <- max(FemaleYearlings[i] + Cows[i] + OtherAdultsFemales[i], 1)
    Calves[i] ~ dbin(Recruitment[Annual[i],PopulationID[i]], AdultsFemales[i])
  }

  for (i in 1:nAnnual) {
    for (k in 1:nPops) {
      mu.R[i,k] <- max(0.01,min(0.99,exp(b0[k] + b1[k]*anthro[i,k] + b2[k]*fire[i,k])))
      sig.R[i,k] <- min(cv.R[k]*mu.R[i,k],0.99*pow(mu.R[i,k]*(1-mu.R[i,k]),0.5)) # Constrain on sig.R to fall within theoretical range
      alpha[i,k] <- ((1-mu.R[i,k])/pow(sig.R[i,k],2) - 1/mu.R[i,k]) * pow(mu.R[i,k],2)
      beta[i,k] <- alpha[i,k] * (1/mu.R[i,k] - 1)
      Recruitment[i,k] ~ dbeta(alpha[i,k],beta[i,k])T(0.01,0.99)
    }
  }

  #adult_female_proportion ~ dbeta(adult_female_proportion_alpha, adult_female_proportion_beta)
  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
    b2[k] ~ dnorm(b2_mu, 1/pow(b2_sd,2))
    cv.R[k]~dunif(sig.R.Prior1,sig.R.Prior2)
  }
}
  
", fill = TRUE)
  sink()
  
  ######## Define data, parameters, initials and settings #####
  # Setting parameters - setting parameters that you want to monitor
  params = c("Recruitment")
  
  # Setting initial values
  inits1 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits2 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits3 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits <- list(inits1,inits2,inits3)
  
  # option2: seed setting
  # inits = parallel.seeds("base::BaseRNG", 2) # For MCMC reproducibility: returns a list of values that may be used to initialize the random number generator of each chain
  
  # MCMC settings - (bboutools default: 1000 MCMC samples from 3 chains, number of )
  nc <- 3      # number of chains
  ni <- 1000   # number of samples for each chain
  nb <- 500    # number of samples to discard as burnin
  nt <- 10     # thinning rate
  
  ########################## Running JAGS ########################################################################
  ### Running JAGS
  # Create a model object - this compiles and initialize the model (if adaptation is required then a prgress bar made of '+' signs will be printed)

  model.fit <- rjags::jags.model(file="recruit_model.txt", data=datal, n.adapt=nb, n.chains = nc)
  
  # To get samples from the posterior distribution of the parameters
  update(model.fit, n.iter=ni)
  model.samples <- rjags::coda.samples(model.fit, params, n.iter=ni, thin = nt) #"bYear",
  rec_data <- subset(rec_data,is.element(Year,disturbance$Year))
  return(list(data=rec_data,samples=model.samples)) 
  
}

