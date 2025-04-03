betaMakeSummaryTable <- function(surv_data, recruit_data, disturbance, ...){
  #TO DO: allow using bboutools to structure the data without fitting the models...
  recruit_fit_in <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, ...)
  recruit_fit <- betaRecruitment(recruit_fit_in,disturbance)

  #TO DO: add survival beta model
  surv_fit <- bboutools::bb_fit_survival(surv_data, multi_pops = TRUE, allow_missing = TRUE, quiet = TRUE, ...)
  return(list(parTab=NULL,surv_fit=surv_fit,recruit_fit=recruit_fit))
}

betaRecruitment <- function(rec_fit, disturbance){
  rec_data <- rec_fit$data
  model_rec_data <- as.data.frame(rec_data)
  model_rec_data <- merge(model_rec_data,disturbance)
  
  model_rec_data$Cows[is.na(model_rec_data$Cows)]=0
  model_rec_data$CowsBulls[is.na(model_rec_data$CowsBull)]=0
  model_rec_data$UnknownAdults[is.na(model_rec_data$UnknownAdults)]=0
  model_rec_data$Yearlings[is.na(model_rec_data$Yearlings)]=0
  
  nObs = nrow(model_rec_data)
  Cows = model_rec_data$Cows
  CowsBulls = model_rec_data$CowsBulls
  UnknownAdults = model_rec_data$UnknownAdults
  Yearlings = model_rec_data$Yearlings
  Calves = model_rec_data$Calves
  Year = model_rec_data$Year
  nYear = length(model_rec_data$Year)
  adult_female_proportion =0.65
  anthro = model_rec_data$Anthro
  fire = model_rec_data$fire_excl_anthro
  
  # Setting priors from caribouMetrics
  priors <- getPriors()

  rec_priors = priors[c("l.R.Prior1","l.R.Prior2","beta.Rec.anthro.Prior1","beta.Rec.anthro.Prior2","beta.Rec.fire.Prior1","beta.Rec.fire.Prior2","sig.R.Prior1","sig.R.Prior2")]
  
  mu.beta0=rec_priors[["l.R.Prior1"]]
  sig.beta0=rec_priors[["l.R.Prior2"]]
  
  mu.beta1=rec_priors[["beta.Rec.anthro.Prior1"]]
  sig.beta1=rec_priors[["beta.Rec.anthro.Prior2"]]
  
  mu.beta2=rec_priors[["beta.Rec.fire.Prior1"]]
  sig.beta2=rec_priors[["beta.Rec.fire.Prior2"]]
  
  sig.R.Prior1=rec_priors[["sig.R.Prior1"]]
  sig.R.Prior2=rec_priors[["sig.R.Prior2"]]
  
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
  Calves[i] ~ dbin(Recruitment[i], AdultsFemales[i])
                }

for(i in 1:nYear) {
  mu.R[i] <- max(0.01, min(0.99, exp(beta0 + beta1*anthro[i] + beta2*fire[i])))
  #sig.R[i] <- min(cv.R*mu.R[i],0.99*pow(mu.R[i]*(1-mu.R[i]),0.5)) # Constrain on sig.R to fall within theoretical range
  sig.R[i] <- cv.R*mu.R[i]
  
  #alpha[i] <- ((1-mu.R[i])/pow(sig.R[i],2) - 1/mu.R[i]) * pow(mu.R[i],2)
  #beta[i] <- alpha[i] * (1/mu.R[i] - 1)
  alpha[i]<-max(0.01,mu.R[i]*(mu.R[i]*(1-mu.R[i])/sigma2-1))
  beta[i]<-max(0.01,(1-mu.R[i])*(mu.R[i]*(1-mu.R[i])/sigma2-1))
  Recruitment[i] ~ dbeta(alpha[i],beta[i])
                }

  #adult_female_proportion ~ dbeta(adult_female_proportion_alpha, adult_female_proportion_beta)
  cv.R~dunif(sig.R.Prior1,sig.R.Prior2)
   
  beta0 ~ dnorm(mu.beta0, tau.beta0)
  beta1 ~ dnorm(mu.beta1, tau.beta1)
  beta2 ~ dnorm(mu.beta2, tau.beta2)

  tau.beta0 <- 1/pow(sig.beta0,2)
  tau.beta1 <- 1/pow(sig.beta1,2)
  tau.beta2 <- 1/pow(sig.beta2,2)
  
  tau~dunif(0.0001,100)
  sigma<-sqrt(1/tau)
  sigma2<-1/tau
        }
  
", fill = TRUE)
  sink()
  
  ######## Define data, parameters, initials and settings #####
  # Setting the data that you want to pass the model objects
  data = list(anthro=anthro,fire=fire,nObs=nObs,nYear=nYear,sex_ratio=0.5,Cows=Cows,CowsBulls=CowsBulls,UnknownAdults=UnknownAdults,Yearlings=Yearlings,Calves=Calves,adult_female_proportion=adult_female_proportion,mu.beta0=mu.beta0,sig.beta0=sig.beta0,mu.beta1=mu.beta1,sig.beta1=sig.beta1,mu.beta2=mu.beta2,sig.beta2=sig.beta2,sig.R.Prior1=sig.R.Prior1,sig.R.Prior2=sig.R.Prior2) # list all the observations ,Year=Year
  
  # Setting parameters - setting parameters that you want to monitor
  params = c("Recruitment")
  
  # Setting initial values
  inits1 <- list(b0=0) #bYear=0)
  inits2 <- list(b0=0) #bYear=0)
  inits3 <- list(b0=0) #bYear=0)
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
  model.fit <- rjags::jags.model(file="recruit_model.txt", data=data, n.adapt=nb, n.chains = nc)
  
  # To get samples from the posterior distribution of the parameters
  update(model.fit, n.iter=ni)
  model.samples <- rjags::coda.samples(model.fit, params, n.iter=ni, thin = nt) #"bYear",
  
  rec_data <- subset(rec_data,is.element(Year,disturbance$Year))
  return(list(data=rec_data,samples=model.samples)) 
  
}

