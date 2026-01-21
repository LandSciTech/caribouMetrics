betaMakeSummaryTable <- function(surv_data, recruit_data, disturbance,priors,nc,nt,ni,nb){
  # if(n_distinct(surv_data$PopulationName) > 1){
  #   testTable(disturbance, req_col_names = "PopulationName", 
  #             req_vals = unique(surv_data$PopulationName))
  # }
  #Note: using bboutools to check and structure the data without fitting the models...0
  surv_fit_in <- bboutools::bb_fit_survival(surv_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)
  surv_fit <- betaSurvival(surv_fit_in,disturbance,priors,nc,nt,ni,nb)
  
  recruit_fit_in <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)
  recruit_fit <- betaRecruitment(recruit_fit_in,disturbance,priors,nc,nt,ni,nb)
  
  
  summaries <- rbind(recruit_fit$summaries,surv_fit$summaries)
  
  parTab <- list(numSteps=1)
  parTab$Recruitment <- subset(summaries,MetricTypeID=="Recruitment")
  parTab$Survival <- subset(summaries,MetricTypeID=="Survival")
  parTab$Rbar <- subset(summaries,MetricTypeID=="Rbar")
  parTab$Sbar <- subset(summaries,MetricTypeID=="Sbar")
  parTab$Siv <- priors#subset(summaries,MetricTypeID=="sig.R")
  parTab$Riv <- priors#subset(summaries,MetricTypeID=="sig.R")
  parTab$type <- "betaTime"

  #numSteps, replicates, N0, R_bar, S_bar, R_sd, S_sd,
  #R_iv_mean, R_iv_shape, S_iv_mean, S_iv_shape,  
  #scn_nm, type = "logistic"
  
  return(list(parTab=parTab,surv_fit=surv_fit,recruit_fit=recruit_fit))
}

betaSurvival <-function(surv_fit,disturbance,priors,nc,nt,ni,nb){
  data <- surv_fit$data
  data <- as.data.frame(data)
  data <- subset(data,is.element(Annual,disturbance$Year))
  data$Annual <- as.factor(as.character(data$Annual))
  disturbance <- merge(disturbance,unique(select(data, any_of(c("Annual", "Year", "PopulationName","PopulationID")))))
  anthro <- spread(subset(disturbance,select=c(Annual,PopulationID,Anthro)), PopulationID, Anthro)
  data <- merge(data,disturbance)
  
  if(any(is.na(anthro))){
    filter(anthro, if_any(-Annual, \(x)is.na(x))) %>% 
      select(Annual, where(\(x) any(is.na(x))))
    
   na_yrs <- anthro %>% summarise(across(-Annual, \(x) {
      paste0(Annual[which(is.na(x))], collapse = ", ")
    })) %>% 
     select(where(\(x) x != ""))
   
   warning("no data for population(s) ", paste0(names(na_yrs), collapse = ", "), " for years ", 
           paste0(unlist(na_yrs), collapse = ", "), 
           ". These years will be removed from the data for all populations")
   
   anthro <- na.omit(anthro)
   data <- subset(data,is.element(Annual, anthro$Annual))
  }
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
                anthro = anthro,
                nMonths = length(unique(data$Month))
  )

  surv_priors <- priors[c("S_b0_mu","S_b0_sd","S_b1_mu","S_b1_sd","S_cv_min","S_cv_max")] 
  names(surv_priors)<- gsub("S_","",names(surv_priors),fixed=T)

  datal <- c(datal,surv_priors)

  ###### Survival Model ######
  surv_mod_fl <- tempfile(pattern = "survival_model_", fileext = ".txt")
  sink(surv_mod_fl)  # assign model file name #
  if(surv_priors$cv_max>0){
    cat("
model {
  for(i in 1:nObs) {
    eSurvival[i] <- AnnualSurvival[Annual[i],PopulationID[i]]
    Mortalities[i] ~ dbin(1 - eSurvival[i], StartTotal[i])
  }

  for (i in 1:nAnnual) {
    for (k in 1:nPops) {
      mu.S[i,k] <- max(0.01,min(0.99,(46*exp(b0[k] + b1[k]*anthro[i,k])-0.5)/45))
      sig.S[i,k] <- min(cv.S[k]*mu.S[i,k],0.99*pow(mu.S[i,k]*(1-mu.S[i,k]),0.5)) \n 
      alpha[i,k] <- ((1-mu.S[i,k])/pow(sig.S[i,k],2) - 1/mu.S[i,k]) * pow(mu.S[i,k],2)
      beta[i,k] <- alpha[i,k] * (1/mu.S[i,k] - 1)
      Sbar[i,k] <- alpha[i,k]/(alpha[i,k]+beta[i,k])
      Survival[i,k]~ dbeta(alpha[i,k],beta[i,k])T(0.01,0.99)
      AnnualSurvival[i,k] <- pow(Survival[i,k],1/nMonths)
    }
  }

  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
    cv.S[k]~dunif(cv_min,cv_max)
  }

}
", fill = TRUE)
  }else{
    datal$cv_min=NULL;datal$cv_max=NULL
    cat("
model {
  for(i in 1:nObs) {
    eSurvival[i] <- AnnualSurvival[Annual[i],PopulationID[i]]
    Mortalities[i] ~ dbin(1 - eSurvival[i], StartTotal[i])
  }

  for (i in 1:nAnnual) {
    for (k in 1:nPops) {
      mu.S[i,k] <- max(0.01,min(0.99,(46*exp(b0[k] + b1[k]*anthro[i,k])-0.5)/45))
      Sbar[i,k] <- mu.S[i,k]
      sig.S <- 0
      Survival[i,k] <- mu.S[i,k]
      AnnualSurvival[i,k] <- pow(Survival[i,k],1/nMonths)
    }
  }

  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
  }

}
", fill = TRUE)
}
    sink()
  
  # Setting parameters - setting parameters that you want to monitor
  params = c("Sbar","Survival")
  
  # Setting initial values - not assigned
  inits1 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits2 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits3 <- list(b0 = rnorm(datal$nPops, 3, 2),b1 = rnorm(datal$nPops, 0, 2)) 
  inits <- list(inits1,inits2,inits3)
  
  return(jagsRunAndSummarize(data,datal,params,fname=surv_mod_fl,inits=inits,nc=nc,ni=ni,nb=nb,nt=nt))

}

betaRecruitment <- function(rec_fit, disturbance,priors,nc,nt,ni,nb){
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
  
  if(any(is.na(anthro))){
    filter(anthro, if_any(-Annual, \(x)is.na(x))) %>% 
      select(Annual, where(\(x) any(is.na(x))))
    
    na_yrs <- anthro %>% summarise(across(-Annual, \(x) {
      paste0(Annual[which(is.na(x))], collapse = ", ")
    })) %>% 
      select(where(\(x) x != ""))
    
    warning("no data for population(s) ", paste0(names(na_yrs), collapse = ", "), " for years ", 
            paste0(unlist(na_yrs), collapse = ", "), 
            ". These years will be removed from the data for all populations")
    
    anthro <- na.omit(anthro)
    data <- subset(data,is.element(Annual, anthro$Annual))
  }
  
  if(any(is.na(fire))){
    filter(fire, if_any(-Annual, \(x)is.na(x))) %>% 
      select(Annual, where(\(x) any(is.na(x))))
    
    na_yrs <- fire %>% summarise(across(-Annual, \(x) {
      paste0(Annual[which(is.na(x))], collapse = ", ")
    })) %>% 
      select(where(\(x) x != ""))
    
    warning("no data for population(s) ", paste0(names(na_yrs), collapse = ", "), " for years ", 
            paste0(unlist(na_yrs), collapse = ", "), 
            ". These years will be removed from the data for all populations")
    
    fire <- na.omit(fire)
    data <- subset(data,is.element(Annual, fire$Annual))
  }
  
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

  rec_priors = priors[c("R_b0_mu","R_b0_sd","R_b1_mu","R_b1_sd","R_b2_mu","R_b2_sd","R_cv_min","R_cv_max")]
  names(rec_priors)<- gsub("R_","",names(rec_priors),fixed=T)
  
  datal <- c(datal,rec_priors)
  
  # Currently set "adult_female_proportion" to be constant - potential to consider to be stochastic (beta distribution)
  # adult_female_proportion_alpha=rec_priors[["adult_female_proportion_alpha"]]
  # adult_female_proportion_beta=rec_priors[["adult_female_proportion_beta"]]
  
  ###### Recruitment Model ######
  rec_mod_fl <- tempfile(pattern = "recruit_model_", fileext = ".txt")
  sink(rec_mod_fl)  # observation model can be connected by CaribouYear[i] 
  
  if(rec_priors$cv_max>0){
  modStr <- "

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
      mu.R[i,k] <- max(0.01,min(0.99,_Rinvlink_(b0[k] + b1[k]*anthro[i,k] + b2[k]*fire[i,k])))
      sig.R[i,k] <- min(cv.R[k]*mu.R[i,k],0.99*pow(mu.R[i,k]*(1-mu.R[i,k]),0.5)) # Constrain on sig.R to fall within theoretical range
      alpha[i,k] <- ((1-mu.R[i,k])/pow(sig.R[i,k],2) - 1/mu.R[i,k]) * pow(mu.R[i,k],2)
      beta[i,k] <- alpha[i,k] * (1/mu.R[i,k] - 1)
      Rbar[i,k] <- alpha[i,k]/(alpha[i,k]+beta[i,k])
      Recruitment[i,k] ~ dbeta(alpha[i,k],beta[i,k])T(0.01,0.99)
    }
  }

  #adult_female_proportion ~ dbeta(adult_female_proportion_alpha, adult_female_proportion_beta)
  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
    b2[k] ~ dnorm(b2_mu, 1/pow(b2_sd,2))
    cv.R[k]~dunif(cv_min,cv_max)
  }
}
  
"
  }else{
    datal$cv_min=NULL;datal$cv_max=NULL
    modStr <- "

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
      mu.R[i,k] <- max(0.01,min(0.99,_Rinvlink_(b0[k] + b1[k]*anthro[i,k] + b2[k]*fire[i,k])))
      sig.R[i,k] <- 0
      Rbar[i,k] <- mu.R[i,k]
      Recruitment[i,k] <- mu.R[i,k]
    }
  }

  #adult_female_proportion ~ dbeta(adult_female_proportion_alpha, adult_female_proportion_beta)
  for (k in 1:nPops) {
    b0[k] ~ dnorm(b0_mu, 1/pow(b0_sd,2))
    b1[k] ~ dnorm(b1_mu, 1/pow(b1_sd,2))
    b2[k] ~ dnorm(b2_mu, 1/pow(b2_sd,2))
  }
}
"
}
  cat(gsub("_Rinvlink_",priors$R_inv_link,modStr,fixed=T), fill = TRUE)
  sink()
  
  ######## Define data, parameters, initials and settings #####
  # Setting parameters - setting parameters that you want to monitor
  params = c("Rbar","Recruitment")

  # Setting initial values
  inits1 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits2 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits3 <- list(b0 = rnorm(datal$nPops,-1, 2)) 
  inits <- list(inits1,inits2,inits3)
  
  return(jagsRunAndSummarize(data,datal,params,fname=rec_mod_fl,inits=inits,nc=nc,ni=ni,nb=nb,nt=nt))
  
}

