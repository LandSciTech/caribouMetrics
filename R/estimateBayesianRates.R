#' Create summary table of demographic rates from survival and recruitment surveys
#'
#' @param surv_data dataframe. Survival data in bboudata format
#' @param recruit_data dataframe. Recruitment data in bboudata format
#' @param N0 dataframe. Optional. Initial population estimates, required columns are PopulationName and N0
#' @param disturbance dataframe. Optional. If provided, fit a Beta model that includes disturbance covariates.
#' @param priors list. Optional. If disturbance is NA, this should be list(priors_survival=c(...),priors_recruitment=c(...)); see `bboutools::bb_priors_survival` and `bboutools::bb_priors_recruitment` for details.
#'               If disturbance is not NA, see `betaNationalPriors` for details.
#' @param shiny_progress logical. Should shiny progress bar be updated. Only set
#'   to TRUE if using in an app.
#' @param return_mcmc boolean. If TRUE return fitted survival and recruitment
#'   models. Default FALSE.
#' @param niters integer. The number of iterations per chain after thinning and burn-in.
#' @param nthin integer. The number of the thinning rate.
#' @param ... Other parameters passed on to `bboutools::bb_fit_survival` and
#'   `bboutools::bb_fit_recruitment`.
#'
#' @return If `return_mcmc` is TRUE then a list with results and fitted models,
#'   if FALSE just the results summaries are returned.
#' @export
#' @family demography
#'
#' @examples
#' s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
#' r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)
#' estimateBayesianRates(s_data, r_data, N0 = 500)

estimateBayesianRates <-function(surv_data, recruit_data, N0=NA, disturbance = NULL, priors = NULL, shiny_progress = FALSE,
                                return_mcmc=FALSE,i18n = NULL, niters = formals(bboutools::bb_fit_survival)$niters, nthin = formals(bboutools::bb_fit_survival)$nthin,...){
  #shiny_progress = FALSE;return_mcmc=FALSE;i18n = NULL
  
  if(length(N0)==1){
    N0= expand.grid(PopulationName=unique(surv_data$PopulationName),N0=N0)
  }
  if(is.null(i18n)){
    i18n <- list(t = function(x)paste0(x))
  }

  # MCMC settings - (bboutools default: 1000 MCMC samples from 3 chains, number of )
  nc <- 3      # number of chains
  ni <- niters * nthin   # number of samples for each chain
  nb <- ni / 2    # number of samples to discard as burnin
  
  if(!is.null(disturbance)){
    if(is.null(priors)){
      priors=betaNationalPriors()
    }
    ret <- betaMakeSummaryTable(surv_data, recruit_data, disturbance, priors, nc,nthin,ni,nb) 
    
    ret$parList$N0 <- merge(N0,disturbance)
    
    if(nrow(unique(subset(disturbance,select=c(-Year))))==1){
      
      rbar <- unique(subset(ret$parList$Rbar,select=c(mean,sd,lower,upper,PopulationName,Anthro,fire_excl_anthro)))
      names(rbar)[1:4] <- c("R_bar","R_sd","R_bar_lower","R_bar_upper")
      sbar <- unique(subset(ret$parList$Sbar,select=c(mean,sd,lower,upper,PopulationName)))
      names(sbar)[1:4] <- c("S_bar","S_sd","S_bar_lower","S_bar_upper")
      ivs <- data.frame(R_iv_mean=(ret$parList$Riv$R_cv_min+ret$parList$Riv$R_cv_max)/2,R_iv_shape=NA,
                        S_iv_mean=(ret$parList$Siv$S_cv_min+ret$parList$Siv$S_cv_max)/2,S_iv_shape=NA)
      ret$parTab <- merge(merge(merge(rbar,sbar),N0),ivs)
    }else{
       ret$parTab <- merge(N0,disturbance)
    }
    return(ret)
  }
  
  if(length(unique(surv_data$Year))<5){
    stop("At least 5 years of survival data are needed to estimate interannual variation using bboutools")
  }
  
  if(length(unique(recruit_data$Year))<5){
    stop("At least 5 years of survival data are needed to estimate interannual variation using bboutools")
  }
  
  if(shiny_progress && !rlang::is_installed("shiny")){
    warning("Package shiny is not installed. Setting shiny_progress to FALSE")
    shiny_progress <- FALSE
  }
  if(shiny_progress) shiny::setProgress(0.2, message = i18n$t("Fitting survival"))

  if(is.element("priors_survival",names(priors))){
    surv_fit <- bboutools::bb_fit_survival(surv_data, priors=priors$priors_survival, multi_pops = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  }else{
    surv_fit <- bboutools::bb_fit_survival(surv_data, multi_pops = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  }
  
  if(shiny_progress) shiny::setProgress(0.4, message = i18n$t("Fitting recruitment"))
  if(is.element("priors_recruitment",names(priors))){
    recruit_fit <- bboutools::bb_fit_recruitment(recruit_data, priors=priors$priors_recruitment, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  }else{
    recruit_fit <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  }
  
  surv_pred_bar <- bboutools::bb_predict_survival(surv_fit, year = FALSE, month = FALSE, conf_level = FALSE)
  rec_pred_bar <- bboutools::bb_predict_calf_cow_ratio(recruit_fit,year = FALSE, conf_level = FALSE)
  
  # summarize model output
  data_sur <- surv_pred_bar$data
  data_rec <- rec_pred_bar$data
  
  # Force matrix b/c if only one it is numeric
  S_samp <- mcmcr::collapse_chains(surv_pred_bar$samples)[, , ] %>% as.matrix()
  R_samp <- mcmcr::collapse_chains(rec_pred_bar$samples)[, , ] %>% as.matrix()
  rownames(S_samp) <- seq(1, nrow(S_samp))
  colnames(S_samp) <- levels(data_sur$PopulationID)
  rownames(R_samp) <- seq(1:nrow(R_samp))
  colnames(R_samp) <- levels(data_rec$PopulationID)
  
  pops <- levels(data_sur$PopulationID)

  #characterize distribution of sAnnual
  x = mcmcr::collapse_chains(surv_fit$samples$sAnnual)[, , ]
  #descdist(x, discrete = FALSE)
  s_dist <- fitdistrplus::fitdist(x, "gamma")
  #plot(gamma_dist)
  S_annual_mean <- s_dist$estimate[1]/s_dist$estimate[2]
  S_annual_shape <- s_dist$estimate[1]
  
  x = mcmcr::collapse_chains(recruit_fit$samples$sAnnual)[, , ]
  #descdist(x, discrete = FALSE)
  r_dist <- fitdistrplus::fitdist(x, "gamma")
  #plot(r_dist)
  R_annual_mean <- r_dist$estimate[1]/r_dist$estimate[2]
  R_annual_shape <- r_dist$estimate[1]
  
  for (i in 1:length(pops)) {
    # i=2
    p <- pops[i]
    R_samp_long <- R_samp[, i]
    S_samp_long <- S_samp[, i]
    R_bar <- inv.logit(mean(logit(R_samp_long)))
    R_sd <- sd(logit(R_samp_long))
    R_qt <- quantile(R_samp_long, probs=c(0.025, 0.975))
    
    S_bar <- inv.logit(mean(logit(S_samp_long)))
    S_sd <- sd(logit(S_samp_long))
    S_qt <- quantile(S_samp_long, probs=c(0.025, 0.975))
    
    if (i == 1) {
      parTab <- data.frame(
        PopulationName = p,
        R_bar = R_bar, R_sd = R_sd, R_iv_mean = R_annual_mean,R_iv_shape = R_annual_shape, R_bar_lower = R_qt[1], R_bar_upper = R_qt[2],
        S_bar = S_bar, S_sd = S_sd, S_iv_mean = S_annual_mean,S_iv_shape = S_annual_shape, S_bar_lower = S_qt[1], S_bar_upper = S_qt[2]
      )
    } else {
      parTab <- rbind(parTab, data.frame(
        PopulationName = p,
        R_bar = R_bar, R_sd = R_sd, R_iv_mean = R_annual_mean,R_iv_shape = R_annual_shape, R_bar_lower = R_qt[1], R_bar_upper = R_qt[2],
        S_bar = S_bar, S_sd = S_sd, S_iv_mean = S_annual_mean,S_iv_shape = S_annual_shape, S_bar_lower = S_qt[1], S_bar_upper = S_qt[2]
      ))
    }
  }

  # dplyr version, not using but might want to some day...
  # R_samp %>% as_tibble(rownames = "id") %>%
  #   # move PopulationName from column name to value
  #   pivot_longer(-id, names_to = "PopulationName", values_to = "R") %>%
  #   # do the same to S and add it
  #   full_join(
  #     S_samp %>% as_tibble(rownames = "id") %>%
  #       pivot_longer(-id, names_to = "PopulationName", values_to = "S"),
  #     by = c("id", "PopulationName")
  #   ) %>%
  #   # calculate mean and sd of recruitment and survival for each pop
  #   group_by(PopulationName) %>%
  #   summarise(across(c(R, S), .fns = list(bar = \(x)inv.logit(mean(logit(x))),
  #                                         sd = \(x)sd(logit(x))))) %>%
  #   # add annual columns
  #   mutate(N0 = N0, R_iv_par = R_Annual,
  #          S_iv_par = S_Annual) %>%
  #   # reorder columns
  #   select(PopulationName, N0, matches("^R"), matches("^S"))
  
  # data amount summary
  surv_data_amt <- surv_data %>% group_by(PopulationName, Year) %>%
    summarise(nCollars = max(StartTotal)) %>%
    summarise(nCollarYears = sum(nCollars),
              nSurvYears = n_distinct(Year))
  
  recruit_data_amt <- recruit_data %>% group_by(PopulationName, Year) %>%
    summarise(nCows = sum(Cows)) %>%
    summarise(nCowsAllYears = sum(nCows),
              nRecruitYears = n_distinct(Year))
  
  data_amt <- merge(surv_data_amt, recruit_data_amt)
  parTab = merge(parTab,N0)
  parTab = merge(parTab, data_amt)

  parList = list()
  parList$Rbar <- subset(parTab,select=c(R_bar,R_sd,R_bar_lower,R_bar_upper,PopulationName))
  names(parList$Rbar)[1:4] = c("mean","sd","lower","upper")
  parList$Sbar <- subset(parTab,select=c(S_bar,S_sd,S_bar_lower,S_bar_upper,PopulationName))
  names(parList$Rbar)[1:4] = c("mean","sd","lower","upper")
  parList$Siv <- subset(parTab,select=c(S_iv_mean,S_iv_shape))
  parList$Riv <- subset(parTab,select=c(R_iv_mean,R_iv_shape))
  parList$type <- "bbou"
  
  if(return_mcmc){
    if(length(unique(surv_fit$data$Month))>1){
      newYr =  surv_fit$data$Year
      newYr[(surv_fit$data$Year==surv_fit$data$Annual)&(as.numeric(as.character(surv_fit$data$Month))<formals(bboutools::bb_fit_survival)$year_start)]=
        newYr[(surv_fit$data$Year==surv_fit$data$Annual)&(as.numeric(as.character(surv_fit$data$Month))<formals(bboutools::bb_fit_survival)$year_start)]+1
      surv_fit$data$Year = newYr
    }
    return(list(parTab=parTab,oarList=parList,surv_fit=surv_fit,recruit_fit=recruit_fit))
  }else{
    return(parTab)
  }
  
}
