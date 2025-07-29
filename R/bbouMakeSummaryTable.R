#' Create summary table of demographic rates from survival and recruitment surveys
#'
#' @param surv_data dataframe. Survival data in bboudata format
#' @param recruit_data dataframe. Recruitment data in bboudata format
#' @param N0 dataframe. Initial population estimates, required columns are
#'   PopulationName and N0
#' @param disturbance dataframe. Optional. If provided, fit a Beta model that includes disturbance covariates.
#' @param priors list. Optional. At present these are only used if disturbance is also provided.
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
#'   if FALSE just the results table is returned.
#' @export
#'
#' @examples
#' s_data <- rbind(bboudata::bbousurv_a, bboudata::bbousurv_b)
#' r_data <- rbind(bboudata::bbourecruit_a, bboudata::bbourecruit_b)
#' bbouMakeSummaryTable(s_data, r_data, 500, FALSE)

bbouMakeSummaryTable <-function(surv_data, recruit_data, N0, disturbance = NULL, priors = NULL, shiny_progress = FALSE,
                                return_mcmc=FALSE,i18n = NULL, niters = formals(bboutools::bb_fit_survival)$niters, nthin = formals(bboutools::bb_fit_survival)$nthin,...){
  #shiny_progress = FALSE;return_mcmc=FALSE;i18n = NULL
  
  if(length(N0)==1){
    N0= expand.grid(PopulationName=unique(surv_data$PopulationName),N0=N0)
  }
  if(is.null(i18n)){
    i18n <- list(t = function(x)paste0(x))
  }
  
  if(length(unique(surv_data$Year))<5){
    stop("At least 5 years of survival data are needed to estimate interannual variation using bboutools")
  }
  
  if(length(unique(recruit_data$Year))<5){
    stop("At least 5 years of survival data are needed to estimate interannual variation using bboutools")
  }

  # MCMC settings - (bboutools default: 1000 MCMC samples from 3 chains, number of )
  nc <- 3      # number of chains
  ni <- niters * nthin * 2   # number of samples for each chain
  nb <- ni / 2    # number of samples to discard as burnin
  
  if(!is.null(disturbance)){
    parTab = N0
    if(is.element("PopulationName",names(parTab))){
      names(parTab)[names(parTab)=="PopulationName"]= "pop_name"  
    }
    ret <- betaMakeSummaryTable(surv_data, recruit_data, disturbance, priors, nc,nthin,ni,nb) 
    ret$parTab <- parTab
    return(ret)
  }
  if(shiny_progress && !rlang::is_installed("shiny")){
    warning("Package shiny is not installed. Setting shiny_progress to FALSE")
    shiny_progress <- FALSE
  }
  if(shiny_progress) shiny::setProgress(0.2, message = i18n$t("Fitting survival"))
  
  surv_fit <- bboutools::bb_fit_survival(surv_data, multi_pops = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  
  if(shiny_progress) shiny::setProgress(0.4, message = i18n$t("Fitting recruitment"))
  
  recruit_fit <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, niters = niters, nthin = nthin, ...)
  
  if(shiny_progress) shiny::setProgress(0.6, message = i18n$t("Predicting survival"))
  surv_pred_bar <- bboutools::bb_predict_survival(surv_fit, year = FALSE, month = FALSE, conf_level = FALSE)
  
  if(shiny_progress) shiny::setProgress(0.8, message = i18n$t("Predicting recruitment"))
  rec_pred_bar <- bboutools::bb_predict_calf_cow_ratio(recruit_fit, year = FALSE, conf_level = FALSE)
  
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
        pop_name = p,
        R_bar = R_bar, R_sd = R_sd, R_iv_mean = R_annual_mean,R_iv_shape = R_annual_shape, R_bar_lower = R_qt[1], R_bar_upper = R_qt[2],
        S_bar = S_bar, S_sd = S_sd, S_iv_mean = S_annual_mean,S_iv_shape = S_annual_shape, S_bar_lower = S_qt[1], S_bar_upper = S_qt[2]
      )
    } else {
      parTab <- rbind(parTab, data.frame(
        pop_name = p,
        R_bar = R_bar, R_sd = R_sd, R_iv_mean = R_annual_mean,R_iv_shape = R_annual_shape, R_bar_lower = R_qt[1], R_bar_upper = R_qt[2],
        S_bar = S_bar, S_sd = S_sd, S_iv_mean = S_annual_mean,S_iv_shape = S_annual_shape, S_bar_lower = S_qt[1], S_bar_upper = S_qt[2]
      ))
    }
  }
  # dplyr version, not using but might want to some day...
  # R_samp %>% as_tibble(rownames = "id") %>%
  #   # move pop_name from column name to value
  #   pivot_longer(-id, names_to = "pop_name", values_to = "R") %>%
  #   # do the same to S and add it
  #   full_join(
  #     S_samp %>% as_tibble(rownames = "id") %>%
  #       pivot_longer(-id, names_to = "pop_name", values_to = "S"),
  #     by = c("id", "pop_name")
  #   ) %>%
  #   # calculate mean and sd of recruitment and survival for each pop
  #   group_by(pop_name) %>%
  #   summarise(across(c(R, S), .fns = list(bar = \(x)inv.logit(mean(logit(x))),
  #                                         sd = \(x)sd(logit(x))))) %>%
  #   # add annual columns
  #   mutate(N0 = N0, R_iv_par = R_Annual,
  #          S_iv_par = S_Annual) %>%
  #   # reorder columns
  #   select(pop_name, N0, matches("^R"), matches("^S"))
  
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
  
  if(is.element("PopulationName",names(N0))){
    parTab = merge(parTab,N0, by.x = "pop_name", by.y = "PopulationName")
  }else{
    parTab = merge(parTab,N0)
  }
  
  parTab = merge(parTab, data_amt, by.x = "pop_name", by.y = "PopulationName")
  
  if(return_mcmc){
    if(length(unique(surv_fit$data$Month))>1){
      newYr =  surv_fit$data$Year
      newYr[(surv_fit$data$Year==surv_fit$data$Annual)&(as.numeric(as.character(surv_fit$data$Month))<formals(bboutools::bb_fit_survival)$year_start)]=
        newYr[(surv_fit$data$Year==surv_fit$data$Annual)&(as.numeric(as.character(surv_fit$data$Month))<formals(bboutools::bb_fit_survival)$year_start)]+1
      surv_fit$data$Year = newYr
    }
    return(list(parTab=parTab,surv_fit=surv_fit,recruit_fit=recruit_fit))
  }else{
    return(parTab)
  }
  
}