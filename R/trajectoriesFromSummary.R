#' Do multiple caribouPopGrowth runs
#'
#' @param addl_params a list of additional parameters for `caribouPopGrowth`
#' @inheritParams caribouPopGrowth 
#' @param replicates 
#' @param R_sd,S_sd standard deviation of R_bar and S_bar
#' @param R_iv_mean,R_iv_shape,S_iv_mean,S_iv_shape define the mean and shape of the interannual variation
#' @param scn_nm Sceanrio name
#' @param type "logistic" or "beta" defines how demographic rates are sampled from the given mean and standard deviation.
#' @param doSummary logical. Default TRUE. If FALSE returns unprocessed outcomes from caribouPopGrowth. 
#'  If TRUE returns summaries and (if returnSamples = T) sample trajectories from prepareTrajectories.
#' @param returnSamples logical. If FALSE returns only summaries. If TRUE
#'   returns example trajectories as well. 
#' @return a data.frame
#' @family demography
#'
#' @export
#' 
#' @examples
#'  outParTab <- trajectoriesFromSummary(
#'    numSteps = 5, replicates = 2, N0 = NA, R_bar = 0.18, S_bar = 0.87,
#'    R_sd = 0.085, S_sd = 0.16,
#'    R_iv_mean = 0.34, S_iv_mean = 0.31,
#'    R_iv_shape = 18, S_iv_shape = 3.3,
#'    scn_nm = "base", addl_params = NULL, type = "logistic"
#'  )
#'  outParTab
#'
trajectoriesFromSummary <- function(numSteps, replicates, N0, R_bar, S_bar, R_sd, S_sd,
                  R_iv_mean, R_iv_shape, S_iv_mean, S_iv_shape,  
                  scn_nm, type = "logistic", addl_params = list(), doSummary = F, returnSamples = T){
  #type="logistic"

  if(type=="beta"){
    varSample <- do.call(caribouPopGrowth,
                          c(list(N0 = rep(NA, replicates),
                                 numSteps = 1,
                                 R_bar = R_bar,
                                 S_bar = S_bar,
                                 interannualVar = list(
                                   R_CV = R_sd/R_bar,
                                   S_CV = S_sd/S_bar
                                 ),
                                 probOption = "continuous",
                                 l_S = 0, h_R = 1),
                            addl_params))
    interannualVar = list(R_CV = R_iv_mean, S_CV = S_iv_mean)
  }
  if(type=="logistic"){
     #convert mean
     R_b0 = rnorm(replicates,logit(R_bar),R_sd)
     S_b0 = rnorm(replicates,logit(S_bar),S_sd)
     varSample= list(R_t = inv.logit(R_b0),
                     S_t = inv.logit(S_b0))
     interannualVar = list(R_annual=rgamma(replicates,R_iv_shape,R_iv_shape/R_iv_mean),S_annual=rgamma(replicates,S_iv_shape,S_iv_shape/S_iv_mean))
     
  }

  mod_samps <- do.call(simPopsOverTime, c(
    list(
      N0 = N0, numSteps = numSteps, R_samp = varSample$R_t, S_samp = varSample$S_t,
      interannualVar = interannualVar,
      # using simplified model version but keeping discrete animals
      l_S = 0, h_R = 1),
    addl_params)) %>%
    mutate(type = "samp", scn = scn_nm)


  mod_mean <- do.call(simPopsOverTime, c(
    list(
      N0 = mean(N0), numSteps = numSteps, R_samp = R_bar, S_samp = S_bar,
      # using simplified model version with no stochasticity
      interannualVar = NA, probOption = "continuous", l_S = 0, h_R = 1),
    addl_params)) %>%
    mutate(type = "mean", scn = scn_nm)

  if(!doSummary){
    return(bind_rows(mod_mean, mod_samps))
  }else{
    mod_samps$Year <- mod_samps$time
    mod_samps$PopulationName <- mod_samps$scn
    simBig <- prepareTrajectories(mod_samps, returnSamples = returnSamples)
  }
}

#' Simulate caribou population over time
#'
#' If `dynamicRates = FALSE` then `R_samp` and `S_samp` are constant rates over
#' time and their length is the number of populations. If `dynamicRates = TRUE`
#' and `R_samp` and `S_samp` are vectors they represent the rate at each time
#' step and there length should be equal to `numSteps`. If If `dynamicRates =
#' TRUE` and `R_samp` and `S_samp` are matrices then rows represent populations
#' and columns represent timesteps, so cell `[i,j]` of the matrix is the rate
#' for population i at timestep j
#' 
#' `numSteps` is the number of timesteps in the simulation while `stepLength` is
#' the `numSteps` in each call to `caribouPopGrowth`
#'
#' @noRd
simPopsOverTime <- function(N0, numSteps, R_samp, S_samp, interannualVar, dynamicRates = FALSE, stepLength = 1, ...) {
  if(dynamicRates){
    if(!is.null(dim(R_samp))){
      stopifnot(ncol(R_samp) == numSteps)
      onePop <- FALSE
    }else {
      stopifnot(length(R_samp) == numSteps)
      onePop <- TRUE
    }
  }
  for (ts in 1:numSteps) {
    if(dynamicRates){
      if(onePop){
        R_use <- data.frame(value = R_samp[ts]) %>% setNames(ts)
        S_use <- data.frame(value = S_samp[ts]) %>% setNames(ts)
      } else {
        R_use <- R_samp[,ts, drop = FALSE]
        S_use <- S_samp[,ts, drop = FALSE]
      } 
    } else {
      R_use <- data.frame(value = R_samp)
      S_use <- data.frame(value = S_samp)
    }
    R_use <- na.omit(R_use)
    S_use <- na.omit(S_use)
    if (ts == 1) {
      if(length(N0) == 1){
        N0 <- rep(N0, nrow(R_use))
      } else if (length(N0) == 2) {
        N0 <- seq(from = N0[1], to = N0[2], by  = 1) %>% round() %>% 
          sample(size = nrow(R_use), replace = TRUE)
      }

      out <- caribouPopGrowth(N0,
                              numSteps = stepLength,
                              interannualVar = interannualVar,
                              R_bar = R_use[,1], S_bar = S_use[,1], ...
      )
      
      if(is.null(rownames(R_use))){
        out$id <- seq(1, nrow(out))
      } else {
        out$id <- rownames(R_use)
      }

      #TO DO: replace temporary fix with something else here.
      out$time <- ts#ifelse(!is.null(colnames(R_use)), colnames(R_use), ts)
      
      outBit <- out
    } else {
      if(length(outBit$N) != length(R_use)){
        if(length(unique(N0)) > 1){
          stop("Range of N0 only supported for static rates when there are NAs in R_samp")
        }
        # Add rows missing in previous rounds and set N to N0
        outBit <- left_join(R_use %>% tibble::rownames_to_column(),
                            outBit %>% select(-N0) %>% mutate(id = as.character(id)),
                            by = join_by(rowname == id)) %>% 
          mutate(N = ifelse(is.na(N), unique(N0), N)) %>% 
          select(any_of(colnames(outBit))) 
      }
      outBit <- caribouPopGrowth(outBit$N,
                                 numSteps = stepLength, interannualVar = interannualVar,
                                 R_bar = R_use[,1], S_bar = S_use[,1], ...
      )
      if(is.null(rownames(R_use))){
        outBit$id <- seq(1, nrow(outBit))
      } else {
        outBit$id <- rownames(R_use)
      }
   
      outBit$time <- ts#ifelse(!is.null(colnames(R_use)), colnames(R_use), ts)
      
      out <- rbind(out, outBit)
    }
  }
  return(out)
}

