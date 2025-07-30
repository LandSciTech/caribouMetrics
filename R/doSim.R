#' Do multiple caribouPopGrowth runs
#'
#' @inheritParams caribouPopGrowth
#' @param addl_params a list of additional parameters for `caribouPopGrowth`
#'
#' @return numeric
#' @family demography
#'
#' @export
doSim <- function(numSteps, numPops, N0, R_bar, S_bar, R_sd, S_sd, R_iv_mean,R_iv_shape,S_iv_mean,S_iv_shape,
                  scn_nm, type="logistic", addl_params){
  #type="logistic"

  if(type=="beta"){
    varSample <- do.call(caribouPopGrowth,
                          c(list(N0 = rep(NA, numPops),
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
     R_b0 = rnorm(numPops,logit(R_bar),R_sd)
     S_b0 = rnorm(numPops,logit(S_bar),S_sd)
     varSample= list(R_t = inv.logit(R_b0),
                     S_t = inv.logit(S_b0))
     interannualVar = list(R_annual=rgamma(numPops,R_iv_shape,R_iv_shape/R_iv_mean),S_annual=rgamma(numPops,S_iv_shape,S_iv_shape/S_iv_mean))
     
  }

  mod_samps <- do.call(caribouPopSim, c(
    list(
      N0 = N0, numSteps = numSteps, R_samp = varSample$R_t, S_samp = varSample$S_t,
      interannualVar = interannualVar,
      # using simplified model version but keeping discrete animals
      l_S = 0, h_R = 1),
    addl_params)) %>%
    mutate(type = "samp", scn = scn_nm)


  mod_mean <- do.call(caribouPopSim, c(
    list(
      N0 = N0, numSteps = numSteps, R_samp = R_bar, S_samp = S_bar,
      # using simplified model version with no stochasticity
      interannualVar = NA, probOption = "continuous", l_S = 0, h_R = 1),
    addl_params)) %>%
    mutate(type = "mean", scn = scn_nm)

  return(bind_rows(mod_mean, mod_samps))
}

caribouPopSim <- function(N0, numSteps, R_samp, S_samp, interannualVar, ...) {
  for (ts in 1:numSteps) {
    if (ts == 1) {
      out <- caribouPopGrowth(rep(N0, length(R_samp)),
                              numSteps = 1,
                              interannualVar = interannualVar,
                              R_bar = R_samp, S_bar = S_samp, ...
      )
      out$id <- seq(1:nrow(out))
      out$time <- ts
      outBit <- out
    } else {
      outBit <- caribouPopGrowth(outBit$N,
                                 numSteps = 1, interannualVar = interannualVar,
                                 R_bar = R_samp, S_bar = S_samp, ...
      )
      outBit$id <- seq(1, nrow(outBit))
      outBit$time <- ts
      out <- rbind(out, outBit)
    }
  }
  return(out)
}

