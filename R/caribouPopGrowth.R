#' National caribou demographic model
#'
#' A two-stage demographic model with density dependence and interannual variability
#' following [Johnson et. al. (2020)](doi:10.1111/1365-2664.13637) with modifications described 
#' in [Dyson et al. (2022)](https://doi.org/10.1101/2022.06.01.494350). Demographic rates vary 
#' with disturbance as estimated by [Johnson et. al. (2020)](doi:10.1111/1365-2664.13637).
#' Default parameter values give the model in Dyson et al. (2022). Set `probOption =
#' "matchJohnson2020"` to reproduce the model used in Johnson et al. 2020. Set
#' `probOption = "continuous"`, `interannualVar = FALSE`, and `K = FALSE` to
#' reproduce the simpler 2-stage demographic model without interannual
#' variability, density dependence, or discrete numbers of animals used by
#' [Stewart et al. (2023)](https://doi.org/10.1002/eap.2816). See `vignette("caribouDemography")` for additional
#' details and examples.
#'
#' Given a population of post-juvenile females at the beginning of year \eqn{t}, \eqn{\dot{N}_t},
#' the number of post-juvenile females that survive from year \eqn{t} to the
#' census \eqn{\dot{W}_t} is binomially distributed with survival probability
#' \eqn{\dot{S}_t}: \eqn{\dot{W}_{t} \sim \text{Binomial}(\dot{N}_t,\dot{S}_t)}.
#' Maximum potential recruitment rate is adjusted for sex ratio, misidentification biases, and (optionally)
#' delayed age at first reproduction
#' \deqn{\dot{X}_t=\frac{\dot{c}\dot{R}_t/2}{1+\dot{c}\dot{R}_t/2}.} Realized recruitment rate
#' varies with population density, and the number of juveniles recruiting to the
#' post-juvenile class at the census is a binomially distributed function of the
#' number of surviving post-juvenile females and the adjusted recruitment rate:
#' \deqn{\dot{J}_{t} \sim
#' \text{Binomial}(\dot{W}_t,\dot{X}_t[p_0-(p_0-p_k)(\frac{\dot{W}_t}{N_0k})^b]\frac{\dot{W}_t}{\dot{W}_t+a}).}
#' Given default parameters, recruitment rate is lowest \eqn{(0.5\dot{X}_t)}
#' when \eqn{\dot{N}_t=1}, approaches a maximum of \eqn{\dot{X}_t} at
#' intermediate population sizes, and declines to \eqn{0.6\dot{X}_t} as the
#' population reaches carrying capacity of \eqn{K=50000}. The post-juvenile female population in the next year
#' includes both survivors and new recruits:
#' \eqn{\dot{N}_{t+1}=\text{min}(\dot{W}_t+\dot{J}_t,r_{max}\dot{N}_t)}.
#'
#' Interannual variation in survival and recruitment is modelled using truncated
#' beta distributions: \eqn{\dot{R}_t
#' \sim \text{TruncatedBeta}(\bar{R}_t,\nu_R,l_R,h_R); \dot{S}_t \sim
#' \text{TruncatedBeta}(\bar{S}_t,\nu_S,l_S,h_S)}. \eqn{(\nu_R,\nu_S)} are coefficients of variation
#' among years and \eqn{l_R,h_R,l_S,h_S} are maximum/minimum values for recruitment and survival.
#'
#' @param N0 Number or vector of numbers. Initial population size for one or
#'   more sample populations.
#' @param numSteps Number. Number of years to project.
#' @param R_bar Number or vector of numbers. Expected recruitment rate (calf:cow
#'   ratio) for one or more sample populations.
#' @param S_bar Number or vector of numbers. Expected adult female survival for
#'   one or more sample populations.
#' @param P_0 Number. Maximum recruitment multiplier.
#' @param P_K Number. Recruitment multiplier at carrying capacity.
#' @param a Number. Density dependence shape parameter.
#' @param b Number. Allee effect parameter.
#' @param K Number. Carrying capacity.
#' @param r_max Number. Maximum population growth rate.
#' @param s Number. Sex ratio.
#' @param l_R Number. Minimum recruitment.
#' @param h_R Number. Maximum recruitment.
#' @param l_S Number. Minimum survival.
#' @param h_S Number. Maximum survival.
#' @param c Number. Bias correction term.
#' @param interannualVar list or logical. List containing interannual
#'   variability parameters. These can be either coefficients of variation
#'   (R_CV, S_CV) or beta precision parameters (R_phi, S_phi). Set to `FALSE` to
#'   ignore interannual variability.
#' @param probOption Character. Choices are "binomial","continuous" or
#'   "matchJohnson2020". See description for details.
#' @param adjustR Logical. Adjust R to account for delayed age at first
#'   reproduction (DeCesare et al. 2012; Eacker et al. 2019). 
#' @param progress Logical. Should progress updates be shown?
#'
#' @return A data.frame of population size (`N`), average growth rate
#'   (`lambda`), apparent annual reproduction rate (`R_t`), adjusted reproduction (`X_t`),
#'   survival (`S_t`), number of recruits (`n_recruits`), and surviving females (`surviving_adFemales`)
#'   for each sample population projected for numSteps years.
#'
#' @references 
#'   Dyson, M., Endicott, S., Simpkins, C., Turner, J. W., Avery-Gomm, S.,
#'   Johnson, C. A., Leblond, M., Neilson, E. W., Rempel, R., Wiebe, P. A.,
#'   Baltzer, J. L., Stewart, F. E. C., & Hughes, J. (2022). Existing
#'   caribou habitat and demographic models need improvement for Ring of Fire
#'   impact assessment: A roadmap for improving the usefulness, transparency,
#'   and availability of models for conservation.
#'   <https://doi.org/10.1101/2022.06.01.494350>
#'   
#'   Johnson, C.A., Sutherland, G.D., Neave, E., Leblond, M., Kirby,
#'   P., Superbie, C. and McLoughlin, P.D., 2020. Science to inform policy:
#'   linking population dynamics to habitat for a threatened species in Canada.
#'   Journal of Applied Ecology, 57(7), pp.1314-1327.
#'   <https://doi.org/10.1111/1365-2664.13637>
#'
#'   Stewart, F.E., Micheletti, T., Cumming, S.G., Barros, C., Chubaty, A.M.,
#'   Dookie, A.L., Duclos, I., Eddy, I., Haché, S., Hodson, J. and Hughes, J.,
#'   2023. Climate‐informed forecasts reveal dramatic local habitat shifts and
#'   population uncertainty for northern boreal caribou. Ecological
#'   Applications, 33(3), p.e2816. <https://doi.org/10.1002/eap.2816>
#' @examples
#' caribouPopGrowth(100, 2, 0.5, 0.7)
#'
#' @family demography
#' @export
caribouPopGrowth <- function(N0,
                             numSteps,
                             R_bar,
                             S_bar,
                             P_0 = 1,
                             P_K = 0.6,
                             a = 1,
                             b = 4,
                             K = 50000,
                             r_max = 1.3,
                             s=0.5,
                             l_R=0,
                             h_R=0.82,
                             l_S=0.61,
                             h_S=1,
                             c=1,
                             interannualVar = list(R_CV=0.46,S_CV=0.08696),
                             probOption="binomial",
                             adjustR=FALSE,
                             progress = interactive()){
  rr=data.frame(N0=N0)
  
  N <- N0

  R_bar[R_bar<0]=0.000001
  S_bar[S_bar<0]=0.000001

  #Warn if S_bar outside of range l_S,h_S, or R_bar outside of range
  #l_R,h_R.
  if(any(S_bar < l_S) || any(S_bar > h_S)){
    warning("Setting expected survival S_bar to be between l_S and h_S.",
            call. = FALSE)
    S_bar = pmax(S_bar,l_S);S_bar=pmin(S_bar,h_S)
  }

  if(any(R_bar < l_R) || any(R_bar > h_R)){
    warning("Setting expected recruitment R_bar to be between l_R and h_R.",
            call. = FALSE)
    R_bar = pmax(R_bar,l_R);R_bar=pmin(R_bar,h_R)
  }

  h_R = s*h_R
  l_R = s*l_R

  if(length(N0) != length(R_bar) && length(R_bar) > 1){
    stop("R_bar must have length = 1 or the same length as N0", call. = FALSE)
  }

  if(length(N0) != length(S_bar) && length(S_bar) > 1){
    stop("S_bar  must have length = 1 or the same length as N0", call. = FALSE)
  }

  if(!is.element("R_phi",names(interannualVar))){
    R_bar=s*R_bar
  }else{
    #Phi is precision of calf cow ratio, not recruitment.
    l_R=l_R/s
    h_R=h_R/s
  }

  if(probOption=="matchJohnson2020"){
    roundDigits=0
    doBinomial=F
    a = a+1e-6
  }else{
    if(probOption=="continuous"){
      roundDigits=100
      doBinomial=F
    }else{
      roundDigits=0
      doBinomial=T
    }
    rK <- K #* N0
  }

  if(a<0){
    stop("a should be greater than or equal to 0")
  }

  for(t in 1:numSteps){
    if(progress){
      message(paste("projecting step ",t))
    }
    if(is.null(interannualVar)||any(is.na(interannualVar))||((length(interannualVar)==1)&&!interannualVar)){
      R_t= R_bar
      S_t = S_bar
    }else{
      if(length(N0)>length(R_bar)){R_bar=rep(R_bar,length(N0))}
      if(length(N0)>length(S_bar)){S_bar=rep(S_bar,length(N0))}
      R_t = addInterannualVar(R_bar,interannualVar,type="R",minV =l_R,maxV=h_R)
      S_t = addInterannualVar(S_bar,interannualVar,type="S",minV =l_S,maxV=h_S)
    }
    if(is.element("R_phi",names(interannualVar))){
      R_t=s*R_t
    }
    
    #adjusting for bias and delayed reproduction
    if(adjustR){
      R_tadj=c*R_t/(1+c*R_t)
    }else{R_tadj=c*R_t}

    Ntm1=N

    if(doBinomial){
      n_deaths <- rbinom(length(N),N,(1 - S_t))
    }else{
      n_deaths <- round(N * (1 - S_t),roundDigits)
    }

    surviving_adFemales <- N - n_deaths

    if(probOption=="matchJohnson2020"){
      rK <- K * N
    }

    n_recruitsUnadjDD <- surviving_adFemales * R_tadj

    if(K){
      adjDDRtProportion <- (P_0 -
                              ((P_0 - P_K) *
                                 (surviving_adFemales/rK)^b)) *
        surviving_adFemales/(surviving_adFemales + a)

      adjDDRtProportion[adjDDRtProportion<0] <- 0
      adjDDRtProportion[adjDDRtProportion>1] <- 1
    }else{
      adjDDRtProportion=1
    }
    if(doBinomial){
      n_recruits <- rbinom(length(N),surviving_adFemales,R_tadj*adjDDRtProportion)
    }else{
      n_recruits <- round(n_recruitsUnadjDD * adjDDRtProportion,roundDigits)
    }
    N <- surviving_adFemales + n_recruits

    if(probOption=="matchJohnson2020"){
      ad = adjustN(N,Ntm1,r_max,denominatorAdjust=1e-06,roundDigits=roundDigits)
    }else{
      ad = adjustN(N,Ntm1,r_max,roundDigits=roundDigits)
    }
    if(sum(is.na(ad$N))>0){stop()}

    N=ad$N
    rr[paste0("lam",t)]=ad$Lambda
  }

  lamBits = names(rr)[grepl("lam",names(rr))]
  rr$lambda=matrixStats::rowMeans2(as.matrix(subset(rr,select=lamBits)),na.rm=T)
  rr=subset(rr,select=setdiff(names(rr),lamBits))
  rr$N=N
  rr$R_t=R_t/s #apparent reproduction
  rr$X_t=R_tadj
  rr$S_t=S_t
  rr$n_recruits = n_recruits
  rr$surviving_adFemales = surviving_adFemales
  return(rr)
}

adjustN<-function(N,Ntm1,r_max=NA,denominatorAdjust=0,roundDigits=50){
  N[N<0]=0
  N[Ntm1==0]=0
  Lambda = N/(Ntm1+ denominatorAdjust)

  if(!is.na(r_max)){
    Lambda[Ntm1==0]=0
    N[Lambda>r_max]=round(Ntm1[Lambda>r_max] * r_max,roundDigits)
    Lambda = N/(Ntm1+ denominatorAdjust)
  }

  if(denominatorAdjust==0){
    Lambda[Ntm1==0]=0
  }
  return(list(N=N,Lambda=Lambda))
}
