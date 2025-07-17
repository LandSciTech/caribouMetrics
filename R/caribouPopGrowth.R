#' Caribou demographic model
#'
#' A two-stage demographic model with density dependence and interannual
#' variability following [Johnson et. al. (2020)](doi:10.1111/1365-2664.13637)
#' with modifications described in 
#' [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095) and 
#' [Dyson et al. (2022)](https://doi.org/10.1101/2022.06.01.494350).
#' Demographic rates vary with disturbance as estimated by [Johnson et. al.
#' (2020)](doi:10.1111/1365-2664.13637). Default parameter values give the model
#' in Dyson et al. (2022). Set `probOption = "matchJohnson2020"` to reproduce
#' the model used in Johnson et al. 2020. Set `probOption = "continuous"`,
#' `interannualVar = FALSE`, and `K = FALSE` to reproduce the simpler 2-stage
#' demographic model without interannual variability, density dependence, or
#' discrete numbers of animals used by [Stewart et al.
#' (2023)](https://doi.org/10.1002/eap.2816). 
#' 
#' 
#' If R_annual and S_annual are provided, interannual variation in survival and
#' recruitment is modelled as in a logistic glmm with random effect of year.
#' 
#' See `vignette("caribouDemography")` 
#' and [Hughes et al. (2025)](https://doi.org/10.1016/j.ecoinf.2025.103095) for 
#' additional details and examples.
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
#'   (R_CV, S_CV), beta precision parameters (R_phi, S_phi), 
#'   or random effects parameters from a logistic glmm (R_annual, S_annual). 
#'   Set to `FALSE` to ignore interannual variability.
#' @param probOption Character. Choices are "binomial","continuous" or
#'   "matchJohnson2020". See description for details.
#' @param progress Logical. Should progress updates be shown?
#'
#' @return A data.frame of population size (`N`), expected growth rate
#'   (`lambda`), true growth rate (`lambdaTrue`), apparent annual reproduction rate (`R_t`), adjusted reproduction (`X_t`),
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
#'   Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
#'   Integration of national demographic-disturbance relationships and local
#'   data can improve caribou population viability projections and inform
#'   monitoring decisions. Ecological Informatics, 87, p.103095.
#'   <https://doi.org/10.1016/j.ecoinf.2025.103095>
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
                             K = 10000,
                             r_max = 1.3,
                             s=0.5,
                             l_R=0,
                             h_R=0.82,
                             l_S=0.61,
                             h_S=1,
                             c=1,
                             interannualVar = list(R_CV=0.46,S_CV=0.08696),
                             probOption="binomial",
                             progress = interactive()){
  if(is.character(interannualVar)){
    interannualVar = eval(parse(text=interannualVar))
  }
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

  lambdaE = S_bar*(1+c*R_bar*s) # To do - integrate adjDDRtProportion into expected mean
  
  h_R = s*h_R
  l_R = s*l_R

  if(length(N0) != length(R_bar) && length(R_bar) > 1){
    stop("R_bar must have length = 1 or the same length as N0", call. = FALSE)
  }

  if(length(N0) != length(S_bar) && length(S_bar) > 1){
    stop("S_bar  must have length = 1 or the same length as N0", call. = FALSE)
  }

  if(is.element("R_CV",names(interannualVar))){
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
    if(!is.element("R_CV",names(interannualVar))){
      R_t=s*R_t
    }
    
    
    #adjusting for composition survey bias
    R_tadj=c*R_t

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
      if(max(R_tadj*adjDDRtProportion) > 1){
        warning("Adjusted recruitment greater than 1.")
        R_tadj[R_tadj>1]=1
      }
      n_recruits <- rbinom(length(N),surviving_adFemales,R_tadj*adjDDRtProportion)
    }else{
      n_recruits <- round(n_recruitsUnadjDD * adjDDRtProportion,roundDigits)
    }
    N <- surviving_adFemales + n_recruits
    if(sum(is.na(N))>0){stop()}

    if(probOption=="matchJohnson2020"){
      ad = adjustN(N,Ntm1,r_max,denominatorAdjust=1e-06,roundDigits=roundDigits)
    }else{
      ad = adjustN(N,Ntm1,r_max,roundDigits=roundDigits)
    }
    if(sum(is.na(ad$N))>0){stop()}

    N=ad$N
    rr[paste0("lam",t)]= ad$Lambda
  }
    
  lamBits = names(rr)[grepl("lam",names(rr))]
  rr$lambdaTrue=matrixStats::rowProds(as.matrix(subset(rr,select=lamBits)),na.rm=T)^(1/length(lamBits)) #geometric mean
  rr=subset(rr,select=setdiff(names(rr),lamBits))
  rr$lambda = lambdaE
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
