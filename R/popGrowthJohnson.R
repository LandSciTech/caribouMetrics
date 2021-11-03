#' Caribou demographic model
#'
#' Given default parameter values, this is an implementation of the 2-stage
#' population growth model described in Dyson et al. (in Prep). Set
#' \code{probOption = "matchJohnson2020"} to reproduce the model used in Johnson
#' et al. 2020. Set \code{probOption = "continuous"}, \code{interannualVar =
#' FALSE}, and \code{K = FALSE} to reproduce the simpler 2-stage demographic
#' model without interannual variability, density dependence, or discrete
#' numbers of animals used by Stewart et al. (in prep). See 
#' \code{vignette("caribouDemography")} for additional details and examples.
#'
#' @param N Number or vector of numbers. Initial population size for one or more
#'   sample populations.
#' @param numSteps Number. Number of years to project.
#' @param R_bar Number or vector of numbers. Expected recruitment rate (calf:cow
#'   ratio) for one or more sample populations.
#' @param S_bar Number or vector of numbers. Expected adult female survival for
#'   one or more sample populations.
#' @param P_0 Number. Maximum recruitment multiplier.
#' @param P_K Number. Recruitment multiplier at carrying capacity.
#' @param a Number. Density dependence shape parameter.
#' @param b Number. Allee effect parameter.
#' @param K Number. Carrying capacity multiplier.
#' @param r_max Number. Maximum population growth rate.
#' @param s Number. Sex ratio.
#' @param l_R Number. Minimum recruitment.
#' @param h_R Number. Maximum recruitment.
#' @param l_S Number. Minimum survival.
#' @param h_S Number. Maximum survival.
#' @param interannualVar list or logical. List containing interannual
#'   variability parameters. These can be either coefficients of variation
#'   (R_CV, S_CV) or beta precision parameters (R_phi, S_phi). Set to
#'   \code{FALSE} ignore interannual variability.
#' @param probOption Character. Choices are "binomial","continuous" or
#'   "matchJohnson2020". See description for details.
#' @param progress Logical. Should progress updates be shown? 
#'
#' @return A data.frame of population size (N) and average growth rate (lambda)
#'   projections for each sample population.
#'
#'
#' @export
popGrowthJohnson <- function(N,
                             numSteps,
                             R_bar,
                             S_bar,
                             P_0 = 1,
                             P_K = 0.6,
                             a = 1,
                             b = 4,
                             K = 100,
                             r_max = 1.3,
                             s=0.5,
                             l_R=0,
                             h_R=0.82,
                             l_S=0.61,
                             h_S=1,
                             interannualVar = list(R_CV=0.46,S_CV=0.08696),
                             probOption="binomial",
                             progress = interactive()){
  rr=data.frame(N=N)
  R_bar[R_bar<0]=0
  S_bar[S_bar<0]=0
  
  #Return error if S_bar outside of range l_S,h_S, or R_bar outside of range
  #l_R,h_R.
  if(any(S_bar < l_S) || any(S_bar > h_S)){
    stop("Expected survival S_bar should be between l_R and h_R.",
         call. = FALSE)
  }
  
  if(any(R_bar < l_R) || any(R_bar > h_R)){
    stop("Expected recruitment R_bar should be between l_R and h_R.", 
         call. = FALSE)
  }
  
  h_R = s*h_R
  l_R = s*l_R
  
  if(length(N) != length(R_bar) && length(R_bar) > 1){
    stop("R_bar must have length = 1 or the same length as N", call. = FALSE)
  }
  
  if(length(N) != length(S_bar) && length(S_bar) > 1){
    stop("S_bar  must have length = 1 or the same length as N", call. = FALSE)
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
    rK <- K * N 
  }
  
  if(a<=0){
    stop("a should be greater than 0")
  }
  
  for(t in 1:numSteps){
    if(progress){
      print(paste("projecting step ",t))
    }
    if(is.null(interannualVar)||any(is.na(interannualVar))||((length(interannualVar)==1)&&!interannualVar)){
      R_t= R_bar
      S_t = S_bar
    }else{
      R_t = addInterannualVar(R_bar,interannualVar,type="R",minV =l_R,maxV=h_R)  
      S_t = addInterannualVar(S_bar,interannualVar,type="S",minV =l_S,maxV=h_S)      
    }
    if(is.element("R_phi",names(interannualVar))){
      R_t=s*R_t 
    }  
    
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

    n_recruitsUnadjDD <- surviving_adFemales * R_t 

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
      n_recruits <- rbinom(length(N),surviving_adFemales,R_t*adjDDRtProportion)
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
