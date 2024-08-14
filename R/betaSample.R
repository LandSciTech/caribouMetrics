#' Take a sample from a beta distribution
#'
#'
#' @noRd
betaSample<-function(x,phi,quantilesToUse=NULL){
  #x=predictedTableSD[1,]
  bShapes = simstudy::betaGetShapes(x,phi) 
  
  if(is.null(quantilesToUse)){
    return(rbeta(length(x),bShapes$shape1,bShapes$shape2))
  }else{
    qq = qbeta(quantilesToUse,bShapes$shape1,bShapes$shape2)
    return(qq)
  }
}

#' Take a sample from a normal distribution
#'
#' 
#' @noRd
#' 
normalSample<-function(x,sd,quantilesToUse=NULL){
  #x=predictedTableSD[1,]
  
  if(is.null(quantilesToUse)){
    return(rnorm(length(x), mean = x, sd = sd))
  }else{
    qq = qnorm(quantilesToUse,mean=x,sd=sd)
    return(qq)
  }
}

#' Take a sample from a lognormal distribution
#'
#' 
#' @noRd
#' 
lnormSample<-function(x,sd,quantilesToUse=NULL){
  if(is.null(quantilesToUse)){
    return(rnorm(length(x),sd))
  }else{
    qq = qlnorm(quantilesToUse,sd)
    return(qq)
  }
}

estBetaParams <- function(mu, sigma){
  
  if(any(mu<0)){
    print("ERROR (estBetaParam): mu must not be less than 0. Returning NULL.")
    return()
  } 
  
  if(any(mu>1)){
    print("ERROR (estBetaParam): mu must not be greater than 1. Returning NULL.")
    return()
  } 
  
  if(any(sigma<0)){
    print("ERROR (estBetaParam): sigma must not be negative. Returning NULL.")
    return()
  } 
  
  alpha <- ((1-mu)/sigma^2 - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(list(alpha=alpha, beta=beta))
  
}
fillNAsWithMean <- function(vector) {
  vector[is.na(vector)] <- mean(vector, na.rm = TRUE)
  vector
}

addInterannualVar<-function(bar,interannualVar,type,minV,maxV){
  if(is.element(paste0(type,"_CV"),names(interannualVar))){
    #reproducing ECCC_CaribouPopnProjection - see line 143 etc of functions.R
    sigma <- (bar * interannualVar[[paste0(type,"_CV")]])
    sigma=fillNAsWithMean(sigma)
    BetaPars  <- estBetaParams(bar, sigma)
    BetaPars$alpha[BetaPars$alpha < 0] <- 0.01
    BetaPars$beta[BetaPars$beta < 0] <- 0.01
    interannualVar[[paste0(type,"_alpha")]]=BetaPars$alpha
    interannualVar[[paste0(type,"_beta")]]=BetaPars$beta
  }
  if(is.element(paste0(type,"_phi"),names(interannualVar))){
    bShapes = simstudy::betaGetShapes(bar,interannualVar[[paste0(type,"_phi")]]) 
    interannualVar[[paste0(type,"_alpha")]]=bShapes$shape1
    interannualVar[[paste0(type,"_beta")]]=bShapes$shape2
  }
  
  bar_t = withCallingHandlers(
    # using calling handler because of problem in truncdist package. Have flagged to developer
    rtrunc(length(bar), 
                 spec="beta", 
                 shape1=interannualVar[[paste0(type,"_alpha")]], 
                 shape2= interannualVar[[paste0(type,"_beta")]],
                 a=minV,
                 b=maxV),
    warning = function(cnd){
      if (startsWith(conditionMessage(cnd), "the condition has length > 1"))
        invokeRestart("muffleWarning")
    }
  )
  
  return(bar_t)
}

# Copyright 2023 Frederick Novomestky, Saralees Nadarajah & Her Majesty the Queen
# in Right of Canada as represented by the Minister of the Environment 
# License GPL-3 
# NOTICE: This function has been modified from the truncdist package
# because because of an issue with if statements with length > 1 in R 4.2
# I modified to add any()
qtrunc <- function (p, spec, a = -Inf, b = Inf, ...) {
  if (a >= b) 
    stop("argument a is greater than or equal to b")
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  G.a <- G(a, ...)
  G.b <- G(b, ...)
  if (any(G.a == G.b)) {
    stop("Trunction interval is not inside the domain of the quantile function")
  }
  result <- pmin(pmax(a, Gin(G(a, ...) + p * (G(b, ...) - G(a, 
                                                            ...)), ...)), b)
  return(result)
}

rtrunc <- function (n, spec, a = -Inf, b = Inf, ...){
  if (any(a >= b)) 
    stop("argument a is greater than or equal to b")
  u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b, ...)
  return(x)
}
