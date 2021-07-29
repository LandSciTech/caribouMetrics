#' Take a sample from a beta distribution
#'
#' @param x
#' @param phi
#' @param useQuantiles
#' 

betaSample<-function(x,phi,useQuantiles=F){
  #x=predictedTableSD[1,]
  bShapes = betaGetShapes(x,phi) 
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rbeta(length(x),bShapes$shape1,bShapes$shape2))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    
    qq = qbeta(q,bShapes$shape1,bShapes$shape2)
    return(qq)
  }
}

#' Take a sample from a normal distribution
#'
#' @param x
#' @param sd
#' @param useQuantiles
#' 

normalSample<-function(x,sd,useQuantiles=F){
  #x=predictedTableSD[1,]
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rnorm(length(x),sd))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    qq = qnorm(q,mean=x,sd=sd)
    return(qq)
  }
}

#' Take a sample from a lognormal distribution
#'
#' @param x
#' @param sd
#' @param useQuantiles
#' 

lnormSample<-function(x,sd,useQuantiles=F){
  #x=predictedTableSD[1,]
  
  if((length(useQuantiles)==1)&&!useQuantiles){
    return(rnorm(length(x),sd))
  }else{
    if(length(useQuantiles)!=length(x)){
      q=getQuantiles(x)
    }else{
      q=useQuantiles
    }
    
    qq = qlnorm(q,sd)
    return(qq)
  }
}
