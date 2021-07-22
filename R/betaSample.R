#' Take a sample from a beta distribution
#'
#' @param x
#' @param phi
#' @param useQuantiles
#' 
#' 
#' @description ...
#' 
#' @return ...

betaSample <- function(x,phi,useQuantiles=F){
  bShapes = betaGetShapes(x,phi) 
  
  if (!useQuantiles) {
    return(rbeta(length(x), bShapes$shape1, bShapes$shape2))
  }
  else {
    q = 0.025 + (seq(0, length(x)) / length(x)) * 0.95
    
    qq = qbeta(q, bShapes$shape1, bShapes$shape2)
    return(qq)
  }
}