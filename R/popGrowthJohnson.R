#' Implement the Johnson 2020 population model
#'
#' @param N
#' @param numSteps
#' @param Rec_bar
#' @param S_bar
#' @param P_0
#' @param P_K
#' @param alpha
#' @param beta
#' @param Kmultiplier
#' @param r_max
#' @param interannualVar
#' 
#' 
#' @description ...
#' 
#' @return ...

popGrowthJohnson <- function(N,
                             numSteps,
                             Rec_bar,
                             S_bar,
                             P_0 = 0.6,
                             P_K = 0.95,
                             alpha = 1,
                             beta = 4,
                             Kmultiplier = 100,
                             r_max = 1.3,
                             interannualVar = list(Rec = 1.62 + 3.44,
                                                   S = 13.98 + 2.51)){
  
  K = N * Kmultiplier  
  rr = data.frame(N = N)
  
  for (t in 1:numSteps) {
    print(paste("projecting step ", t))
    
    if (is.null(interannualVar)) {
      Rec_t = Rec_bar
      S_t = S_bar
    }
    else {
      Rec_phi = unique(precisionRec)
      if (length(Rec_phi) > 1) {
        stop("handle vector of recruitment precision parameters")
      }
      Rec_t = betaSample(Rec_bar,Rec_phi)
      
      S_phi = unique(precisionS)
      if (length(S_phi) > 1) {
        stop("handle vector of survival precision parameters")
      }
      S_t = betaSample(S_bar, S_phi)
      
    }
    
    adjustedRec = (P_0 * (1 - N/K) ^ beta + (P_K * N/K) ^ beta) * N / 
      (alpha + N)    
    
    N = N - N * (1 - S_t) + N * S_t * adjustedRec
    
  }
  
  rr$lambda = (N/rr$N) / numSteps
  rr$N = N
  return(rr)
}