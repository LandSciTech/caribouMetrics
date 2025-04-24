#' Calculate bias correction term for calf:cow composition survey.
#'
#' When composition surveys are conducted there is a possibility of bias in calf
#' cow ratios due to misidentifying young bulls as adult females and vice versa
#' or missing calves. Here we address this gap with a bias term derived from a simple
#' model of the recruitment survey observation process. See
#' [Hughes et al. (2025) Section 2.2](https://doi.org/10.1016/j.ecoinf.2025.103095) 
#' for a detailed description of the model.
#' 
#' 
#' @param w number. The apparent number of adult females per collared animal in
#'   composition survey.
#' @param q number in 0, 1. Ratio of bulls to cows in composition survey groups.
#' @param u number in 0, 1. Probability of misidentifying young bulls as adult
#'   females and vice versa in composition survey.
#' @param z number in 0, <1. Probability of missing calves in composition
#'   survey.
#' @param approx logical. If TRUE approximate the uncertainty about the value of
#'   the composition bias correction value (c) with the log-normal distribution
#'   of c given all the supplied values of `q`, `u`, and `z`. If FALSE the
#'   composition bias correction value (c) is returned for each value of `q`,
#'   `u`, and `z`
#'   
#' @return number or tibble. If `approx = FALSE` a vector of composition bias
#'   correction values (c) of the same length as `q`, `u`, and `z`. If `approx =
#'   TRUE` a tibble with on row per unique value of `w` and columns `w`, `m`,
#'   `v`, `sig2`, `mu` representing `w`, mean `c`, variance of `c`, and parameters for a
#'   log-normal approximation of the distribution of `c`.
#' @export
#' @family demography
#' 
#' 
#' @references    
#'   Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
#'   Integration of national demographic-disturbance relationships and local
#'   data can improve caribou population viability projections and inform
#'   monitoring decisions. Ecological Informatics, 87, p.103095.
#'   <https://doi.org/10.1016/j.ecoinf.2025.103095>
#'   
#' @examples
#' # number or reps
#' nr <- 10
#'
#' compositionBiasCorrection(w = 6,
#'                           q = runif(nr, 0, 0.6),
#'                           u = runif(nr, 0, 0.2),
#'                           z = runif(nr, 0, 0.2),
#'                           approx = FALSE)
#'
#' compositionBiasCorrection(w = 6,
#'                           q = runif(nr, 0, 0.6),
#'                           u = runif(nr, 0, 0.2),
#'                           z = runif(nr, 0, 0.2),
#'                           approx = TRUE)
#'
#' 

compositionBiasCorrection<-function(w,q,u,z,approx=F){
  #q=runif(nr,0,0.6);w=cr$w;u=runif(nr,0,0.2);z=runif(nr,0,0.2);approx=T
  w = pmax(w,1) 
  if(max(z)>=1){
    stop("Composition bias correct term is undefined when the probability of missing calves z >= 1.")
  }
  c = w*(q*u+1-u)/((q*u+w-u)*(1-z))
  
  if(approx){
    #get mean and sd as fn of w
    cr = data.frame(w=w,c=c)
    cs <- cr %>%
      group_by(w) %>%
      summarise(m = mean(c), v = var(c))

    #now find lognormal parameters from mean and variance
    #https://www.johndcook.com/blog/2022/02/24/find-log-normal-parameters/
    cs$sig2 = NA
    cs$sig2[(cs$v==0)|(cs$m==0)]=0
    cs$sig2[!((cs$v==0)|(cs$m==0))] = log(1+cs$v/cs$m^2)
    cs$mu = log(cs$m)-cs$sig2/2
    cs$mu[cs$m==0]=NA
    
    return(cs)
  }else{
    return(c)
  }
}