#' Calculate bias correction term for calf:cow composition survey.
#' TO DO: documentation borrowing from manuscript.
#' 
#' @param w number. The apparent number of adult females per collared animal in composition survey.
#' @param q number in 0, 1. Ratio of bulls to cows in composition survey groups.
#' @param u number in 0, 1. Probability of misidentifying young bulls as adult females and vice versa in composition survey.
#' @param z number in 0, <1. Probability of missing calves in composition survey.
#' @return number. Composition bias correction value.
#' @export
#' @family demography
#' @examples
#' 
compositionBiasCorrection<-function(w,q,u,z,approx=F){
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
    cs$v[cs$v==0]=0.000001
    
    #now find lognormal parameters from mean and variance
    #https://www.johndcook.com/blog/2022/02/24/find-log-normal-parameters/
    cs$sig2 = log(1+cs$v/cs$m^2)
    cs$mu = log(cs$m)-cs$sig2/2
    return(cs)
  }else{
    return(c)
  }
}