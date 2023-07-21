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
compositionBiasCorrection<-function(w,q,u,z){
  if(max(z)>=1){
    stop("Composition bias correct term is undefined when the probability of missing calves z >= 1.")
  }
  c = w*(q*u+1-u)/((q*u+w-u)*(1-z))
  return(c)
}