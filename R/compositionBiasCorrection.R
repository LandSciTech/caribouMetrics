#' Calculate bias correction term for calf:cow composition survey.
#'
#' When composition surveys are conducted there is a possibility of bias in calf
#' cow ratios due to misidentifying young bulls as adult females and vice versa
#' or missing calves.
#' 
#' # Model of bias in recruitment estimates from calf:cow surveys
#'
#' We assume each group of animals in a calf:cow composition survey contains one
#' or more collared adult females (\eqn{T}), and may also include: uncollared adult
#' females misidentified as young bulls or unknown sex (\eqn{U}); correctly
#' identified uncollared adult females (\eqn{V}); young bulls correctly identified
#' as male or unknown sex (\eqn{O}); young bulls misidentified as uncollared adult
#' females (\eqn{P}); observed calves (\eqn{J}); and unobserved calves (\eqn{K}). The
#' apparent number of adult females in the group is \eqn{T+V+P=Tw}, where \eqn{w} is a
#' multiplier that defines the apparent number of adult females as a function of
#' the number of collared animals. The ratio of young bulls to uncollared adult
#' females in the group is: \deqn{q = \frac{P+O}{U+V}}. Assuming an equal
#' probability \eqn{u} of misidentifying young bulls as adult females and vice
#' versa, we get \eqn{V=(U+V)(1-u)} and \eqn{P=(O+P)u}. Given a probability \eqn{z} of
#' missing calves, we get \eqn{J=(J+K)(1-z)}.
#'
#' Our objective is to model the sex and bias-corrected recruitment rate
#' \eqn{X=\frac{J+K}{2(T+U+V)}} as a function of the observed calf:cow ratio
#' \eqn{R=J/(T+V+P)}, the cow multiplier \eqn{w}, the ratio of young bulls to adult
#' females \eqn{q}, and the misidentification probabilities \eqn{u} and \eqn{z}. We start by
#' solving for \eqn{T+U+V} as a function of \eqn{q,w,u} and \eqn{T}. Recognize that
#' \eqn{P=Tw-T-V}, \eqn{U+V=V/(1-u)}, and \eqn{P+O=P/u} to write \eqn{q} as
#' \deqn{q=\frac{Tw-T-V}{uV/(1-u)}.} Rearrange to get
#' \deqn{V=\frac{T(w-1)(1-u)}{qu+1-u}.} Recognize that \eqn{U=Vu/(1-u)} to write
#' \eqn{T+U+V} as a function of \eqn{q,w,u} and \eqn{T}: \deqn{T+U+V=T\frac{qu+w-u}{qu+1-u}.}
#' Recognize that the number of observed calves \eqn{J} is the product of the
#' apparent recruitment rate and the apparent number of adult females \eqn{J=RTw},
#' and that therefore \eqn{J+K=RTw/(1-z)} to rewrite the bias corrected recruitment
#' rate \eqn{X=\frac{J+K}{2(T+U+V)}} as a function of \eqn{w,u,z} and \eqn{R}:
#' \deqn{X=R\frac{w(1+qu-u)}{2(w+qu-u)(1-z)}.} For simplicity, we write \eqn{X} as a
#' function of a bias correction term \eqn{c}: \deqn{c=\frac{w(1+qu-u)}{(w+qu-u)(1-z)};
#' X=cR/2.}
#' If we also adjust for delayed age at first reproduction (DeCesare et al.
#' 2012; Eacker et al. 2019), the adjusted recruitment rate becomes
#' \deqn{X=\frac{cR/2}{1+cR/2}.}
#' 
#' Uncertainty about the value of the bias correction term \eqn{c} can be
#' approximated with a Log-normal distribution. Given the apparent number of
#' adult females per collared animal \eqn{w} the mean and standard deviation of
#' \eqn{\log{c}} can be calculated for samples from the expected range of values
#' of \eqn{q}, \eqn{u} and \eqn{z}.
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
#' @return number or tibble. If `approx = FALSE` a vector of composition bias
#'   correction values (c) of the same length as `q`, `u`, and `z`. If `approx =
#'   TRUE` a tibble with on row per unique value of `w` and columns `w`, `m`,
#'   `v`, `sig2`, `mu` representing `w`, mean c, variance of c, and variance and
#'   mean of the approximated log-normal distribution of c.
#' @export
#' @family demography
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