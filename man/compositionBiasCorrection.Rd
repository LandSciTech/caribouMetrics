% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositionBiasCorrection.R
\name{compositionBiasCorrection}
\alias{compositionBiasCorrection}
\title{Calculate bias correction term for calf:cow composition survey.}
\usage{
compositionBiasCorrection(w, q, u, z, approx = F)
}
\arguments{
\item{w}{number. The apparent number of adult females per collared animal in
composition survey.}

\item{q}{number in 0, 1. Ratio of bulls to cows in composition survey groups.}

\item{u}{number in 0, 1. Probability of misidentifying young bulls as adult
females and vice versa in composition survey.}

\item{z}{number in 0, <1. Probability of missing calves in composition
survey.}

\item{approx}{logical. If TRUE approximate the uncertainty about the value of
the composition bias correction value (c) with the log-normal distribution
of c given all the supplied values of \code{q}, \code{u}, and \code{z}. If FALSE the
composition bias correction value (c) is returned for each value of \code{q},
\code{u}, and \code{z}}
}
\value{
number or tibble. If \code{approx = FALSE} a vector of composition bias
correction values (c) of the same length as \code{q}, \code{u}, and \code{z}. If \code{approx = TRUE} a tibble with on row per unique value of \code{w} and columns \code{w}, \code{m},
\code{v}, \code{sig2}, \code{mu} representing \code{w}, mean \code{c}, variance of \code{c}, and parameters for a
log-normal approximation of the distribution of \code{c}.
}
\description{
When composition surveys are conducted there is a possibility of bias in calf
cow ratios due to misidentifying young bulls as adult females and vice versa
or missing calves. Here we address this gap with a bias term derived from a simple
model of the recruitment survey observation process. See
\href{https://doi.org/10.1016/j.ecoinf.2025.103095}{Hughes et al. (2025) Section 2.2}
for a detailed description of the model.
}
\examples{
# number or reps
nr <- 10

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = FALSE)

compositionBiasCorrection(w = 6,
                          q = runif(nr, 0, 0.6),
                          u = runif(nr, 0, 0.2),
                          z = runif(nr, 0, 0.2),
                          approx = TRUE)


}
\references{
Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
Integration of national demographic-disturbance relationships and local
data can improve caribou population viability projections and inform
monitoring decisions. Ecological Informatics, 87, p.103095.
\url{https://doi.org/10.1016/j.ecoinf.2025.103095}
}
\seealso{
Caribou demography functions:
\code{\link{caribouBayesianPM}()},
\code{\link{caribouPopGrowth}()},
\code{\link{demographicCoefficients}()},
\code{\link{demographicProjectionApp}()},
\code{\link{demographicRates}()},
\code{\link{getOutputTables}()},
\code{\link{getPriors}()},
\code{\link{getScenarioDefaults}()},
\code{\link{getSimsNational}()},
\code{\link{plotRes}()},
\code{\link{popGrowthTableJohnsonECCC}},
\code{\link{runScnSet}()},
\code{\link{simulateObservations}()}
}
\concept{demography}
