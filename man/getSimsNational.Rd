% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSimsNational.R
\name{getSimsNational}
\alias{getSimsNational}
\title{Get a set of simulation results from the national demographic model}
\usage{
getSimsNational(
  replicates = 1000,
  N0 = 1000,
  Anthro = seq(0, 100, by = 1),
  fire_excl_anthro = 0,
  useQuantiles = NULL,
  populationGrowthTable = NULL,
  cPars = getScenarioDefaults(),
  forceUpdate = F,
  interannualVar = eval(formals(caribouPopGrowth)$interannualVar)
)
}
\arguments{
\item{replicates}{integer. Number of replicate populations.}

\item{N0}{initial population size}

\item{Anthro, fire_excl_anthro}{numeric. A vector of numbers between 0 and 100
representing the percentage of the landscape covered by anthropogenic
disturbance buffered by 500 m, and the percentage covered by fire that does
not overlap anthropogenic disturbance. The two vectors will be combined
with \code{expand.grid()} to give the set of scenarios simulated.}

\item{useQuantiles}{logical or numeric. If it is a numeric vector it must be
length 2 and give the low and high limits of the quantiles to use. If
\code{useQuantiles != FALSE}, each replicate population is assigned to a
quantile of the distribution of variation around the expected values, and
remains in that quantile as covariates change. If \code{useQuantiles = TRUE},
replicate populations will be assigned to quantiles in the default range of
0.025 and 0.975.}

\item{populationGrowthTable}{data.frame.\link{popGrowthTableJohnsonECCC} is
included in the package and should be used in most cases. A custom table of
model coefficients and standard errors or confidence intervals can be
provided but it must match the column names of \link{popGrowthTableJohnsonECCC}.
If the table does not contain the standard error it is calculated from the
confidence interval.}

\item{cPars}{optional. Parameters for calculating composition survey bias term.}

\item{forceUpdate}{logical. If the default inputs are used the result is
cached. Set \code{forceUpdate} to TRUE to ensure the simulations are re-run.}

\item{interannualVar}{list or logical. List containing interannual
variability parameters. These can be either coefficients of variation
(R_CV, S_CV), beta precision parameters (R_phi, S_phi),
or random effects parameters from a logistic glmm (R_annual, S_annual).
Set to \code{FALSE} to ignore interannual variability.}
}
\value{
a list with two elements:
\itemize{
\item summary: a tibble with a summary of parameter values for each scenario.
Column names are Anthro, Mean, lower, upper, Parameter.
\item samples: a tibble with parameter values for each scenario and replicate
4 rows per replicate \* scenario. Column names are Anthro, Parameter and Value
}
}
\description{
Get a set of simulation results from the national demographic model
}
\examples{
getSimsNational()
}
\seealso{
Caribou demography functions:
\code{\link{caribouBayesianPM}()},
\code{\link{caribouPopGrowth}()},
\code{\link{compositionBiasCorrection}()},
\code{\link{demographicCoefficients}()},
\code{\link{demographicProjectionApp}()},
\code{\link{demographicRates}()},
\code{\link{getOutputTables}()},
\code{\link{getPriors}()},
\code{\link{getScenarioDefaults}()},
\code{\link{plotRes}()},
\code{\link{popGrowthTableJohnsonECCC}},
\code{\link{runScnSet}()},
\code{\link{simulateObservations}()}
}
\concept{demography}
