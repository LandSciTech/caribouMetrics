# Comparing the caribouMetrics Beta model to bboutools models

## 0.1 Similarities and differences between the caribouMetrics Beta model and the bboutools model

Hughes et al. ([2025](#ref-hughes_integration_2025)) described a Beta
model with disturbance covariates and informative priors for analysis of
a single boreal caribou population with survival data aggregated by year
and only cows and calves reported in the composition survey. To align
with and allow comparison to the models implemented in bboutools
([Dalgarno et al. 2025](#ref-dalgarno_bbousuite_2025)) we extended the
Hughes et al. ([2025](#ref-hughes_integration_2025)) models as follows:

- Allow analysis of composition survey data in bboutools form that
  includes Yearlings, Bulls, and UnknownAdults. See [Analytic Methods
  for Estimation of Boreal Caribou Survival, Recruitment and Population
  Growth](https://poissonconsulting.github.io/bboutools/articles/methods.html)
  for details.
- Allow analysis of monthly survival data in bboutools form
  (e.g. [`bboudata::bbousurv_a`](https://poissonconsulting.github.io/bboudata/reference/bbousurv_a.html)).
- Allow for calculation of population growth rate without initial
  population size information. See [Analytic Methods for Estimation of
  Boreal Caribou Survival, Recruitment and Population
  Growth](https://poissonconsulting.github.io/bboutools/articles/methods.html)
  for details.
- Allow for analysis of groups of populations with shared interannual
  variation. See [Multi-Population Analysis and Other
  Extensions](https://poissonconsulting.github.io/bboutools/articles/extensions.html)
  for details.

Note that bboutools includes a number of different model variants.
Throughout this document we compare to the Bayesian model with a random
effect of year and no annual trend. For brevity, we refer to that as the
“bbou” model.

Key remaining differences between the caribouMetrics Beta model and the
bbou model include:

- The Beta model includes anthropogenic and fire disturbance covariates.
- The Beta model intercept and disturbance covariate slope parameters
  are informed by national demographic-disturbance relationships.
- The distribution of interannual variation differs between the models.
- The bbou model allows for variation in survival among months.

## 0.2 Example data

We compare the Beta and bbou models using three examples derived from
the
[`bboudata::bbousurv_a`](https://poissonconsulting.github.io/bboudata/reference/bbousurv_a.html)
and
[`bboudata::bbourecruit_a`](https://poissonconsulting.github.io/bboudata/reference/bbourecruit_a.html)
data sets:

- No local data: used to examine and compare prior predictions from the
  models.
- Informative local data: Observations from 2010 to 2016.
- Limited local data: Observations from 2014 to 2016.

To examine and compare predictions for unobserved years (2018 to 2021)
we add missing data.

``` r
library(caribouMetrics)
# use local version on local and installed on GH
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all()
library(bboudata)
library(bboutools)
library(dplyr)
library(ggplot2)
library(patchwork)

figWidth <- 8
figHeight <- 10

useSaved <- T # option to skip slow step of fitting bboutools model
bbouInformativeFile <- here::here("results/vignetteBbbouExample.rds")
bbouLimitedFile <- here::here("results/vignetteBbbouExample2.rds")
bbouPriorFile <- here::here("results/vignetteBbbouExample1.rds")
bbouPriorNationalFile <- here::here("results/vignetteBbbouExample1p.rds")
betaAnthroFile <- here::here("results/betaAnthroExample.rds")

surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
surv_data_add <- expand.grid(Year = seq(2017, 2022), Month = seq(1:12),
                             PopulationName = unique(surv_data$PopulationName))
surv_data <- merge(surv_data, surv_data_add, all.x = TRUE, all.y = TRUE)

recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
recruit_data_add <- expand.grid(Year = seq(2017, 2022), PopulationName = unique(recruit_data$PopulationName))
recruit_data <- merge(recruit_data, recruit_data_add, all.x = TRUE, all.y = TRUE)

surv_dataNone <- surv_data %>% filter(Year > 2017)
recruit_dataNone <- recruit_data %>% filter(Year > 2017)

surv_dataLimited <- surv_data %>% filter(Year > 2014)
recruit_dataLimited <- recruit_data %>% filter(Year > 2014)
```

## 0.3 Comparison of the Beta and national models

As shown by Hughes et al. ([2025](#ref-hughes_integration_2025)) the
prior means and 95% prior predictive intervals from the Beta model are
similar to the means and ranges between the 2.5% and 97.5% quantiles of
3000 simulated survival and recruitment trajectories from the national
model (Figure [1](#fig:fig-plot1)).

``` r
disturbance <- data.frame(Year = unique(surv_data$Year), Anthro = 5,
                          Fire_excl_anthro = 1)
betaPrior <- bayesianTrajectoryWorkflow(surv_dataNone, recruit_dataNone, disturbance)
simNational <- trajectoriesFromNational(disturbance = disturbance)
out_tbls <- compareTrajectories(betaPrior, simInitial = simNational)
typeLabs <- c("Beta", "National")
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs)
surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs)
lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs,
               lowBound = 0.5, highBound = 1.5)

rec / surv / lam
```

![95% prior predictive intervals from Beta model with 5% anthropogenic
disturbance, default priors, and no local data compared to simulations
from the national
model.](compare-bayesian-models_files/figure-html/fig-plot1-1.png)

Figure 1: 95% prior predictive intervals from Beta model with 5%
anthropogenic disturbance, default priors, and no local data compared to
simulations from the national model.

## 0.4 Comparison of Beta and bbou prior predictions

The bbou model default priors are less informative than the Beta model
default priors (Figure [2](#fig:fig-plotPriorsDefault)). The bbou model
assumes that little is known about a population with no local data,
whereas when there is no local data the Beta model assumes that local
recruitment and survival rates are within the range of variation among
boreal caribou populations there were included in the national model.

``` r
if (useSaved & file.exists(bbouPriorFile)) {
  bbouPrior <- readRDS(bbouPriorFile)
} else {
  bbouPrior <- estimateBayesianRates(surv_dataNone, recruit_dataNone, 
                                    return_mcmc = TRUE)
  if (dir.exists(dirname(bbouPriorFile))) {
    saveRDS(bbouPrior, bbouPriorFile)
  }
}
simBbouPrior <- trajectoriesFromBayesian(bbouPrior)
out_tbls <- compareTrajectories(betaPrior, simInitial = simBbouPrior)
typeLabs <- c("Beta", "Bbou")
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs)

surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs)

lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs, 
               lowBound = 0, highBound = 2)

rec / surv / lam
```

![Prior means and 95% predictive intervals from Beta and bbou models
with 5% anthropogenic disturbance and default
priors.](compare-bayesian-models_files/figure-html/fig-plotPriorsDefault-1.png)

Figure 2: Prior means and 95% predictive intervals from Beta and bbou
models with 5% anthropogenic disturbance and default priors.

The bbou priors can be set to (approximately) match the Beta model
priors. When disturbance is constant, the models with informative priors
make similar prior predictions of expected demographic rates. The bbou
model assumes more prior uncertainty about interannual variation, so
prior predictions of the distribution of observed demographic rates in
each year continue to differ. TODO explain what expected demographic
parameters are and how they are different from the plots shown on the
right.

``` r
b0Priors <- bbouNationalPriors(anthro = unique(disturbance$Anthro), fire_excl_anthro = unique(disturbance$Fire_excl_anthro), month = TRUE)
if (useSaved & file.exists(bbouPriorNationalFile)) {
  bbouPriorNational <- readRDS(bbouPriorNationalFile)
} else {
  bbouPriorNational <- estimateBayesianRates(surv_dataNone, recruit_dataNone,
                                            return_mcmc = TRUE, priors = b0Priors)
  if (dir.exists(dirname(bbouPriorNationalFile))) {
    saveRDS(bbouPriorNational, bbouPriorNationalFile)
  }
}
simBbouPriorNational <- trajectoriesFromBayesian(bbouPriorNational)
out_tbls <- compareTrajectories(betaPrior, simInitial = simBbouPriorNational)
typeLabs <- c("Beta", "Bbou National")
```

``` r
recBar <- plotCompareTrajectories(out_tbls, "Expected recruitment", typeLabels = typeLabs)

rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs)

survBar <- plotCompareTrajectories(out_tbls, "Expected survival", typeLabels = typeLabs)

surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs)

lamBar <- plotCompareTrajectories(out_tbls, "Expected growth rate", typeLabels = typeLabs,
                  lowBound = 0, highBound = 2)

lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs, 
               lowBound = 0, highBound = 2)
(recBar + rec) / (survBar + surv) / (lamBar + lam)
```

![Prior means and 95% predictive intervals from Beta and bbou models
with 5% anthropogenic disturbance and priors informed by national
demographic-disturbance
relationships.](compare-bayesian-models_files/figure-html/fig-plotPriorsNational-1.png)

Figure 3: Prior means and 95% predictive intervals from Beta and bbou
models with 5% anthropogenic disturbance and priors informed by national
demographic-disturbance relationships.

## 0.5 Comparison of Beta and bbou models with informative local data

With more informative data, the Beta model with constant disturbance
covariates and informative priors often gives results that are
comparable to the bbou model with default (less informative) priors.
When there is enough local data available differences in the priors are
less important.

``` r
betaInformative <- bayesianTrajectoryWorkflow(surv_data, recruit_data, disturbance)
if (useSaved & file.exists(bbouInformativeFile)) {
  bbouInformative <- readRDS(bbouInformativeFile)
} else {
  bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
                                          return_mcmc = TRUE)
  if (dir.exists(dirname(bbouInformativeFile))) {
    saveRDS(bbouInformative, bbouInformativeFile)
  }
}
simBbouInformative <- trajectoriesFromBayesian(bbouInformative)
out_tbls <- compareTrajectories(betaInformative, simInitial = simBbouInformative)
typeLabs <- c("Beta", "Bbou")
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs)

surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs)

lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs,
               lowBound = 0.5, highBound = 1.5)
rec / surv / lam
```

![Posterior means and 95% posterior predictive intervals from Beta and
bbou models with 5% anthropogenic disturbance, default priors, and
informative local
data.](compare-bayesian-models_files/figure-html/fig-plotInformativeLocal-1.png)

Figure 4: Posterior means and 95% posterior predictive intervals from
Beta and bbou models with 5% anthropogenic disturbance, default priors,
and informative local data.

The bbou model does not include disturbance covariates, so we expect the
models to project different outcomes in a case where disturbance changes
over time.

``` r
disturbance <- data.frame(Year = seq(2011, 2051),
                          Anthro = seq(20, 2 * 40 + 20, length.out = 41),
                          Fire_excl_anthro = 1)
if (useSaved & file.exists(betaAnthroFile)){
  betaAnthroChange <- readRDS(betaAnthroFile)
} else {
  betaAnthroChange <- bayesianTrajectoryWorkflow(surv_data, recruit_data, disturbance)  
  if(dir.exists(dirname(betaAnthroFile))){
    saveRDS(betaAnthroChange, betaAnthroFile)
  }
}

if (useSaved & file.exists(bbouInformativeFile)) {
  bbouInformative <- readRDS(bbouInformativeFile)
} else {
  bbouInformative <- estimateBayesianRates(surv_data, recruit_data,
                                          return_mcmc = TRUE)
  if (dir.exists(dirname(bbouInformativeFile))) {
    saveRDS(bbouInformative, bbouInformativeFile)
  }
}
simBbouInformative <- trajectoriesFromBayesian(bbouInformative)
out_tbls <- compareTrajectories(betaAnthroChange, simInitial = simBbouInformative)
typeLabs <- c("Beta", "Bbou")
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs,
               breakInterval = 5)

surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs, 
                breakInterval = 5)

lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs, 
               lowBound = 0.5, highBound = 1.5, breakInterval = 5)
rec / surv / lam 
```

![Posterior means and 95% posterior predictive intervals from Beta and
bbou models with default priors, informative local data, and
anthropogenic disturbance increasing from 20 to 100% over 40
years.](compare-bayesian-models_files/figure-html/fig-plotDisturbanceChange-1.png)

Figure 5: Posterior means and 95% posterior predictive intervals from
Beta and bbou models with default priors, informative local data, and
anthropogenic disturbance increasing from 20 to 100% over 40 years.

## 0.6 Comparison of Beta and bbou models with limited local data

If local data is limited the bbou model with default (less informative)
priors and the Beta model with disturbance covariates and informative
priors may also give different results, particularly if observations are
inconsistent with expectations from the observed distribution of
outcomes across the country. For example (Figure
[6](#fig:fig-plotLimitedLocal)) when disturbance is high the prior
expectation from the Beta model is that recruitment and survival will be
low. In this example there are a couple of observations of survival and
recruitment that are higher than expected given knowledge of
disturbance. Both models acknowledge high uncertainty about the true
state of the population given limited local data. The Beta model
predicts lower growth than the bbou model in this case because it is
combining the limited local information with a prior expectation that
demographic rates will be low when disturbance is high.

``` r
disturbance <- data.frame(Year = unique(surv_data$Year), Anthro = 90, Fire_excl_anthro = 5)
betaLimited <- bayesianTrajectoryWorkflow(surv_dataLimited, recruit_dataLimited, disturbance)
if (useSaved & file.exists(bbouLimitedFile)) {
  bbouLimited <- readRDS(bbouLimitedFile)
} else {
  bbouLimited <- estimateBayesianRates(surv_dataLimited, recruit_dataLimited,
                                      return_mcmc = TRUE)
  if (dir.exists(dirname(bbouLimitedFile))) {
    saveRDS(bbouLimited, bbouLimitedFile)
  }
}
simBbouLimited <- trajectoriesFromBayesian(bbouLimited)
out_tbls <- compareTrajectories(betaLimited, simInitial = simBbouLimited)
typeLabs <- c("Beta", "Bbou")
```

``` r
rec <- plotCompareTrajectories(out_tbls, "Recruitment", typeLabels = typeLabs)

surv <- plotCompareTrajectories(out_tbls, "Adult female survival", typeLabels = typeLabs)

lam <- plotCompareTrajectories(out_tbls, "Population growth rate", typeLabels = typeLabs,
               lowBound = 0.4, highBound = 1.6)
rec / surv / lam
```

![Posterior means and 95% predictive intervals from Beta and bbou models
with 90% anthropogenic disturbance, default priors, and limited local
data.](compare-bayesian-models_files/figure-html/fig-plotLimitedLocal-1.png)

Figure 6: Posterior means and 95% predictive intervals from Beta and
bbou models with 90% anthropogenic disturbance, default priors, and
limited local data.

## References

Dalgarno, Seb, John Boulanger, Ayla Pearson, Joe Thorley, Troy Hegel,
Barry Nobert, and Dave Hervieux. 2025. “Bbousuite: A Set of R Packages
to Facilitate Analysis of Boreal Caribou Survival and Recruitment Data.”
*Journal of Open Source Software* 10 (109): 7997.
<https://doi.org/10.21105/joss.07997>.

Hughes, Josie, Sarah Endicott, Anna M. Calvert, and Cheryl A. Johnson.
2025. “Integration of National Demographic-Disturbance Relationships and
Local Data Can Improve Caribou Population Viability Projections and
Inform Monitoring Decisions.” *Ecological Informatics* 87 (July):
103095. <https://doi.org/10.1016/j.ecoinf.2025.103095>.
