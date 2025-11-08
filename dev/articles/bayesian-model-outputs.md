# Options for using Bayesian model results

``` r
library(caribouMetrics)
# use local version on local and installed on GH
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all()
library(bboudata)
library(bboutools)
library(dplyr)
library(ggplot2)
library(patchwork)

resFile <- here::here("results/vignetteBbbouExample.rds")
useSaved <- TRUE # option to skip slow step of fitting bboutools model
```

## 0.1 Demographic trajectories and summaries from bboutools models

[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md)
returns fitted bboutools survival and recruitment models (assuming no
trends over time) and summaries of those models. In this example, we add
5 missing years that will be predicted by the model.

``` r
if (useSaved & file.exists(resFile)) {
  mod <- readRDS(resFile)
} else {
  surv_data <- bboudata::bbousurv_a %>% filter(Year > 2010)
  surv_data_add <- expand.grid(
    Year = seq(2017, 2022), 
    Month = seq(1:12),
    PopulationName = unique(surv_data$PopulationName)
  )
  surv_data <- merge(surv_data, surv_data_add, all.x = T, all.y = T)
  surv_data$StartTotal[is.na(surv_data$StartTotal)] <- 1

  recruit_data <- bboudata::bbourecruit_a %>% filter(Year > 2010)
  recruit_data_add <- expand.grid(
    Year = seq(2017, 2022), 
    PopulationName = unique(recruit_data$PopulationName)
  )
  recruit_data <- merge(recruit_data, recruit_data_add, all.x = T, all.y = T)
  recruit_data$Month[is.na(recruit_data$Month)] <- 3
  recruit_data$Day[is.na(recruit_data$Day)] <- 15

  mod <- estimateBayesianRates(surv_data,
    recruit_data,
    N0 = NA, return_mcmc = TRUE
  )
  if (dir.exists(dirname(resFile))) {
    saveRDS(mod, resFile)
  }
}
```

``` r
names(mod)
#> [1] "parTab"      "surv_fit"    "recruit_fit"
mod$parTab
#>   pop_name     R_bar      R_sd R_iv_mean R_iv_shape R_bar_lower R_bar_upper
#> 1        A 0.1927174 0.1974911 0.3554032   1.912094   0.1359539   0.2575156
#>       S_bar      S_sd S_iv_mean S_iv_shape S_bar_lower S_bar_upper N0
#> 1 0.9415866 0.5854533  0.607396   1.547401   0.8530804   0.9833124 NA
#>   nCollarYears nSurvYears nCowsAllYears nRecruitYears
#> 1          185         12            NA            12
```

Use
[`simulateTrajectoriesFromPosterior()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md)
to predict demographic trajectories from the mcmc samples. Use
[`summarizeTrajectories()`](https://landscitech.github.io/caribouMetrics/dev/reference/simulateTrajectoriesFromPosterior.md)
to get 95% prediction intervals from the demographic trajectories.

``` r
outmcmc <- simulateTrajectoriesFromPosterior(popInfo = NA, mod$recruit_fit, mod$surv_fit)
names(outmcmc)
#>  [1] "N0"                  "lambda"              "lambdaE"            
#>  [4] "N"                   "R_t"                 "X_t"                
#>  [7] "S_t"                 "n_recruits"          "surviving_adFemales"
#> [10] "lab"                 "Year"                "PopulationName"     
#> [13] "id"

# get 95% prediction intervals from demographic trajectories
PImcmc <- summarizeTrajectories(convertTrajectories(outmcmc))
PImcmc$summary$id <- NA
```

``` r

rec_traj_plt <- ggplot(subset(outmcmc, id <= 35), 
                       aes(x = Year, y = R_t, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "recruitment"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)


surv_traj_plt <- ggplot(subset(outmcmc, id <= 35), 
                        aes(x = Year, y = S_t, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "survival"),
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)


lambda_traj_plt <- ggplot(subset(outmcmc, id <= 35),
                          aes(x = Year, y = lambda, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "lambda"),
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

rec_traj_plt / surv_traj_plt / lambda_traj_plt
```

![Example trajectories and 95% prediction intervals from mcmc
samples.](bayesian-model-outputs_files/figure-html/plot1-1.png)

Figure 1: Example trajectories and 95% prediction intervals from mcmc
samples.

Note that the 95% prediction intervals from the mcmc samples match the
results returned by bboutools.

``` r

bbou_rec_plt <- bb_plot_year_calf_cow_ratio(mod$recruit_fit) + 
  theme_bw() + 
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "recruitment"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

surv_pred <- bb_predict_survival(mod$surv_fit)

bbou_surv_plt <- bb_plot_year_survival(surv_pred) + 
  theme_bw() + 
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "survival"),
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

bbou_lambda <- bb_predict_growth(survival = mod$surv_fit, 
                                 recruitment = mod$recruit_fit)

bbou_lambda_plt <- bb_plot_year_growth(bbou_lambda) + 
  theme_bw() + 
  geom_ribbon(data = subset(PImcmc$summary, MetricTypeID == "lambda"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

bbou_rec_plt / bbou_surv_plt / bbou_lambda_plt
```

![Comparison of 95% prediction intervals from mcmc samples (grey ribbon)
with results returned by bboutools (points and error
bars).](bayesian-model-outputs_files/figure-html/plot2-1.png)

Figure 2: Comparison of 95% prediction intervals from mcmc samples (grey
ribbon) with results returned by bboutools (points and error bars).

The model summary returned by
[`estimateBayesianRates()`](https://landscitech.github.io/caribouMetrics/dev/reference/estimateBayesianRates.md)
can also be used to simulate demographic trajectories. This simulates
trajectories using the modeled mean and standard deviation of
demographic rates and their interannual variation across the years in
the data. Therefore the distribution across years stays similar. **TODO
is there any difference or benefit to creating projections for the
future by adding missing years before fitting the bboutools model vs
using trajectoriesFromSummary?**

``` r
yrs <- unique(subset(PImcmc$summary, select = c(PopulationName, Year)))
yrs$time <- seq(1:nrow(yrs))
outParTab <- trajectoriesFromSummary(
  numSteps = nrow(yrs), replicates = 10000, N0 = NA, R_bar = mod$parTab$R_bar,
  S_bar = mod$parTab$S_bar,
  R_sd = mod$parTab$R_sd, S_sd = mod$parTab$S_sd,
  R_iv_mean = mod$parTab$R_iv_mean, S_iv_mean = mod$parTab$S_iv_mean,
  R_iv_shape = mod$parTab$R_iv_shape, S_iv_shape = mod$parTab$S_iv_shape,
  scn_nm = "base", addl_params = NULL, type = "logistic"
)
outParTab <- merge(outParTab, yrs)
PIParTab <- summarizeTrajectories(convertTrajectories(subset(outParTab, type == "samp")))
PIParTab$summary$id <- NA
PIParTab$summary$scn <- unique(outParTab$scn)
```

``` r
rec_sum_traj_plt <- ggplot(subset(outParTab, id <= 35), aes(x = Year, y = R_t, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "recruitment"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)


surv_sum_traj_plt <- ggplot(subset(outParTab, id <= 35), aes(x = Year, y = S_t, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "survival"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

lambda_sum_traj_plt <- ggplot(subset(outParTab, id <= 35), 
                              aes(x = Year, y = lambda, group = id)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "lambda"),
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

rec_sum_traj_plt / surv_sum_traj_plt / lambda_sum_traj_plt
```

![Example trajectories and 95% prediction intervals from the model
summary.](bayesian-model-outputs_files/figure-html/plot3-1.png)

Figure 3: Example trajectories and 95% prediction intervals from the
model summary.

Note that the 95% prediction intervals derived from the model summary
only match the results returned by bboutools in unobserved years. This
is because year-specific information is not included in the model
summary.

``` r
bbou_sum_rec_plt <- bb_plot_year_calf_cow_ratio(mod$recruit_fit) + 
  theme_bw() + 
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "recruitment"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

surv_pred <- bb_predict_survival(mod$surv_fit)

bbou_sum_surv_plt <- bb_plot_year_survival(surv_pred) + 
  theme_bw() + 
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "survival"), 
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

bbou_lambda <- bb_predict_growth(survival = mod$surv_fit, 
                                 recruitment = mod$recruit_fit)
bbou_sum_lambda_plt <- bb_plot_year_growth(bbou_lambda) + 
  theme_bw() +
  geom_ribbon(data = subset(PIParTab$summary, MetricTypeID == "lambda"),
              aes(x = Year, y = Mean, ymin = lower, ymax = upper), alpha = 0.1)

bbou_sum_rec_plt / bbou_sum_surv_plt / bbou_sum_lambda_plt
```

![Comparison of 95% prediction intervals derived from the model summary
with results returned by
bboutools.](bayesian-model-outputs_files/figure-html/plot4-1.png)

Figure 4: Comparison of 95% prediction intervals derived from the model
summary with results returned by bboutools.

## 0.2 Scenarios informed by bboutools models

We can explore the implications of changing demographic model parameters
by modifying the model summary table. In this example, we increase the
expected calf:cow ratio by 30%. The [Boreal Demographic Projection
Explorer app](https://github.com/LandSciTech/CaribouDemographyBasicApp)
provides an R Shiny interface for fitting bboutools models and exploring
scenarios informed by those models, with more thorough explanations.

``` r
parTabHighRec <- mod$parTab
parTabHighRec$R_bar <- parTabHighRec$R_bar * 1.3

outParTabHR <- trajectoriesFromSummary(
  numSteps = nrow(yrs), replicates = 10000, N0 = NA, R_bar = parTabHighRec$R_bar, 
  S_bar = parTabHighRec$S_bar,
  R_sd = parTabHighRec$R_sd, S_sd = parTabHighRec$S_sd,
  R_iv_mean = parTabHighRec$R_iv_mean, S_iv_mean = parTabHighRec$S_iv_mean,
  R_iv_shape = parTabHighRec$R_iv_shape, S_iv_shape = parTabHighRec$S_iv_shape,
  scn_nm = "highRec", addl_params = NULL, type = "logistic"
)
outParTabHR <- merge(outParTabHR, yrs)
PIParTabHR <- summarizeTrajectories(convertTrajectories(subset(outParTabHR, 
                                                                type == "samp")))
PIParTabHR$summary$id <- NA
PIParTabHR$summary$scn <- unique(outParTabHR$scn)

comp <- rbind(PIParTab$summary, PIParTabHR$summary)
comp <- subset(comp, !is.na(Mean))
```

``` r
comp_scn_plt <- ggplot(comp, 
                       aes(x = Year, y = Mean, ymin = lower, ymax = upper,
                           group = scn, colour = scn, fill = scn)) +
  geom_line(alpha = 0.3) +
  geom_ribbon(alpha = 0.1) +
  theme_bw() +
  facet_wrap(~Parameter)

comp_scn_plt
```

![Comparison of 95% prediction intervals derived from the model summary
with a scenario in which expected calf:cow ratio increases by
30%.](bayesian-model-outputs_files/figure-html/plot5-1.png)

Figure 5: Comparison of 95% prediction intervals derived from the model
summary with a scenario in which expected calf:cow ratio increases by
30%.

## 0.3 Other demographic modeling options and outputs

TO DO: show & discuss expected values of parameters - without
interannual variation.

TO DO: show effects of other demographic modeling options & outputs
i.e. changes in population size with and without carrying capacity.

TO DO: Add lots more text and explanation of what the functions are for,
what the different rates are and why you would want to calculate them
one way vs another. As well as general context within the package and
relative to bboutools.
