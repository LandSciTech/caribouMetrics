# Instructions for using the app

## Overview

This tool is designed to allow users to explore the impacts of various factors that affect projections of caribou population dynamics using the Johnson et al. (2020) national models as a starting point.

A Bayesian model is fit that uses the national model as the priors and is updated based on a simulated observed data set. 


## Inputs

The model can be run with a simulated disturbance scenario and default values for all parameters by clicking **Run model**

### Scenario

The simulation scenario can be defined using either a csv file or several parameters that are used to generate a scenario. 

#### Scenario csv file

1. Choose the **Last year of observations**. Years before this will be included in the simulated observations and years after will be used for projections
2. Click **Select File** to select a csv file using a dialog box. The file must contain the columns "Year", "Anthro", and "fire_excl_anthro". The file must have one row per year that will be included in the observation and projection periods.
3. **Number of years of observation/projections** will appear after the file is loaded and can be used to reduce the number of years used for either period (described in detail below).

#### Simulate a scenario

1. Choose the **Last year of observations**. Years before this will be included in the simulated observations and years after will be used for projections
2. Choose a **Number of years of observations**. This value determines how much data the model will use in addition to the priors set from the national model. If it is low the projections will follow the national model. If it is higher and the observations are different from the national model the projections will be different from the national model. 
3. Choose a **Number of years of projections**. This is the number of years after the last observation that projections should be made for. 
4. Define the amount of anthropogenic disturbance in the scenario by selecting an **Initial % anthropogenic disturbance** and the **% increase in anthropogenic disturbance per year in observation/projection period**. The amount of disturbance caused by fire is assumed to be constant and can be set as the **Initial % natural disturbance**

### True population parameters

These parameters define the population that simulated observations will be sampled from relative to the national model. 

1. The **Initial population size** is the number of individuals at the start of the observations
2. The **Recruitment and Survival quantiles** determine the quantile of the distribution in the national model that the recruitment and survival coefficients are taken from. A quantile of 0.5 will give a population that follows the expected behaviour from the national model, while a quantile of 0.025 would follow the lower bound of the confidence interval from the national model. 
3. The **Multipliers for the effect of disturbance on recruitment/survival** are used to adjust the coefficients from the national model. Setting a value of 0 assumes that disturbance has no effect on recruitment/survival, while a value of 2 would double the effect of disturbance compared to the national model.


### Observation data parameters

These parameters define the assumptions made about the sampling strategies used to collect observations. These affect the amount of data behind survival and recruitment observations for each year and therefore the amount of weight they are given in the model.

1. **Target number of collars**
2. **Number of years between collar deployments**
3. **Number of cows per collared cow in aerial surveys for calf:cow ratio each year**
2. **Number of years until collar falls off**
3. **Month that collars are deployed**. A number from 1 (January) to 12 (December)
4. **Month that collars fall off**. A number from 1 (January) to 12 (December)

### Model priors 

The first choice is which **version of the national model** to use. Johnson et al. 2020 includes 5 different recruitment and survival models the defaults are the best supported models but you can choose a different model number if desired. Note only models that use anthropogenic disturbance and/or fire excluding anthropogenic disturbance are supported at this time. Another option is to provide a csv file with custom coefficients where the column names are "responseVariable", "Coefficient", "Value", "lowerCI", and "upperCI". "responseVariable" must be one of "recruitment" or "femaleSurvival". And "Coefficient" must be one of "Intercept", "Anthro", "fire_excl_anthro" or "Precision". 
The remaining parameters in this section affect the level of uncertainty around the coefficients from the national model. If uncertainty in the national model is higher then the projections can deviate more from the national model if the observed data does. 

1. **Multiplier for uncertainty of effect of disturbance on survival**
2. **Multiplier for uncertainty of effect of disturbance on recruitment**
3. **Multiplier for uncertainty of survival intercept**
4. **Multiplier for uncertainty of recruitment intercept**
5. **Interannual coefficient of variation for survival**
6. **Uncertainty around interannual coefficient of variation for survival**
7. **Interannual coefficient of variation for recruitment**
8. **Uncertainty around interannual coefficient of variation for recruitment**

### Baysian model parameters

These parameters affect the fitting of the Bayesian model used to make predictions.

1. **Number of chains**
2. **Number of iterations**
3. **Length of burn-in**
4. **Thinning rate**

## Outputs

Running the model can take some time, see the R console for details on the model progress. 

Once the model is finished the results graphs and tables tabs will be populated. You can also view the MCMC diagnostic plots in a new window by clicking **Open MCMC diagnostic plots**. This will also take some time with progress shown in the R console. 

Click **Save results to csv or rds** to save the outputs. The csv file will contain all the parameters used and the projected population metrics. If you save the results as a .rds file the whole model object will be saved and can be loaded into an R session using `readRDS()`




















