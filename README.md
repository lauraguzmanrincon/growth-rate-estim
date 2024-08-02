
Source code to use the model for estimating and short-term projecting cases of infectious diseases. 

Full manuscript:
*Bayesian estimation of real-time epidemic growth rates using Gaussian processes*. https://academic.oup.com/jrsssc/article/72/5/1413/7209012

## 1. Set up workspace
Load required packages (main and optional packages) and source code.

```R
# Main packages
library(data.table)
library(dplyr)
library(INLA)
library(rstan)

# Optional packages (for figures)
library(ggplot2)
library(ggnewscale)
library(scales)
library(egg)
library(mvtnorm)

# Load source code
source("functions.R")
source("MCMC_functions.R")
source("functions_Forecast.R")
```

## 2. Set up model parameters
Set up prior distributions for the parameters of the model.
```R
parametersModel <- setParametersFn(modelType = "positives",
                                   interceptPrior = list(prior = "gaussian", mean = 0, prec = 1/100),
                                   dispersionPrior = list(prior = "gaussian", mean = 0, prec = 0.5),
                                   unitTime = "day",
                                   randomEffect = "weekday",
                                   randomEffectPrior = list(prior = "loggamma", a = 1, b = 0.01),
                                   GPcovariance = list(list(type = "rationalQ",
                                                            param = c(sigma = T, range = T, alpha = T),
                                                            theta0 = c(120, 1, 1))),
                                   GPHyperparamPrior = list(prior = "multivariateGaussian",
                                                            mean = c(120, 1, 1),
                                                            B = diag(3)),
                                   sizeSample = 1000,
                                   derivativeFromGP = F)
```

## 3. Run model
Run model and produces output object. Requires data frame containing the columns 'date', 'positiveResults' and 'numberTest'. (The column 'numberTest' is optional for the positives model).
```R
modelObject <- runModelGrowthRate(countTable = countTable,
                                  parametersModel = parametersModel,
                                  minDate = as.Date("2021-12-23"),
                                  maxDate = as.Date("2022-03-22"))
modelOutput <- processINLAOutput(objectInla = modelObject,
                                 parametersModel = parametersModel,
                                 saveSamples = T)
```

## 4. Compute projections
Compute projections based on the model output. The optional data frame 'futureTable' contains the columns 'date', 'positiveResults' and 'numberTest' for the dates being projected.
```R
projectionObject <- getProjectionGP2(parametersModel = parametersModel,
                               outputModel = modelOutput,
                               currentDate = as.Date("2022-03-22"),
                               sizePrediction = 60,
                               daysForPrediction = 90)
projectionTable <- createTableProjection(projectionSamplesGP = projectionObject,
                                    outputModel = outputModelAll,
                                    parametersModel = parametersModel,
                                    nameProjection = "Example",
                                    testingVector = NULL)
projectionPlot <- plotProjection(listProjection = projectionTable,
                                 outputModel = modelOutput,
                                 parametersModel = parametersModel,
                                 futureTable = futureTable,
                                 minDate = NULL,
                                 plotGP = F)
```
## 5. Plot output
Plots the output from the model and the projection.

```R
## Model
egg::ggarrange(plotFitting(outputModel = modelOutput,
                           parametersModel),
               plotGR(outputModel = modelOutput),
               ncol = 1)
## Projection
egg::ggarrange(plots = projectionPlot,
               ncol = 1)


```
