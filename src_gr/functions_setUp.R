
#' Set parameters of model
#'
#' `setParametersFn()` sets the parameter priors of the model, sample size
#'
#' @param modelType 
#' @param interceptPrior 
#' @param dispersionPrior 
#' @param unitTime 
#' @param randomEffect 
#' @param randomEffectPrior 
#' @param GPcovariance 
#' @param GPHyperparamPrior 
setModelParameters <- function(modelType,
                               interceptPrior = NULL,
                               dispersionPrior = NULL,
                               unitTime = "day", # day week integer
                               randomEffect = "weekday", # weekday all none # ?? if day, adds day effect automatically # old dayWeek
                               randomEffectPrior = NULL,
                               GPcovariance = list(list(type = "matern12", param = c(sigma = T, range = T), theta0 = c(1, 100)),
                                                   list(type = "periodic", param = c(sigma = T, range = T, period = T), theta0 = c(1, 1, 10)),
                                                   list(type = "rationalQ", param = c(sigma = T, range = T, alpha = T), theta0 = c(1, 10, 1))),
                               GPHyperparamPrior = list(prior = "multivariateGaussian", B = NULL)){
  # Checks
  modelTypeValues <- c("positives", "proportions")
  if(!modelType %in% modelTypeValues) stop(paste0("Parameter 'modelType' can only take value from '", paste(modelTypeValues, collapse = "','"),"'."))
  randomEffectValues <- c("weekday", "all", "none")
  if(!randomEffect %in% randomEffectValues) stop(paste0("Parameter 'randomEffect' can only take value from '", paste(randomEffectValues, collapse = "','"),"'."))
  covarianceValues <- c("matern12", "periodic", "rationalQ")
  covarianceNotValid <- sapply(GPcovariance, function(cov) !cov$type %in% covarianceValues)
  if(sum(covarianceNotValid) != 0) stop(paste0("Type of 'GPcovariance' can only take value from '", paste(covarianceValues, collapse = "','"),"'."))
  if(randomEffect == "weekday" & unitTime != "day") stop("randomEffect cannot be 'weekday' if unitTime is not 'day'.")
  
  # New variables
  linkType <- ifelse(modelType == "positives", "NB", "BB")
  
  # Default
  if(is.null(interceptPrior)) interceptPrior <- list(prior = "gaussian", mean = 0, prec = 1) # TODO which default should we choose?
  if(is.null(dispersionPrior)) dispersionPrior <- list(prior = "gaussian", mean = 0, prec = ifelse(linkType == "NB", 0.5, 0.01))
  if(is.null(randomEffectPrior)) randomEffectPrior <- list(prior = "loggamma", a = 1, b = 0.01)
  if(is.null(GPcovariance)) GPcovariance <- list(list(type = "matern12", param = c(sigma = T, range = T), theta0 = c(1, 100)))
  GPcovarianceList <- list(type = sapply(GPcovariance, function(cov) cov$type),
                           indexCov = do.call(c, lapply(1:length(GPcovariance), function(i) rep(i, length(GPcovariance[[i]]$param)))),
                           indexPar = do.call(c, lapply(1:length(GPcovariance), function(i) 1:length(GPcovariance[[i]]$param))),
                           theta0 = do.call(c, lapply(GPcovariance, function(cov) cov$theta0)))
  if(is.null(GPHyperparamPrior)) GPHyperparamPrior <- list(prior = "multivariateGaussian",
                                                           mean = GPcovarianceList$theta0,
                                                           B = diag(length(GPcovarianceList$theta0)))
  if(!is.null(GPHyperparamPrior) & is.null(GPHyperparamPrior$mean)) GPHyperparamPrior$mean <- GPcovarianceList$theta0
  if(!is.null(GPHyperparamPrior) & is.null(GPHyperparamPrior$B)) GPHyperparamPrior$B <- diag(length(GPcovarianceList$theta0))
  # TODO what if not updating one of the parameters?
  
  # Output
  parameters <- list(
    # Model type
    modelType = modelType,
    linkType = linkType,
    
    # Overdispersion
    interceptPrior = interceptPrior,
    dispersionPrior = dispersionPrior,
    
    # Random effect
    unitTime = unitTime,
    randomEffect = randomEffect,
    randomEffectPrior = randomEffectPrior,
    
    # GP hyperparameters
    numKernels = length(GPcovariance),
    GPcovarianceList = GPcovarianceList,
    GPHyperparamPrior = GPHyperparamPrior
  )
  return(parameters)
}



#' Title
#'
#' @param inferenceType 
#' @param derivativeFromGP 
#' @param numSamples (LA-INLA)
#' @param numIterations (MCMC-STAN)
#' @param warmUp (MCMC-STAN)
#' @param thinning (MCMC-STAN)
#' @param sampleFile (MCMC-STAN)
#'
#' @return
#' @export
#'
#' @examples
setInferenceSettings <- function(inferenceType = "LA-INLA",
                                 derivativeFromGP = F,
                                 seed = sample.int(.Machine$integer.max, 1),
                                 numSamples = 1000,
                                 numIterations = 1000,
                                 warmUp = 100,
                                 thinning = 1,
                                 numChains = 1,
                                 sampleFile = NULL
                                 ){
  # Checks
  typeValues <- c("LA-INLA", "MCMC-STAN", "MCMC-Int")
  if(!inferenceType %in% typeValues) stop(paste0("Parameter 'inferenceType' can only take value from '", paste(typeValues, collapse = "','"),"'."))
  if(warmUp >= numIterations) stop("Parameter 'warmUp' should be an integer less than 'numIterations'")
  
  config <- list(
    inferenceType = inferenceType,
    derivativeFromGP = derivativeFromGP,
    seed = seed
  )
  
  if(inferenceType == "LA-INLA"){
    inferenceParam <- list(numSamples = numSamples)
  }else{
    inferenceParam <- list(numChains = numChains,
                           numIterations = numIterations,
                           warmUp = warmUp,
                           thinning = thinning,
                           sampleFile = sampleFile)
  }
  
  return(c(config, inferenceParam))
}

#' Title
#' 
#' Equally spaced time steps between minDate and maxDate
#'
#' @param countTable 
#' @param unitTime "day" "week" "integer"
#' @param minDate (integer or date)
#' @param maxDate (integer or date)
#'
#' @return
#' @export
#'
#' @examples
constructDateList <- function(countTable, unitTime, minDate, maxDate){
  if(unitTime == "integer"){
    dateSeq <- seq(from = minDate, to = maxDate, by = unitTime)
  }else{
    dateSeq <- seq.Date(from = minDate, to = maxDate, by = unitTime)
  }
  
  minDay <- 1
  maxDay <- length(dateSeq)
  numDays <- maxDay
  dateTable <- data.table(dayId = minDay:maxDay,
                          date = dateSeq)
  if(nrow(dateTable) <= 1) stop("There must be at least 2 time units between minDate and maxDate")
  
  dateList <- list(dateTable = dateTable,
                   numDays = numDays,
                   minDay = minDay,
                   maxDay = maxDay,
                   minDate = minDate,
                   maxDate = maxDate,
                   unitTime = unitTime)
  
  return(dateList)
}

#' Title
#' 
#' Create data with all days, including the ones with missing data
#'
#' @param countTable 
#' @param dateList 
#'
#' @return
#' @export
#'
#' @examples
constructInputDataTable <- function(countTable, dateList){
  dataForModel <- dateList$dateTable[order(dayId), .(dayId,
                                                     date,
                                                     numberTest = as.integer(NA),
                                                     positiveResults = as.integer(NA))]
  setkey(dataForModel, date)
  setkey(countTable, date)
  dataForModel[countTable, ":="(numberTest = i.numberTest, positiveResults = i.positiveResults)]
  dataForModel[numberTest == 0, ":="(numberTest = NA, positiveResults = NA)]
  if(dateList$unitTime == "day"){
    dataForModel[, dayWeek := weekdays(date)]
  }else{
    dataForModel[, dayWeek := NA]
  }
  
  return(dataForModel)
}

#' (Private)
#'
#' @return
#' @export
#'
#' @examples
getInternalSettings <- function(){
  internal <- list(
    levelsWeek = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  )
  return(internal)
}

