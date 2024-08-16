

#' Title
#' 
#' Wrapper
#'
#' @param countTable data table with one row per day (of days with available data) and these columns:
#' - numberTest: number of test on the day (>= 0)
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - date: date of count in R date format
#' @param parametersModel output of setModelParameters
#' @param inferenceSettings output of setInferenceSettings
#' @param minDate (optional) minimum date to include in the model
#' @param maxDate (optional) maximum date to include in the model
#' @param dateList (optional, advanced) structure as from constructDateList... personalised
#' @param saveModelObject (optional, advanced) return output from INLA/STAN as output$objectInla/objectStan
#'
#' @return ...returns matrixSampleGP and matrixSampleGPDerivative,
#'              two matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' @export
#'
#' @examples
runModelGrowthRate <- function(countTable, parametersModel, inferenceSettings, minDate = NULL, maxDate = NULL, dateList = NULL, saveSamples = T, saveModelObject = F){
  # TODO check if min dates in countTable are aligned as in unitTime
  # TODO include cases with no date
  
  # Checks
  if(is.null(minDate)) minDate <- min(countTable$date)
  if(is.null(maxDate)) maxDate <- max(countTable$date)
  if(!is.null(dateList)){
    # Advanced
    if(dateList$unitTime != parametersModel$unitTime) stop("dateList$unitTime is not the same as parametersModel$unitTime")
  }
  
  # ---------------------------------------------------- #
  #                      SHAPE DATA                      #
  # ---------------------------------------------------- #
  
  if(is.null(dateList)){
    # Create date table with time points to make inference
    dateList <- constructDateList(countTable = countTable,
                                  unitTime = parametersModel$unitTime,
                                  minDate = minDate,
                                  maxDate = maxDate)
  }
  
  # Create data table with all days in dateList, including the ones with missing data
  dataForModel <- constructInputDataTable(countTable = countTable,
                                          dateList = dateList)
  
  # ---------------------------------------------------- #
  #                      FIT MODEL                       #
  # ---------------------------------------------------- #
  
  if(inferenceSettings$inferenceType == "LA-INLA"){
    outputInla <- runModelGrowthRate_INLA(dataForModel = dataForModel,
                                          dateList = dateList,
                                          parametersModel = parametersModel,
                                          inferenceSettings = inferenceSettings)
    output <- processINLAOutput(objectInla = outputInla,
                                parametersModel = parametersModel,
                                inferenceSettings = settingsINLA,
                                saveSamples = saveSamples,
                                saveInlaObject = saveModelObject)
  }else if(inferenceSettings$inferenceType == "MCMC-STAN"){
    outputStan <- runModelGrowthRate_STAN(dataForModel = dataForModel,
                                          dateList = dateList,
                                          parametersModel = parametersModel,
                                          inferenceSettings = inferenceSettings)
    output <- processSTANOutput(objectStan = outputStan,
                                parametersModel = parametersModel,
                                inferenceSettings = inferenceSettings,
                                saveSamples = saveSamples,
                                saveStanObject = saveModelObject)
  }else if(inferenceSettings$inferenceType == "MCMC-Int"){
    # TODO
    output <- NULL
  }else{
    stop("Invalid inferenceType")
  }
  
  return(output)
}


#' (Private) Title
#'
#' @param matrixSampleGP matrix [dayId x samples]
#' @param matrixSampleGPDerivative matrix [dayId x samples]
#' @param matrixSampleHyperparameters matrix [c("overdispersion", "precision", "logParam1", "logParam2") x samples]
#' @param matrixSampleNu matrix [dayId x samples]
#' @param matrixSampleRandomEffect matrix [dayId OR dayWeek OR other x samples]
#' @param matrixSampleIntercept vector [samples]
#' @param matrixTestData vector [dayId]
#' @param parametersModel 
#'
#' @return
#' @export
#'
#' @examples
computeSummaryPosteriors <- function(matrixSampleGP, matrixSampleGPDerivative, matrixSampleHyperparameters,
                                     matrixSampleNu, matrixSampleRandomEffect, matrixSampleIntercept,
                                     matrixTestData, parametersModel){
  
  numDays <- nrow(matrixSampleGP)
  sizeSample <- ncol(matrixSampleGP)
  internalSettings <- getInternalSettings()
  
  # ---------------------------------------------------- #
  #                POSTERIOR GROWTH RATE                 #
  # ---------------------------------------------------- #
  # Compute log (or inv. logit) of posterior of Gaussian process derivative in transformed space
  cat("Computing posterior of growth rate ... ")
  tempList <- getSamplesGrowthRate(matrixSampleGP, matrixSampleGPDerivative, parametersModel)
  tempDoubling <- abs(log(2)/tempList)
  posteriorGrowth <- data.table(dayId = 1:numDays,
                                mean = rowMeans(tempList),
                                sd = apply(tempList, 1, sd),
                                median = apply(tempList, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025 = apply(tempList, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975 = apply(tempList, 1, quantile, probs = 0.975, na.rm = T),
                                q0.25 = apply(tempList, 1, quantile, probs = 0.250, na.rm = T),
                                q0.75 = apply(tempList, 1, quantile, probs = 0.750, na.rm = T),
                                prob0 = apply(tempList, 1, function(x) sum(x >= 0)/sizeSample),
                                medianDoubT = apply(tempDoubling, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025DoubT = apply(tempDoubling, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975DoubT = apply(tempDoubling, 1, quantile, probs = 0.975, na.rm = T))
  
  # ---------------------------------------------------- #
  #                 POSTERIOR INCIDENCE                  #
  # ---------------------------------------------------- #
  # Compute posterior of Gaussian process in real space (incidence or positivity)
  cat("Computing posterior of incidence... ")
  
  # (this version is slow and applies only to INLA, and no constant)
  # depurar?
  #samplesGP <- objectInla$marginals.random$day[objectInla$nonBoundaryIndices]
  #if(parametersModel$linkType %in% c("NB")){
  #  transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(exp, x))
  #}else if(parametersModel$linkType == "BB"){
  #  transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(function(x) exp(x)/(1 + exp(x)), x))
  #}
  #posteriorTransfGP <- data.table(dayId = 1:numDays,
  #                                t(sapply(transformedSamples, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))),
  #                                t(sapply(samplesGP, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))))
  #setnames(posteriorTransfGP,
  #         c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
  #         c("median_transConsGP", "q0.025_transConsGP", "q0.975_transConsGP", "q0.25_transConsGP", "q0.75_transConsGP",
  #           "median_GP", "q0.025_GP", "q0.975_GP", "q0.25_GP", "q0.75_GP")) # NEW 19.09.2023
  
  # NEW 09.01.2023
  repSampleIntercept <- matrix(c(matrixSampleIntercept), nrow = numDays, ncol = ncol(matrixSampleGP), byrow = T)
  if(parametersModel$linkType %in% c("NB")){
    transformedSamples <- exp(matrixSampleGP + repSampleIntercept)
  }else if(parametersModel$linkType == "BB"){
    transformedSamples <- exp(matrixSampleGP + repSampleIntercept)/(1 + exp(matrixSampleGP + repSampleIntercept))
  }
  posteriorTransfGP <- data.table(dayId = 1:numDays,
                                  median_transConsGP = apply(transformedSamples, 1, quantile, probs = 0.5, na.rm = T),
                                  q0.025_transConsGP = apply(transformedSamples, 1, quantile, probs = 0.025, na.rm = T),
                                  q0.975_transConsGP = apply(transformedSamples, 1, quantile, probs = 0.975, na.rm = T),
                                  q0.25_transConsGP = apply(transformedSamples, 1, quantile, probs = 0.250, na.rm = T),
                                  q0.75_transConsGP = apply(transformedSamples, 1, quantile, probs = 0.750, na.rm = T),
                                  median_GP = apply(matrixSampleGP, 1, quantile, probs = 0.5, na.rm = T),
                                  q0.025_GP = apply(matrixSampleGP, 1, quantile, probs = 0.025, na.rm = T),
                                  q0.975_GP = apply(matrixSampleGP, 1, quantile, probs = 0.975, na.rm = T),
                                  q0.25_GP = apply(matrixSampleGP, 1, quantile, probs = 0.250, na.rm = T),
                                  q0.75_GP = apply(matrixSampleGP, 1, quantile, probs = 0.750, na.rm = T))
  
  # ---------------------------------------------------- #
  #                 POSTERIOR MODEL FIT                  #
  # ---------------------------------------------------- #
  # Compute the model posterior (as in R31.R)
  matrixSampleOverdisp <- matrix(matrixSampleHyperparameters[c("overdispersion"),], nrow = numDays, ncol = sizeSample, byrow = T)
  if(parametersModel$linkType == "NB"){
    samplesFit <- matrix(rnbinom(n = matrix(1, nrow = numDays, ncol = sizeSample),
                                 size = matrixSampleOverdisp,
                                 mu = exp(matrixSampleNu)),
                         nrow = numDays, ncol = sizeSample, byrow = F)
  }else{
    samplesMu <- exp(matrixSampleNu)/(1 + exp(matrixSampleNu))
    matrixAlpha <- samplesMu*(1 - matrixSampleOverdisp)/matrixSampleOverdisp
    matrixBeta <- (1 - samplesMu)*(1 - matrixSampleOverdisp)/matrixSampleOverdisp
    
    matrixP <- matrix(rbeta(matrix(1, nrow = numDays, ncol = sizeSample), shape1 = matrixAlpha, shape2 = matrixBeta),
                      nrow = numDays, ncol = sizeSample, byrow = F)
    matrixTests <- matrix(matrixTestData, nrow = numDays, ncol = sizeSample, byrow = F)
    samplesFit <- matrix(rbinom(n = matrix(1, nrow = numDays, ncol = sizeSample),
                                size = matrixTests,
                                prob = matrixP),
                         nrow = numDays, ncol = sizeSample, byrow = F)/matrixTests
  }
  posteriorTransfGP[order(dayId), ":="(medianFT = apply(samplesFit, 1, quantile, 0.5),
                                       q0.975FT = apply(samplesFit, 1, quantile, 0.975),
                                       q0.025FT = apply(samplesFit, 1, quantile, 0.025),
                                       q0.75FT = apply(samplesFit, 1, quantile, 0.75),
                                       q0.25FT = apply(samplesFit, 1, quantile, 0.25))]
  
  # ---------------------------------------------------- #
  #                 POSTERIOR DAY EFFECT                 #
  # ---------------------------------------------------- #
  if(parametersModel$randomEffect != "none"){
    # NEW 19.09.2023
    posteriorRandomEffect <- data.table(index = 1:nrow(matrixSampleRandomEffect),
                                        median = apply(matrixSampleRandomEffect, 1, quantile, probs = 0.5, na.rm = T),
                                        q0.025 = apply(matrixSampleRandomEffect, 1, quantile, probs = 0.025, na.rm = T),
                                        q0.975 = apply(matrixSampleRandomEffect, 1, quantile, probs = 0.975, na.rm = T),
                                        q0.25 = apply(matrixSampleRandomEffect, 1, quantile, probs = 0.250, na.rm = T),
                                        q0.75 = apply(matrixSampleRandomEffect, 1, quantile, probs = 0.750, na.rm = T))
    if(parametersModel$randomEffect == "weekday") posteriorRandomEffect[order(index), weekDay := internalSettings$levelsWeek]
  }else{
    posteriorRandomEffect <- NULL
  }
  
  # ---------------------------------------------------- #
  #                 POSTERIOR INTERCEPT                  #
  # ---------------------------------------------------- #
  posteriorIntercept <- data.table(median = quantile(matrixSampleIntercept, probs = 0.5, na.rm = T),
                                   q0.025 = quantile(matrixSampleIntercept, probs = 0.025, na.rm = T),
                                   q0.975 = quantile(matrixSampleIntercept, probs = 0.975, na.rm = T),
                                   q0.25 = quantile(matrixSampleIntercept, probs = 0.250, na.rm = T),
                                   q0.75 = quantile(matrixSampleIntercept, probs = 0.750, na.rm = T))
  
  # ---------------------------------------------------- #
  #              POSTERIOR HYPERPARAMETERS               #
  # ---------------------------------------------------- #
  # NEW 12.08.2024
  posteriorHyperparameters <- data.table(index = 1:nrow(matrixSampleHyperparameters),
                                         parameter = rownames(matrixSampleHyperparameters),
                                         median = apply(matrixSampleHyperparameters, 1, quantile, probs = 0.5, na.rm = T),
                                         q0.025 = apply(matrixSampleHyperparameters, 1, quantile, probs = 0.025, na.rm = T),
                                         q0.975 = apply(matrixSampleHyperparameters, 1, quantile, probs = 0.975, na.rm = T),
                                         q0.25 = apply(matrixSampleHyperparameters, 1, quantile, probs = 0.250, na.rm = T),
                                         q0.75 = apply(matrixSampleHyperparameters, 1, quantile, probs = 0.750, na.rm = T))
  
  return(list(posteriorGrowth = posteriorGrowth,
              posteriorTransfGP = posteriorTransfGP,
              posteriorRandomEffect = posteriorRandomEffect,
              posteriorIntercept = posteriorIntercept,
              posteriorHyperparameters = posteriorHyperparameters))
}


#' (Private) Title
#'
#' @param matrixSampleGP 
#'
#' @return
#' @export
#'
#' @examples
getSamplesGPDerivative_approx <- function(matrixSampleGP){
  # Compute approximate derivative using windows (+-3)
  # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
  numDays <- nrow(matrixSampleGP)
  matrixSampleGPDerivative <- rbind(matrix(NA, nrow = 1, ncol = ncol(matrixSampleGP)),
                             ((matrixSampleGP[3,] - matrixSampleGP[1,])/2),
                             ((matrixSampleGP[5,] - matrixSampleGP[1,])/4),
                             (matrixSampleGP[7:numDays,] - matrixSampleGP[1:(numDays - 7 + 1),])/6,
                             ((matrixSampleGP[numDays,] - matrixSampleGP[numDays - 5 + 1,])/4),
                             ((matrixSampleGP[numDays,] - matrixSampleGP[numDays - 3 + 1,])/2),
                             matrix(NA, nrow = 1, ncol = ncol(matrixSampleGP)))
  return(matrixSampleGPDerivative)
}

#' Title
#' (not to trust from INLA?... I don't trust INLA sampling of hyperparameters)
#'
#' @param matrixSampleGP days x samples
#' @param samplesHyperparam hyperparams X samples. sigma = sigma0 exp(logParam1), range = range0 exp(logParam2) ??
#' @param sigma0 
#' @param range0 
#'
#' @return
#' @export
#'
#' @examples
getSamplesGPDerivative_GPmatern12 <- function(matrixSampleGP, samplesHyperparam, sigma0, range0){
  numDays <- nrow(matrixSampleGP)
  sizeSample <- ncol(matrixSampleGP)
  
  # Compute distance matrix for ordered days
  distanceMatrix <- sapply(1:numDays, function(nd) abs(nd - (1:numDays)))
  auxRelativeDistanceMatrix <- outer(X = 1:numDays, Y = 1:numDays, FUN = "-")
  
  # Compute auxiliar vectors, with vGP = 3/2
  sig2Vector <- (sigma0*exp(samplesHyperparam["logParam1",]))^2
  kappaVector <- sqrt(12)/(range0*exp(samplesHyperparam["logParam2",]))
  
  # Loop per sample of (f1, ..., fn, log.tau, log.kappa)
  matrixSampleGPDerivative <- matrix(0, nrow = sizeSample, ncol = numDays)
  for(indexSample in 1:sizeSample){
    sig2Value <- sig2Vector[indexSample]
    kappaVal <- kappaVector[indexSample]
    
    expMatrix <- exp(-kappaVal*distanceMatrix)
    deltaMatrix <- sig2Value*(1 + kappaVal*distanceMatrix)*expMatrix # K(X,X)
    invDeltaMatrix <- chol2inv(chol(deltaMatrix)) # solve vs. chol2inv system.time(31700*system.time(solve(deltaMatrix))/60)
    fVector <- matrixSampleGP[, indexSample]
    
    # Compute derivative matrices of f:
    KXpX <- -sig2Value*kappaVal^2*auxRelativeDistanceMatrix*expMatrix # K(X*,X) - first partial derivative
    KXpXp <- sig2Value*kappaVal^2*expMatrix*(1 - kappaVal*distanceMatrix) # K(X*,X*) - second partial derivative
    #KXpXp <- sig2Value*kappaVal^2*diag(expMatrix)*(1 - kappaVal*diag(distanceMatrix)) # here we only compute the diagonal
    
    # Draw samples
    meanMVN <- KXpX%*%invDeltaMatrix%*%fVector
    iSample <- MASS::mvrnorm(n = 1, mu = KXpX%*%invDeltaMatrix%*%fVector, Sigma = KXpXp - KXpX%*%invDeltaMatrix%*%t(KXpX))
    #iSample <- sapply(1:numDays, function(x) rnorm(n = 1, mean = meanMVN[x], sd = sqrt( KXpXp[x] - KXpX[x,]%*%invDeltaMatrix%*%KXpX[x,] )))
    
    # Store
    matrixSampleGPDerivative[indexSample,] <- iSample
  }
  return(t(matrixSampleGPDerivative))
}

#' Title
#' 
#' Gives a matrix of dayId x samples of samples of the growth rate
#'
#' @param matrixSampleGP 
#' @param matrixSampleGPDerivative 
#' @param parametersModel 
#'
#' @return
#' @export
#'
#' @examples
getSamplesGrowthRate <- function(matrixSampleGP, matrixSampleGPDerivative, parametersModel){
  # Transform samples from GP derivative to GR
  if(parametersModel$linkType == "BB"){
    tempExpGP <- exp(matrixSampleGP)
    tempLogit <- tempExpGP/(1 + tempExpGP)
    samplesGR <- matrixSampleGPDerivative/(1 + tempLogit)
  }else if(parametersModel$linkType == "NB"){
    samplesGR <- matrixSampleGPDerivative
  }
  rownames(samplesGR) <- NULL
  return(samplesGR)
}


#' Title
#'
#' (22 sep 2023)
#' Gives a matrix of dayId x samples of samples of the random effect
#' TODO when do I use it? needed?
#'
#' @param outputModel 
#' @param parametersModel 
#'
#' @return
#' @export
#'
#' @examples
getSamplesRandomEffect <- function(outputModel, parametersModel){
  internalSettings <- getInternalSettings()
  if(parametersModel$randomEffect == "weekday"){
    predictionWeekday <- outputModel$dateList$dateTable[order(dayId), weekdays(date)]
    matrixRandomEffect <- outputModel$matrixSampleRandomEffect[match(predictionWeekday, internalSettings$levelsWeek),]
  }else if(parametersModel$randomEffect == "all"){
    # TODO
    matrixRandomEffect <- matrix(0, nrow = outputModel$dateList$numDays, ncol = outputModel$sizeSample)
  }else{
    matrixRandomEffect <- matrix(0, nrow = outputModel$dateList$numDays, ncol = outputModel$sizeSample)
  }
  return(matrixRandomEffect)
}

#' Title
#'
#' (22 sep 2023)
#' Gives a matrix of dayId x samples of samples of the random effect
#' (called getSamplesGPConsTrans)
#' TODO when do I use it? needed?
#'
#' @param outputModel 
#' @param parametersModel 
#'
#' @return
#' @export
#'
#' @examples
getSamplesGPPlusConstantTransformed <- function(outputModel, parametersModel){
  repSampleIntercept <- matrix(c(outputModel$matrixSampleIntercept), nrow = outputModel$dateList$numDays, ncol = outputModel$sizeSample, byrow = T)
  if(parametersModel$linkType %in% c("NB")){
    matrixGPConsTrans <- exp(outputModel$matrixSampleGP + repSampleIntercept)
  }else if(parametersModel$linkType == "BB"){
    matrixGPConsTrans <- exp(outputModel$matrixSampleGP + repSampleIntercept)/(1 + exp(outputModel$matrixSampleGP + repSampleIntercept))
  }
  return(matrixGPConsTrans)
}
