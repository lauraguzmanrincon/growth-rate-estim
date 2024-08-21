
# Compute projection ----

# TODO write function to create empty projection object (to compare other non-GP models)
# OR create simpler object that accommodates all types of projections (non-GP ones)

# Definition
# currentDate: last day of observations (<= maxDate in model)
# sizePredicition: days into the future
# daysForPrediction:
# testingVector: tests in the future days to project

#' Title
#' 
#' Wrapper for projection using different GPs
#'
#' @param outputModel 
#' @param currentDate Last observed day to compute the projection. Recommended to be outputModel$dateList$maxDay
#' @param sizePrediction 
#' @param daysForPrediction 
#'
#' @return
#' @export
#'
#' @examples
getProjectionGP <- function(outputModel, currentDate, sizePrediction, daysForPrediction){
  
  if(outputModel$inferenceSettings$inferenceType == "LA-INLA"){
    # Only matern12 is available for INLA
    output <- getProjectionGP_matern12(outputModel = outputModel,
                                       currentDate = currentDate,
                                       sizePrediction = sizePrediction,
                                       daysForPrediction = daysForPrediction)
  }else if(outputModel$inferenceSettings$inferenceType == "MCMC-STAN"){
    covarianceFunction <- rationalQuadraticMatrixFn # Only one available for now
    output <- getProjectionGP_covInput(outputModel = outputModel,
                                           currentDate = currentDate,
                                           sizePrediction = sizePrediction,
                                           daysForPrediction = daysForPrediction,
                                           covarianceFunction = covarianceFunction)
  }else if(outputModel$inferenceSettings$inferenceType == "MCMC-Int"){
    # TODO
    output <- NULL
  }else{
    stop("Invalid inferenceType")
  }
  
  return(output)
}


#' Title
#' 
#'
#' @param outputModel 
#' @param currentDate Last observed day to compute the projection. Recommended to be outputModel$dateList$maxDay
#' @param sizePrediction 
#' @param daysForPrediction 
#'
#' @return
#' @export
#'
#' @examples
getProjectionGP_matern12 <- function(outputModel, currentDate, sizePrediction, daysForPrediction){
  
  if(outputModel$dateList$numDays < daysForPrediction) warning("daysForPrediction is longer than available observations.")
  
  # ---------------------------------------------------- #
  #                       SET-UP                         #
  # ---------------------------------------------------- #
  
  # Get basis configuration
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  sizeSample <- outputModel$numSamples
  numDaysPast <- min(outputModel$dateList$numDays, daysForPrediction)
  matrixSampleGP <- outputModel$matrixSampleGP[(currentDayId - daysForPrediction + 1):currentDayId,]
  
  # Compute distance matrix for ordered days
  pastInterval <- 1:numDaysPast
  predInterval <- (numDaysPast + 1):(numDaysPast + sizePrediction)
  relativeDistanceMatrixObsObs <- outer(X = pastInterval, Y = pastInterval, FUN = "-")
  relativeDistanceMatrixPredObs <- outer(X = predInterval, Y = pastInterval, FUN = "-")
  relativeDistanceMatrixPredPred <- outer(X = predInterval, Y = predInterval, FUN = "-")
  
  # Get basis parameters
  vGP <- 2 - 1/2
  sigma0 <- outputModel$parametersModel$GPcovarianceList$theta0[1]
  range0 <- 2*outputModel$parametersModel$GPcovarianceList$theta0[2]
  
  # Compute auxiliar vectors
  sig2Vector <- (sigma0*exp(outputModel$matrixSampleHyperparameters["logParam1",]))^2
  kappaVector <- sqrt(8*vGP)/(range0*exp(outputModel$matrixSampleHyperparameters["logParam2",]))
  
  # ---------------------------------------------------- #
  #                  PRODUCE SAMPLES                     #
  # ---------------------------------------------------- #
  
  # Loop per sample of (f1, ..., fn, log.tau, log.kappa)
  sampleProjections <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  sampleDerivatives <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  errorIterations <- rep(0, sizeSample)
  for(indexSample in 1:sizeSample){
    sig2Value <- sig2Vector[indexSample]
    kappaVal <- kappaVector[indexSample]
    
    expMatrix <- exp(-kappaVal*abs(relativeDistanceMatrixObsObs))
    deltaMatrix <- sig2Value*(1 + kappaVal*abs(relativeDistanceMatrixObsObs))*expMatrix # K(X,X)
    invDeltaMatrix <- chol2inv(chol(deltaMatrix)) # solve vs. chol2inv system.time(31700*system.time(solve(deltaMatrix))/60)
    fVector <- matrixSampleGP[, indexSample]
    
    expMatrixPredObs <- exp(-kappaVal*abs(relativeDistanceMatrixPredObs))
    expMatrixPredPred <- exp(-kappaVal*abs(relativeDistanceMatrixPredPred))
    
    # Compute matrices for Gaussian conditionals for GP prediction:
    KPredObs <- rbind(sig2Value*(1 + kappaVal*abs(relativeDistanceMatrixPredObs))*expMatrixPredObs, # K
                      -sig2Value*kappaVal^2*relativeDistanceMatrixPredObs*expMatrixPredObs) # K'
    KPredPred <- rbind(cbind(sig2Value*(1 + kappaVal*abs(relativeDistanceMatrixPredPred))*expMatrixPredPred,# K
                             -sig2Value*kappaVal^2*relativeDistanceMatrixPredPred*expMatrixPredPred),
                       cbind(t(-sig2Value*kappaVal^2*relativeDistanceMatrixPredPred*expMatrixPredPred), # symmetric?
                             sig2Value*kappaVal^2*expMatrixPredPred*(1 - kappaVal*abs(relativeDistanceMatrixPredPred))))
    # TODO is this slow?
    
    # Draw samples
    meanMVN <- KPredObs%*%invDeltaMatrix%*%fVector
    sigmaMVN <- KPredPred - KPredObs%*%invDeltaMatrix%*%t(KPredObs)
    ###
    iSampleGP <- MASS::mvrnorm(n = 1, mu = meanMVN[1:sizePrediction], Sigma = sigmaMVN[1:sizePrediction,1:sizePrediction])
    iSampleGR <- MASS::mvrnorm(n = 1, mu = meanMVN[(sizePrediction + 1):(2*sizePrediction)],
                               Sigma = sigmaMVN[(sizePrediction + 1):(2*sizePrediction),(sizePrediction + 1):(2*sizePrediction)])
    tryCatch({
      iSampleGP <- MASS::mvrnorm(n = 1, mu = meanMVN[1:sizePrediction], Sigma = sigmaMVN[1:sizePrediction,1:sizePrediction])
      iSampleGR <- MASS::mvrnorm(n = 1, mu = meanMVN[(sizePrediction + 1):(2*sizePrediction)],
                                 Sigma = sigmaMVN[(sizePrediction + 1):(2*sizePrediction),(sizePrediction + 1):(2*sizePrediction)])
    }, error = function(e) {
      errorIterations[indexSample] <- 1
      print(indexSample)
      iSample <- rep(NA, 2*sizePrediction)
    })
    
    # Save
    sampleProjections[,indexSample] <- iSampleGP
    sampleDerivatives[,indexSample] <- iSampleGR
  }
  
  if(sum(errorIterations) > 0) warning("Sigma not positive definite for ", sum(errorIterations), " samples out of ", sizeSample)
  
  # ---------------------------------------------------- #
  #                SAMPLES DERIVATIVE/GR                 #
  # ---------------------------------------------------- #
  
  # sampleDerivatives already computed from the GP above
  # We do it again if derivative is calculated from approximation
  if(outputModel$inferenceSettings$derivativeFromGP == F){
    sampleDerivatives <- getSamplesGPDerivative_approx(matrixSampleGP = rbind(matrixSampleGP, sampleProjections))[predInterval,]
  }
  
  # Transform derivatives to GR
  samplesGR <- getSamplesGrowthRate(matrixSampleGP = sampleProjections,
                                    matrixSampleGPDerivative = sampleDerivatives,
                                    parametersModel = outputModel$parametersModel)
  
  # ---------------------------------------------------- #
  #                    SAMPLES BOUNDARY                  #
  # ---------------------------------------------------- #
  
  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
  if(outputModel$inferenceSettings$derivativeFromGP == F){
    projectionGR_onBoundary <- rep(0, sizeSample) #projectionGR_onBoundary <- c(tempDerivative[daysForPrediction,])
    # TODO compute GR on boundary (derivative + transformation)
  }else{
    projectionGR_onBoundary <- NA
  }
  
  return(list(projectionGP = sampleProjections, projectionGR = samplesGR,
              projectionGR_onBoundary = projectionGR_onBoundary,
              currentDate = outputModel$dateList$dateTable[dayId == outputModel$dateList$maxDay, date],
              sizePrediction = sizePrediction))
}

#' Title
#'
#' @param outputModel 
#' @param currentDate Last observed day to compute the projection. Recommended to be outputModel$dateList$maxDay
#' @param sizePrediction 
#' @param daysForPrediction 
#' @param covarianceFunction Function with three inputs: params, inputMatrix, sqrInputMatrix
#'
#' @return
#' @export
#'
#' @examples
getProjectionGP_covInput <- function(outputModel, currentDate, sizePrediction, daysForPrediction, covarianceFunction){
  
  if(outputModel$dateList$numDays < daysForPrediction) warning("daysForPrediction is longer than available observations.")
  
  # ---------------------------------------------------- #
  #                       SET-UP                         #
  # ---------------------------------------------------- #
  
  # Get basis configuration
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  sizeSample <- outputModel$numSamples
  numDaysPast <- min(outputModel$dateList$numDays, daysForPrediction)
  matrixSampleGP <- outputModel$matrixSampleGP[(currentDayId - daysForPrediction + 1):currentDayId,]
  
  # Compute distance matrix for ordered days
  pastInterval <- 1:numDaysPast
  predInterval <- (numDaysPast + 1):(numDaysPast + sizePrediction)
  relativeDistanceMatrixObsObs <- outer(X = pastInterval, Y = pastInterval, FUN = "-")
  relativeDistanceMatrixPredObs <- outer(X = predInterval, Y = pastInterval, FUN = "-")
  relativeDistanceMatrixPredPred <- outer(X = predInterval, Y = predInterval, FUN = "-")
  sqrDistanceMatrixObsObs <- relativeDistanceMatrixObsObs^2
  sqrDistanceMatrixPredObs <- relativeDistanceMatrixPredObs^2
  sqrDistanceMatrixPredPred <- relativeDistanceMatrixPredPred^2
  
  numParameters <- nrow(outputModel$matrixSampleHyperparameters) - 2
  
  # ---------------------------------------------------- #
  #                  PRODUCE SAMPLES                     #
  # ---------------------------------------------------- #
  
  sampleProjections <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  #sampleDerivatives <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  errorIterations <- rep(0, sizeSample)
  for(indexSample in 1:sizeSample){
    paramSample <- sapply(X = paste("logParam", 1:numParameters, sep = ""), FUN = function(colname) exp(outputModel$matrixSampleHyperparameters[colname, indexSample]))
    
    # Can we done in one but maybe saves space/time?
    covarianceObsObs <- covarianceFunction(params = paramSample, inputMatrix = relativeDistanceMatrixObsObs, sqrInputMatrix = sqrDistanceMatrixObsObs) + 0.00001*diag(numDaysPast)
    covariancePredObs <- covarianceFunction(params = paramSample, inputMatrix = relativeDistanceMatrixObsObs, sqrInputMatrix = sqrDistanceMatrixPredObs)
    covariancePredPred <- covarianceFunction(params = paramSample, inputMatrix = relativeDistanceMatrixObsObs, sqrInputMatrix = sqrDistanceMatrixPredPred)
    observation <- matrixSampleGP[, indexSample]
    
    importance <- covariancePredObs%*%chol2inv(chol(covarianceObsObs))
    meanGP <- importance%*%observation
    covGP <- covariancePredPred - importance%*%t(covariancePredObs)
    tryCatch({
      iSampleGP <- MASS::mvrnorm(n = 1, mu = meanGP, Sigma = covGP) # plot(c(observation, samplesGP))
    }, error = function(e) {
      errorIterations[indexSample] <- 1
      iSampleGP <- rep(NA, 2*sizePrediction)
    })
    
    sampleProjections[,indexSample] <- iSampleGP
    #sampleDerivatives[,indexSample] <- iSampleGR
  }
  
  if(sum(errorIterations) > 0) warning("Sigma not positive definite for ", sum(errorIterations), " samples out of ", sizeSample)
  
  # ---------------------------------------------------- #
  #                SAMPLES DERIVATIVE/GR                 #
  # ---------------------------------------------------- #
  
  # We are not sampling from GR in this case, but using finite differences with the sampled GP
  sampleDerivatives <- getSamplesGPDerivative_approx(matrixSampleGP = rbind(matrixSampleGP, sampleProjections))[predInterval,]
  
  # Transform derivatives to GR
  samplesGR <- getSamplesGrowthRate(matrixSampleGP = sampleProjections,
                                    matrixSampleGPDerivative = sampleDerivatives,
                                    parametersModel = outputModel$parametersModel)
  
  # ---------------------------------------------------- #
  #                    SAMPLES BOUNDARY                  #
  # ---------------------------------------------------- #
  
  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
  if(outputModel$inferenceSettings$derivativeFromGP == F){
    projectionGR_onBoundary <- NA
    #projectionGR_onBoundary <- rep(0, sizeSample) #projectionGR_onBoundary <- c(tempDerivative[daysForPrediction,])
    # TODO compute GR on boundary (derivative + transformation)
  }else{
    projectionGR_onBoundary <- NA
  }
  
  return(list(projectionGP = sampleProjections, projectionGR = samplesGR,
              projectionGR_onBoundary = projectionGR_onBoundary,
              currentDate = outputModel$dateList$dateTable[dayId == outputModel$dateList$maxDay, date],
              sizePrediction = sizePrediction))
}


#' Title
#' 
#' Function for projections that provide a GP and a GR !!!
#' Creates a table with samples for relevant parameters plus median, 95% CI for model posterior
#' testingVector: vector of length projectionSamplesGP$sizePrediction with testing (for proportions model)
#' TODO this is if we give GP and GR as inputs... maybe we'll need something more general?
#'
#' @param projectionSamplesGP 
#' @param outputModel 
#' @param nameProjection 
#' @param testingVector 
#'
#' @return TableProjection object
#' @export
#'
#' @examples
createTableProjection <- function(projectionSamplesGP, outputModel, nameProjection, testingVector = NULL){
  # Check input
  if(outputModel$parametersModel$linkType == "BB" & is.null(testingVector)) stop("testingVector is required for the proportions model")
  
  internalConstants <- getInternalSettings()
  
  # Auxiliar variables
  sizeSample <- outputModel$numSamples
  sizePrediction <- nrow(projectionSamplesGP$projectionGP)
  currentDayId <- outputModel$dateList$dateTable[date == projectionSamplesGP$currentDate, dayId]
  
  # ---------------------------------------------------- #
  #                CREATE SAMPLES TABLE                  #
  # ---------------------------------------------------- #
  
  # Create matrices of parameters
  # TODO same as in computePosteriors()... Unify!
  if(outputModel$parametersModel$randomEffect == "weekday"){
    predictionWeekday <- weekdays(projectionSamplesGP$currentDate + (1:sizePrediction))
    matrixRandomEffect <- outputModel$matrixSampleRandomEffect[match(predictionWeekday, internalConstants$levelsWeek),]
  }else if(outputModel$parametersModel$randomEffect == "all"){
    matrixRandomEffect <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  }else{
    matrixRandomEffect <- matrix(0, nrow = sizePrediction, ncol = sizeSample)
  }
  matrixSampleOverdisp <- outputModel$matrixSampleHyper[c("overdispersion"),]
  matrixIntercept <- matrix(outputModel$matrixSampleIntercept,
                            nrow = sizePrediction, ncol = sizeSample, byrow = T)
  etaSample <- matrixIntercept + projectionSamplesGP$projectionGP + matrixRandomEffect
  rhoSample <- matrix(matrixSampleOverdisp, nrow = sizePrediction, ncol = sizeSample, byrow = T)
  
  # Create projection table with samples of relevant model parameters
  if(outputModel$parametersModel$linkType == "NB"){
    mupSample <- exp(etaSample)
    gpTransSample <- exp(matrixIntercept + projectionSamplesGP$projectionGP)
    tableSamples <- data.table(lastDayId = currentDayId,
                               dayId = currentDayId + (1:sizePrediction),
                               sample = rep(1:sizeSample, each = sizePrediction),
                               gp = c(projectionSamplesGP$projectionGP),
                               gr = c(projectionSamplesGP$projectionGR),
                               gpConsTrans = c(gpTransSample),
                               eta = c(etaSample),
                               mu = c(mupSample),
                               rho = c(rhoSample))
  }else{
    mupSample <- exp(etaSample)/(1 + exp(etaSample))
    gpTransSample <- exp(matrixIntercept + projectionSamplesGP$projectionGP)/
      (1 + exp(matrixIntercept + projectionSamplesGP$projectionGP))
    
    # Recall: mup = a/(a+b), rho = 1/(a+b+1). See Other_files/.../NLADoc_betabinomial.pdf or inla.doc("betabinomial")
    aS <- mupSample*(1 - rhoSample)/rhoSample
    bS <- (mupSample*rhoSample - mupSample - rhoSample + 1)/rhoSample
    pS <- matrix(mapply(function(a,b) rbeta(1, a, b), aS, bS), nrow = sizePrediction, ncol = sizeSample, byrow = F)
    
    tableSamples <- data.table(lastDayId = currentDayId,
                               dayId = currentDayId + (1:sizePrediction),
                               sample = rep(1:sizeSample, each = sizePrediction),
                               gp = c(projectionSamplesGP$projectionGP),
                               gr = c(projectionSamplesGP$projectionGR),
                               gpConsTrans = c(gpTransSample),
                               eta = c(etaSample),
                               mu = c(mupSample),
                               rho = c(rhoSample),
                               alpha = c(aS),
                               beta = c(bS),
                               pp = c(pS))
  }
  
  # Add dates to the samples table
  setkey(tableSamples, lastDayId)
  setkey(outputModel$dateList$dateTable, dayId)
  tableSamples[outputModel$dateList$dateTable, dateLastDay := date]
  tableSamples[, date := dateLastDay + dayId - lastDayId]
  
  # Add tests to the samples table (relevant for BB)
  if(outputModel$parametersModel$linkType == "NB"){
    tableSamples[, numberTest := NA]
  }else{
    tableSamples[order(dayId), numberTest := testingVector[dayId - lastDayId]]
  }
  
  # Compute median and 95% CI of model posterior and add it to samples table
  if(outputModel$parametersModel$linkType == "NB"){
    tableSamples[, yS := rnbinom(n = nrow(tableSamples), size = rho, mu = mu)]
  }else{
    tableSamples[, yS := rbinom(n = nrow(tableSamples), size = numberTest, prob = pp)] # eg. rbinom(n = 21, size = c(2,30), c(0.1,0.9))
  }
  
  # ---------------------------------------------------- #
  #                CREATE SUMMARY TABLE                  #
  # ---------------------------------------------------- #
  # Create tableProjections with quantiles of every relevant parameter
  tableProjections <- tableSamples[, .(gp_q025 = quantile(gp, 0.025),
                                       gp_q25 = quantile(gp, 0.25),
                                       gp_median = quantile(gp, 0.5),
                                       gp_q75 = quantile(gp, 0.75),
                                       gp_q975 = quantile(gp, 0.975),
                                       #
                                       gr_q025 = quantile(gr, 0.025, na.rm = T),
                                       gr_q25 = quantile(gr, 0.25, na.rm = T),
                                       gr_median = quantile(gr, 0.5, na.rm = T),
                                       gr_q75 = quantile(gr, 0.75, na.rm = T),
                                       gr_q975 = quantile(gr, 0.975, na.rm = T),
                                       #
                                       gpConsTrans_q025 = quantile(gpConsTrans, 0.025),
                                       gpConsTrans_q25 = quantile(gpConsTrans, 0.25),
                                       gpConsTrans_median = quantile(gpConsTrans, 0.5),
                                       gpConsTrans_q75 = quantile(gpConsTrans, 0.75),
                                       gpConsTrans_q975 = quantile(gpConsTrans, 0.975),
                                       #
                                       eta_q025 = quantile(eta, 0.025),
                                       eta_q25 = quantile(eta, 0.25),
                                       eta_median = quantile(eta, 0.5),
                                       eta_q75 = quantile(eta, 0.75),
                                       eta_q975 = quantile(eta, 0.975),
                                       #
                                       mu_q025 = quantile(mu, 0.025),
                                       mu_q25 = quantile(mu, 0.25),
                                       mu_median = quantile(mu, 0.5),
                                       mu_q75 = quantile(mu, 0.75),
                                       mu_q975 = quantile(mu, 0.975),
                                       #
                                       yS_q025 = quantile(yS, 0.025, na.rm = T),
                                       yS_q25 = quantile(yS, 0.25, na.rm = T),
                                       yS_median = quantile(yS, 0.5, na.rm = T),
                                       yS_q75 = quantile(yS, 0.75, na.rm = T),
                                       yS_q975 = quantile(yS, 0.975, na.rm = T)),
                                   .(dayId, date, lastDayId, dateLastDay, numberTest)]
  
  # Add GR on last day of observation if not provided by model (when GP from finite differences)
  if(outputModel$inferenceSettings$derivativeFromGP == F){
    tableSamples_boundary <- data.table(lastDayId = currentDayId,
                                        dayId = currentDayId,
                                        sample = 1:sizeSample,
                                        gr = projectionSamplesGP$projectionGR_onBoundary)
    tableProjections_boundary <- tableSamples_boundary[, .(gr_q025 = quantile(gr, 0.025, na.rm = T),
                                                           gr_q25 = quantile(gr, 0.25, na.rm = T),
                                                           gr_median = quantile(gr, 0.5, na.rm = T),
                                                           gr_q75 = quantile(gr, 0.75, na.rm = T),
                                                           gr_q975 = quantile(gr, 0.975, na.rm = T)),
                                                       .(dayId, lastDayId)]
    # Add extras to the projection table
    setkey(tableProjections_boundary, lastDayId)
    setkey(outputModel$dateList$dateTable, dayId)
    tableProjections_boundary[outputModel$dateList$dateTable, dateLastDay := date]
    tableProjections_boundary[, date := dateLastDay + dayId - lastDayId]
    
    setkey(tableProjections_boundary, dayId)
    setkey(outputModel$posteriorTransfGP, dayId)
    tableProjections_boundary[outputModel$posteriorTransfGP, ":="(numberTest = i.numberTest,
                                                     gp_q025 = i.q0.025_GP,
                                                     gp_q25 = i.q0.25_GP,
                                                     gp_median = i.median_GP,
                                                     gp_q75 = i.q0.75_GP,
                                                     gp_q975 = i.q0.975_GP,
                                                     gpConsTrans_q025 = i.q0.025_transConsGP,
                                                     gpConsTrans_q25 = i.q0.25_transConsGP,
                                                     gpConsTrans_median = i.median_transConsGP,
                                                     gpConsTrans_q75 = i.q0.75_transConsGP,
                                                     gpConsTrans_q975 = i.q0.975_transConsGP)]
    # TODO Double-check the above correct!!!!
  }else{
    tableSamples_boundary <- NA
    tableProjections_boundary <- NA
  }
    
  return(list(tableProjections = tableProjections,
              tableSamples = tableSamples,
              tableProjections_boundary = tableProjections_boundary,
              tableSamples_boundary = tableSamples_boundary,
              currentDate = projectionSamplesGP$currentDate,
              sizePrediction = sizePrediction,
              nameProjection = nameProjection))
}



#' Title
#' 
#' Creates projection under null model (future observations are in 95% CI of normal around latest observations/proportion avg and sd = sdPast)
#' daysForPrediction: days before currentDate defined as 'latest observations'/'tests'
#' TODO caveat: quantiles from samples could be calculated with the qnorm function instead of using tableSample, but it's more practical to provide samples
#'
#' @param outputModel 
#' @param currentDate 
#' @param sizePrediction 
#' @param daysForPrediction 
#' @param nameProjection 
#' @param testingVector 
#'
#' @return TableProjection object
#' @export
#'
#' @examples
createTableProjectionNull <- function(outputModel, currentDate,
                                      sizePrediction, daysForPrediction = 7, nameProjection, testingVector = NULL){
  # Check input
  if(outputModel$parametersModel$linkType == "BB" & is.null(testingVector)) stop("testingVector is required for the proportions model")
  if(is.null(testingVector)) testingVector <- rep(NA, sizePrediction)
  
  sizeSample <- outputModel$numSamples
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  casesBeforePrediction <- outputModel$dataForModel[date > currentDate - daysForPrediction & date <= currentDate][order(date), positiveResults]
  testsBeforePrediction <- outputModel$dataForModel[date > currentDate - daysForPrediction & date <= currentDate][order(date), numberTest]
  
  # Table projection output
  overviewNull <- data.table(dayId = currentDayId + 1:sizePrediction,
                             date = currentDate + 1:sizePrediction,
                             lastDayId = currentDayId,
                             dateLastDay = currentDate,
                             numberTest = testingVector)
  overviewNull[, ":="(gp_q025 = as.numeric(NA),
                      gp_q25 = as.numeric(NA),
                      gp_median = as.numeric(NA),
                      gp_q75 = as.numeric(NA),
                      gp_q975 = as.numeric(NA),
                      #
                      gr_q025 = 0,
                      gr_q25 = 0,
                      gr_median = 0,
                      gr_q75 = 0,
                      gr_q975 = 0,
                      #
                      gpConsTrans_q025 = as.numeric(NA),
                      gpConsTrans_q25 = as.numeric(NA),
                      gpConsTrans_median = as.numeric(NA),
                      gpConsTrans_q75 = as.numeric(NA),
                      gpConsTrans_q975 = as.numeric(NA)
                      )]
  
  # Compute yS
  if(outputModel$parametersModel$linkType == "NB"){
    #alphaIS <- 0.05
    #positiveResultsLastDay <- casesBeforePrediction[daysForPrediction] # around current observation
    positiveResultsLastDay <- round(mean(casesBeforePrediction)) # around latest observations avg
    sdPast <- mean(abs(diff(casesBeforePrediction)))
    
    overviewNull[, ":="(yS_q025 = pmax(0, qnorm(0.025, mean = positiveResultsLastDay, sd = sdPast)),
                        yS_q25 = pmax(0, qnorm(0.25, mean = positiveResultsLastDay, sd = sdPast)),
                        yS_median = positiveResultsLastDay,
                        yS_q75 = pmax(0, qnorm(0.75, mean = positiveResultsLastDay, sd = sdPast)),
                        yS_q975 = pmax(0, qnorm(0.975, mean = positiveResultsLastDay, sd = sdPast)))]
  }else{
    proportionsLastDay <- mean(casesBeforePrediction/testsBeforePrediction) # around latest proportion avg
    sdPast <- mean(abs(diff(casesBeforePrediction/testsBeforePrediction)))
    overviewNull[, ":="(yS_q025 = testingVector*pmax(0, pmin(1, qnorm(0.025, mean = proportionsLastDay, sd = sdPast))),
                        yS_q25 = testingVector*pmax(0, pmin(1, qnorm(0.25, mean = proportionsLastDay, sd = sdPast))),
                        yS_median = testingVector*proportionsLastDay,
                        yS_q75 = testingVector*pmax(0, pmin(1, qnorm(0.75, mean = proportionsLastDay, sd = sdPast))),
                        yS_q975 = testingVector*pmax(0, pmin(1, qnorm(0.975, mean = proportionsLastDay, sd = sdPast))))]
  }
  
  # Copy y into gpConsTrans for visualisation purposes
  overviewNull[, ":="(gpConsTrans_median = yS_median,
                      gpConsTrans_q025 = yS_q025,
                      gpConsTrans_q975 = yS_q975,
                      gpConsTrans_q25 = yS_q25,
                      gpConsTrans_q75 = yS_q75)]
  
  # Table samples output - artificial samples that can be useful when using scoring functions (see caveat in description)
  tableSamples <- CJ(dayId = currentDayId + 1:sizePrediction, sample = 1:sizeSample)
  setkey(tableSamples, dayId)
  setkey(overviewNull, dayId)
  tableSamples[overviewNull, ":="(lastDayId = lastDayId, date = date, dateLastDay = dateLastDay, numberTest = numberTest)]
  
  # Produce samples
  if(outputModel$parametersModel$linkType == "NB"){
    # TODO this might produce a warning if casesBeforePrediction are unknown
    tableSamples[, yS := rnorm(n = sizePrediction*sizeSample, mean = positiveResultsLastDay, sd = sdPast)]
  }else{
    tableSamples[, yS := numberTest*pmax(0, pmin(1, rnorm(n = sizePrediction*sizeSample, mean = proportionsLastDay, sd = sdPast)))]
  }
  
  return(list(tableProjections = overviewNull,
              tableSamples = tableSamples,
              currentDate = currentDate,
              sizePrediction = sizePrediction,
              nameProjection = nameProjection))
}

# Scoring ----

#' Title
#'
#' Function for comparing scores
#' projectionsList: list of listProjection objects
#' scoringList: list of tableScore objects
#' positivesVector: positives observed in the future
#' TODO add list of optional arguments?
#' 
#' @param projectionsList List of TableProjection objects (e.g. output of createTableProjection)
#' @param scoringList List of TableScore objects
#' @param positivesVector Numeric verctor of "true positives" observed in the future
#'
#' @return
#' @export
#'
#' @examples
compareProjections <- function(projectionsList, scoringList, positivesVector){
  outputList <- vector("list", length(projectionsList)*length(projectionsList))
  for(p in 1:length(projectionsList)){
    for(s in 1:length(scoringList)){
      scoreTable <- scoringList[[s]](listProjection = projectionsList[[p]], positivesVector = positivesVector)
      scoreTable[, nameProjection := projectionsList[[p]]$nameProjection]
      outputList[[(p - 1)*length(scoringList) + s]] <- scoreTable
    }
  }
  return(do.call("rbind", outputList))
}

#' Title
#'
#' @param listProjection TableProjection object (output of createTableProjection)
#' @param positivesVector Numeric verctor of "true positives" observed in the future
#' @param alphaIS 
#'
#' @return TableScore object
#' @export
#'
#' @examples
scoreIS <- function(listProjection, positivesVector, alphaIS = 0.5){
  if(length(positivesVector) != listProjection$sizePrediction) stop("positivesVector should have length equal to the prediction size in listProjection")
  
  tableScore <- listProjection$tableSamples[, .(l = quantile(yS, probs = alphaIS/2, na.rm = T),
                                                u = quantile(yS, probs = 1 - (alphaIS/2), na.rm = T)),
                                            .(lastDayId, dayId, dateLastDay, date)]
  tableScore[order(dayId), positiveResults := positivesVector]
  tableScore[, scoreIS := (u - l) + 2*(l - positiveResults)*(positiveResults < l)/alphaIS + 2*(positiveResults - u)*(positiveResults > u)/alphaIS]
  tableScore[, ":="(score = log(scoreIS),
                    scoreName = "IS score")]
  return(tableScore[, .(lastDayId, dayId, dateLastDay, date, positiveResults, score, scoreName)])
}

#' Title
#'
#' "A lower DSS indicates a better forecast."
#' 
#' @param listProjection TableProjection object (e.g. output of createTableProjection)
#' @param positivesVector Numeric verctor of "true positives" observed in the future
#'
#' @return TableScore object
#' @export
#'
#' @examples
scoreDSS <- function(listProjection, positivesVector){
  if(length(positivesVector) != listProjection$sizePrediction) stop("positivesVector should have length equal to the prediction size in listProjection")
  
  tableScore <- listProjection$tableSamples[, .(mean = mean(yS, na.rm = T),
                                                sd = sd(yS, na.rm = T)),
                                            .(lastDayId, dayId, dateLastDay, date)]
  tableScore[order(dayId), positiveResults := positivesVector]
  tableScore[, ":="(score = ((positiveResults - mean)/sd)^2 + 2+log(sd),
                    scoreName = "DSS score")]
  return(tableScore[, .(lastDayId, dayId, dateLastDay, date, positiveResults, score, scoreName)])
}







