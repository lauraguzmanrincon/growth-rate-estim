# From W01.R

# Definition
# currentDate: last day of observations (<= maxDate in model)
# sizePredicition: days into the future
# daysForPrediction:
# testingVector: tests in the future days to project

#' currentDayId: Last observed day to compute the projection. Recommended to be outputModel$dateList$maxDay
#' DEPRICATED 22.09.2023
#getProjectionLinear <- function(parametersModel, outputModel, currentDate, sizePrediction, daysForPrediction){
#  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
#  
#  # Prediction - linear
#  lastDay <- currentDayId - 1
#  xVal <- (currentDayId - daysForPrediction):(currentDayId - 1)
#  xxVal <- currentDayId + (0:sizePrediction) # + (0:(sizePrediction - 1)) # *** !!!!!!!!! ????
#  yMat <- outputModel$matrixSampleDays[xVal,]
#  dMat <- outputModel$sampleDerivatives[xVal,]
#  xxValDF <- data.frame(xVal = xxVal)
#  length(intersect(xVal, xxVal)) == 0
#  
#  # Create matrix with projection
#  tempProjectionGPLinear <- matrix(0, sizePrediction + 1, sizeSample) # ***
#  tempProjectionGRLinear <- matrix(0, sizePrediction + 1, sizeSample) # ***
#  for(jjj in 1:sizeSample){
#    fitProj <- lm(dMat[,jjj] ~ poly(xVal, degree = 1, raw = TRUE))
#    tempProjectionGPLinear[,jjj] <- 0.5*fitProj$coefficient[2]*(xxVal^2 - lastDay^2) + fitProj$coefficient[1]*(xxVal - lastDay) +
#      outputModel$matrixSampleDays[lastDay, jjj] # see notes 19.05.2021
#    tempProjectionGRLinear[,jjj] <- predict(fitProj, newdata = xxValDF)
#    # e.g.  plot(c(dMat[,jjj], tempProjectionGRLinear[,jjj]))
#    #       plot(c(yMat[,jjj], tempProjectionGPLinear[,jjj])) # ***
#    # TODO check: if x is too large there might be a numerical error? If so, move the points to 0 before fitting
#  }
#  projectionGPLinear <- tempProjectionGPLinear[2:(sizePrediction + 1),]
#  projectionGRLinear <- tempProjectionGRLinear[2:(sizePrediction + 1),]
#  
#  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
#  if(parametersModel$config$derivativeFromGP == F)
#    projectionGR_onBoundary <- c(tempProjectionGRLinear[1,])
#  else
#    projectionGR_onBoundary <- NA
#  
#  return(list(projectionGP = projectionGPLinear, projectionGR = projectionGRLinear, projectionGR_onBoundary = projectionGR_onBoundary,
#              currentDate = currentDate, sizePrediction = sizePrediction))
#}

#' Create projection matrices using GP, computed in INLA
#' projectionGP: matrix sizePrediction x samples of samples of GP
#' projectionGR: matrix sizePrediction x samples of samples of GR
#' projectionGR_onBoundary: vector samples of GR in last days of observation if using finite differences (since it is NA from model)
#' DEPRICATED 22.09.2023
#getProjectionGP <- function(parametersModel, outputModel){
#  daysForPrediction <- 14
#  #firstPredictionPoint <- outputModel$dateList$maxDay
#  sizePrediction <- parametersModel$config$sizeGPProjection # >7
#  
#  numDays <- outputModel$dateList$numDays
#  matrixSampleDays <- outputModel$matrixSampleDays
#  #BORRARmatrixDerivatives <- outputModel$sampleDerivatives
#  
#  # Prediction - GP
#  projectionGPInla <- outputModel$projectionGP
#  
#  matrixSampleAndPred <- rbind(matrixSampleDays[(numDays - daysForPrediction + 1):numDays,], projectionGPInla)
#  tempDerivative <- getGrowthFromSamples(matrixSampleAndPred)
#  projectionGRInla <- tempDerivative[daysForPrediction + (1:sizePrediction),]
#  
#  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
#  if(parametersModel$config$derivativeFromGP == F)
#    projectionGR_onBoundary <- c(tempDerivative[daysForPrediction,])
#  else
#    projectionGR_onBoundary <- NA
#  
#  return(list(projectionGP = projectionGPInla[1:sizePrediction,], projectionGR = projectionGRInla[1:sizePrediction,],
#              projectionGR_onBoundary = projectionGR_onBoundary,
#              currentDate = outputModel$dateList$dateTable[dayId == outputModel$dateList$maxDay, date], sizePrediction = sizePrediction))
#}

#' 20.09.2023 - code adapted from getGrowthFromSamples_GP() ... Key: find d1Matrix and dMatrixAll
#'  currentDayId: Last observed day to compute the projection. Recommended to be outputModel$dateList$maxDay
getProjectionGP2 <- function(parametersModel, outputModel, currentDate, sizePrediction, daysForPrediction){
  
  if(outputModel$dateList$numDays < daysForPrediction) warning("daysForPrediction is longer than available observations.")
  
  # Get basis configuration
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  sizeSample <- parametersModel$config$sizeSample
  numDaysPast <- min(outputModel$dateList$numDays, daysForPrediction)
  matrixSampleGP <- outputModel$matrixSampleDays[(currentDayId - daysForPrediction + 1):currentDayId,]
  
  # Compute distance matrix for ordered days
  pastInterval <- 1:numDaysPast
  predInterval <- (numDaysPast + 1):(numDaysPast + sizePrediction)
  relativeDistanceMatrixObsObs <- matrix(data = pastInterval, nrow = numDaysPast, ncol = numDaysPast, byrow = F) -
    matrix(data = pastInterval, nrow = numDaysPast, ncol = numDaysPast, byrow = T)
  relativeDistanceMatrixPredObs <- matrix(data = predInterval, nrow = sizePrediction, ncol = numDaysPast, byrow = F) -
    matrix(data = pastInterval, nrow = sizePrediction, ncol = numDaysPast, byrow = T)
  relativeDistanceMatrixPredPred <- matrix(data = predInterval, nrow = sizePrediction, ncol = sizePrediction, byrow = F) -
    matrix(data = predInterval, nrow = sizePrediction, ncol = sizePrediction, byrow = T)
  
  # Get basis parameters
  vGP <- 2 - 1/2
  sigma0 <- parametersModel$params$prior2.sigma0
  range0 <- parametersModel$params$prior2.range0
  
  # Compute auxiliar vectors
  sig2Vector <- (sigma0*exp(outputModel$matrixSampleHyperAll["theta1",]))^2
  kappaVector <- sqrt(8*vGP)/(range0*exp(outputModel$matrixSampleHyperAll["theta2",]))
  
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
    tryCatch({
      #iSample <- MASS::mvrnorm(n = 1, mu = meanMVN, Sigma = sigmaMVN + 0.1*diag(2*sizePrediction))
      iSampleGP <- MASS::mvrnorm(n = 1, mu = meanMVN[1:10], Sigma = sigmaMVN[1:10,1:10])
      iSampleGR <- MASS::mvrnorm(n = 1, mu = meanMVN[11:20], Sigma = sigmaMVN[11:20,11:20])
    }, error = function(e) {
      errorIterations[indexSample] <- 1
      iSample <- rep(NA, 2*sizePrediction)
    })
    
    # Save
    #sampleProjections[,indexSample] <- iSample[1:sizePrediction]
    #sampleDerivatives[,indexSample] <- iSample[(sizePrediction + 1):(sizePrediction)]
    sampleProjections[,indexSample] <- iSampleGP
    sampleDerivatives[,indexSample] <- iSampleGR
  }
  
  if(sum(errorIterations) > 0) warning("Sigma not positive definite for ", sum(errorIterations), " samples out of ", sizeSample)
  
  # Transform derivatives to GR
  samplesGR <- getSamplesGR(matrixSampleDays = sampleProjections,
                            sampleDerivatives = sampleDerivatives,
                            parametersModel = parametersModel)
  #?samplesGR <- sampleDerivatives
  
  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
  if(parametersModel$config$derivativeFromGP == F){
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

#' Function for projections that provide a GP and a GR !!!
#' Creates a table with samples for relevant parameters plus median, 95% CI for model posterior
#' testingVector: vector of length projectionSamplesGP$sizePrediction with testing (for proportions model)
#' TODO this is if we give GP and GR as inputs... maybe we'll need something more general?
createTableProjection <- function(projectionSamplesGP, outputModel, parametersModel, nameProjection, testingVector = NULL){
  # Check input
  if(parametersModel$params$linkType == "BB" & is.null(testingVector)) stop("testingVector is required for the proportions model")
  
  # Auxiliar variables
  sizePrediction <- nrow(projectionSamplesGP$projectionGP)
  levelsWeek <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  currentDayId <- outputModel$dateList$dateTable[date == projectionSamplesGP$currentDate, dayId]
  
  # Create matrices of parameters
  # TODO same as in computePosteriors()... Unify!
  if(parametersModel$params$randomEffect == "weekday"){
    predictionWeekday <- weekdays(projectionSamplesGP$currentDate + (1:sizePrediction))
    matrixRandomEffect <- outputModel$matrixSampleRandomEffect[match(predictionWeekday, levelsWeek),]
  }else if(parametersModel$params$randomEffect == "all"){
    matrixRandomEffect <- matrix(0, nrow = sizePrediction, ncol = parametersModel$config$sizeSample)
  }else{
    matrixRandomEffect <- matrix(0, nrow = sizePrediction, ncol = parametersModel$config$sizeSample)
  }
  matrixSampleOverdisp <- outputModel$matrixSampleHyper[c("overdispersion"),]
  matrixIntercept <- matrix(outputModel$matrixSampleIntercept,
                            nrow = sizePrediction, ncol = parametersModel$config$sizeSample, byrow = T)
  etaSample <- matrixIntercept + projectionSamplesGP$projectionGP + matrixRandomEffect
  rhoSample <- matrix(matrixSampleOverdisp, nrow = sizePrediction, ncol = parametersModel$config$sizeSample, byrow = T)
  
  #?matrixIntercept <- 0
  
  # Create projection table with samples of relevant model parameters
  if(parametersRun$params$linkType == "NB"){
    mupSample <- exp(etaSample)
    gpTransSample <- exp(matrixIntercept + projectionSamplesGP$projectionGP)
    tableSamples <- data.table(lastDayId = currentDayId,
                               dayId = currentDayId + (1:sizePrediction),
                               sample = rep(1:parametersModel$config$sizeSample, each = sizePrediction),
                               #dayId = rep(currentDayId + (1:sizePrediction), each = parametersModel$config$sizeSample), # BUG!!!
                               #sample = 1:parametersModel$config$sizeSample, # BUG!!!
                               gp = c(projectionSamplesGP$projectionGP),
                               gr = c(projectionSamplesGP$projectionGR),
                               gpTrans = c(gpTransSample),
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
    pS <- matrix(mapply(function(a,b) rbeta(1, a, b), aS, bS), nrow = sizePrediction, ncol = parametersModel$config$sizeSample, byrow = F)
    
    tableSamples <- data.table(lastDayId = currentDayId,
                               dayId = currentDayId + (1:sizePrediction),
                               sample = rep(1:parametersModel$config$sizeSample, each = sizePrediction),
                               gp = c(projectionSamplesGP$projectionGP),
                               gr = c(projectionSamplesGP$projectionGR),
                               gpTrans = c(gpTransSample),
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
  if(parametersRun$params$linkType == "NB")
    tableSamples[, numberTest := NA]
  else
    tableSamples[order(dayId), numberTest := testingVector[dayId - lastDayId]]
  
  # Compute median and 95% CI of model posterior and add it to samples table
  if(parametersRun$params$linkType == "NB"){
    tableSamples[, yS := rnbinom(n = nrow(tableSamples), size = rho, mu = mu)]
  }else{
    tableSamples[, yS := rbinom(n = nrow(tableSamples), size = numberTest, prob = pp)] # eg. rbinom(n = 21, size = c(2,30), c(0.1,0.9))
  }
  
  # Create tableProjections with quantiles of every relevant parameter
  tableProjections <- tableSamples[, .(gp_q025 = quantile(gp, 0.025), gp_q25 = quantile(gp, 0.25), gp_median = quantile(gp, 0.5),
                                       gp_q75 = quantile(gp, 0.75), gp_q975 = quantile(gp, 0.975),
                                       #
                                       gr_q025 = quantile(gr, 0.025), gr_q25 = quantile(gr, 0.25), gr_median = quantile(gr, 0.5),
                                       gr_q75 = quantile(gr, 0.75), gr_q975 = quantile(gr, 0.975),
                                       #
                                       gpConsTrans_q025 = quantile(gpTrans, 0.025), gpConsTrans_q25 = quantile(gpTrans, 0.25),
                                       gpConsTrans_median = quantile(gpTrans, 0.5),
                                       gpConsTrans_q75 = quantile(gpTrans, 0.75), gpConsTrans_q975 = quantile(gpTrans, 0.975),
                                       #
                                       eta_q025 = quantile(eta, 0.025), eta_q25 = quantile(eta, 0.25), eta_median = quantile(eta, 0.5),
                                       eta_q75 = quantile(eta, 0.75), eta_q975 = quantile(eta, 0.975),
                                       #
                                       mu_q025 = quantile(mu, 0.025), mu_q25 = quantile(mu, 0.25),
                                       mu_median = quantile(mu, 0.5), mu_q75 = quantile(mu, 0.75), mu_q975 = quantile(mu, 0.975),
                                       #
                                       yS_q025 = quantile(yS, 0.025, na.rm = T), yS_q25 = quantile(yS, 0.25, na.rm = T),
                                       yS_median = quantile(yS, 0.5, na.rm = T),
                                       yS_q75 = quantile(yS, 0.75, na.rm = T), yS_q975 = quantile(yS, 0.975, na.rm = T)),
                                   .(dayId, date, lastDayId, dateLastDay, numberTest)]
  
  # Add GR on last day of observation if not provided by model (when GP from finite differences)
  if(parametersModel$config$derivativeFromGP == F){
    tableSamples_boundary <- data.table(lastDayId = currentDayId,
                                        dayId = currentDayId,
                                        sample = 1:parametersModel$config$sizeSample,
                                        gr = projectionSamplesGP$projectionGR_onBoundary)
    tableProjections_boundary <- tableSamples_boundary[, .(gr_q025 = quantile(gr, 0.025), gr_q25 = quantile(gr, 0.25), gr_median = quantile(gr, 0.5),
                                                           gr_q75 = quantile(gr, 0.75), gr_q975 = quantile(gr, 0.975)),
                                                       .(dayId, lastDayId)]
    # Add extras to the projection table
    setkey(tableProjections_boundary, lastDayId)
    setkey(outputModel$dateList$dateTable, dayId)
    tableProjections_boundary[outputModel$dateList$dateTable, dateLastDay := date]
    tableProjections_boundary[, date := dateLastDay + dayId - lastDayId]
    
    setkey(tableProjections_boundary, dayId)
    setkey(outputModel$posteriorTransfGP, dayId)
    tableProjections_boundary[outputModel$posteriorTransfGP, ":="(numberTest = i.numberTest,
                                                     gp_q025 = i.q0.025_GP, gp_q25 = i.q0.25_GP, gp_median = i.median_GP, gp_q75 = i.q0.75_GP, gp_q975 = i.q0.975_GP,
                                                     #gpTrans_q025 = i.q0.025FT, gpTrans_q25 = i.q0.25FT, gpTrans_median = i.medianFT,
                                                     #gpTrans_q75 = i.q0.75FT, gpTrans_q975 = i.q0.975FT) # removed 19.09.2023
                                                     gpConsTrans_q025 = i.q0.025_transConsGP, gpConsTrans_q25 = i.q0.25_transConsGP,
                                                     gpConsTrans_median = i.median_transConsGP,
                                                     gpConsTrans_q75 = i.q0.75_transConsGP, gpConsTrans_q975 = i.q0.975_transConsGP)]
    # TODO Double-check the above correct!!!!
  }else{
    tableSamples_boundary <- NA
    tableProjections_boundary <- NA
  }
    
  return(list(tableProjections = tableProjections, tableSamples = tableSamples,
              tableProjections_boundary = tableProjections_boundary, tableSamples_boundary = tableSamples_boundary,
              currentDate = projectionSamplesGP$currentDate, sizePrediction = sizePrediction, nameProjection = nameProjection))
}


#' Creates projection under null model (future observations are in 95% CI of normal around latest observations/proportion avg and sd = sdPast)
#' daysForPrediction: days before currentDate defined as 'latest observations'/'tests'
#' TODO caveat: quantiles from samples could be calculated with the qnorm function instead of using tableSample, but it's more practical to provide samples
createTableProjectionNull <- function(outputModel = outputModel, parametersModel = parametersRun, currentDate = currentDate,
                                      sizePrediction, daysForPrediction = 7, testingVector = NULL){
  # Check input
  if(parametersModel$params$linkType == "BB" & is.null(testingVector)) stop("testingVector is required for the proportions model")
  
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  casesBeforePrediction <- outputModel$dataForModel[date > currentDate - daysForPrediction & date <= currentDate][order(date), positiveResults]
  testsBeforePrediction <- outputModel$dataForModel[date > currentDate - daysForPrediction & date <= currentDate][order(date), numberTest]
  
  # Table projection output
  overviewNull <- data.table(dayId = currentDayId + 1:sizePrediction,
                             date = currentDate + 1:sizePrediction,
                             lastDayId = currentDayId,
                             dateLastDay = currentDate,
                             numberTest = testingVector)
  overviewNull[, ":="(gp_q025 = as.numeric(NA), gp_q25 = as.numeric(NA), gp_median = as.numeric(NA), gp_q75 = as.numeric(NA), gp_q975 = as.numeric(NA),
                      #
                      gr_q025 = 0, gr_q25 = 0, gr_median = 0, gr_q75 = 0, gr_q975 = 0,
                      #
                      gpConsTrans_q025 = as.numeric(NA), gpConsTrans_q25 = as.numeric(NA), gpConsTrans_median = as.numeric(NA),
                      gpConsTrans_q75 = as.numeric(NA), gpConsTrans_q975 = as.numeric(NA)
                      #gpTrans_q025 = pmax(0, qnorm(0.025, mean = positiveResultsLastDay, sd = sdPast)),
                      #gpTrans_q25 = pmax(0, qnorm(0.25, mean = positiveResultsLastDay, sd = sdPast)),
                      #gpTrans_median = positiveResultsLastDay,
                      #gpTrans_q75 = pmax(0, qnorm(0.75, mean = positiveResultsLastDay, sd = sdPast)),
                      #gpTrans_q975 = pmax(0, qnorm(0.975, mean = positiveResultsLastDay, sd = sdPast)),
                      )]
  
  # Compute yS
  if(parametersModel$params$linkType == "NB"){
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
  
  # Table samples output - artificial samples that can be useful when using scoring functions (see caveat in description)
  # Number of samples = parametersModel$config$sizeSample
  tableSamples <- CJ(dayId = currentDayId + 1:sizePrediction, sample = 1:parametersModel$config$sizeSample)
  setkey(tableSamples, dayId)
  setkey(overviewNull, dayId)
  tableSamples[overviewNull, ":="(lastDayId = lastDayId, date = date, dateLastDay = dateLastDay, numberTest = numberTest)]
  
  # Produce samples
  if(parametersModel$params$linkType == "NB")
    tableSamples[, yS := rnorm(n = parametersModel$config$sizeSample, mean = proportionsLastDay, sd = sdPast)]
  else
    tableSamples[, yS := numberTest*pmax(0, pmin(1, rnorm(n = parametersModel$config$sizeSample, mean = proportionsLastDay, sd = sdPast)))]
  
  return(list(tableProjections = overviewNull, tableSamples = tableSamples,
              currentDate = currentDate, sizePrediction = sizePrediction, nameProjection = "null"))
}

#' Plots posterior of GP and model fitting before the last observation date (30 days before as default) and after the size of the projection
#' listProjection: output of createTableProjection...()
#' outputModel:
#' parametersModel:
#' minDate (optional): min date for plot
#' futureTable (optiona): data from the future for comparison (date, positiveResults, numberTest)
#' TODO requires additional package: ggtree
plotProjection <- function(listProjection, outputModel, parametersModel, futureTable = NULL, minDate = NULL, plotGP = F){
  # TODO make better plot
  if(is.null(minDate)) minDate <- listProjection$currentDate - 30
  maxDate <- listProjection$currentDate + sizePrediction
  
  # Create auxiliar tables
  countTableAll <- outputModel$dataForModel[date >= minDate & date <= maxDate, .(date, positiveResults, numberTest)]
  if(!is.null(futureTable)) countTableAll <- rbind(countTableAll, futureTable)
  countTableAll[, ratio := pmin(1, pmax(0, positiveResults/numberTest))]
  
  currentDate <- listProjection$currentDate
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  sizePrediction <- listProjection$sizePrediction
  
  dataToPlotBeforeGP <- outputModel$posteriorTransfGP[date >= minDate & date <= maxDate] #[dayId %in% (currentDayId - (0:30))]
  dataToPlotBeforeGR <- outputModel$posteriorGrowth[date >= minDate & date <= maxDate] #[dayId %in% (currentDayId - (0:30))]
  dataToPlotAfter <- listProjection$tableProjections
  
  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
  if(parametersModel$config$derivativeFromGP == F & !is.null(listProjection$tableProjections_boundary)){
    setkey(dataToPlotBeforeGR, dayId)
    setkey(listProjection$tableProjections_boundary, dayId)
    dataToPlotBeforeGR[listProjection$tableProjections_boundary,
                                ":="(median = i.gr_median, q0.025 = i.gr_q025, q0.975 = i.gr_q975, q0.25 = i.gr_q25, q0.75 = i.gr_q75)]
  }
  
  # Compute model posterior
  if(parametersModel$params$linkType == "NB"){
    dataToPlotAfter[, ":="(yTrans_median = yS_median, yTrans_q025 = yS_q025, yTrans_q975 = yS_q975, yTrans_q25 = yS_q25, yTrans_q75 = yS_q75)]
  }else{
    dataToPlotAfter[, ":="(yTrans_median = yS_median/numberTest, yTrans_q025 = yS_q025/numberTest, yTrans_q975 = yS_q975/numberTest,
                           yTrans_q25 = yS_q25/numberTest, yTrans_q75 = yS_q75/numberTest)]
  }
  
  # Plot
  #jointDatesTrick <- rbind(outputModel$dateList$dateTable, dataToPlotAfter[, .(dayId, date)])
  ggP1 <- ggplot(dataToPlotBeforeGP, aes(x = date)) + theme_laura() +
    #geom_vline(xintercept = lockdownDates, linetype = 2, colour = "gray50") +
    #geom_line(data = allPosterior[order(sample, date)], aes(y = value, group = sample))
    geom_ribbon(aes(ymin = q0.025FT, ymax = q0.975FT), fill = "gray70", alpha = 0.5) +
    geom_ribbon(aes(ymin = q0.25FT, ymax = q0.75FT), fill = "gray40", alpha = 0.5) +
    geom_line(aes(y = medianFT)) +
    geom_ribbon(data = dataToPlotAfter, aes(x = date, ymin = yTrans_q025, ymax = yTrans_q975), fill = "#A8D0DB", alpha = 0.5) +
    geom_ribbon(data = dataToPlotAfter, aes(ymin = yTrans_q25, ymax = yTrans_q75), fill = "#0C7489", alpha = 0.5) +
    geom_line(data = dataToPlotAfter, aes(y = yTrans_median), colour = "#13505B") +
    ##scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
    #if(linkType == "BB") gg9.1 <- gg9.1p +
    geom_point(data = countTableAll,
               aes(y = ratio, colour = weekdays(date) %in% c("Saturday", "Sunday"), shape = weekdays(date) %in% c("Saturday", "Sunday")), show.legend = F, size = 1) +
    scale_colour_manual(values = c("black", "#D41B19")) +
    #labs(title = paste0(partitionToPlot, " - Gaussian process (positive cases/number of tests)"), x = "day", y = "relative effect of GP (positive cases/number of tests)")
    labs(x = "day", y = "posterior fitting", title = listProjection$nameProjection)
  ggP2 <- ggplot(dataToPlotBeforeGR, aes(x = date)) + theme_laura() +
    #geom_vline(xintercept = lockdownDates, linetype = 2, colour = "gray50") +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "gray70", alpha = 0.5) +
    geom_ribbon(aes(ymin = q0.25, ymax = q0.75), fill = "gray40", alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_ribbon(data = dataToPlotAfter, aes(ymin = gr_q025, ymax = gr_q975), fill = "#EFDABD", alpha = 0.5) +
    geom_ribbon(data = dataToPlotAfter, aes(ymin = gr_q25, ymax = gr_q75), fill = "#D7A35B", alpha = 0.5) +
    geom_line(data = dataToPlotAfter, aes(y = gr_median), colour = "#31220C") +
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "growth rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if(plotGP == T){
    ggP1 <- ggP1 +
      geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "#EFDABD", alpha = 0.5) +
      geom_ribbon(aes(ymin = q0.25, ymax = q0.75), fill = "#D7A35B", alpha = 0.5) +
      geom_line(aes(y = median), colour = "#31220C") +
      geom_ribbon(data = dataToPlotAfter, aes(ymin = gpConsTrans_q025, ymax = gpConsTrans_q975), fill = "#EFDABD", alpha = 0.5) +
      geom_ribbon(data = dataToPlotAfter, aes(ymin = gpConsTrans_q25, ymax = gpConsTrans_q75), fill = "#D7A35B", alpha = 0.5) +
      geom_line(data = dataToPlotAfter, aes(y = gpConsTrans_median), colour = "#31220C")
  }
  
  #return(suppressWarnings(multiplot(ggP1, ggP2, cols = 1)))
  return(list(ggP1, ggP2))
}

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

#' Function for comparing scores
#' projectionsList: list of listProjection objects
#' scoringList: list of tableScore objects
#' positivesVector: positives observed in the future
#' TODO add list of optional arguments?
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

