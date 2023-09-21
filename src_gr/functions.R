
# Functions for running model ----

# Auxiliar function
#' Loads parameters
#' linkType:
#' prior2.sigma0: base-line value for 1 / precision of covariance function
#' prior2.lengthscale0: base-line value for the length scale of covariance function
#' theta.prior2.mean: mean of theta prior
#' theta.prior2.prec: precision matrix of theta prior
#' BB.prior.rho: prior for log precision if beta-binomial case
#' NB.prior.rho: prior for log precision if negative binomial case
#' dayWeek: if T, the day-of-the-week effect is used. If F, the noise term is used
#' dayWeek.prior.prec: prior for day-of-the-week effect or noise term
#' sizeSample: samples to estimate posterior of parameters
#' derivativeFromGP: if growth rate distribution is estimated from samples of GP derivative.
#'   If F, dataset must have at least 7 time points between minDate and maxDate. If not, it will automatically use the numerical version.
setParametersFn <- function(linkType,
                            prior2.sigma0 = 1,
                            prior2.lengthscale0 = 100,
                            theta.prior2.mean = c(0,0),
                            theta.prior2.prec = diag(2),
                            BB.prior.rho = list(overdispersion = list(prior = "gaussian", param = c(0, 0.5))),
                            NB.prior.rho = list(size = list(prior = "gaussian", param = c(0, 0.01))),
                            unitTime = "day", # day week integer
                            randomEffect = "weekday", # weekday all none # if day, adds day effect automatically # old dayWeek
                            dayWeek.prior.prec = list(theta = list(prior = 'loggamma', param = c(1, 0.01))), # TODO fix name
                            sizeSample = sizeSample,
                            derivativeFromGP = F,
                            computeGPProjection = F,
                            sizeGPProjection = 10){
  # Checks
  if(!randomEffect %in% c("weekday", "all", "none")) stop("Parameters randomEffect can only take a values from 'weekday', 'all', 'none'.")
  
  parameters <- list(
    params = list(
      linkType = linkType,
      prior2.sigma0 = prior2.sigma0,
      prior2.range0 = 2*prior2.lengthscale0,
      theta.prior2.mean = theta.prior2.mean,
      theta.prior2.prec = theta.prior2.prec,
      BB.prior.rho = BB.prior.rho,
      NB.prior.rho = NB.prior.rho,
      unitTime = unitTime,
      randomEffect = randomEffect,
      dayWeek.prior.prec = dayWeek.prior.prec
    ),
    config = list(
      sizeSample = sizeSample,
      derivativeFromGP = derivativeFromGP,
      computeGPProjection = computeGPProjection,
      sizeGPProjection = sizeGPProjection
      #computeTaylorProjection = F,
      #computeLinearProjection = T
      )
    )
    return(parameters)
}

# TODO create function to check countTable before using it. It should have numberTest if BB, or create an empty one if non existant and NB.
# It should have positives <= number tests

#' countTable: data table with one row per day (of days with available data) and these columns:
#' - numberTest: number of test on the day (>= 0)
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - date: date of count in R date format
#' parametersModel: output of setParametersFn()
#' saveSamples: Default F. If T, returns matrixSampleDays and sampleDerivatives,
#'              two matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' minDate (optional): minimum date to include in the model
#' maxDate (optional): maximum date to include in the model
runModelGrowthRate <- function(countTable, parametersModel, minDate = NULL, maxDate = NULL){
  # TODO check and complete content of parametersModel from function (i.e. if NULL then default)
  
  # Load parameters into function environment and add auxiliar variables
  # TODO check if min dates in countTable are aligned as in unitTime
  # TODO include cases with no date
  # TODO remove this:
  #list2env(parametersModel$params, envir = environment())
  #list2env(parametersModel$config, envir = environment())
  levelsWeek <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
  
  # Create auxiliar table with dates
  if(is.null(minDate)) minDate <- min(countTable$date)
  if(is.null(maxDate)) maxDate <- max(countTable$date)
  minDay <- 1
  maxDay <- length(seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime)) #as.integer(maxDate - minDate + 1)
  numDays <- maxDay #maxDay - minDay + 1
  dateTable <- data.table(dayId = minDay:maxDay,
                          date = seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime))
  if(nrow(dateTable) <= 1) stop("There must be at least 2 time units between minDate and maxDate")
  
  # ---------------------------------------------------- #
  #                      FIT MODEL                       #
  # ---------------------------------------------------- #
  # Create data with all days, including the ones with missing data
  dataForModel <- data.table(dayId = minDay:maxDay,
                             date = seq.Date(from = minDate, to = maxDate, by = parametersModel$params$unitTime),
                             numberTest = as.integer(NA),
                             positiveResults = as.integer(NA))
  setkey(dataForModel, date)
  setkey(countTable, date)
  dataForModel[countTable, ":="(numberTest = i.numberTest, positiveResults = i.positiveResults)]
  dataForModel[numberTest == 0, ":="(numberTest = NA, positiveResults = NA)]
  dataForModel[, dayWeek := weekdays(date)]
  
  # Define grid for finite differences
  boundaryVertices <- 200 # Points on the left and right
  boundaryStepSize <- 10 # Distance between boundary points
  boundaryPoints <- c(seq(min(dataForModel$dayId) - boundaryStepSize*boundaryVertices, min(dataForModel$dayId) - boundaryStepSize, by = boundaryStepSize),
                      seq(min(dataForModel$dayId), max(dataForModel$dayId), by = 1), 
                      seq(max(dataForModel$dayId) + boundaryStepSize, max(dataForModel$dayId) + boundaryStepSize*boundaryVertices, by = boundaryStepSize))
  nonBoundaryIndices <- (boundaryVertices + 1):(length(boundaryPoints) - boundaryVertices)
  
  # Create mesh in 1d (mesh1d), projection matrix (A1), model (spde1), and indexes (spde1.idx)
  
  # - Mesh and projection matrix
  mesh1d <- inla.mesh.1d(boundaryPoints, boundary = "free") # NEW 09.02.2021 boundary issues
  A1 <- inla.spde.make.A(mesh1d, dataForModel[order(dayId), dayId])
  
  # - Set priors
  # Prior for day of the week
  priorDayWeek <- parametersModel$params$dayWeek.prior.prec
  # Prior for Gaussian process parameters, following (Lindgren, 2015, v63i19)
  vGP <- 2 - 1/2
  sigma0 <- parametersModel$params$prior2.sigma0
  range0 <- parametersModel$params$prior2.range0
  # Convert into tau and kappa:
  kappa0 <- sqrt(8*vGP)/range0
  tau0 <- sqrt(gamma(vGP)/(gamma(2)*sqrt(4*pi)))/(kappa0^(vGP)*sigma0)
  basis.prior2.tau <- c(-1, 3/2)
  basis.prior2.kappa <- c(0, -1)
  spde1 <- inla.spde2.matern(mesh1d,
                             B.tau = cbind(log(tau0), basis.prior2.tau[1], basis.prior2.tau[2]),
                             B.kappa = cbind(log(kappa0), basis.prior2.kappa[1], basis.prior2.kappa[2]),
                             theta.prior.mean = parametersModel$params$theta.prior2.mean,
                             theta.prior.prec = parametersModel$params$theta.prior2.prec)
  
  # - Create stack data
  spde1.idx <- inla.spde.make.index("day", n.spde = spde1$n.spde)
  # Create stacked data with two datasets (observations and prediction). For paper, I'm ignoring prediction.
  stackD <- inla.stack(data = list(positiveResults = dataForModel[order(dayId), positiveResults]),
                       A = list(A1, 1, 1),
                       effects = list(c(list(Intercept = 1), spde1.idx),
                                      list(dayWeek = factor(dataForModel[order(dayId), dayWeek], levels = levelsWeek)),
                                      list(dayId2 = dataForModel[order(dayId), dayId])),
                       tag = "est")
  # Create prediction data for the next parametersModel$config$sizeGPProjection days
  daysProjection <- parametersModel$config$sizeGPProjection + 1
  predictionGPWeekday <- dateTable[dayId == max(dataForModel$dayId), weekdays(date + 1:daysProjection)]
  xx <- seq(max(dataForModel$dayId) + 1, max(dataForModel$dayId) + daysProjection, by = 1)
  A.xx <- inla.spde.make.A(mesh1d, xx)
  stack.pred <- inla.stack(data = list(positiveResults = NA),
                           A = list(A.xx),
                           effects = list(c(list(Intercept = 1), spde1.idx)),
                           tag = "pred")
  joint.stack <- inla.stack(stackD, stack.pred)
  
  # Fit model (using INLA)
  cat("Fitting model ... ")
  if(parametersModel$params$randomEffect == "weekday"){
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayWeek, model = 'iid', hyper = priorDayWeek, constr = T)
  }else if(parametersModel$params$randomEffect == "all"){
    formula <- positiveResults ~ -1 + f(day, model = spde1) + f(dayId2, model = 'iid', hyper = priorDayWeek, constr = T)
  }else{
    formula <- positiveResults ~ -1 + f(day, model = spde1)
  }
  if(parametersModel$params$linkType == "BB"){
    # Betabinomial
    objectInla <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "betabinomial",
                   Ntrials =  c(dataForModel[order(dayId), numberTest], rep(NA, length(xx))),
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = parametersModel$params$BB.prior.rho))
  }else if(parametersModel$params$linkType == "NB"){
    # Negative binomial
    objectInla <- inla(formula = formula, data = inla.stack.data(joint.stack),
                   control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                   family = "nbinomial",
                   control.compute = list(config = TRUE),
                   control.family = list(hyper = parametersModel$params$NB.prior.rho))
  }
  
  objectInla$dateList <- list(dateTable = dateTable, numDays = numDays, minDay = minDay, maxDay = maxDay)
  objectInla$dataForModel <- dataForModel
  objectInla$nonBoundaryIndices <- nonBoundaryIndices # TODO better way?
  return(objectInla)
}

processINLAOutput <- function(objectInla, parametersModel, saveSamples = F){
  # TODO make it faster
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  cat("Estimating growth rate ... ")
  
  # Get samples of posterior distribution of parameters
  sampleList <- inla.posterior.sample(parametersModel$config$sizeSample, objectInla)
  
  # Extract relevant indexes from output
  rownms <- rownames(sampleList[[1]]$latent)
  set1 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "APredictor"))
  #set2 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "Predictor"))
  set3 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "day"))
  set4 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == case_when(parametersModel$params$randomEffect == "weekday" ~ "dayWeek",
                                                                                    parametersModel$params$randomEffect == "all" ~ "dayId2",
                                                                                    parametersModel$params$randomEffect == "none" ~ "NA")))
  dayIndexInSample <- set3[objectInla$nonBoundaryIndices]
  
  # Create matrix of samples (and hyperparameters/weekday effect if needed)
  matrixSampleDays <- sapply(sampleList, function(x) x$latent[dayIndexInSample,1])
  matrixSampleWeekday <- sapply(sampleList, function(x) x$latent[set4,1])
  matrixSampleEta <- sapply(sampleList, function(x) x$latent[set1[1:objectInla$dateList$numDays],1])
  
  nameDispersionINLA <- ifelse(parametersModel$params$linkType == "NB",
                               "size for the nbinomial observations (1/overdispersion)",
                               "overdispersion for the betabinomial observations")
  namePrecisionINLA <- case_when(parametersModel$params$randomEffect == "weekday" ~ "Precision for dayWeek",
                                 parametersModel$params$randomEffect == "all" ~ "Precision for dayId2",
                                 parametersModel$params$randomEffect == "none" ~ "NA")
  #matrixSampleHyper <- sapply(sampleList, function(x) exp(x$hyperpar[c("Theta1 for day", "Theta2 for day")]))
  #matrixSampleHyperAll <- sapply(sampleList, function(x) x$hyperpar) # overdisp log(theta1) log(theta2) precision
  #matrixSampleHyperAll <- rbind(sapply(sampleList, function(x) log(x$hyperpar[c(nameDispersionINLA, namePrecisionINLA)])), ...
  matrixSampleHyperAll <- sapply(sampleList, function(x) x$hyperpar[c("Theta1 for day", "Theta2 for day",
                                                                      nameDispersionINLA, namePrecisionINLA)]) # NEW 08.09.2023
  rownames(matrixSampleHyperAll) <- c("theta1", "theta2", "overdispersion", "precision")
  
  # Create matrix with GP projection (sizeGPProjection + 1 as we lose one in the GR estimation)
  daysProjection <- parametersModel$config$sizeGPProjection + 1
  predicIndexInSample <- set1[objectInla$dateList$numDays + (1:daysProjection)]
  projectionGP <- sapply(sampleList, function(x) x$latent[predicIndexInSample,1]) # previous projectionGPInla
  
  # Compute derivative
  if(parametersModel$config$derivativeFromGP == T | nrow(objectInla$dateList$dateTable) < 7){
    # Sample from derivative
    sampleDerivatives <- getGrowthFromSamples_GP(matrixSampleGP = matrixSampleDays,
                                                 samplesHyperparam = matrixSampleHyperAll,
                                                 sigma0 = parametersModel$params$prior2.sigma0,
                                                 range0 = parametersModel$params$prior2.range0)
  }else{
    # Compute approximate derivative using windows (+-3)
    # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
    sampleDerivatives <- getGrowthFromSamples(matrixSampleDays = matrixSampleDays)
  }
  
  # ---------------------------------------------------- #
  #                  PRODUCE OUTPUT                      #
  # ---------------------------------------------------- #
  listPosteriors <- computePosteriors(matrixSampleDays, sampleDerivatives, matrixSampleHyperAll, matrixSampleEta, matrixSampleWeekday, objectInla, parametersModel)
  
  setkey(listPosteriors$posteriorGrowth, dayId)
  setkey(objectInla$dateList$dateTable, dayId)
  listPosteriors$posteriorGrowth[objectInla$dateList$dateTable, ":="(date = i.date)]
  
  setkey(listPosteriors$posteriorTransfGP, dayId)
  setkey(objectInla$dateList$dateTable, dayId)
  listPosteriors$posteriorTransfGP[objectInla$dateList$dateTable, ":="(date = i.date)]
  setkey(listPosteriors$posteriorTransfGP, date)
  setkey(objectInla$dataForModel, date)
  listPosteriors$posteriorTransfGP[objectInla$dataForModel, ":="(positiveResults = i.positiveResults, numberTest = i.numberTest)]
  
  # Output
  output_main <- list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP,
                      posteriorRandomEffect = listPosteriors$posteriorRandomEffect,
                      dateList = objectInla$dateList, dataForModel = objectInla$dataForModel)
  output_samples <- list(matrixSampleDays = matrixSampleDays, sampleDerivatives = sampleDerivatives, matrixSampleHyperAll = matrixSampleHyperAll)
  output_projection <- list(matrixSampleWeekday = matrixSampleWeekday, projectionGP = projectionGP)
  
  output <- output_main
  if(saveSamples == T) output <- c(output, output_samples)
  if(parametersModel$config$computeGPProjection == T) output <- c(output, output_projection)
  return(output)
}

# INTERNAL
# NOT AVAILABLE ANYMORE as I don't trust INLA sampling of hyperparameters
getGrowthFromSamples <- function(matrixSampleDays){
  # Compute approximate derivative using windows (+-3)
  # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
  numDays <- nrow(matrixSampleDays)
  sampleDerivatives <- rbind(matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays)),
                             ((matrixSampleDays[3,] - matrixSampleDays[1,])/2),
                             ((matrixSampleDays[5,] - matrixSampleDays[1,])/4),
                             (matrixSampleDays[7:numDays,] - matrixSampleDays[1:(numDays - 7 + 1),])/6,
                             ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 5 + 1,])/4),
                             ((matrixSampleDays[numDays,] - matrixSampleDays[numDays - 3 + 1,])/2),
                             matrix(NA, nrow = 1, ncol = ncol(matrixSampleDays)))
  return(sampleDerivatives)
}

# INTERNAL
# TODO objectInla should not be arg? should be undo into dataForModel as arg only and objectInla optional?
computePosteriors <- function(matrixSampleDays, sampleDerivatives, matrixSampleHyper, matrixSampleEta, matrixSampleWeekday, objectInla, parametersModel,
                              ifINLAMarginal = FALSE){
  numDays <- nrow(matrixSampleDays)
  
  # ---------------------------------------------------- #
  #                POSTERIOR GROWTH RATE                 #
  # ---------------------------------------------------- #
  # Compute log (or inv. logit) of posterior of Gaussian process derivative in transformed space
  cat("Computing posterior of growth rate ... ")
  if(parametersModel$params$linkType == "BB"){
    tempExpGP <- exp(matrixSampleDays)
    tempLogit <- tempExpGP/(1 + tempExpGP)
    tempList <- sampleDerivatives/(1 + tempLogit)
  }else if(parametersModel$params$linkType == "NB"){
    tempList <- sampleDerivatives
  }
  tempDoubling <- abs(log(2)/tempList)
  posteriorGrowth <- data.table(dayId = 1:numDays,
                                mean = rowMeans(tempList),
                                sd = apply(tempList, 1, sd),
                                median = apply(tempList, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025 = apply(tempList, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975 = apply(tempList, 1, quantile, probs = 0.975, na.rm = T),
                                q0.25 = apply(tempList, 1, quantile, probs = 0.250, na.rm = T),
                                q0.75 = apply(tempList, 1, quantile, probs = 0.750, na.rm = T),
                                prob0 = apply(tempList, 1, function(x) sum(x >= 0)/parametersModel$config$sizeSample),
                                medianDoubT = apply(tempDoubling, 1, quantile, probs = 0.5, na.rm = T),
                                q0.025DoubT = apply(tempDoubling, 1, quantile, probs = 0.025, na.rm = T),
                                q0.975DoubT = apply(tempDoubling, 1, quantile, probs = 0.975, na.rm = T))
  
  # ---------------------------------------------------- #
  #                 POSTERIOR INCIDENCE                  #
  # ---------------------------------------------------- #
  # Compute posterior of Gaussian process in real space (incidence or positivity)
  cat("Computing posterior of incidence... ")
  if(ifINLAMarginal){
    # (this version is slow and applies only to INLA)
    samplesGP <- objectInla$marginals.random$day[objectInla$nonBoundaryIndices]
    if(parametersModel$params$linkType %in% c("NB")){
      transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(exp, x))
    }else if(parametersModel$params$linkType == "BB"){
      transformedSamples <- lapply(samplesGP, function(x) inla.tmarginal(function(x) exp(x)/(1 + exp(x)), x))
    }
    posteriorTransfGP <- data.table(dayId = 1:numDays,
                                    t(sapply(transformedSamples, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))),
                                    t(sapply(samplesGP, function(x) inla.qmarginal(c(0.5, 0.025, 0.975, 0.25, 0.75), x))))
    setnames(posteriorTransfGP,
             c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
             c("median_transGP", "q0.025_transGP", "q0.975_transGP", "q0.25_transGP", "q0.75_transGP",
               "median_GP", "q0.025_GP", "q0.975_GP", "q0.25_GP", "q0.75_GP"))
    # TODO median_GP
  }else{
    # NEW 09.01.2023
    if(parametersModel$params$linkType %in% c("NB")){
      transformedSamples <- exp(matrixSampleDays)
    }else if(parametersModel$params$linkType == "BB"){
      transformedSamples <- exp(matrixSampleDays)/(1 + exp(matrixSampleDays))
    }
    posteriorTransfGP <- data.table(dayId = 1:numDays,
                                  median_transGP = apply(transformedSamples, 1, quantile, probs = 0.5, na.rm = T),
                                  q0.025_transGP = apply(transformedSamples, 1, quantile, probs = 0.025, na.rm = T),
                                  q0.975_transGP = apply(transformedSamples, 1, quantile, probs = 0.975, na.rm = T),
                                  q0.25_transGP = apply(transformedSamples, 1, quantile, probs = 0.250, na.rm = T),
                                  q0.75_transGP = apply(transformedSamples, 1, quantile, probs = 0.750, na.rm = T),
                                  median_GP = apply(matrixSampleDays, 1, quantile, probs = 0.5, na.rm = T),
                                  q0.025_GP = apply(matrixSampleDays, 1, quantile, probs = 0.025, na.rm = T),
                                  q0.975_GP = apply(matrixSampleDays, 1, quantile, probs = 0.975, na.rm = T),
                                  q0.25_GP = apply(matrixSampleDays, 1, quantile, probs = 0.250, na.rm = T),
                                  q0.75_GP = apply(matrixSampleDays, 1, quantile, probs = 0.750, na.rm = T))
    
  }
  
  # ---------------------------------------------------- #
  #                 POSTERIOR MODEL FIT                  #
  # ---------------------------------------------------- #
  # Compute the model posterior (as in R31.R)
  sizeSample <- parametersModel$config$sizeSample
  sampleOverdisp <- matrixSampleHyper[c("overdispersion"),]
  #predictionWeekday <- ?
  #matrixSampleEta <- matrixSampleDays + matrixSampleWeekday[match(predictionWeekday, levelsWeek),]
  if(parametersModel$params$linkType == "NB"){
    samplesFit <- matrix(rnbinom(n = matrix(1, nrow = numDays, ncol = sizeSample),
                                size = matrix(sampleOverdisp, nrow = numDays, ncol = sizeSample, byrow = T),
                                mu = exp(matrixSampleEta)),
                        nrow = numDays, ncol = sizeSample, byrow = F)
  }else{
    testData <- objectInla$dataForModel[order(dayId), numberTest]
    
    matrixSampleOverdisp <- matrix(sampleOverdisp, nrow = numDays, ncol = sizeSample, byrow = T)
    samplesMu <- exp(matrixSampleEta)/(1 + exp(matrixSampleEta))
    matrixAlpha <- samplesMu*(1 - matrixSampleOverdisp)/matrixSampleOverdisp
    matrixBeta <- (1 - samplesMu)*(1 - matrixSampleOverdisp)/matrixSampleOverdisp
    
    matrixP <- matrix(rbeta(matrix(1, nrow = numDays, ncol = sizeSample), shape1 = matrixAlpha, shape2 = matrixBeta),
                      nrow = numDays, ncol = sizeSample, byrow = F)
    matrixTests <- matrix(testData, nrow = numDays, ncol = sizeSample, byrow = F)
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
  if(parametersModel$params$randomEffect != "none"){
    posteriorRandomEffect <- data.table(index = 1:nrow(matrixSampleWeekday),
                                        median = apply(matrixSampleWeekday, 1, quantile, probs = 0.5, na.rm = T),
                                        q0.025 = apply(matrixSampleWeekday, 1, quantile, probs = 0.025, na.rm = T),
                                        q0.975 = apply(matrixSampleWeekday, 1, quantile, probs = 0.975, na.rm = T),
                                        q0.25 = apply(matrixSampleWeekday, 1, quantile, probs = 0.250, na.rm = T),
                                        q0.75 = apply(matrixSampleWeekday, 1, quantile, probs = 0.750, na.rm = T))
  }else{
    posteriorRandomEffect <- NULL
  }
  
  return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP, posteriorRandomEffect = posteriorRandomEffect))
}

#' matrixSampleDays: days x samples
#' samplesHyperparam: OLD (exp(theta 1), exp(theta2)) X samples. sigma = sigma0 exp(theta1), range = range0 exp(theta2)
#' samplesHyperparam: hyperparams X samples. sigma = sigma0 exp(theta1), range = range0 exp(theta2)
#' Copied FROM R47/Sandbox_deerivative.R
getGrowthFromSamples_GP <- function(matrixSampleGP, samplesHyperparam, sigma0, range0){
  numDays <- nrow(matrixSampleGP)
  sizeSample <- ncol(matrixSampleGP)
  
  # Compute distance matrix for ordered days
  distanceMatrix <- sapply(1:numDays, function(nd) abs(nd - (1:numDays)))
  auxRelativeDistanceMatrix <- matrix(data = 1:numDays, nrow = numDays, ncol = numDays, byrow = F) -
    matrix(data = 1:numDays, nrow = numDays, ncol = numDays, byrow = T)
  
  # Compute auxiliar vectors, with vGP = 3/2
  sig2Vector <- (sigma0*exp(samplesHyperparam["theta1",]))^2
  kappaVector <- sqrt(12)/(range0*exp(samplesHyperparam["theta2",]))
  
  # Loop per sample of (f1, ..., fn, log.tau, log.kappa)
  sampleDerivatives <- matrix(0, nrow = sizeSample, ncol = numDays)
  for(indexSample in 1:sizeSample){
    sig2Value <- sig2Vector[indexSample]
    kappaVal <- kappaVector[indexSample]
    
    expMatrix <- exp(-kappaVal*distanceMatrix)
    deltaMatrix <- sig2Value*(1 + kappaVal*distanceMatrix)*expMatrix
    invDeltaMatrix <- chol2inv(chol(deltaMatrix)) # solve vs. chol2inv system.time(31700*system.time(solve(deltaMatrix))/60)
    fVector <- matrixSampleGP[, indexSample]
    
    # Compute derivative matrices of f:
    # (OLD) d1Matrix -> KXpX, dMatrixAll -> KXpXp
    KXpX <- -sig2Value*kappaVal^2*auxRelativeDistanceMatrix*expMatrix # K(X*,X) - first partial derivative
    KXpXp <- sig2Value*kappaVal^2*expMatrix*(1 - kappaVal*distanceMatrix) # K(X*,X*) - second partial derivative
    #KXpXp <- sig2Value*kappaVal^2*diag(expMatrix)*(1 - kappaVal*diag(distanceMatrix)) # here we only compute the diagonal
    
    # Draw samples
    meanMVN <- KXpX%*%invDeltaMatrix%*%fVector
    iSample <- MASS::mvrnorm(n = 1, mu = KXpX%*%invDeltaMatrix%*%fVector, Sigma = KXpXp - KXpX%*%invDeltaMatrix%*%t(KXpX))
    #iSample <- sapply(1:numDays, function(x) rnorm(n = 1, mean = meanMVN[x], sd = sqrt( KXpXp[x] - KXpX[x,]%*%invDeltaMatrix%*%KXpX[x,] )))
    
    # Store
    sampleDerivatives[indexSample,] <- iSample
  }
  return(t(sampleDerivatives))
}

# Functions for running multiple groups ----

# TODO update functions in this section

#' Wrapper for runModelGrowthRate(), to run the model in multiple /groupsareas/partitions/locations if same settings/priors
#' countTableAll: data table with one row per day per group and these columns:
#' - labelPartition: name of group or 'partition'
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - numberTest: number of test on the day (>= 0)
#' - date: date of count in R date format
#' partitionTable
#' - labelPartition: name of group or 'partition'
#' - idPartition: unique identifier of each 'partition', numeric
#' - minDate: min date for which the model is going to be fitted for that partition
#' - maxDate: max date for which the model is going to be fitted for that partition
#' parametersModel: output of setParametersFn()
#' saveSamples: Default F. If T, returns matrixSampleDays_list and sampleDerivatives_list,
#'              a list of matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' runSubsetId: vector ids as in partitionTable$idPartition. Default: all ids ordered
runModelMultipleGroups <- function(countTableAll, partitionTable, parametersModel, saveSamples = F,
                                   runSubsetId = sort(partitionTable$idPartition)){
  # TODO note that try catch does not exist for unique runs
  # TODO replace "saveSamples" and other inputs for "..."
  # Run model for all areas
  output <- vector("list", length(runSubsetId))
  for(ii in 1:length(runSubsetId)){
    tryCatch({
      labelPartition_ii <- partitionTable[idPartition == runSubsetId[ii], labelPartition]
      minDateModel_ii <- partitionTable[idPartition == runSubsetId[ii], minDate]
      maxDateModel_ii <- partitionTable[idPartition == runSubsetId[ii], maxDate]
      cat(labelPartition_ii, ": ", sep = "")
      
      countTable <- countTableAll[labelPartition == labelPartition_ii & date >= minDateModel_ii & date <= maxDateModel_ii,
                                  .(date, positiveResults, numberTest)]
      inlaObject <- runModelGrowthRate(countTable = countTable,
                                         parametersModel = parametersModel,
                                         minDate = minDateModel_ii, maxDate = maxDateModel_ii)
      output[[ii]] <- processINLAOutput(objectInla = inlaObject, parametersModel = parametersModel, saveSamples = saveSamples)
      
      cat("\n")
    }, error = function(e) {cat("SKIPPED ERROR for labelPartition", labelPartition_ii, ":", conditionMessage(e), "\n")})
  }
  return(output)
}

#' Stack output of runModelMultipleGroups() into a data table.
#' Intput:
#' - output: output of runModelMultipleGroups()
#' Output:
#' - posteriorGrowth
stackOutputMultipleGroups <- function(output, partitionTable){
  # Stack output
  posteriorGrowth <- do.call("rbind", lapply(1:length(output), function(ii) output[[ii]]$posteriorGrowth[, c(.SD, idPartition = ii)]))
  posteriorTransfGP <- do.call("rbind", lapply(1:length(output), function(ii) output[[ii]]$posteriorTransfGP[, c(.SD, idPartition = ii)]))
  
  #if(saveSamples == T){
  #  matrixSampleDays_list <- lapply(1:length(output), function(ii) output[[ii]]$matrixSampleDays)
  #  sampleDerivatives_list <- lapply(1:length(output), function(ii) output[[ii]]$sampleDerivatives)
  #}
  
  setkeyv(posteriorGrowth, c("date", "idPartition"))
  setkeyv(posteriorTransfGP, c("date", "idPartition"))
  posteriorGrowth[posteriorTransfGP, ":="(median_transGP = i.median_transGP, q0.025_transGP = i.q0.025_transGP, q0.975_transGP = i.q0.975_transGP,
                                          q0.25_transGP = i.q0.25_transGP, q0.75_transGP = i.q0.75_transGP,
                                          median_GP = i.median_GP, q0.025_GP = i.q0.025_GP, q0.975_GP = i.q0.975_GP,
                                          q0.25_GP = i.q0.25_GP, q0.75_GP = i.q0.75_GP,
                                          positiveResults = i.positiveResults)]
  
  setkey(posteriorGrowth, idPartition)
  setkey(partitionTable, idPartition)
  posteriorGrowth[partitionTable, labelPartition := i.labelPartition]
  
  # Calculate probability growth rate is greater than 0. Add to posteriorGrowth
  #hasSamples <- sapply(1:length(output), function(i) !is.null(output[[i]]$matrixSampleDays))
  #if(sum(hasSamples > 0)){
  #  prob0Table <- do.call("rbind",
  #                        lapply(which(hasSamples), function(i) data.table(idPartition = i,
  #                                                                         dayId = 1:nrow((output[[i]]$sampleDerivatives)),
  #                                                                         prob0 = apply(output[[i]]$sampleDerivatives, 1, function(x) sum(x >= 0)/length(x)))))
  #  setkeyv(posteriorGrowth, c("idPartition", "dayId"))
  #  setkeyv(prob0Table, c("idPartition", "dayId"))
  #  posteriorGrowth[prob0Table, prob0 := i.prob0]
  #}
  
  return(posteriorGrowth)
  
  # Add extras to output
  #posteriorGrowth[, dayId := NULL]
  #setkey(posteriorGrowth, idPartition)
  #setkey(partitionTable, idPartition)
  #posteriorGrowth[partitionTable, ":="(labelPartition = i.labelPartition, pop = i.pop)]
  #posteriorGrowth[partitionTable, ":="(labelPartition = i.labelPartition)]
  #
  #posteriorGrowth[, ":="(median_incidencePer100TH = median_incidence*100000/pop,
  #                       q0.025_incidencePer100TH = q0.025_incidence*100000/pop,
  #                       q0.975_incidencePer100TH = q0.975_incidence*100000/pop)]
  #posteriorGrowth[, positivePer100TH := positiveResults*100000/pop]
  
  #setkey(posteriorTransfGP, idPartition)
  #setkey(partitionTable, idPartition)
  #posteriorTransfGP[partitionTable, ":="(labelPartition = i.labelPartition)]
  #posteriorTransfGP[, ratio := NA]
  
  #if(saveSamples == F){
  #return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP))
  #}else{
  # TODO
  #return(list(posteriorGrowth = posteriorGrowth, posteriorTransfGP = posteriorTransfGP,
  #            matrixSampleDays = matrixSampleDays, sampleDerivatives = t(sampleDerivatives)))
  #}
}

# Figures ----
# TODO test they are alright, they are ok
# TODO sort out axis bit, sort out week/weekend dots for non daily functions
# TODO recycle code if needed (e.g. stacking up data for plots)
# TODO choose better names
# TODO create basic format function for all plots

#' Requires: library(ggnewscale), library(scales)
#' outputModel$posteriorTransfGP, parametersModel$params$linkType
plotFitting <- function(outputModel, parametersModel){
  if(!"ggnewscale" %in% (.packages())) stop("This action requires package ggnewscale. Try library(ggnewscale)")
  if(!"scales" %in% (.packages())) stop("This action requires package scales. Try library(scales)")
  
  colourDots <- c("black", "black") # c("#D41B19", "black")
  colourFit <- c("gray", "#F5CC14")
  
  dataToPlotPosterior <- outputModel$posteriorTransfGP
  dataToPlotPosterior[, ratio := pmin(1, pmax(0, positiveResults/numberTest))]
  dataToPlotPosterior[, isWeekend := weekdays(date) %in% c("Saturday", "Sunday")]
  dataToPlotPosterior[, weekendLabel := factor(ifelse(isWeekend, "weekend", "weekday"), levels = c("weekend", "weekday"))]
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio, positiveResults, median = NA, q0.025 = NA, q0.975 = NA, type = "points")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA, median = medianFT, q0.025 = q0.025FT, q0.975 = q0.975FT, type = "model fit")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA, median = median_transGP, q0.025 = q0.025_transGP, q0.975 = q0.975_transGP, type = "Gaussian process")])
  dataToPlot[, typeLevel := factor(type, levels = c("model fit", "Gaussian process"))]
  p01 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_ribbon(data = dataToPlot[typeLevel != "points"], aes(ymin = q0.025, ymax = q0.975, fill = typeLevel), alpha = 0.5) +
    geom_line(data = dataToPlot[typeLevel != "points"], aes(y = median, colour = typeLevel, linetype = typeLevel)) +
    scale_colour_manual(name = "posterior (median, 95% CI)", values = colourFit) +
    scale_fill_manual(name = "posterior (median, 95% CI)", values = colourFit) +
    scale_linetype_manual(name = "posterior (median, 95% CI)", values = c(2, 1)) +
    #scale_x_date(labels = scales::label_date(formatBreaks), breaks = dateBreaks2) +
    scale_x_date(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(),
          legend.key = element_blank(), panel.grid.major.x = element_line(linetype = 2, colour = "gray90")) # 0.6
  if(parametersModel$params$linkType == "BB"){
    p02 <- p01 +
      new_scale_color() +
      scale_colour_manual(name = "proportion of positives", values = colourDots) +
      scale_shape_discrete(name = "proportion of positives") +
      geom_point(aes(y = ratio, colour = weekendLabel, shape = weekendLabel)) +
      labs(x = "day", y = "model posterior\nand proportion of positives") +
      guides(shape = guide_legend(order = 1), colour = guide_legend(order = 1))
  }else{
    p02 <- p01 +
      new_scale_color() +
      scale_colour_manual(name = "count of positives", values = colourDots) +
      scale_shape_discrete(name = "count of positives") +
      geom_point(aes(y = positiveResults, colour = weekendLabel, shape = weekendLabel)) +
      scale_y_continuous(labels = unit_format(unit = "", scale = 1e-3)) + #unit = "K"
      labs(x = "day", y = "model posterior\nand count of positives (thousands)") +
      guides(shape = guide_legend(order = 1), colour = guide_legend(order = 1))
  }
  
  return(p02)
}

#' outputModel$posteriorGrowth
plotGR <- function(outputModel){
  dataToPlotPosterior <- outputModel$posteriorGrowth
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, median = median, qlow = NA, qhigh = NA, type = "median")],
                       dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.025, qhigh = q0.975, type = "95% CI")],
                       dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.25, qhigh = q0.75, type = "50% CI")])
  dataToPlot[, typeLevel := factor(type, levels = c("median", "50% CI", "95% CI"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_ribbon(aes(ymin = qlow, ymax = qhigh, fill = typeLevel), alpha = 0.5) +
    geom_line(aes(y = median, colour = typeLevel)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
    scale_colour_manual(name = "posterior", values = c("black", "gray40", "gray70")) +
    scale_fill_manual(name = "posterior", values = c("black", "gray40", "gray70"),
                      guide = guide_legend(override.aes = list(colour = c("black", NA, NA), fill = c(NA, "gray40", "gray70")))) + # trick for legend
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "growth rate posterior") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank())
  return(p0)
}

#' INTERNAL
#' Requires library(mvtnorm)
#' outputModel$matrixSampleHyperAll parametersModel$params (theta.prior2.mean theta.prior2.prec prior2.sigma0 prior2.range0)
#' TODO remove week/weekend colour if week
plotHyperparametersGP <- function(outputModel, parametersModel, seed = 123){
  if(!"mvtnorm" %in% (.packages())) stop("This action requires package mvtnorm. Try library(mvtnorm)")
  
  matrixSampleHyper <- outputModel$matrixSampleHyperAll
  
  # Get random samples for visualisation
  set.seed(seed)
  samplesPrior <- mvtnorm::rmvnorm(n = 100,
                                   mean = parametersModel$params$theta.prior2.mean,
                                   sigma = solve(parametersModel$params$theta.prior2.prec))
  samplesPosterior <- data.table(sigma = parametersModel$params$prior2.sigma0*exp(matrixSampleHyper["theta1",]),
                                 range = parametersModel$params$prior2.range0*exp(matrixSampleHyper["theta2",]))
  
  plotHyperparameters <- rbind(data.table(samplesPosterior[,.N,.(sigma, range)],
                                          type = "posterior"),
                               data.table(sigma = parametersModel$params$prior2.sigma0*exp(samplesPrior[,1]),
                                          range = parametersModel$params$prior2.range0*exp(samplesPrior[,2]),
                                          N = 1,
                                          type = "prior"))
  plotHyperparameters[, typeLevel := factor(type, levels = c("prior", "posterior"))]
  p0 <- ggplot(plotHyperparameters, aes(x = sigma, y = range/2, colour = typeLevel, shape = typeLevel)) + theme_laura() + geom_point() +
    scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') + stat_ellipse(type = "norm", linewidth = 0.5) +
    scale_colour_manual(name = "", values = c("black", "red")) + scale_shape_manual(name = "", values = c(16, 17)) +
    theme(legend.title = element_blank(), legend.key = element_blank(), #legend.position = c(0.82, 0.07),
          legend.background = element_rect(fill = NA),
          legend.key.height = unit(0.1, "cm")) +
    labs(x = expression(standard~deviation~sigma), y = expression(length~scale~plain(l)~(days)))
  return(p0)
}

#' INTERNAL
plotHyperparametersDispersion <- function(outputModel, parametersModel, log = T){
  if(parametersModel$params$linkType == "NB"){
    samplesPrior <- rnorm(n = 300,
                          mean = parametersModel$params$NB.prior.rho$size$param[1],
                          sd = sqrt(1/parametersModel$params$NB.prior.rho$size$param[2]))
    samplesPriorTransf <- exp(samplesPrior)
  }else{
    samplesPrior <- rnorm(n = 300,
                          mean = parametersModel$params$NB.prior.rho$overdispersion$param[1],
                          sd = sqrt(1/parametersModel$params$BB.prior.rho$overdispersion$param[2]))
    samplesPriorTransf <- exp(samplesPrior)/(1 + exp(samplesPrior))
  }
  
  dataToPlot <- rbind(data.table(samples = outputModel$matrixSampleHyper[c("overdispersion"),], type = "posterior"),
                      data.table(samples = samplesPriorTransf, type = "prior"))
  
  p0 <- ggplot(dataToPlot, aes(x = samples, group = type)) + theme_laura() +
    geom_density(aes(fill = type), alpha = 0.6) +
    #geom_point(aes(x = samples, y = 0, colour = type)) +
    scale_fill_manual(values = c("#305D72", "#C36577")) +
    labs(x = ifelse(parametersModel$params$linkType == "NB", expression(1/overdispersion~eta), expression(overdispersion~rho)), y = "density") +
    scale_x_continuous(trans = ifelse(!log, "identity", "log10"))
  return(p0)
}

#' INTERNAL
plotHyperparametersPrecision <- function(outputModel, parametersModel, log = T){
  samplesPrior <- rgamma(n = 300,
                         shape = parametersModel$params$dayWeek.prior.prec$theta$param[1],
                         rate = parametersModel$params$dayWeek.prior.prec$theta$param[2])
  dataToPlot <- rbind(data.table(samples = outputModel$matrixSampleHyper[c("precision"),], type = "posterior"),
                      data.table(samples = samplesPrior, type = "prior"))
  
  p0 <- ggplot(dataToPlot, aes(x = samples, group = type)) + theme_laura() +
    geom_density(aes(fill = type), alpha = 0.6) +
    #geom_point(aes(x = samples, y = 0, colour = type)) +
    scale_fill_manual(values = c("#305D72", "#C36577")) +
    labs(x = expression(precision~tau[w]), y = "density") +
    scale_x_continuous(trans = ifelse(!log, "identity", "log10"))
  
  #dataToPlot <- data.table(samples = outputModel$marginalsHyper$`Precision for dayWeek`)
  #p0 <- ggplot(dataToPlot, aes(x = samples)) + theme_laura() +
  #  geom_density(alpha = 0.6, fill = "gray90") +
  #  geom_point(aes(x = samples, y = 0)) +
  #  labs(x = expression(precision~tau[w]), y = "density")
  #coord_cartesian(xlim = c(0, 200)) +
  #theme(legend.position = c(0.85,0.85), legend.title = element_blank())
  
  return(p0)
}

plotHyperparametersAll <- function(outputModel, parametersModel){
  p1 <- plotHyperparametersGP(outputModel, parametersModel)
  p2 <- plotHyperparametersDispersion(outputModel, parametersModel)
  p3 <- plotHyperparametersPrecision(outputModel, parametersModel)
  pAll <- ggarrange(p1, p2, p3, ncol = 3)
  return(pAll)
}

printHyperparametersSummary <- function(outputModel, parametersModel){
  # TODO ?
  summary(parametersModel$params$prior2.sigma0*exp(outputModel$matrixSampleHyperAll[c("theta1"),]))
  summary(parametersModel$params$prior2.range0*exp(outputModel$matrixSampleHyperAll[c("theta2"),]))
  summary(parametersModel$params$prior2.sigma0*outputModel$matrixSampleHyperAll[c("overdispersion"),])
}

plotFittingSamples <- function(outputModel, parametersModel){
  # TODO !!!!! check
  # TODO do dots for BB
  samples <- outputModel$matrixSampleDays
  rownames(samples) <- NULL
  dataSamplesGP <- data.table(suppressWarnings(melt(samples, varnames = c("dayId", "sample"))))
  
  setkey(dataSamplesGP, dayId)
  setkey(outputModel$dataForModel, dayId)
  dataSamplesGP[outputModel$dataForModel, date := i.date]
  
  if(parametersModel$params$linkType == "BB"){
    dataSamplesGP[, valueTr := exp(value)/(1 + exp(value))]
  }else if(parametersModel$params$linkType == "NB"){
    dataSamplesGP[, valueTr := exp(value)]
  }
  randomCols <- sample(ncol(outputModel$matrixSampleDays)) # for colour of samples
  dataSamplesGP[, colourId := randomCols[sample]]
  
  dataToPlotPosterior <- outputModel$posteriorTransfGP
  dataToPlotPosterior[, ratio := pmin(1, pmax(0, positiveResults/numberTest))]
  dataToPlotPosterior[, isWeekend := weekdays(date) %in% c("Saturday", "Sunday")]
  dataToPlotPosterior[, weekendLabel := factor(ifelse(isWeekend, "weekend", "weekday"), levels = c("weekend", "weekday"))]
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio, positiveResults, median = NA, q0.025 = NA, q0.975 = NA, type = "points")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA, median = medianFT, q0.025 = q0.025FT, q0.975 = q0.975FT, type = "model fit")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA, median = median_transGP, q0.025 = q0.025_transGP, q0.975 = q0.975_transGP, type = "Gaussian process")])
  dataToPlot[, typeLevel := factor(type, levels = c("model fit", "Gaussian process"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_line(data = dataSamplesGP, aes(y = valueTr, group = sample, colour = colourId), alpha = 0.5) + # colour = "gray50"
    #geom_ribbon(data = dataToPlot[typeLevel == "Gaussian process"], aes(ymin = q0.025, ymax = q0.975, fill = typeLevel), alpha = 0.1) +
    geom_line(data = dataToPlot[typeLevel == "Gaussian process"], aes(y = median), linetype = 2, colour = "#F5CC14") +
    geom_point(aes(y = positiveResults), colour = "red", size = 0.5) +
    scale_colour_gradient(guide = "none", low = "gray90", high = "gray10") +
    #scale_colour_manual(name = "posterior (median, 95% CI)", values = "#F5CC14") +
    #scale_fill_manual(name = "posterior (median, 95% CI)", values = "#F5CC14") +
    #scale_x_date(labels = scales::label_date(formatBreaks), breaks = dateBreaks2) +
    scale_x_date(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(),
          legend.key = element_blank(), panel.grid.major.x = element_line(linetype = 2, colour = "gray90")) # 0.6
  
  return(p0)
}


plotFittingSamplesGR <- function(outputModel, parametersModel){
  # TODO !!!!! check
  # TODO do dots for BB
  samples <- outputModel$sampleDerivatives
  rownames(samples) <- NULL
  dataSamplesGR <- data.table(suppressWarnings(melt(samples, varnames = c("dayId", "sample"))))
  
  setkey(dataSamplesGR, dayId)
  setkey(outputModel$dataForModel, dayId)
  dataSamplesGR[outputModel$dataForModel, date := i.date]
  
  
  
  dataToPlotPosterior <- outputModel$posteriorGrowth
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, median = median, qlow = NA, qhigh = NA, type = "median")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.025, qhigh = q0.975, type = "95% CI")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.25, qhigh = q0.75, type = "50% CI")])
  dataToPlot[, typeLevel := factor(type, levels = c("median", "50% CI", "95% CI"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    #geom_ribbon(aes(ymin = qlow, ymax = qhigh, fill = typeLevel), alpha = 0.5) +
    geom_line(data = dataSamplesGR, aes(y = value, group = sample), colour = "gray50", alpha = 0.5) +
    geom_line(aes(y = median, colour = typeLevel)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
    scale_colour_manual(name = "posterior", values = c("black", "gray40", "gray70")) +
    scale_fill_manual(name = "posterior", values = c("black", "gray40", "gray70"),
                      guide = guide_legend(override.aes = list(colour = c("black", NA, NA), fill = c(NA, "gray40", "gray70")))) + # trick for legend
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "growth rate posterior") +
    scale_x_date(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank(), panel.grid.major.x = element_line(linetype = 2, colour = "gray90"))
  return(p0)
}

#' outputModel$
plotLatent <- function(outputModel){
  dataToPlotPosterior <- outputModel$posteriorTransfGP
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, median = median_GP, qlow = NA, qhigh = NA, type = "median")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.025_GP, qhigh = q0.975_GP, type = "95% CI")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.25_GP, qhigh = q0.75_GP, type = "50% CI")])
  dataToPlot[, typeLevel := factor(type, levels = c("median", "50% CI", "95% CI"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_ribbon(aes(ymin = qlow, ymax = qhigh, fill = typeLevel), alpha = 0.5) +
    geom_line(aes(y = median, colour = typeLevel)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
    scale_colour_manual(name = "posterior", values = c("black", "gray40", "gray70")) +
    scale_fill_manual(name = "posterior", values = c("black", "gray40", "gray70"),
                      guide = guide_legend(override.aes = list(colour = c("black", NA, NA), fill = c(NA, "gray40", "gray70")))) + # trick for legend
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    #scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "GP in latent space") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank())
  return(p0)
}

plotLatent <- function(outputModel){
  dataToPlotPosterior <- outputModel$posteriorTransfGP
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, median = median_GP, qlow = NA, qhigh = NA, type = "median")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.025_GP, qhigh = q0.975_GP, type = "95% CI")],
                      dataToPlotPosterior[, .(dayId, date, median = NA, qlow = q0.25_GP, qhigh = q0.75_GP, type = "50% CI")])
  dataToPlot[, typeLevel := factor(type, levels = c("median", "50% CI", "95% CI"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_ribbon(aes(ymin = qlow, ymax = qhigh, fill = typeLevel), alpha = 0.5) +
    geom_line(aes(y = median, colour = typeLevel)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
    scale_colour_manual(name = "posterior", values = c("black", "gray40", "gray70")) +
    scale_fill_manual(name = "posterior", values = c("black", "gray40", "gray70"),
                      guide = guide_legend(override.aes = list(colour = c("black", NA, NA), fill = c(NA, "gray40", "gray70")))) + # trick for legend
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    #scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "GP in latent space") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank())
  return(p0)
}

plotRandomEffect <- function(outputModel, parametersModel){
  dataToPlotPosterior <- outputModel$posteriorRandomEffect
  dataToPlot <- rbind(dataToPlotPosterior[, .(index, median = median, qlow = NA, qhigh = NA, type = "median")],
                      dataToPlotPosterior[, .(index, median = NA, qlow = q0.025, qhigh = q0.975, type = "95% CI")],
                      dataToPlotPosterior[, .(index, median = NA, qlow = q0.25, qhigh = q0.75, type = "50% CI")])
  dataToPlot[, typeLevel := factor(type, levels = c("median", "50% CI", "95% CI"))]
  p0 <- ggplot(dataToPlot, aes(x = index)) + theme_laura() +
    geom_segment(aes(xend = index, y = qlow, yend = qhigh, colour = typeLevel), linetype = 1) +
    geom_point(aes(y = median, colour = typeLevel)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
    scale_colour_manual(name = "posterior", values = c("black", "gray40", "gray70")) +
    scale_fill_manual(name = "posterior", values = c("black", "gray40", "gray70"),
                      guide = guide_legend(override.aes = list(colour = c("black", NA, NA), fill = c(NA, "gray40", "gray70")))) + # trick for legend
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))), breaks = dateBreaks) +
    #scale_y_continuous(labels = scales::percent_format(accuracy = NULL)) +
    labs(x = "day", y = "random effect") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank())
  if(parametersModel$params$randomEffect == "weekday"){
    levelsWeek <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
    p0 <- p0 + scale_x_continuous(breaks = 1:7, labels = levelsWeek)
  }
  return(p0)
}
