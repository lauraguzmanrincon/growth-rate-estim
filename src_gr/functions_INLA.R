
#' Title
#'
#' @param dataForModel 
#' @param dateList 
#' @param parametersModel output of setModelParameters
#' @param inferenceSettings output of setInferenceSettings
#'
#' @return objectInla
#' @export
#'
#' @examples
runModelGrowthRate_INLA <- function(dataForModel, dateList, parametersModel, inferenceSettings){
  
  internalConstants <- getInternalSettings()
  
  # ---------------------------------------------------- #
  #                      FIT MODEL                       #
  # ---------------------------------------------------- #
  
  # Define grid for finite differences
  # TODO generalise boundary construction?
  boundaryVertices <- 200 # Points on the left and right
  boundaryStepSize <- 10 # Distance between boundary points
  boundaryPoints <- c(seq(dateList$minDay - boundaryStepSize*boundaryVertices, dateList$minDay - boundaryStepSize, by = boundaryStepSize),
                      dateList$dateTable[order(dayId), dayId], 
                      seq(dateList$maxDay + boundaryStepSize, dateList$maxDay + boundaryStepSize*boundaryVertices, by = boundaryStepSize))
  nonBoundaryIndices <- (boundaryVertices + 1):(length(boundaryPoints) - boundaryVertices)
  
  # Create mesh in 1d (mesh1d), projection matrix (A1), model (spde1), and indexes (spde1.idx)
  
  # - Mesh and projection matrix
  mesh1d <- inla.mesh.1d(boundaryPoints, boundary = "free") # NEW 09.02.2021 boundary issues
  A1 <- inla.spde.make.A(mesh1d, dataForModel[order(dayId), dayId])
  
  # - Set priors
  # Prior for day of the week
  priorRandomEffect <- list(theta = list(prior = parametersModel$randomEffectPrior$prior,
                                         param = c(parametersModel$randomEffectPrior$a, parametersModel$randomEffectPrior$b)))
  
  # Prior for Gaussian process parameters, following (Lindgren, 2015, v63i19)
  vGP <- 2 - 1/2
  sigma0 <- parametersModel$GPcovarianceList$theta0[1]
  range0 <- 2*parametersModel$GPcovarianceList$theta0[2]
  # Convert into tau and kappa:
  kappa0 <- sqrt(8*vGP)/range0
  tau0 <- sqrt(gamma(vGP)/(gamma(2)*sqrt(4*pi)))/(kappa0^(vGP)*sigma0)
  basisPriorTau <- c(-1, 3/2)
  basisPriorKappa <- c(0, -1)
  spde1 <- inla.spde2.matern(mesh1d,
                             B.tau = cbind(log(tau0), basisPriorTau[1], basisPriorTau[2]),
                             B.kappa = cbind(log(kappa0), basisPriorKappa[1], basisPriorKappa[2]),
                             theta.prior.mean = c(0,0),
                             theta.prior.prec = parametersModel$GPHyperparamPrior$B,
                             constr = T) # NEW 21.09.2023
  
  # - Create stack data
  spde1.idx <- inla.spde.make.index("day", n.spde = spde1$n.spde)
  # Create stacked data with two datasets (observations and prediction). For paper, I'm ignoring prediction.
  if(parametersModel$randomEffect == "weekday"){ # NEW 22.09.2023' - correct extra information done at the beginning
    A <- list(A1, 1)
    effectsList <- list(spde1.idx,
                        list(randomEffect = factor(dataForModel[order(dayId), dayWeek], levels = internalConstants$levelsWeek)))
  }else if(parametersModel$randomEffect == "all"){
    A <- list(A1, 1)
    effectsList <- list(spde1.idx,
                        list(randomEffect = dataForModel[order(dayId), dayId]))
  }else{
    A <- list(A1)
    effectsList <- list(spde1.idx)
  }
  stackD <- inla.stack(data = list(positiveResults = dataForModel[order(dayId), positiveResults]),
                       A = A,
                       effects = effectsList,
                       tag = "est")
  # Create prediction data for the next parametersModel$config$sizeGPProjection days
  sizeGPProjection <- 10
  daysProjection <- sizeGPProjection + 1
  predictionGPWeekday <- dateList$dateTable[dayId == max(dataForModel$dayId), weekdays(date + 1:daysProjection)]
  xx <- seq(max(dataForModel$dayId) + 1, max(dataForModel$dayId) + daysProjection, by = 1)
  A.xx <- inla.spde.make.A(mesh1d, xx)
  stack.pred <- inla.stack(data = list(positiveResults = NA),
                           A = list(A.xx),
                           effects = list(spde1.idx), # NEW 22.09.2023' - correct extra information done at the beginning
                           tag = "pred")
  joint.stack <- inla.stack(stackD, stack.pred)
  
  # Fit model (using INLA)
  cat("Fitting model ... ")
  constantTerm <- 1
  if(parametersModel$randomEffect %in% c("weekday", "all")){
    formulaText <- paste0("positiveResults ~ ", constantTerm, " + f(day, model = spde1) + f(randomEffect, model = 'iid', hyper = priorRandomEffect, constr = T)")
  }else{
    formulaText <- paste0("positiveResults ~ ", constantTerm, " + f(day, model = spde1)")
  }
  formula <- as.formula(formulaText)
  
  if(parametersModel$linkType == "BB"){
    # Betabinomial
    priorOverdispersion <- list(overdispersion = list(prior = parametersModel$dispersionPrior$prior,
                                                      param = c(parametersModel$dispersionPrior$mean,
                                                                parametersModel$dispersionPrior$prec)))
    objectInla <- inla(formula = formula, data = inla.stack.data(joint.stack),
                       control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                       family = "betabinomial",
                       Ntrials =  c(dataForModel[order(dayId), numberTest], rep(NA, length(xx))),
                       control.compute = list(config = TRUE),
                       control.family = list(hyper = priorOverdispersion))
  }else if(parametersModel$linkType == "NB"){
    # Negative binomial
    priorOverdispersion <- list(size = list(prior = parametersModel$dispersionPrior$prior,
                                            param = c(parametersModel$dispersionPrior$mean,
                                                      parametersModel$dispersionPrior$prec)))
    objectInla <- inla(formula = formula, data = inla.stack.data(joint.stack),
                       control.predictor = list(A = inla.stack.A(joint.stack), compute = TRUE, link = 1),
                       family = "nbinomial",
                       control.compute = list(config = TRUE),
                       control.family = list(hyper = priorOverdispersion))
  }
  
  objectInla$dateList <- dateList
  objectInla$dataForModel <- dataForModel
  objectInla$nonBoundaryIndices <- nonBoundaryIndices # TODO better way?
  return(objectInla)
}

#' Title
#'
#' @param objectInla 
#' @param parametersModel 
#' @param inferenceSettings 
#' @param saveSamples (T) only of internal use?
#'
#' @return
#' @export
#'
#' @examples
processINLAOutput <- function(objectInla, parametersModel, inferenceSettings, saveSamples = T, saveInlaObject = F){
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  cat("Estimating growth rate ... ")
  
  computeGPProjection <- F
  hasConstant <- T
  sizeGPProjection <- 10
  
  # Get samples of posterior distribution of parameters
  sampleList <- inla.posterior.sample(inferenceSettings$numSamples, objectInla)
  
  # Extract relevant indexes from output (see unique(sapply(strsplit(rownms, ":"), function(x) x[[1]])))
  rownms <- rownames(sampleList[[1]]$latent)
  set1 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "APredictor"))
  #set2 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "Predictor"))
  set3 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "day"))
  set4 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "randomEffect"))
  if(hasConstant) set5 <- sort(which(sapply(strsplit(rownms, ":"), function(x) x[[1]]) == "(Intercept)"))
  dayIndexInSample <- set3[objectInla$nonBoundaryIndices]
  
  # Create matrix of samples (and hyperparameters/weekday effect if needed)
  matrixSampleGP <- sapply(sampleList, function(x) x$latent[dayIndexInSample,1])
  if(parametersModel$randomEffect != "none"){
    matrixSampleRandomEffect <- sapply(sampleList, function(x) x$latent[set4,1])
  }else{
    matrixSampleRandomEffect <- NULL # TODO ??
  }
  matrixSampleNu <- sapply(sampleList, function(x) x$latent[set1[1:objectInla$dateList$numDays],1])
  if(hasConstant){
    matrixSampleIntercept <- (1 + (parametersModel$randomEffect != "none"))*sapply(sampleList, function(x) x$latent[set5,1]) # TODO why x2?
  }else{
    matrixSampleIntercept <- rep(0, inferenceSettings$numSamples)
  }
  
  nameDispersionINLA <- ifelse(parametersModel$linkType == "NB",
                               "size for the nbinomial observations (1/overdispersion)",
                               "overdispersion for the betabinomial observations")
  namePrecisionINLA <- ifelse(parametersModel$randomEffect %in% c("weekday", "all"), "Precision for randomEffect", "NA")
  matrixSampleHyperparameters <- sapply(sampleList, function(x) x$hyperpar[c(nameDispersionINLA, namePrecisionINLA,
                                                                      "Theta1 for day", "Theta2 for day")]) # NEW 08.09.2023
  rownames(matrixSampleHyperparameters) <- c("overdispersion", "precision", "logParam1", "logParam2")
  
  # Create matrix with GP projection (sizeGPProjection + 1 as we lose one in the GR estimation)
  daysProjection <- sizeGPProjection + 1
  predicIndexInSample <- set1[objectInla$dateList$numDays + (1:daysProjection)]
  projectionGP <- sapply(sampleList, function(x) x$latent[predicIndexInSample,1]) # previous projectionGPInla
  
  # Compute derivative
  if(inferenceSettings$derivativeFromGP == T | nrow(objectInla$dateList$dateTable) < 7){
    # Sample from derivative
    matrixSampleGPDerivative <- getSamplesGPDerivative_GPmatern12(matrixSampleGP = matrixSampleGP,
                                                         samplesHyperparam = matrixSampleHyperparameters,
                                                         sigma0 = parametersModel$GPcovarianceList$theta0[1],
                                                         range0 = 2*parametersModel$GPcovarianceList$theta0[2])
  }else{
    # Compute approximate derivative using windows (+-3)
    # Compute +- 3 days window derivative (+-2 and +-1 for the third and second last point)
    matrixSampleGPDerivative <- getSamplesGPDerivative_approx(matrixSampleGP = matrixSampleGP)
  }
  
  vectorTestData <- objectInla$dataForModel[order(dayId), numberTest] # 27.11.2023
  
  # ---------------------------------------------------- #
  #                  PRODUCE OUTPUT                      #
  # ---------------------------------------------------- #
  listPosteriors <- computeSummaryPosteriors(matrixSampleGP, matrixSampleGPDerivative, matrixSampleHyperparameters,
                                      matrixSampleNu, matrixSampleRandomEffect, matrixSampleIntercept,
                                      vectorTestData, parametersModel)
  
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
  output_main <- list(posteriorGrowth = listPosteriors$posteriorGrowth,
                      posteriorTransfGP = listPosteriors$posteriorTransfGP,
                      posteriorRandomEffect = listPosteriors$posteriorRandomEffect,
                      posteriorIntercept = listPosteriors$posteriorIntercept,
                      posteriorHyperparameters = listPosteriors$posteriorHyperparameters,
                      dateList = objectInla$dateList,
                      dataForModel = objectInla$dataForModel,
                      parametersModel = parametersModel,
                      inferenceSettings = inferenceSettings,
                      objectInla = objectInla)
  output_samples <- list(numSamples = inferenceSettings$numSamples,
                         matrixSampleGP = matrixSampleGP,
                         matrixSampleGPDerivative = matrixSampleGPDerivative,
                         matrixSampleRandomEffect = matrixSampleRandomEffect,
                         matrixSampleIntercept = matrixSampleIntercept,
                         matrixSampleHyperparameters = matrixSampleHyperparameters)
  #output_projection <- list(projectionGP = projectionGP)
  
  # Remove redundant data
  output_main$objectInla$dateList <- NULL
  output_main$objectInla$dataForModel <- NULL
  if(saveInlaObject == F) output_main$objectInla <- NULL
  
  if(saveSamples == F){
    output <- output_main
  }else{
    output <- c(output_main, output_samples)
  }
  #if(computeGPProjection == T) output <- c(output, output_projection)
  return(output)
}
