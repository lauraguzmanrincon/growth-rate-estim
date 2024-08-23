
#' Title
#'
#' @param dataForModel 
#' @param dateList 
#' @param parametersModel output of setModelParameters
#' @param inferenceSettings output of setInferenceSettings
#'
#' @return ...returns objectStan (includes dataForModel)
#'         ... dataForModel: this data includes all days in range. Recall STAN does not recieve the ones without positives and tests
#' @export
#'
#' @examples
runModelGrowthRate_STAN <- function(dataForModel, dateList, parametersModel, inferenceSettings){
  
  internalConstants <- getInternalSettings()
  
  # ---------------------------------------------------- #
  #                   SHAPE DATA FOR STAN                #
  # ---------------------------------------------------- #
  
  if(parametersModel$linkType == "BB"){
    # TODO what if NA in numberTest for BB???
    dataForModel[, includedInStan := !is.na(positiveResults) & !is.na(numberTest)]
    dataForModelNoNa <- dataForModel[includedInStan == T]
  }else{
    dataForModel[, includedInStan := !is.na(positiveResults)]
    dataForModelNoNa <- dataForModel[includedInStan == T]
    dataForModelNoNa[, numberTest := 0]
  }
  
  # ---------------------------------------------------- #
  
  # ---------------------------------------------------- #
  #                        STAN                          #
  # ---------------------------------------------------- #
  
  # Load executable Stan model
  # TODO how to set dirSource inside package?
  # TODO can we do this without relying on saving the exec in Globals?
  dirSource <- "/Users/lauraguzmanrincon/Documents/GitHub/growth-rate-estim/src_gr"
  constructStanExec(linkType = parametersModel$linkType, dirSource = dirSource)
  if(parametersModel$linkType == "BB"){
    modelExec <- modelExec_BB
  }
  else{
    modelExec <- modelExec_NB
  }
  
  # Data to Stan
  # TODO generalise
  # We assume dayId runs from 1 in steps of unitTime
  modelData <- list(
    num_days = nrow(dataForModelNoNa), #as.integer(numDays),
    num_groups = as.integer(7),
    # data
    t = dataForModelNoNa[order(dayId), dayId], #1:numDays,
    pos = dataForModelNoNa[order(dayId), positiveResults],
    tot = dataForModelNoNa[order(dayId), as.integer(numberTest)],
    day_to_group = as.integer(dataForModelNoNa[order(dayId), factor(dayWeek, levels = internalConstants$levelsWeek)]),
    # prior intercept
    m_int = parametersModel$interceptPrior$mean,
    sig_int = 1/sqrt(parametersModel$interceptPrior$prec),
    # prior dispersion
    m_eta = parametersModel$dispersionPrior$mean, # positives
    sig_eta = 1/sqrt(parametersModel$dispersionPrior$prec), # positives
    m_rho = parametersModel$dispersionPrior$mean, # proportions
    sig_rho = 1/sqrt(parametersModel$dispersionPrior$prec), # proportions
    # prior day-of-the-week effect
    a_w = parametersModel$randomEffectPrior$a,
    b_w = parametersModel$randomEffectPrior$b,
    # prior GP
    log_theta_0 = 0 + log(parametersModel$GPHyperparamPrior$mean),
    B = parametersModel$GPHyperparamPrior$B)
    
  
  # Initial values
  set.seed(inferenceSettings$seed)
  if(parametersModel$linkType == "NB"){
    x_t_notTrasnform <- log(modelData$pos)
    modelInits <- replicate(inferenceSettings$numChains, list(intercept = mean(log(modelData$pos)),
                                                        eta = exp(rnorm(n = 1, mean = modelData$m_eta, sd = modelData$sig_eta)),
                                                        x_t = x_t_notTrasnform - mean(x_t_notTrasnform),
                                                        #TODO
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(parametersModel$GPHyperparamPrior$mean) + 0,
                                                                                       sigma = parametersModel$GPHyperparamPrior$B)[1,],
                                                        w_d = rep(0, length(internalConstants$levelsWeek)),
                                                        tau_w = rgamma(1, shape = modelData$a_w, rate = modelData$b_w)), simplify = F)
  }else if(parametersModel$linkType == "BB"){
    logit <- function(p) log(p/(1 - p))
    inv.logit <- function(x) exp(x)/(1 + exp(x))
    x_t_notTrasnform <- logit(pmin(pmax(modelData$pos, 1), modelData$tot - 1)/modelData$tot)
    modelInits <- replicate(inferenceSettings$numChains, list(intercept = rnorm(n = 1, mean = modelData$m_int, sd = modelData$sig_int),
                                                        rho = inv.logit(rnorm(n = 1, mean = modelData$m_rho, sd = modelData$sig_rho)), #modelData$m_rho OR -2?
                                                        x_t = x_t_notTrasnform - mean(x_t_notTrasnform),
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(parametersModel$GPHyperparamPrior$mean) + 0,
                                                                                       sigma = parametersModel$GPHyperparamPrior$B)[1,],
                                                        w_d = rep(0, length(internalConstants$levelsWeek)),
                                                        tau_w = rgamma(1, shape = modelData$a_w, rate = modelData$b_w)), simplify = F)
  }
  
  # Run model
  cat("Running model ... \n")
  modelFit <- sampling(modelExec,
                       data = modelData,
                       init = modelInits,
                       chains = inferenceSettings$numChains,
                       iter = inferenceSettings$numIterations,
                       warmup = inferenceSettings$warmUp,
                       thin = inferenceSettings$thinning,
                       seed = c(inferenceSettings$seed + 1),
                       sample_file = inferenceSettings$sampleFile,
                       include = T,
                       pars = c("x_transf", "w_d", "intercept", ifelse(parametersModel$linkType == "BB", "rho", "eta"), "tau_w", "log_theta_x"))
  # TODO what does STAN do with NA values?
  
  # Create output
  modelStanc <- modelExec@model_code
  objectStan <- list(modelFit = modelFit,
                     modelStanc = modelStanc,
                     dataForModel = dataForModel,
                     dateList = dateList)
  
  #cat("Saving objectStan in ", inferenceSettings$sampleFile, ".RData\n", sep = "")
  #save(objectStan, file = paste0(parametersStan$sampleFile, ".RData")) # 30.11.2023
  
  return(objectStan)
}


#' Title
#'
#' @param objectStan 
#' @param parametersModel 
#' @param inferenceSettings 
#' @param saveSamples 
#'
#' @return
#' @export
#'
#' @examples
processSTANOutput <- function(objectStan, parametersModel, inferenceSettings, saveSamples = T, saveStanObject = F){
  # TODO process less samples than iterations??
  
  set.seed(12345)
  samplesToExtract <- 1000 # TODO input? # getting random subset
  internalConstants <- getInternalSettings()
  
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  cat("Estimating growth rate ... ")
  
  # Samples GP
  #outputFit <- data.table(summary(objectStan$modelFit)$summary, keep.rownames = T)
  samplesFit <- extract(objectStan$modelFit)
  subsetSamples <- sample(x = length(samplesFit$intercept), size = min(length(samplesFit$intercept), samplesToExtract))
  matrixSampleGP <- t(samplesFit$x_transf)
  
  # Samples derivative
  # TODO Condition this calculation on the value of the parameter derivativeFromGP
  matrixSampleGPDerivative <- getSamplesGPDerivative_approx(matrixSampleGP = matrixSampleGP)
  
  # Samples parameters
  matrixSampleRandomEffect <- t(samplesFit$w_d)
  if(parametersModel$linkType == "BB"){
    matrixSampleOverdisp <- samplesFit$rho
  }else{
    matrixSampleOverdisp <- samplesFit$eta
  }
  matrixSampleIntercept <- samplesFit$intercept
  orderDayWeek <- as.integer(objectStan$dataForModel[includedInStan == T][order(dayId), sapply(dayWeek, function(dw) which(dw == internalConstants$levelsWeek))])
  matrixSampleNu <- t(c(matrixSampleIntercept) + t(matrixSampleGP) + t(matrixSampleRandomEffect)[, orderDayWeek])
  vectorTestData <- objectStan$dataForModel[includedInStan == T][order(dayId), numberTest]
  matrixSampleHyperparameters <- rbind(matrixSampleOverdisp,
                                       samplesFit$tau_w,
                                       t(samplesFit$log_theta_x))
  rownames(matrixSampleHyperparameters) <- c("overdispersion", "precision", paste("logParam", 1:ncol(samplesFit$log_theta_x), sep = ""))
  
  # ---------------------------------------------------- #
  #                  PRODUCE OUTPUT                      #
  # ---------------------------------------------------- #
  listPosteriors <- computeSummaryPosteriors(matrixSampleGP, matrixSampleGPDerivative, matrixSampleHyperparameters,
                                             matrixSampleNu, matrixSampleRandomEffect, matrixSampleIntercept,
                                             vectorTestData, parametersModel)
  
  setkey(listPosteriors$posteriorGrowth, dayId)
  setkey(objectStan$dateList$dateTable, dayId)
  listPosteriors$posteriorGrowth[objectStan$dateList$dateTable, ":="(date = i.date)]
  
  setkey(listPosteriors$posteriorTransfGP, dayId)
  setkey(objectStan$dateList$dateTable, dayId)
  listPosteriors$posteriorTransfGP[objectStan$dateList$dateTable, ":="(date = i.date)]
  setkey(listPosteriors$posteriorTransfGP, date)
  setkey(objectStan$dataForModel, date)
  listPosteriors$posteriorTransfGP[objectStan$dataForModel, ":="(positiveResults = i.positiveResults, numberTest = i.numberTest)]
  
  # Output
  output_main <- list(posteriorGrowth = listPosteriors$posteriorGrowth,
                      posteriorTransfGP = listPosteriors$posteriorTransfGP,
                      posteriorRandomEffect = listPosteriors$posteriorRandomEffect,
                      posteriorIntercept = listPosteriors$posteriorIntercept,
                      posteriorHyperparameters = listPosteriors$posteriorHyperparameters,
                      dateList = objectStan$dateList,
                      dataForModel = objectStan$dataForModel,
                      parametersModel = parametersModel,
                      inferenceSettings = inferenceSettings,
                      objectStan = objectStan)
  output_samples <- list(numSamples = ncol(matrixSampleGP), # TODO check
                         matrixSampleGP = matrixSampleGP,
                         matrixSampleGPDerivative = matrixSampleGPDerivative,
                         matrixSampleRandomEffect = matrixSampleRandomEffect,
                         matrixSampleIntercept = matrixSampleIntercept,
                         matrixSampleHyperparameters = matrixSampleHyperparameters)
  
  # Remove redundant data
  output_main$objectInla$dateList <- NULL
  output_main$objectInla$dataForModel <- NULL
  if(saveStanObject == F) output_main$objectInla <- NULL
  
  if(saveSamples == F){
    output <- output_main
  }else{
    output <- c(output_main, output_samples)
  }
  
  return(output)
}

#' Creates and stores executable in global environment
#' TODO how to do this inside package?
#' TODO do it before running? Annoying error when using Github desktop
constructStanExec <- function(linkType, dirSource){
  # TODO test if exists. Add indicator to variable to show (model type + covariance types).
  # TODO write code dynamically
  
  # Temporal code
  if(1){
    if(!exists("modelExec_NB", envir = .GlobalEnv)){
      cat("Creating temporary executable for the positives model...\n")
      
      # Translate Stan Code to C++ with stanc
      modelStanc <- stanc(file = paste0("/Users/lauraguzmanrincon/Downloads/", "/Stan_fit_temp.stan"))
      
      # Make an Executable Stan Model with stan model
      modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
      #readLines("src/Stan_fit_NB.stan")[readLines("src/Stan_fit_NB.stan") != ""][67]
      assign("modelExec_NB", modelExec, envir = .GlobalEnv)
    }else{
      cat("Reusing executable for the positives model...\n")
    }
  }
  
  
  if(0){
    if(linkType == "BB"){
      if(!exists("modelExec_BB", envir = .GlobalEnv)){
        cat("Creating executable...\n")
        
        # Translate Stan Code to C++ with stanc
        modelStanc <- stanc(file = paste0(dirSource, "/Stan_fit_BB.stan"))
        
        # Make an Executable Stan Model with stan model
        modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
        assign("modelExec_BB", modelExec, envir = .GlobalEnv)
      }
    }else{
      if(!exists("modelExec_NB", envir = .GlobalEnv)){
        cat("Creating executable...\n")
        
        # Translate Stan Code to C++ with stanc
        modelStanc <- stanc(file = paste0(dirSource, "/Stan_fit_NB_period.stan"))
        
        # Make an Executable Stan Model with stan model
        modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
        #readLines("src/Stan_fit_NB.stan")[readLines("src/Stan_fit_NB.stan") != ""][67]
        assign("modelExec_NB", modelExec, envir = .GlobalEnv)
      }
    }
  }
}

