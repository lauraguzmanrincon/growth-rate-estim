
# Functions for running model with STAN (MCMC) ----
#' countTable: data table with one row per day (of days with available data) and these columns:
#' - numberTest: number of test on the day (>= 0)
#' - positiveResults: number of positive tests (<= numberTest and >= 0)
#' - date: date of count in R date format
#' parametersModel: output of setParametersFn()
#' saveSamples: Default F. If T, returns matrixSampleDays and sampleDerivatives,
#'              two matrices of size [days, num. samples] containing samples of the posterior of the GP and GP derivative respectively.
#' minDate (optional): minimum date to include in the model
#' maxDate (optional): maximum date to include in the model
#' 
#' Output: dataForModel: this data includes all days in range. Recall STAN does not recieve the ones without positives and tests
#' 
#' Example:
#' parametersStan <- list(sampleFile = "R_storage/R47_Output/Test", chains = 1, iter = 10, warmup = 2, thin = 1)
#' modelFit <- runModelGrowthRate_STAN(countTable, parametersRun, minDate = NULL, maxDate = NULL, parametersStan = parametersStan)
#' outputStan <- processSTANOutput(modelFit, parametersRun, saveSamples = F)
#' runModelGrowthRate_STAN saves this in [parametersStan$sampleFile].RData: modelFit, modelStanc, modelData, parametersStan
#' TODO minDate, maxDate lost. FIX dateTable!
#' TODO using minDayInla and maxDayInla from outside!
runModelGrowthRate_STAN <- function(countTable, parametersModel, minDate = NULL, maxDate = NULL, seed = sample(1000, 1),
                               parametersStan = list(sampleFile = "MCMC_samples", chains = 1, iter = 10000, warmup = 5000, thin = 1)){
  # Load parameters into function environment and add auxiliar variables
  #list2env(parametersModel$params, envir = environment())
  #list2env(parametersModel$config, envir = environment())
  
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
  
  if(parametersModel$param$linkType == "BB"){
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
  dirSource <- "/Users/lauraguzmanrincon/Documents/GitHub/growth-rate-estim/src_gr" # TODO how to do in package?
  constructStanExec(linkType = parametersModel$param$linkType, dirSource = dirSource)
  if(parametersModel$param$linkType == "BB"){
    modelExec <- modelExec_BB
  }
  else{
    modelExec <- modelExec_NB
  }
  # TODO today... write other file
  
  # Data to Stan
  # We assume dayId runs from 1 in steps of unitTime
  # TODO works only for day of the week effect
  modelData <- list(
    # dimensions
    num_days = nrow(dataForModelNoNa), #as.integer(numDays),
    num_groups = as.integer(7),
    # data
    t = dataForModelNoNa[order(dayId), dayId], #1:numDays,
    pos = dataForModelNoNa[order(dayId), positiveResults],
    tot = dataForModelNoNa[order(dayId), as.integer(numberTest)],
    day_to_group = as.integer(dataForModelNoNa[order(dayId), factor(dayWeek, levels = parametersModel$internal$levelsWeek)]),
    # prior intercept
    m_int = parametersModel$params$interceptPrior$mean,
    sig_int = 1/sqrt(parametersModel$params$interceptPrior$prec),
    # prior dispersion
    m_eta = parametersModel$params$dispersionPrior$mean, # positives
    sig_eta = 1/sqrt(parametersModel$params$dispersionPrior$prec), # positives
    m_rho = parametersModel$params$dispersionPrior$mean, # proportions
    sig_rho = 1/sqrt(parametersModel$params$dispersionPrior$prec), # proportions
    # prior day-of-the-week effect
    a_w = parametersModel$params$randomEffectPrior$a,
    b_w = parametersModel$params$randomEffectPrior$b,
    # prior GP
    log_theta_0 = 0 + log(parametersModel$params$GPHyperparamPrior$mean),
    #OLDparametersModel$params$theta.prior2.mean[2:1] + log(c(parametersModel$params$prior2.range0, parametersModel$params$prior2.sigma0)),
    B = parametersModel$params$GPHyperparamPrior$B)
  
  # Initial values
  set.seed(seed)
  if(parametersModel$param$linkType == "NB"){
    x_t_notTrasnform <- log(modelData$pos)
    modelInits <- replicate(parametersStan$chains, list(intercept = mean(log(modelData$pos)),
                                                        eta = exp(rnorm(n = 1, mean = modelData$m_eta, sd = modelData$sig_eta)),
                                                        x_t = x_t_notTrasnform - mean(x_t_notTrasnform),
                                                        #TODO
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(parametersModel$params$GPHyperparamPrior$mean) + 0,
                                                                                       sigma = parametersModel$params$GPHyperparamPrior$B)[1,],
                                                        w_d = rep(0, length(parametersModel$internal$levelsWeek)),
                                                        tau_w = rgamma(1, shape = modelData$a_w, rate = modelData$b_w)), simplify = F)
  }else if(parametersModel$param$linkType == "BB"){
    logit <- function(p) log(p/(1 - p))
    inv.logit <- function(x) exp(x)/(1 + exp(x))
    x_t_notTrasnform <- logit(pmin(pmax(modelData$pos, 1), modelData$tot - 1)/modelData$tot)
    modelInits <- replicate(parametersStan$chains, list(intercept = rnorm(n = 1, mean = modelData$m_int, sd = modelData$sig_int),
                                                        rho = inv.logit(rnorm(n = 1, mean = modelData$m_rho, sd = modelData$sig_rho)), #modelData$m_rho OR -2?
                                                        x_t = x_t_notTrasnform - mean(x_t_notTrasnform),
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(parametersModel$params$GPHyperparamPrior$mean) + 0,
                                                                                       sigma = parametersModel$params$GPHyperparamPrior$B)[1,],
                                                        w_d = rep(0, length(levelsWeek)),
                                                        tau_w = rgamma(1, shape = modelData$a_w, rate = modelData$b_w)), simplify = F)
  }
  
  # Run model
  cat("Running model\n")
  modelFit <- sampling(modelExec,
                       data = modelData,
                       init = modelInits,
                       chains = parametersStan$chains,
                       iter = parametersStan$iter,
                       warmup = parametersStan$warmup,
                       thin = parametersStan$thin,
                       seed = c(12321),
                       sample_file = parametersStan$sampleFile)
  # TODO what does STAN do with NA values?
  
  # Create output
  modelStanc <- modelExec@model_code
  objectStan <- list(modelFit = modelFit,
                     modelStanc = modelStanc,
                     dataForModel = dataForModel,
                     dateList = list(dateTable = dateTable, numDays = numDays, minDay = minDay, maxDay = maxDay))
  
  cat("Saving modelFit, modelStanc, modelData, parametersStan in ", parametersStan$sampleFile, ".RData\n", sep = "")
  #save(modelFit, modelStanc, dataForModel, parametersStan, file = paste0(parametersStan$sampleFile, ".RData")) # 30.11.2023
  
  return(objectStan)
}

processSTANOutput <- function(objectStan, parametersModel, saveSamples = F){
  # ---------------------------------------------------- #
  #                      GET SAMPLES                     #
  # ---------------------------------------------------- #
  
  # Samples GP
  outputFit <- data.table(summary(objectStan$modelFit)$summary, keep.rownames = T)
  samplesFit <- extract(objectStan$modelFit)
  matrixSampleDays <- t(samplesFit$x_t)
  
  # Samples derivative
  # TODO Condition on parametersModel$config$derivativeFromGP !! can we use getGrowthFromSamples_GP here?
  sampleDerivatives <- getGrowthFromSamples(matrixSampleDays)
  
  # Samples parameters
  # TODO
  #matrixSampleRandomEffect matrix [ x samples]
  #matrixSampleHyperAll matrix [c("theta1", "theta2", "overdispersion", "precision") x samples]
  #matrixSampleIntercept vector [samples]
  #matrixSampleNu
  #vectorTestData
  
  # ---------------------------------------------------- #
  #                  PRODUCE OUTPUT                      #
  # ---------------------------------------------------- #
  listPosteriors <- computePosteriors(matrixSampleDays, sampleDerivatives, matrixSampleOverdisp,
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
  output_main <- list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP,
                      posteriorRandomEffect = listPosteriors$posteriorRandomEffect,
                      dateList = objectStan$dateList, dataForModel = objectStan$dataForModel)
  output_samples <- list(matrixSampleDays = matrixSampleDays, sampleDerivatives = sampleDerivatives,
                         matrixSampleRandomEffect = matrixSampleRandomEffect, matrixSampleHyperAll = matrixSampleHyperAll, matrixSampleIntercept = matrixSampleIntercept)
  
  output <- output_main
  if(saveSamples == T) output <- c(output, output_samples, list(outputType = "MCMC_Stan"))
  return(output)
}

#' Creates and stores executable in global environment
#' TODO how to do in package?
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
        cat("Creating executable for the proportions model...\n")
        
        # Translate Stan Code to C++ with stanc
        modelStanc <- stanc(file = paste0(dirSource, "/Stan_fit_BB.stan"))
        
        # Make an Executable Stan Model with stan model
        modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
        assign("modelExec_BB", modelExec, envir = .GlobalEnv)
      }
    }else{
      if(!exists("modelExec_NB", envir = .GlobalEnv)){
        cat("Creating executable for the positives model...\n")
        
        # Translate Stan Code to C++ with stanc
        modelStanc <- stanc(file = paste0(dirSource, "/Stan_fit_NB.stan"))
        
        # Make an Executable Stan Model with stan model
        modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
        #readLines("src/Stan_fit_NB.stan")[readLines("src/Stan_fit_NB.stan") != ""][67]
        assign("modelExec_NB", modelExec, envir = .GlobalEnv)
      }
    }
  }
}

