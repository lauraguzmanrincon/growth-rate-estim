
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
runModelGrowthRate_STAN <- function(countTable, parametersModel, minDate = NULL, maxDate = NULL,
                               parametersStan = list(sampleFile = "MCMC_samples", chains = 1, iter = 10000, warmup = 5000, thin = 1)){
  # Load parameters into function environment and add auxiliar variables
  list2env(parametersModel$params, envir = environment())
  list2env(parametersModel$config, envir = environment())
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
  
  # Data to Stan
  # We assume dayId runs from 1 in steps of unitTime
  modelData <- list(
    # dimensions
    num_days = nrow(dataForModelNoNa), #as.integer(numDays),
    num_groups = as.integer(7),
    # data
    t = dataForModelNoNa[order(dayId), dayId], #1:numDays,
    pos = dataForModelNoNa[order(dayId), positiveResults],
    tot = dataForModelNoNa[order(dayId), as.integer(numberTest)],
    day_to_group = as.integer(dataForModelNoNa[order(dayId), factor(dayWeek, levels = levelsWeek)]),
    # prior dispersion
    m_eta = parametersModel$params$NB.prior.rho$size$param[1], # NB
    sig_eta = 1/sqrt(parametersModel$params$NB.prior.rho$size$param[2]), # NB
    m_rho = parametersModel$params$BB.prior.rho$overdispersion$param[1], # BB
    sig_rho = 1/sqrt(parametersModel$params$BB.prior.rho$overdispersion$param[2]), # BB
    # prior day-of-the-week effect
    a_w = parametersModel$params$dayWeek.prior.prec$theta$param[1],
    b_w = parametersModel$params$dayWeek.prior.prec$theta$param[2],
    # prior GP
    log_theta_0 = parametersModel$params$theta.prior2.mean[2:1] + log(c(parametersModel$params$prior2.range0, parametersModel$params$prior2.sigma0)),
    B = solve(parametersModel$params$theta.prior2.prec))
  
  # Initial values
  set.seed(9812)
  if(parametersModel$param$linkType == "NB"){
    modelInits <- replicate(parametersStan$chains, list(eta = exp(rnorm(n = 1, mean = modelData$m_eta, sd = modelData$sig_eta)),
                                                        x_t = log(modelData$pos),
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(c(parametersModel$params$prior2.range0, parametersModel$params$prior2.sigma0)) +
                                                                                         parametersModel$params$theta.prior2.mean,
                                                                                       sigma = solve(parametersModel$params$theta.prior2.prec))[1,],
                                                        w_d = rep(0, length(parametersModel$internal$levelsWeek)),
                                                        tau_w = rgamma(1, shape = modelData$a_w, rate = modelData$b_w)), simplify = F)
  }else if(parametersModel$param$linkType == "BB"){
    logit <- function(p) log(p/(1 - p))
    inv.logit <- function(x) exp(x)/(1 + exp(x))
    modelInits <- replicate(parametersStan$chains, list(rho = inv.logit(rnorm(n = 1, mean = -2, sd = modelData$sig_rho)), #modelData$m_rho
                                                        x_t = logit(pmin(pmax(modelData$pos, 1), modelData$tot - 1)/modelData$tot),
                                                        log_theta_x = mvtnorm::rmvnorm(n = 1,
                                                                                       mean = log(c(parametersModel$params$prior2.range0, parametersModel$params$prior2.sigma0)) +
                                                                                         parametersModel$params$theta.prior2.mean,
                                                                                       sigma = solve(parametersModel$params$theta.prior2.prec))[1,],
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
                     dataForModel = dataForModel)
  
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
  listPosteriors <- computePosteriors(matrixSampleDays, sampleDerivatives, matrixSampleHyperAll,
                                      matrixSampleNu, matrixSampleRandomEffect, matrixSampleIntercept,
                                      vectorTestData, parametersModel)
  
  setkey(listPosteriors$posteriorGrowth, dayId)
  setkey(dateTable, dayId)
  listPosteriors$posteriorGrowth[dateTable, ":="(date = i.date)]
  
  setkey(listPosteriors$posteriorTransfGP, dayId)
  setkey(dateTable, dayId)
  listPosteriors$posteriorTransfGP[dateTable, ":="(date = i.date)]
  setkey(listPosteriors$posteriorTransfGP, date)
  setkey(countTable, date)
  listPosteriors$posteriorTransfGP[countTable, ":="(positiveResults = i.positiveResults)]
  
  if(saveSamples == F){
    return(list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP))
  }else{
    return(list(posteriorGrowth = listPosteriors$posteriorGrowth, posteriorTransfGP = listPosteriors$posteriorTransfGP,
                matrixSampleDays = matrixSampleDays, sampleDerivatives = sampleDerivatives,
                matrixSampleRandomEffect = matrixSampleRandomEffect, matrixSampleHyperAll = matrixSampleHyperAll, matrixSampleIntercept = matrixSampleIntercept))
  }
}

#' Creates and stores executable in global environment
#' TODO how to do in package?
#' TODO do it before running? Annoying error when using Github desktop
constructStanExec <- function(linkType, dirSource){
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
      modelStanc <- stanc(file = paste0(dirSource, "/Stan_fit_NB.stan"))
      
      # Make an Executable Stan Model with stan model
      modelExec <- stan_model(stanc_ret = modelStanc) # ~30s
      #readLines("src/Stan_fit_NB.stan")[readLines("src/Stan_fit_NB.stan") != ""][67]
      assign("modelExec_NB", modelExec, envir = .GlobalEnv)
    }
  }
}

