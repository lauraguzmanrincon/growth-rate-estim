
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
#' NOdayWeek: if T, the day-of-the-week effect is used. If F, the noise term is used
#' NOdayWeek.prior.prec: prior for day-of-the-week effect or noise term
#' ... TODO
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
                            hasConstant = F, # 21.09.2023
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
      hasConstant = hasConstant,
      dayWeek.prior.prec = dayWeek.prior.prec
    ),
    config = list(
      sizeSample = sizeSample,
      derivativeFromGP = derivativeFromGP,
      computeGPProjection = computeGPProjection,
      sizeGPProjection = sizeGPProjection
      #computeTaylorProjection = F,
      #computeLinearProjection = T
      ),
    internal = list(
      levelsWeek = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
    ))
    return(parameters)
}

# TODO create function to check countTable before using it. It should have numberTest if BB, or create an empty one if non existant and NB.
# It should have positives <= number tests










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
  posteriorGrowth[posteriorTransfGP, ":="(median_GP = i.median_GP, q0.025_GP = i.q0.025_GP, q0.975_GP = i.q0.975_GP,
                                          q0.25_GP = i.q0.25_GP, q0.75_GP = i.q0.75_GP,
                                          #median_transGP = i.median_transGP, q0.025_transGP = i.q0.025_transGP, q0.975_transGP = i.q0.975_transGP,
                                          #q0.25_transGP = i.q0.25_transGP, q0.75_transGP = i.q0.75_transGP,
                                          median_transConsGP = i.median_transConsGP, q0.025_transConsGP = i.q0.025_transConsGP, q0.975_transConsGP = i.q0.975_transConsGP,
                                          q0.25_transConsGP = i.q0.25_transConsGP, q0.75_transConsGP = i.q0.75_transConsGP,
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







# Auxiliar post-run ----
#' Function that updates outputs from before 19 Sep 2023 to the 19.09.2023 version
#' To obtain the 19.09.2023 version use [ git checkout 'master@{2023-09-19}'" ]
updateOutputVersion18Sep2023 <- function(outputModel, parametersModel){
  parametersModel$params$randomEffect <- ifelse(parametersModel$params$dayWeek == "day", "weekday", "all")
  
  rownames(outputModel$matrixSampleHyperAll) <- c("overdispersion", "theta1", "theta2", "precision")
  outputModel$posteriorTransfGP[, ":="(median_transConsGP = median, q0.025_transConsGP = q0.025, q0.975_transConsGP = q0.975,
                                       q0.25_transConsGP = q0.25, q0.75_transConsGP = q0.75)]
  outputModel$posteriorTransfGP[order(dayId), ":="(median_GP = apply(outputModel$matrixSampleDays, 1, quantile, probs = 0.5, na.rm = T),
                                                   q0.025_GP = apply(outputModel$matrixSampleDays, 1, quantile, probs = 0.025, na.rm = T),
                                                   q0.975_GP = apply(outputModel$matrixSampleDays, 1, quantile, probs = 0.975, na.rm = T),
                                                   q0.25_GP = apply(outputModel$matrixSampleDays, 1, quantile, probs = 0.250, na.rm = T),
                                                   q0.75_GP = apply(outputModel$matrixSampleDays, 1, quantile, probs = 0.750, na.rm = T))]
  
  outputModel$matrixSampleRandomEffect <- outputModel$matrixSampleWeekday
  outputModel$posteriorRandomEffect <- data.table(index = 1:nrow(outputModel$matrixSampleRandomEffect),
                                                  median = apply(outputModel$matrixSampleRandomEffect, 1, quantile, probs = 0.5, na.rm = T),
                                                  q0.025 = apply(outputModel$matrixSampleRandomEffect, 1, quantile, probs = 0.025, na.rm = T),
                                                  q0.975 = apply(outputModel$matrixSampleRandomEffect, 1, quantile, probs = 0.975, na.rm = T),
                                                  q0.25 = apply(outputModel$matrixSampleRandomEffect, 1, quantile, probs = 0.250, na.rm = T),
                                                  q0.75 = apply(outputModel$matrixSampleRandomEffect, 1, quantile, probs = 0.750, na.rm = T))
  
  return(list(outputModel = outputModel, parametersModel = parametersModel))
}

# Figures ----


