
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
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA,
                                              median = median_transConsGP, q0.025 = q0.025_transConsGP, q0.975 = q0.975_transConsGP, type = "Gaussian process")])
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
                          mean = parametersModel$params$BB.prior.rho$overdispersion$param[1],
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
  samples <- outputModel$matrixSampleDays + # TODO ok? # 22.09.2023
    matrix(outputModel$matrixSampleIntercept, nrow = outputModel$dateList$numDays, ncol = parametersModel$config$sizeSample, byrow = T)
  rownames(samples) <- NULL
  dataSamplesGPCons <- data.table(suppressWarnings(melt(samples, varnames = c("dayId", "sample"))))
  
  setkey(dataSamplesGPCons, dayId)
  setkey(outputModel$dataForModel, dayId)
  dataSamplesGPCons[outputModel$dataForModel, date := i.date]
  
  if(parametersModel$params$linkType == "BB"){
    dataSamplesGPCons[, valueTr := exp(value)/(1 + exp(value))]
  }else if(parametersModel$params$linkType == "NB"){
    dataSamplesGPCons[, valueTr := exp(value)]
  }
  randomCols <- sample(ncol(outputModel$matrixSampleDays)) # for colour of samples
  dataSamplesGPCons[, colourId := randomCols[sample]]
  
  dataToPlotPosterior <- outputModel$posteriorTransfGP
  dataToPlotPosterior[, ratio := pmin(1, pmax(0, positiveResults/numberTest))]
  dataToPlotPosterior[, isWeekend := weekdays(date) %in% c("Saturday", "Sunday")]
  dataToPlotPosterior[, weekendLabel := factor(ifelse(isWeekend, "weekend", "weekday"), levels = c("weekend", "weekday"))]
  dataToPlot <- rbind(dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio, positiveResults, median = NA, q0.025 = NA, q0.975 = NA, type = "points")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA, median = medianFT, q0.025 = q0.025FT, q0.975 = q0.975FT, type = "model fit")],
                      dataToPlotPosterior[, .(dayId, date, weekendLabel, ratio = NA, positiveResults = NA,
                                              median = median_transConsGP, q0.025 = q0.025_transConsGP, q0.975 = q0.975_transConsGP, type = "Gaussian process")])
  if(parametersModel$params$linkType == "NB"){
    dataToPlot[, points := positiveResults]
  }else{
    dataToPlot[, points := ratio]
  }
  dataToPlot[, typeLevel := factor(type, levels = c("model fit", "Gaussian process"))]
  p0 <- ggplot(dataToPlot, aes(x = date)) + theme_laura() +
    geom_line(data = dataSamplesGPCons, aes(y = valueTr, group = sample, colour = colourId), alpha = 0.5) + # colour = "gray50"
    #geom_ribbon(data = dataToPlot[typeLevel == "Gaussian process"], aes(ymin = q0.025, ymax = q0.975, fill = typeLevel), alpha = 0.1) +
    geom_line(data = dataToPlot[typeLevel == "Gaussian process"], aes(y = median), linetype = 2, colour = "black") + # #F5CC14
    geom_point(aes(y = points), colour = "red", size = 0.5) +
    scale_colour_gradient(guide = "none", low = "gray90", high = "gray10") +
    #scale_colour_manual(name = "posterior (median, 95% CI)", values = "#F5CC14") +
    #scale_fill_manual(name = "posterior (median, 95% CI)", values = "#F5CC14") +
    #scale_x_date(labels = scales::label_date(formatBreaks), breaks = dateBreaks2) +
    scale_x_date(expand = c(0,0)) +
    scale_y_continuous(labels = unit_format(unit = "", scale = ifelse(parametersModel$params$linkType == "NB", 1, 1e-3))) + #unit = "K"
    labs(x = "day", y = ifelse(parametersModel$params$linkType == "NB",
                               "model posterior samples\nand proportion of positives",
                               "model posterior samples\nand count of positives (thousands)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(),
          legend.key = element_blank(), panel.grid.major.x = element_line(linetype = 2, colour = "gray90")) # 0.6
  
  return(p0)
}


plotFittingSamplesGR <- function(outputModel, parametersModel){
  # TODO !!!!! check
  
  # Transform samples from GP derivative to GR
  samplesGR <- getSamplesGR(outputModel$matrixSampleDays, outputModel$sampleDerivatives, parametersModel)
  
  # Create data table samples
  dataSamplesGR <- data.table(suppressWarnings(melt(samplesGR, varnames = c("dayId", "sample"))))
  setkey(dataSamplesGR, dayId)
  setkey(outputModel$dataForModel, dayId)
  dataSamplesGR[outputModel$dataForModel, date := i.date]
  
  # Create plot
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
    labs(x = "day", y = "growth rate posterior samples") +
    scale_x_date(expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), #legend.position = c(0.62, 0.75),
          legend.key = element_blank(), panel.grid.major.x = element_line(linetype = 2, colour = "gray90"))
  return(p0)
}

#' .
plotLatent <- function(outputModel){
  # TODO plot GP or plot latent (GP + constant)?
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
    p0 <- p0 + scale_x_continuous(breaks = 1:7, labels = parametersModel$internal$levelsWeek)
  }
  return(p0)
}

plotIntercept <- function(outputModel, parametersModel){
  dataToPlot <- data.table(samples = outputModel$matrixSampleIntercept)
  p0 <- ggplot(dataToPlot, aes(x = samples)) + theme_laura() +
    geom_density(alpha = 0.6, fill = "#87A696") +
    #geom_point(aes(x = samples, y = 0, colour = type)) +
    labs(x = "intercept", y = "density")
  return(p0)
}

#' Copied from C01Fun.R
plotHyperparametersGPStackFn <- function(outputModelList, parametersModel, listLabels, listColours){
  set.seed(159)
  samplesPrior <- mvtnorm::rmvnorm(100,
                                   mean = parametersRun$params$theta.prior2.mean,
                                   sigma = solve(parametersModel$params$theta.prior2.prec))
  samplesPosteriorList <- vector("list", length(outputModelList))
  for(iip in 1:length(outputModelList)){
    samplesPosteriorList[[iip]] <- data.table(sigma = parametersModel$params$prior2.sigma0*exp(outputModelList[[iip]]$matrixSampleHyperAll["theta1",]),
                                              range = parametersModel$params$prior2.range0*exp(outputModelList[[iip]]$matrixSampleHyperAll["theta2",]),
                                              type = paste0("posterior ", listLabels[iip]))[,.N,.(sigma, range, type)]
  }
  plotHyperparameters <- rbind(do.call("rbind", samplesPosteriorList),
                               data.table(sigma = parametersModel$params$prior2.sigma0*exp(samplesPrior[,1]),
                                          range = parametersModel$params$prior2.range0*exp(samplesPrior[,2]),
                                          N = 1,
                                          type = "prior"))
  plotHyperparameters[, typeLevel := factor(type, levels = c("prior", paste0("posterior ", listLabels)))]
  reviewHyperparameters <- plotHyperparameters[, .(medianSigma = median(sigma), medianRange = median(range)), .(type, typeLevel)]
  gg <- ggplot(plotHyperparameters, aes(x = sigma, y = range/2, colour = typeLevel, shape = typeLevel)) + theme_laura() +
    geom_point(data = reviewHyperparameters, aes(x = medianSigma, y = medianRange/2)) + #geom_point(aes(size = N))
    scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') + stat_ellipse(type = "norm", linewidth = 0.5) +
    scale_colour_manual(name = "", values = c("black", listColours)) +
    scale_shape_manual(name = "", values = rep(c(16, 17), length(outputModelList))) +
    theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = c(0.82, 0.07), legend.background = element_rect(fill = NA),
          legend.key.height = unit(0.1, "cm")) +
    labs(x = expression(standard~deviation~sigma), y = expression(length~scale~plain(l)~(days)))
  return(gg)
}

