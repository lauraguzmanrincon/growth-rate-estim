
# Model fitting ----

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
  if(parametersModel$linkType == "BB"){
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

# Projections ----

#' Plots posterior of GP and model fitting before the last observation date (30 days before as default) and after the size of the projection
#' listProjection: output of createTableProjection...()
#' outputModel:
#' parametersModel:
#' minDate (optional): min date for plot
#' futureTable (optiona): data from the future for comparison (date, positiveResults, numberTest)
#' TODO requires additional package: ggtree
plotProjection <- function(listProjection, outputModel, futureTable = NULL, minDate = NULL, plotGP = F){
  # TODO make better plot
  if(is.null(minDate)) minDate <- listProjection$currentDate - 30
  sizePrediction <- listProjection$sizePrediction
  maxDate <- listProjection$currentDate + sizePrediction
  
  # Create auxiliar tables
  countTableAll <- outputModel$dataForModel[date >= minDate & date <= maxDate, .(date, positiveResults, numberTest)]
  if(!is.null(futureTable)) countTableAll <- rbind(countTableAll, futureTable)
  countTableAll[, ratio := pmin(1, pmax(0, positiveResults/numberTest))]
  
  currentDate <- listProjection$currentDate
  currentDayId <- outputModel$dateList$dateTable[date == currentDate, dayId]
  
  dataToPlotBeforeGP <- outputModel$posteriorTransfGP[date >= minDate & date <= maxDate] #[dayId %in% (currentDayId - (0:30))]
  dataToPlotBeforeGR <- outputModel$posteriorGrowth[date >= minDate & date <= maxDate] #[dayId %in% (currentDayId - (0:30))]
  dataToPlotAfter <- listProjection$tableProjections
  
  # Recover GR on last day of observation if not provided by model (when GP from finite differences)
  if(outputModel$inferenceSettings$derivativeFromGP == F & !is.null(listProjection$tableProjections_boundary)){
    setkey(dataToPlotBeforeGR, dayId)
    setkey(listProjection$tableProjections_boundary, dayId)
    dataToPlotBeforeGR[listProjection$tableProjections_boundary,
                       ":="(median = i.gr_median, q0.025 = i.gr_q025, q0.975 = i.gr_q975, q0.25 = i.gr_q25, q0.75 = i.gr_q75)]
  }
  
  # Compute model posterior
  if(outputModel$parametersModel$linkType == "NB"){
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
               aes(y = positiveResults, colour = weekdays(date) %in% c("Saturday", "Sunday"), shape = weekdays(date) %in% c("Saturday", "Sunday")), show.legend = F, size = 1) +
    scale_colour_manual(values = c("black", "#D41B19")) +
    #labs(title = paste0(partitionToPlot, " - Gaussian process (positive cases/number of tests)"), x = "day", y = "relative effect of GP (positive cases/number of tests)")
    labs(x = "day", y = "posterior fitting", title = listProjection$nameProjection)
  ggP2 <- ggplot(dataToPlotBeforeGR, aes(x = date)) + theme_laura() +
    #geom_vline(xintercept = lockdownDates, linetype = 2, colour = "gray50") +
    geom_hline(yintercept = 0, linetype = 2, colour = "gray50") +
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
      geom_ribbon(aes(ymin = q0.025_GP, ymax = q0.975_GP), fill = "#EFDABD", alpha = 0.5) +
      geom_ribbon(aes(ymin = q0.25_GP, ymax = q0.75_GP), fill = "#D7A35B", alpha = 0.5) +
      geom_line(aes(y = median_GP), colour = "#31220C") +
      geom_ribbon(data = dataToPlotAfter, aes(ymin = gpConsTrans_q025, ymax = gpConsTrans_q975), fill = "#EFDABD", alpha = 0.5) +
      geom_ribbon(data = dataToPlotAfter, aes(ymin = gpConsTrans_q25, ymax = gpConsTrans_q75), fill = "#D7A35B", alpha = 0.5) +
      geom_line(data = dataToPlotAfter, aes(y = gpConsTrans_median), colour = "#31220C")
  }
  
  #return(suppressWarnings(multiplot(ggP1, ggP2, cols = 1)))
  return(list(ggP1, ggP2))
}

plotProjectionComparison <- function(comparisonTable){
  dataToPlot <- comparisonTable
  dataToPlot[, projectionLevel := factor(nameProjection)] # , levels = nameMethods, labels = nameMethods
  dataToPlot[, scoreLevel := factor(scoreName)] # levels = c("IS score", "DSS score", "RE score"), labels = c("IS score", "DSS score", "RE score")
  #dataToPlot[date %in% dataToPlot[scoreName == "DSS score" & score > 100, date], score := NA] # manual adjustment ro remove too high scores
  ggSc <- ggplot(dataToPlot, aes(x = date, y = score)) + theme_laura() + facet_grid(scoreLevel ~ ., scales = "free_y") +
    geom_line(aes(colour = projectionLevel), alpha = 0.3) + geom_point(aes(colour = projectionLevel)) +
    #scale_x_date(limits = range(outputModel$posteriorTransfGP$date), expand = c(0,0)) +
    #scale_x_date(labels = scales::label_date(c("%d %b\n(%Y)", rep("%d %b", length(dateBreaks) - 1))),
    #             breaks = dateBreaks, limits = c(as.Date(c("2020-09-01", "2021-07-31")))) +
    scale_colour_manual(values = c("#1D1E18", "#E3B505", "gray")) +
    #theme(legend.position = c(0.7, 0.8), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
    #coord_cartesian(ylim = c(0, 2000)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.grid.major.x = element_line(linetype = 2, colour = "gray90")) +
    labs(x = "current date", y = "score", colour = "method")
  return(ggSc)
}

