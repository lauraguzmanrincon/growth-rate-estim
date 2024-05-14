loadFile <- function(fileName, fileType){
  # TODO load xlsx files
  if(fileType == "text/csv"){
    options_csv <- list(header = T, sep = ",", quote = "\"", dec = ".")
    dataFile <- data.table(read.csv(file = fileName,
                                    header = options_csv$header,
                                    sep = options_csv$sep,
                                    quote = options_csv$quote,
                                    dec = options_csv$dec))
    return(list(data = dataFile, fileName = fileName, type = "csv", options = options_csv))
  }else{
    stop("Wrong type of file (", fileType, ")")
  }
}

#' Find suitable date format if not chosen, or test format if already chosen
#' Return list of best/chosen format and message list(format, newDates)
findBestDateFormat <- function(dateColumn, dataTable, listFormats, defaultFormat){ # formatColumn,
  tryCatch({
    dateList <- dataTable[order(id), get(dateColumn)]
    output <- list(format = defaultFormat,
                   newDates = NULL)
    
    if(is.numeric(dateList[1])){
      # check if first entry is numeric, assuming first entry is representative
      output <- list(format = listFormats[1],
                     newDates = dateList)
    }else{
      # try date formats. Give the one that gives less NA's
      dateFormats <- listFormats[2:length(listFormats)]
      #chosenFormat <- dateFormats[which.max(sapply(dateFormats, function(x) sum(is.na(as.Date(dateList, format = x))) == 0))] # No NA's
      chosenFormat <- dateFormats[which.max(sapply(dateFormats, function(x) sum(!is.na(as.Date(dateList, format = x)))))] # less NA's
      output <- list(format = chosenFormat,
                     newDates = as.Date(dateList, format = chosenFormat))
    }
    
    return(output)
    # BORRAR
    #dateFormatR <- readDateFormat()
    #dd <- as.Date(dateInput, format = formatInput, tryFormats = listFormats)
    #print(range(dd, na.rm = T))
    #lapply(dateRFormat, function(x) as.Date(rep("23.06.01", 3), format = x))
    #dateRFormat[2:37][which.max(sapply(dateRFormat[2:37], function(x) sum(is.na(as.Date(rep("23.06.01", 3), format = x))) == 0))]
    
    #return(list(format = defaultFormat,
    #            message = ""))
  }, error = function(e) {
    print("error")
    print(e)
    return(list(format = defaultFormat,
                newDates = NULL))
  })
}

#readDateFormat <- function(inputFormat){
#  dateFormatR <- ""
#  if(inputFormat == "YYYY-MM-DD") dateFormatR <- "%Y-%m-%d"
#  return(dateFormatR)
#}

createNewData <- function(){
  data <- list(fileName = "",
               dataTable = data.table(),
               fileType = "",
               columnDate = "",
               columnPositives = "",
               columnTests = "",
               columnCategories = "",
               isRead = F,
               isChecked = F)
}

createNewModel <- function(){
  model <- list()
}

createNewProject <- function(){
  project <- list(code = paste0("proj", sample(x = 1:9000, size = 1) + 999),
                  name = "New project",
                  creationDate = Sys.Date(),
                  date = createNewData(),
                  model = createNewModel())
}