## 06/08/2016: v1.0

## Functions to plot results of tests for Stan Pmetrics

Moustache <- function(testName, parameterNames, resultDir, exclude = NULL){
  result <- read.csv(file.path(resultDir,paste(testName,".summary.csv",sep="")))
  result$X <- NULL ## elimiate useless column X 
  N <- nrow(result) ## number of simulations (expected: 100)
  nParameters <- ncol(result)/3  ## number of parameters in tested model
  
  NamesCol <- c(paste0(parameterNames,"_diff (%)"), 
                paste0(parameterNames,"_neff"),
                paste0(parameterNames,"_time"))
  
  colnames(result) <- NamesCol
  
  ## Find indexes of parameters we want to plot
  excludeIndexes <- c()
  for(i in 1:length(exclude)){
    excludeIndexes <- append(excludeIndexes, which(parameterNames ==  exclude[i]))
  } 
  Indexes <- 1:nParameters
  Indexes <- Indexes[! Indexes %in% excludeIndexes]
  plotParameters <- parameterNames[Indexes]
  nPlotParameters <- length(plotParameters)
  
  ## difference (non-absolute) (fractional and in percentage)
  diffData <- c()
  for(i in Indexes) diffData <- append(diffData, result[[i]])

  BoxData <- data.frame(diffData,
                        rep(plotParameters, each=N),
                        rep(testName, N*nPlotParameters))
  colnames(BoxData)[1] <- "FracDiff"
  colnames(BoxData)[2] <- "parameter"
  colnames(BoxData)[3] <- "Model"

  ## effective number of samples
  neffData <- c()
  for(i in (Indexes + nParameters)) neffData <- append(neffData, result[[i]])
  BoxData$neff <- neffData
  
  ## computation time (in units of time per iteration * 1000)
  timeData <- c()
  for(i in (Indexes + 2*nParameters)) timeData <- append(timeData, result[[i]])
  BoxData$time <- timeData
  
  ## return data frame for box plots 
  BoxData
}

## function for estimating the mode
## (source: http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode)
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

