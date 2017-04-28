## Analyze result from performance test

## Names of two performance tests we want to compare
testName <- "fTwoCpt"
testName2 <- "fTwoCpt_mixed"

## Relative paths
scriptDir <- getwd()
resultDir <- file.path(scriptDir, "test", "deliv")
figDir <- file.path(scriptDir, "test", "deliv", testName)
compFigDir <- file.path(scriptDir, "test", "deliv", paste0(testName, "To", testName2))

dir.create(resultDir)

library(ggplot2)
source("test/testPlotTools.R")

###############################################################################

parameters <- c("CL", "Q", "VC", "VP", "ka", "sigma", "mtt", "circ0", "alpha",
                "gamma", "sigmaNeut")

###############################################################################
## Generate plots for only one model
BoxData1 <- Moustache(testName = testName, 
                      parameterNames = parameters,
                      resultDir = resultDir)

## Can exclude some parameters (for readibility purposes)
# BoxData1 <- Moustache(testName, parameters, resultDir, exclude = c("alpha"))

dir.create(figDir)
## opem graphics device
pdf(file = file.path(figDir, paste(testName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

## BoxPlots
ggplot(BoxData1, aes(parameter, FracDiff)) + geom_boxplot() + ggtitle(testName)
ggplot(BoxData1, aes(parameter, abs(FracDiff))) + geom_boxplot() + ggtitle(testName)
ggplot(BoxData1, aes(parameter, neff)) + geom_boxplot() + ggtitle(testName)
ggplot(BoxData1, aes(parameter, time)) + geom_boxplot() + ggtitle(testName)

dev.off()

###############################################################################
## Comparison with model 2

BoxData2 <- Moustache(testName2, parameters, resultDir)
# BoxData2 <- Moustache(testName2, parameters, resultDir, exclude = c("alpha"))

BoxDataComp <- rbind(BoxData1, BoxData2)

dir.create(compFigDir)
pdf(file = file.path(compFigDir, paste(testName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

compPlotDiff <- ggplot(BoxDataComp, aes(parameter, FracDiff, color = Model)) + geom_boxplot()
compPlotDiff + ggtitle(paste("Comparison between", testName, "and", testName2))

compPlotNeff <- ggplot(BoxDataComp, aes(parameter, neff, color=Model)) + geom_boxplot() 
compPlotNeff + ggtitle(paste("Comparison between", testName, "and", testName2))

compPlotTime <- ggplot(BoxDataComp, aes(parameter, time, color=Model)) + geom_boxplot() 
compPlotTime + ggtitle(paste("Comparison between", testName, "and", testName2))

dev.off()


## Additional PLots of interest
# Compare log times to reduce effects of outliers
compPlotTime <- ggplot(BoxDataComp, aes(parameter, log(time), color=Model)) + geom_boxplot() 
compPlotTime + ggtitle(paste("Comparison between", testName, "and", testName2))


