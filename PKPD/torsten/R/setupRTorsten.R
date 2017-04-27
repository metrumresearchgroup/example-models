## Install RStan, with Torsten built inside of it.
## This file won't installs the RStan's dependencies.

## Adjust directories to your setting.
scriptDir <- getwd()
lib <- file.path(scriptDir, "lib")
.libPaths(lib)

# library(devtools)

## Download and edit StanHeaders (version 2.14)
install.packages("https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.14.0-1.tar.gz", 
                 repos = NULL)
setwd(lib)
system("git clone https://github.com/charlesm93/stan.git")
setwd("stan")
system("git checkout torsten-master")
setwd(lib)
system("rm -rf StanHeaders/include/src/stan")
system("mv stan/src/stan StanHeaders/include/src/stan")
system("rm -rf stan")
system("git clone https://github.com/charlesm93/math.git")
setwd("math")
system("git checkout torsten-master")
setwd(lib)
system("rm -rf StanHeaders/include/stan")
system("mv math/stan StanHeaders/include/stan")
system("rm -rf math")

## Download rstan 2.14.1 without the dependencies.
install.packages("https://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.14.2.tar.gz", repos = NULL)

setwd(scriptDir)
