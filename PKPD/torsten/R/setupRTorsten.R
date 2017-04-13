## Install RStan, with Torsten built inside of it.

## Adjust directories to your setting.
scriptDir <- getwd()
lib <- file.path(scriptDir, "lib")
.libPaths(lib)

## Download and edit StanHeaders
install.packages("StanHeaders")
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

install.packages("rstan", type = "source")

setwd(scriptDir)
