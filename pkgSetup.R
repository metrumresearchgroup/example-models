pkgs <- c('ggplot2', 'bayesplot', 'devtools',
          'plyr', 'dplyr', 'tidyr',
          'KernSmooth', 'xml2', 'rvest')

parentScriptDir <- getwd()  ## You may need to mod this to be in the top level of scriptDir
pkgDir <- file.path(parentScriptDir, "pkg")
dir.create(pkgDir)
libDir <- file.path(parentScriptDir, "lib")
dir.create(libDir)

.libPaths(libDir)

install.packages(pkgs,
                 lib=libDir,
                 contriburl=c(contrib.url("http://r-forge.r-project.org","source"),
                              contrib.url("https://cran.rstudio.com/","source")),
                 destdir=pkgDir,
                 type="source",
                 dependencies = c("Depends", "Imports", "LinkingTo"),
                 INSTALL_opts="--no-multiarch")
 
library(devtools)

devtools::load_all('../TorstenHeaders')
##install_torsten(StanHeaders_version = "2.18.0",
##                rstan_version = "2.18.2")
install_torsten(rstan_version = "2.18.2")
 
