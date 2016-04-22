#! /usr/bin/env Rscript

options(repos=structure(c(CRAN="http://cran-mirror.cs.uu.nl/")))
packages <- c("Rcpp", "lpSolve", "mvtnorm", "stringr", "matrixcalc", "Matrix");
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])];
if(length(new.packages))
  install.packages(new.packages);
packages_github <- c("Karel-Kroeze/MultiGHQuad");
new.packages_github <- packages_github[!(packages_github %in% installed.packages()[,"Package"])];
if(length(new.packages_github)) {
  library(devtools)
  install_github(new.packages_github); 
}
update.packages(lib.loc=Sys.getenv("R_LIBS_USER"), ask=FALSE);
