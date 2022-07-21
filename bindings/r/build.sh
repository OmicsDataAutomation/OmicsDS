#!/bin/bash

./cleanup

R -e 'library(Rcpp); compileAttributes(".")'
R CMD build . &&
R CMD INSTALL --preclean omicsds_0.0.1.tar.gz &&
R CMD check --no-manual omicsds_0.0.1.tar.gz # &&
# R -e 'library(testthat); test()'
