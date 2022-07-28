/**
 * @file omicsds_query.cpp
 *
 * @section LICENSE
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2022 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * Rcpp code to interface with the omicds library
 *
 **/

#include <RcppCommon.h>
#include "omicsds_types.h"

#include <Rcpp.h>
#include <exception>
#include <iostream>

// [[Rcpp::export]]
void version() {
  Rcpp::Rcout << "OmicsDS Version=" << OmicsDS::version()<< std::endl;
}

// [[Rcpp::export]]
size_t connect(const std::string& workspace, const std::string& array) {
  Rcpp::Rcout << "Got OmicsDS" << std::endl;
  return OmicsDS::connect(workspace, array);;
}

// [[Rcpp::export]]
void disconnect(size_t handle) {
  OmicsDS::disconnect(handle);
}

// [[Rcpp::export]]
void query(size_t handle) {
  std::array<int64_t, 2> sample_range = {0, std::numeric_limits<int64_t>::max()};
  std::array<int64_t, 2> position_range = {0, std::numeric_limits<int64_t>::max()};
  OmicsDS::query(handle, sample_range, position_range);
}
