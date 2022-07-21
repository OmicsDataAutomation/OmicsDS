/**
 * src/main/cpp/api/omicsds.cc
 *
 * The MIT License (MIT)
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
 * OmicsDS API Implementation
 */


#include "omicsds.h"
#include "omicsds_export.h"

#include <mutex>
#include <map>



std::string OmicsDS::version() {
  return "0.0.1";
}

std::map<OmicsDSHandle, std::shared_ptr<OmicsExporter>> omicsds_instances {};
std::mutex omicsds_mutex;  // protects omicsds_instances
OmicsDSHandle current_handle = 0ul;

OmicsDSHandle OmicsDS::connect(const std::string& workspace, const std::string& array) {
  const std::lock_guard<std::mutex> lock(omicsds_mutex);
  std::shared_ptr<OmicsExporter> instance = std::make_shared<OmicsExporter>(workspace, array);
  omicsds_instances.emplace(std::piecewise_construct, std::make_tuple(current_handle),
                            std::make_tuple(instance));
  return current_handle++;
}

void OmicsDS::disconnect(OmicsDSHandle handle) {
  const std::lock_guard<std::mutex> lock(omicsds_mutex);
  omicsds_instances.erase(handle);
}

void OmicsDS::query(OmicsDSHandle handle, std::array<int64_t, 2>& sample_range, std::array<int64_t, 2>& position_range) {
  std::shared_ptr<OmicsExporter> instance;
  {
    const std::lock_guard<std::mutex> lock(omicsds_mutex);
    instance = omicsds_instances.at(handle);
  }
  instance->query(sample_range, position_range);
}

