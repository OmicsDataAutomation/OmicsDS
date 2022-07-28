#
# CMakeLists.txt
#
# The MIT License
#
# Copyright (c) 2022 Omics Data Automation, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#


include_directories(Catch INTERFACE Catch2/single_include)

set(CPP_TEST_SOURCES
        test_omicsds_loader.cc
        test_omics_field_data.cc
        test_file_utility.cc
        test_logger.cc)
add_executable(ctests ${CPP_TEST_SOURCES})
target_include_directories(ctests PRIVATE ${CMAKE_SOURCE_DIR}/src/main/cpp/loader ${CMAKE_SOURCE_DIR}/src/main/cpp/utils)
target_compile_definitions(ctests PRIVATE -DOMICSDS_TEST_INPUTS="${CMAKE_CURRENT_SOURCE_DIR}/../inputs/")
target_link_libraries(ctests omicsds_static ${OMICSDS_DEPENDENCIES} ${CMAKE_DL_LIBS})
add_test(ctests ctests -d yes)
