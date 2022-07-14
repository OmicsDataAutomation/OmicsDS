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

set(OMICSDS_CPP "${CMAKE_CURRENT_SOURCE_DIR}/cpp")

set(OMICSDB_INCLUDE
  ${OMICSDS_CPP}/loader/omicsdb_loader.h
  )

set(OMICSDS_SOURCES
  ${OMICSDS_CPP}/loader/omicsds_loader.cc
  ${OMICSDS_CPP}/loader/omicsds_export.cc
  )

# Use PIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)
add_library(omicsds_objs OBJECT ${OMICSDS_SOURCES})

target_include_directories(omicsds_objs
  PUBLIC ${CMAKE_BINARY_DIR}/src
  PRIVATE ${OMICSDS_CPP} ${TILEDB_INCLUDE_DIR}
  )

# Create static library
add_library(omicsds_static STATIC $<TARGET_OBJECTS:omicsds_objs>)
set_target_properties(omicsds_static PROPERTIES OUTPUT_NAME "omicsds" VERSION ${OMICSDS_VERSION})
add_dependencies(omicsds_static htslib)
target_link_libraries(omicsds_static
  tiledb_static
  ${OMICSDS_DEPENDENCIES}
  ${CMAKE_DL_LIBS} # add -ldl
  ) 

# Create shared library
add_library(omicsds_shared SHARED $<TARGET_OBJECTS:omicsds_objs>)
set_target_properties(omicsds_shared PROPERTIES OUTPUT_NAME "omicsds" VERSION ${OMICSDS_VERSION} SOVERSION 1)
add_dependencies(omicsds_shared htslib)
target_link_libraries(omicsds_shared 
  tiledb_static
  ${OMICSDS_DEPENDENCIES}
  pthread
  ${CMAKE_DL_LIBS} # add -ldl
  )

# Install libraries
install(
  TARGETS omicsds_static omicsds_shared
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

install(
  FILES ${OMICSDS_INCLUDE} ${EXPORTS_HEADER} DESTINATION include
  )
