/**
 * src/main/cpp/loader/omicsds_schema.cc
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
 * Header files for OmicsDS schema.
 */

#pragma once

#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "tiledb_constants.h"
#include "tiledb_utils.h"

#define CHECK_RC(...)                                      \
do {                                                       \
  int rc = __VA_ARGS__;                                    \
  if (rc) {                                                \
    printf("%s", &tiledb_errmsg[0]);                       \
    printf("[Examples::%s] Runtime Error.\n", __FILE__);   \
    return rc;                                             \
  }                                                        \
} while (false)

// for reading/writing local/cloud files using TileDBUtils api
struct FileUtility {
  // only need to construct for reading, write functionality is static
  FileUtility(const std::string& filename): filename(filename) {
    if(!TileDBUtils::is_file(filename)) {
      std::cerr << "Note: file " << filename << " does not exist" << std::endl;
    }
    buffer = new char[buffer_size];
    file_size = TileDBUtils::file_size(filename);
    //      m_file_size = 0;
  }
  ~FileUtility() {
    delete[] buffer;
  }

  std::string filename;
  ssize_t file_size = 0;
  ssize_t chars_read = 0;
  const int buffer_size = 512;
  char* buffer;
  std::string str_buffer;

  // returns true if line was read
  // provides similar functionality to std::getline but also supports cloud files
  bool generalized_getline(std::string& retval);
  // read specified number of bytes from file, 
  // should work with generalized_getline but not tested
  // returns tiledb return code
  int read_file(void* buffer, size_t chars_to_read);

  // write string to file
  // returns tiledb return code
  static int write_file(std::string filename, const std::string& str, const bool overwrite=false) {
    auto rcode = TileDBUtils::write_file(filename, str.c_str(), str.size(), overwrite);
    CHECK_RC(rcode);
    return rcode;
  }
  // write buffer to file
  // returns tiledb return code
  static int write_file(std::string filename, const void* buffer, size_t length, const bool overwrite=false) {
    auto rcode = TileDBUtils::write_file(filename, buffer, length, overwrite);
    CHECK_RC(rcode);
    return rcode;
  }
};

// datastructure that keeps contigs sorted by name and position
class GenomicMap {
public:
  GenomicMap() {}
  // construct from file at path mapping_file
  GenomicMap(const std::string& mapping_file);
  // construct from extant FileUtility
  // used to deserialize omics schema after attributes are read from same file
  GenomicMap(std::shared_ptr<FileUtility> mapping_reader);
  // map from contig_name and offset to single coordinate for use in tiledb
  uint64_t flatten(std::string contig_name, uint64_t offset);
  // reverse of flatten
  std::pair<std::string, uint64_t> unflatten(uint64_t position);
  // human readably serialize contig information
  // used as subset of serialized schema
  void serialize(std::string path);

  // struct to represent contigs in memory
  struct contig {
    std::string name;
    uint64_t length;
    uint64_t starting_index;

    contig(const std::string& name, uint64_t length, uint64_t starting_index): name(name), length(length), starting_index(starting_index) {}
    // serialize individual contig
    // used by GenomicMap::serialize
    void serialize(std::string path) {
      std::string str = name + "\t" + std::to_string(length) + "\t" + std::to_string(starting_index) + "\n";
      FileUtility::write_file(path, str);
    }
  };

private:
  std::shared_ptr<FileUtility> m_mapping_reader;
  std::vector<contig> contigs;
  // indices from GenomicMap::contigs sorted by contig name
  // used by flatten
  std::vector<size_t> idxs_name;
  // indices from GenomicMap::contigs sorted by starting offeset
  // used by unflatten
  std::vector<size_t> idxs_position;
};

// struct to store type/length information about fields
// also provides utility fucntions to go to/from string/tiledb types
struct OmicsFieldInfo {
  enum OmicsFieldType { omics_char, omics_uint8_t, omics_int8_t,
                        omics_uint16_t, omics_int16_t, omics_uint32_t,
                        omics_int32_t, omics_uint64_t, omics_int64_t, omics_float_t };

  // OmicsFieldInfo treats negative lengths as variable
  OmicsFieldInfo(OmicsFieldType type, int _length) : type(type) {
    if(_length < 0) {
      length = TILEDB_VAR_NUM;
    }
    else {
      length = _length;
    }
  }
  
  // string constructor, used when deserializing human readable schema
  OmicsFieldInfo(const std::string& stype, int _length) {
    if(_length < 0) {
      length = TILEDB_VAR_NUM;
    }
    else {
      length = _length;
    }
    if(stype == "omics_char") { type = omics_char; return; }
    if(stype == "omics_uint8_t") { type = omics_uint8_t; return; }
    if(stype == "omics_int8_t") { type = omics_int8_t; return; }
    if(stype == "omics_uint16_t") { type = omics_uint16_t; return; }
    if(stype == "omics_int16_t") { type = omics_int16_t; return; }
    if(stype == "omics_uint32_t") { type = omics_uint32_t; return; }
    if(stype == "omics_int32_t") { type = omics_int32_t; return; }
    if(stype == "omics_uint64_t") { type = omics_uint64_t; return; }
    if(stype == "omics_int64_t") { type = omics_int64_t; return; }
    if(stype == "omics_float_t") { type = omics_float_t; return; }
    type = omics_uint8_t;
    return;
  }

  OmicsFieldType type;
  int length; // number of elements, -1 encodes variable

  int tiledb_type() const {
    switch(type) {
      case omics_char:     return TILEDB_CHAR;
      case omics_uint8_t:  return TILEDB_UINT8;
      case omics_int8_t:   return TILEDB_INT8;
      case omics_uint16_t: return TILEDB_UINT16;
      case omics_int16_t:  return TILEDB_INT16;
      case omics_uint32_t: return TILEDB_UINT32;
      case omics_int32_t:  return TILEDB_INT32;
      case omics_uint64_t: return TILEDB_UINT64;
      case omics_int64_t:  return TILEDB_INT64;
      case omics_float_t:  return TILEDB_FLOAT32;
    }
    return TILEDB_CHAR;
  }

  // used to serialize schema
  std::string type_to_string() const {
    switch(type) {
      case omics_char:     return "omics_char";
      case omics_uint8_t:  return "omics_uint8_t";
      case omics_int8_t:   return "omics_int8_t";
      case omics_uint16_t: return "omics_uint16_t";
      case omics_int16_t:  return "omics_int16_t";
      case omics_uint32_t: return "omics_uint32_t";
      case omics_int32_t:  return "omics_int32_t";
      case omics_uint64_t: return "omics_uint64_t";
      case omics_int64_t:  return "omics_int64_t";
      case omics_float_t:  return "omics_float_t";
    }
    return "unknown_type";
  }

  std::string length_to_string() const {
    if(length == TILEDB_VAR_NUM) {
      return "variable";
    }
    return std::to_string(length);
  }

  int element_size() {
    switch(type) {
      case omics_char:     return 1;
      case omics_uint8_t:  return 1;
      case omics_int8_t:   return 1;
      case omics_uint16_t: return 2;
      case omics_int16_t:  return 2;
      case omics_uint32_t: return 4;
      case omics_int32_t:  return 4;
      case omics_uint64_t: return 8;
      case omics_int64_t:  return 8;
      case omics_float_t:  return 4;
    }
    return 1;
  }

  bool is_variable() {
    return length == TILEDB_VAR_NUM;
  }

  bool operator==(const OmicsFieldInfo& o) {
    return type == o.type && length == o.length;
  }
};

// schema for array
// contains attributes and a GenomicMap
struct OmicsSchema {
  enum OmicsStorageOrder { POSITION_MAJOR, SAMPLE_MAJOR };
  OmicsStorageOrder order;

  OmicsSchema() {}
  OmicsSchema(const std::string& mapping_file, OmicsStorageOrder order = POSITION_MAJOR): genomic_map(mapping_file), order(order) {}
  OmicsSchema(const std::string& mapping_file, bool position_major = true): genomic_map(mapping_file) {
    order = position_major ? POSITION_MAJOR : SAMPLE_MAJOR;
  }
  bool create_from_file(const std::string& path);
  bool position_major() const {
    return order == POSITION_MAJOR;
  }
  // swaps between standard and schema order
  // standard order is SAMPLE, POSITION, will swap if position major
  template<class T, size_t U>
  std::array<T, U> swap_order(const std::array<T, U>& coords) {
    std::array<T, U> retval = coords;
    if(position_major()) {
      std::swap(retval[0], retval[1]);
    }
    return retval;
  }
  std::map<std::string, OmicsFieldInfo> attributes; // implies canonical order (lexicographically sorted by name)
  GenomicMap genomic_map;
  void serialize(std::string path);
  // get index of attribute by name
  // useful because fields in OmicsCell are in the same order as in OmicsSchema
  int index_of_attribute(const std::string& name);
};

bool equivalent_schema(const OmicsSchema& l, const OmicsSchema& r);

// stores field data backed by vector
struct OmicsFieldData {
  std::vector<uint8_t> data;
  size_t size() const {
    return data.size();
  }
  uint8_t operator[](size_t idx) {
    return data[idx];
  }
  // templated utilty function to insert elements of potentially more than one byte into vector
  template<class T>
  static void push_back(std::vector<uint8_t>& v, const T& elem) {
    auto size = v.size();
    v.resize(v.size() + sizeof(elem));
    T* ptr = reinterpret_cast<T*>(v.data() + size);
    *ptr = elem;
  }

  template<class T>
  void push_back(const T& elem) {
    push_back(data, elem);
  }

  // templated utilty function to insert elements of potentially more than one byte into vector via pointer
  // useful for c strings/vectors
  template<class T>
  void push_pointer_back(const T* elem_ptr, int n) {
    size_t size = data.size();
    data.resize(data.size() + sizeof(T)*n);
    T* ptr = reinterpret_cast<T*>(data.data() + size);

    for(int i = 0; i < n; i++) {
      ptr[i] = elem_ptr[i];
    }
  }
  
  template<class T>
  const T* get_ptr() const {
    return (T*)data.data();
  }
  template<class T>
  T get(int idx = 0) const { // FIXME check bounds?
    return ((T*)data.data())[idx];
  }
  template<class T>
  int typed_size() const {
    return data.size() / sizeof(T);
  }
};




