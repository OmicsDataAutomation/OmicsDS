/**
 * src/main/cpp/loader/omicsds_loader.h
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
 * Specification for a generic SAM reader
 */

#pragma once

#include <array>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <memory>
#include <cmath>
#include <string>
#include <deque>
#include <queue>
#include <functional>
#include <sstream>
#include <regex>
#include <htslib/sam.h>
#include "tiledb_utils.h"

void read_sam_file(std::string filename);

struct OmicsFieldInfo {
  enum OmicsFieldType { omics_char, omics_uint8_t, omics_int8_t,
                                omics_uint16_t, omics_int16_t, omics_uint32_t,
                                omics_int32_t, omics_uint64_t, omics_int64_t };

  OmicsFieldInfo(OmicsFieldType type) : type(type) {}
  OmicsFieldType type;
  bool operator==(const OmicsFieldInfo& o) {
    return type == o.type;
  }
};

struct OmicsFieldData {
  std::vector<uint8_t> data;
  size_t size() const {
    return data.size();
  }
  uint8_t operator[](size_t idx) {
    return data[idx];
  }
  template<class T>
  void push_back(const T& elem) {
    auto size = data.size();
    data.resize(data.size() + sizeof(elem));
    T* ptr = reinterpret_cast<T*>(data.data() + size);
    *ptr = elem;
  }

  template<class T>
  void push_pointer_back(T* elem_ptr, int n) {
    size_t size = data.size();
    data.resize(data.size() + sizeof(T)*n);
    T* ptr = reinterpret_cast<T*>(data.data() + size);

    for(int i = 0; i < n; i++) {
      ptr[i] = elem_ptr[i];
    }
  }
};

typedef std::map<std::string, OmicsFieldInfo> OmicsSchema; // implies canonical order
bool equivalent_schema(const OmicsSchema& l, const OmicsSchema& r);

struct OmicsCell {
  std::vector<OmicsFieldData> fields; // must be in schema order
  std::shared_ptr<OmicsSchema> schema;
  int file_idx = -1;

  OmicsCell(std::shared_ptr<OmicsSchema> schema, int file_idx) : schema(schema), file_idx(file_idx), fields(std::vector<OmicsFieldData>(schema->size())) {}
  OmicsCell(OmicsCell&& o) : schema(o.schema), fields(std::move(o.fields)), file_idx(o.file_idx) {}
  OmicsCell(const OmicsCell& o) : schema(o.schema), fields(o.fields), file_idx(o.file_idx) {}
  OmicsCell& operator=(const OmicsCell& o) {
    schema = o.schema;
    fields = o.fields;
    file_idx = o.file_idx;
    return *this;
  }
  // return value indicates if successfully added
  template<class T>
  bool add_field(std::string name, T elem) { // helper function to place field in correct position, if in schema
    if(!schema) return false;

    auto it = schema->find(name);
    if(it == schema->end()) return false;

    size_t idx = std::distance(schema->begin(), it);
    fields[idx].push_back(elem);
    return true;
  }

  template<class T>
  bool add_field_ptr(std::string name, T* elem_ptr, int n) {
    if(!schema) return false;

    auto it = schema->find(name);
    if(it == schema->end()) return false;

    size_t idx = std::distance(schema->begin(), it);
    fields[idx].push_pointer_back(elem_ptr, n);
    return true;
  }
};

struct OmicsMultiCell {
  std::vector<uint8_t> as_cell();

  std::shared_ptr<OmicsSchema> schema;  
  std::array<int64_t, 2> coords;
  std::vector<OmicsCell> subcells;

  OmicsMultiCell() {}
  OmicsMultiCell(std::array<int64_t, 2> coords, std::shared_ptr<OmicsSchema> schema) : coords(coords), schema(schema) {}

  bool validate_cell(const OmicsCell& cell) {
    return cell.fields.size() == schema->size();
  }

  size_t push_back(OmicsCell&& elem) {
    if(!validate_cell(elem)) {
      return -1;
    }
    subcells.emplace_back(elem);
    return subcells.size() - 1;
  }

  size_t push_back(const OmicsCell& elem) {
    if(!validate_cell(elem)) {
      return -1;
    }
    subcells.push_back(elem);
    return subcells.size() - 1;
  }

  size_t push_empty_cell(int file_idx = -1) {
    subcells.emplace_back(schema, file_idx);
    return subcells.size() - 1;
  }

  size_t size() {
    return subcells.size();
  }

  bool merge(const OmicsMultiCell& o) {
    if(!equivalent_schema(*schema, *(o.schema))) {
      return false;
    }
    subcells.insert(subcells.begin(), o.subcells.begin(), o.subcells.end());
    return true;
  }

  static OmicsMultiCell create_invalid_cell();

  std::vector<std::vector<uint8_t>> operator[](size_t idx) {
    std::vector<std::vector<uint8_t>> rv;
    for(auto& oc : subcells) {
      rv.push_back(oc.fields[idx].data);
    }
    return rv;
  }

  std::string to_string() {
    std::stringstream ss;
    ss << "beginning of cell {" << coords[0] << ", " << coords[1] << "}" << std::endl;

    for(auto& sc : subcells) {
      ss << "\tsubcell" << std::endl;
      for(auto& f : sc.fields) {
        ss << "\t\tfield" << std::endl << "\t\t\t";
        for(auto& e : f.data) {
          ss << (int)e << "=" << (char)e << " ";
        }
        ss << std::endl;
      }
    }

    return ss.str();
  }

  std::string coords_to_string() {
    return "{" + std::to_string(coords[0]) + ", " + std::to_string(coords[1]) + "}";
  }

  std::vector<int> get_file_idxs() { // returns file indices that are not -1 (not associated to a file/end cell)
    std::vector<int> rv;
    for(auto& sc : subcells) {
      if(sc.file_idx >= 0) {
        rv.push_back(sc.file_idx);
      }
    }
    return rv;
  }
};

class OmicsFileReader {
  public:
    OmicsFileReader(std::string filename, std::shared_ptr<OmicsSchema> schema, int file_idx) : /*m_file(filename),*/ m_filename(filename), m_schema(schema), m_file_idx(file_idx) {
      m_buffer = new char[m_buffer_size];
      //m_file_size = TileDBUtils::file_size(filename);
      m_file_size = 0;
      // FIXME get file size
      std::cout << "REMOVE FIXME get file size" << std::endl;
    }
    ~OmicsFileReader() {
      delete[] m_buffer;
    }

    /*const std::array<int64_t, 2> peek() {
      if(manage_buffer()) {
        return {-1, -1};
      }

      return m_cell_buffer.front().coords;
    }*/

    /*OmicsMultiCell pop() {
      OmicsMultiCell rv;

      if(!manage_buffer()) {
        rv.coords = {-1, -1};
      }
      else {
        rv = m_cell_buffer.front();
        m_cell_buffer.pop_front();
      }

      return rv;
    }*/

    const std::string& get_filename() {
      return m_filename;
    }

    //OmicsMultiCell next_cell_info();
    virtual std::vector<OmicsMultiCell> get_next_cells() = 0;

  protected:
    std::shared_ptr<OmicsSchema> m_schema;
    int m_file_idx;
    std::string m_filename;
    //std::deque<OmicsMultiCell> m_cell_buffer;
    //bool manage_buffer();
    //std::ifstream m_file;
    ssize_t m_file_size = 0;
    ssize_t m_chars_read = 0;
    bool generalized_getline(std::string& retval);

    const int m_buffer_size = 512;
    char* m_buffer;
    std::string m_str_buffer;
};

class SamReader : public OmicsFileReader {
  public:
    SamReader(std::string filename, std::shared_ptr<OmicsSchema> schema, int file_idx);
    ~SamReader();
    std::vector<OmicsMultiCell> get_next_cells() override;

  protected:
    samFile* m_fp; // file pointer
    bam_hdr_t* m_hdr; // header
    bam1_t* m_align; // alignment
};

// intervals are represented as start and end cells, but have no special query
// support as intervals cannot always be prevented from overlapping
class OmicsLoader {
  public:
    enum OmicsStorageOrder { COLUMN_MAJOR, ROW_MAJOR };
    OmicsLoader(
      const std::string& file_list,
      OmicsStorageOrder order = COLUMN_MAJOR,
      const bool superimpose = false
    );
    void import();// import data from callsets
    //virtual void query() = 0; // query
  protected:
    bool m_superimpose; // whether to contain data for multiple logical cells within one cell
    std::shared_ptr<OmicsSchema> m_schema;
    //std::vector<std::string> m_filenames; // list of input file names
    std::vector<std::shared_ptr<OmicsFileReader>> m_files;
    typedef std::shared_ptr<OmicsFileReader> omics_fptr;
    //std::priority_queue<omics_fptr, std::vector<omics_fptr>, std::function<bool(omics_fptr, omics_fptr)>> m_pq;
    std::priority_queue<OmicsMultiCell, std::vector<OmicsMultiCell>, std::function<bool(OmicsMultiCell, OmicsMultiCell)>> m_pq;
    int m_idx;
    int m_array_descriptor;
    OmicsStorageOrder m_order;
    static bool m_column_major_comparitor(OmicsMultiCell _l, OmicsMultiCell _r) {
      std::cout << "REMOVE *********************************** m_column_major_comparitor" << std::endl;
      std::cout << _l.coords_to_string() << " > " << _r.coords_to_string() << std::endl;
      auto l = _l.coords;
      auto r = _r.coords;
      std::cout << "REMOVE " << ((l[1] > r[1]) || (l[1] == r[1] && l[0] > r[0])) << std::endl;
      return (l[1] > r[1]) || (l[1] == r[1] && l[0] > r[0]);
    };
    static bool m_row_major_comparitor(OmicsMultiCell _l, OmicsMultiCell _r) {
      std::cout << "REMOVE *********************************** m_row_major_comparitor" << std::endl;
      std::cout << _l.coords_to_string() << " > " << _r.coords_to_string() << std::endl;
      auto l = _l.coords;
      auto r = _r.coords;
      std::cout << "REMOVE " << ((l[0] > r[0]) || (l[0] == r[0] && l[1] > r[1])) << std::endl;
      return (l[0] > r[0]) || (l[0] == r[0] && l[1] > r[1]);
    };
    bool less_than(const std::array<int64_t, 2>& l, const std::array<int64_t, 2>& r) {
      if(m_order == COLUMN_MAJOR) {
        return (l[1] > r[1]) || (l[1] == r[1] && l[0] > r[0]);
      }
      else return (l[0] > r[0]) || (l[0] == r[0] && l[1] > r[1]);
    }
    bool less_than(const OmicsMultiCell& l, const OmicsMultiCell& r) {
      return l.coords < r.coords;
    }
    // bool info_to_cell(const transcriptomics_cell& tc);
    //cell_info next_cell_info(std::ifstream& file, int type, int ind);
    // std::map<std::string, std::pair<long, long>> transcript_map;
    // void read_uncompressed_gtf(std::istream& input, std::string format);
    // void read_compressed_gtf(std::string filename);
    // void serialize_transcript_map(std::ostream& output);
    // void deserialize_transcript_map(std::istream&);
};
