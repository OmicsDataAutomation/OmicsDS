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

#include "omicsds_schema.h"
#include "omicsds_storage.h"

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
#include <set>
#include <functional>
#include <sstream>
#include <regex>
#include <limits>
#include <cassert>
#include <htslib/sam.h>
#include <numeric>
#include <algorithm>
#include <utility>

void read_sam_file(std::string filename);

// split str into tokens by sep
// similar to java/python split
std::vector<std::string> split(std::string str, std::string sep);

struct contig {
  std::string name;
  uint64_t length;
  uint64_t starting_index;

  contig(const std::string& name, uint64_t length, uint64_t starting_index): name(name), length(length), starting_index(starting_index) {}
  void serialize(std::string path) {
    std::string str = name + "\t" + std::to_string(length) + "\t" + std::to_string(starting_index) + "\n";
    FileUtility::write_file(path, str);
  }
};

// maps from sample name to (logical) row in OmicsDS
struct SampleMap {
  std::map<std::string, size_t> map;
  // file specifying SampleMap should be tab seperated with lines
  // consisting of a single sample name and row index
  SampleMap(const std::string& sample_map);
  size_t& operator[](const std::string& name) {
    return map[name];
  }
  size_t count(const std::string& name) const {
    return map.count(name);
  }
  size_t size() const {
    return map.size();
  }
};

// Maps from gene name to chrom, start, end
// Used in transcriptomics import
struct GeneIdMap {
  struct Gene {
    Gene(const std::string& chrom, uint64_t start, uint64_t end, uint64_t flattened_start, uint64_t flattened_end): chrom(chrom), 
                                                                                                                    start(start), 
                                                                                                                    end(end), 
                                                                                                                    flattened_start(flattened_start), 
                                                                                                                    flattened_end(flattened_end) {
    }

    Gene() = default;

    std::string chrom;
    uint64_t start;
    uint64_t end;
    uint64_t flattened_start;
    uint64_t flattened_end;
  };

  size_t count(const std::string& name) {
    return map.count(name);
  }
  size_t size() {
    return map.size();
  }
  Gene& operator[](const std::string& name) {
    return map[name];
  }

  std::map<std::string, Gene> map;
  std::shared_ptr<OmicsSchema> schema;

  // gene_map is a path to a gtf/gff/gi/gbed file
  // use transcript indicates whether to use transcript_id (true) or gene_id (false)
  // drop_version removes trailing version number (e.g. ENS001.22 -> ENS001)
  GeneIdMap(const std::string& gene_map, std::shared_ptr<OmicsSchema> schema, bool use_transcript = true, bool drop_version = true);
  // standard gtf, looks for genes in "transcript" lines
  // gi/gbed files are preferred as they are more compact and
  // there is less ambiguity around where to look for information
  void create_from_gtf(const std::string& gene_map, bool use_transcript = true, bool drop_version = true);
  // gi is a bespoke binary file format that contains only the gene/transcript names, chromosomes, and offsets
  void create_from_gi(const std::string& gene_map);
  // tab separated file where each line consists of gene name, contig, starting offset in contig, ending offset in contig
  // no header
  void create_from_gbed(const std::string& gene_map);
  void export_as_gi(const std::string& filename);
};

// contains cell information before it is written to disk
struct OmicsCell {
  std::array<int64_t, 2> coords; // sample index, position--does not change with schema order
  std::vector<OmicsFieldData> fields; // must be in schema order (lexicographically sorted by attribute name)
  std::shared_ptr<OmicsSchema> schema;
  int file_idx = -1; // index of file in OmicsLoader::m_files, -1 encodes that OmicsLoader should not try to obtain more cells from any file after processing this (probably end cell)

  OmicsCell() {}
  OmicsCell(std::array<int64_t, 2> coords, std::shared_ptr<OmicsSchema> schema, int file_idx) : coords(coords), schema(schema), file_idx(file_idx), fields(std::vector<OmicsFieldData>(schema->attributes.size())) {}
  OmicsCell(OmicsCell&& o) : coords(o.coords), schema(o.schema), fields(std::move(o.fields)), file_idx(o.file_idx) {}
  OmicsCell(const OmicsCell& o) : coords(o.coords), schema(o.schema), fields(o.fields), file_idx(o.file_idx) {}
  OmicsCell& operator=(const OmicsCell& o) {
    coords = o.coords;
    schema = o.schema;
    fields = o.fields;
    file_idx = o.file_idx;
    return *this;
  }
  // helper function to place field in correct position, if in schema
  // return value indicates successfully added
  template<class T>
  bool add_field(std::string name, T elem) {
    if(!schema) return false;

    auto it = schema->attributes.find(name);
    if(it == schema->attributes.end()) return false;

    size_t idx = std::distance(schema->attributes.begin(), it);
    fields[idx].push_back(elem);
    return true;
  }

  // helper function to place field in correct position, if in schema
  // return value indicates successfully added
  template<class T>
  bool add_field_ptr(std::string name, T* elem_ptr, int n) {
    if(!schema) return false;

    auto it = schema->attributes.find(name);
    if(it == schema->attributes.end()) return false;

    size_t idx = std::distance(schema->attributes.begin(), it);
    fields[idx].push_pointer_back(elem_ptr, n);
    return true;
  }

  bool validate() {
    return fields.size() == schema->attributes.size();
  }

  static OmicsCell create_invalid_cell(); // possibly deprecated because OmicsFileReader::get_next_cells now returns a vector, so can return empty instead of invalid to indicate end of file

  static bool is_invalid_cell(const OmicsCell& cell);

  std::string to_string() {
    std::stringstream ss;
    ss << "beginning of cell {" << coords[0] << ", " << coords[1] << "}" << std::endl;
    auto fiter = fields.begin();
    auto aiter = schema->attributes.begin();
    for(; fiter != fields.end() && aiter != schema->attributes.end(); fiter++, aiter++) {
      ss << "\t\t" << aiter->first << std::endl << "\t\t\t\t";
      for(auto& e : fiter->data) {
        ss << (int)e << "=" << (char)e << " ";
      }
      ss << std::endl;
    }
    return ss.str();
  }

  std::string coords_to_string() {
    return "{" + std::to_string(coords[0]) + ", " + std::to_string(coords[1]) + "}";
  }
};

// utility function to print containers that support foreach loops
template<class T>
std::string container_to_string(const T& c) {
  std::stringstream ss;

  for(auto& e : c) {
    ss << e << " ";
  }

  return ss.str();
}

// class used to ingest information from files
// to support new file formats:
// * derive from OmicsFileReader
// * customize contructor if relevant (e.g. to process file header)
// * override get_next_cells
class OmicsFileReader { // FIXME encode whether cells are end cells (can be deduced by checking position against flattened end, but cumbersome)
  public:
    OmicsFileReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, int file_idx) : /*m_file(filename),*/ m_reader_util(std::make_shared<FileUtility>(filename)), m_schema(schema), m_sample_map(sample_map), m_file_idx(file_idx) {
    }
    virtual ~OmicsFileReader() = default;

    const std::string& get_filename() {
      return m_reader_util->filename;
    }

    // returns cells for use in OmicsLoader::import
    // empty return vector or vector with invalid first element indicate end of file
    // standard order for coords in returned cells is SAMPLE, POSITION regardless of order, will be transformed by loader
    // subsequent calls must return cells that are either at the same coordinates as previously or later
    // i.e. the coords of subsequent cells cannot compare OmicsLoader::less_than previous cells
    virtual std::vector<OmicsCell> get_next_cells() = 0;

  protected:
    std::shared_ptr<OmicsSchema> m_schema;
    std::shared_ptr<SampleMap> m_sample_map;
    int m_file_idx;
    std::shared_ptr<FileUtility> m_reader_util;
};

// uses htslib to read SAM files (must have .sam extension)
// uses file name (without path) as sample name
// TODO (for both SamReader and SamExporter)
// * support BAM format (should only require minimal additions if htslib has good support)
// * potentially transform POS to 0 based and back for query (1 based in SAM)
// * look into RNEXT field, htslib seems to transform =,* into 0,-1, but might also need to write something to header when exporting
// * figure out why htslib is making PNEXT 1 based
// * figure out what htslib is doing to QUAL (goes from ascii to some kind of binary format with mixed ascii (?). Unsure)
class SamReader : public OmicsFileReader {
  public:
    SamReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, int file_idx);
    ~SamReader();
    std::vector<OmicsCell> get_next_cells() override;

  protected:
    uint64_t m_row_idx; // row corresponding to sample
    samFile* m_fp; // file pointer
    bam_hdr_t* m_hdr; // header
    bam1_t* m_align; // alignment
};

// reads ucsc bed files (must have .bed extension)
// looks for sample name in description field of header
class BedReader : public OmicsFileReader {
  public:
    BedReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, int file_idx);
    std::vector<OmicsCell> get_next_cells() override;
  protected:
    std::string m_sample_name;
    uint64_t m_row_idx; // row corresponding to sample
};

// reads matrix files (must have .resort extension)
// expected format is somewhat nonstandard (tab separated)
// either:
// SAMPLE      [sample name] [sample name]
// [gene name] [score]       [score]
//
// or
//
// GENE          [gene name]   [gene name]
// [sample name] [score]       [score]
// 
// Currently matrix file can be either sample or position major, but must match schema order
// TODO: support transposed orders (potentially by reading entire file and transposing in memory)
class MatrixReader : public OmicsFileReader {
  public:
    MatrixReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, std::shared_ptr<GeneIdMap> gene_id_map, int file_idx);
    std::vector<OmicsCell> get_next_cells() override;
  protected:
    std::shared_ptr<GeneIdMap> m_gene_id_map;
    std::vector<std::string> m_columns; // can be samples or genes depending on m_position_major
    bool m_position_major;
    std::vector<float> m_row_scores; // buffer of scores in current row
    size_t m_column_idx = 0; // current column position in matrix
    std::string m_current_token; // can be sample or gene depending on m_position_major
    bool parse_next(std::string& sample, std::string& gene, float& score);
};

// used to ingest information into OmicsDS
// intervals are represented as start and end cells
// to support new file types:
// * derive from OmicsLoader
// * customize constructor if required
// * override create_schema
// * override add_reader
// TODO:
// sharding?
// when some specified pq size is reached stop reading from files and write everything out to fragment
// (continued) create new fragment and resume reading from files
// This also helps with the worst case scenario where the end cells returned from file readers are very far from the starts
// in this case the pq size would increase linearly with reads, rather than staying mostly constant (when end cells are not retained for a very long time in memory)
class OmicsLoader : public OmicsModule {
  public:
    OmicsLoader(
      const std::string& workspace,
      const std::string& array,
      const std::string& file_list, // file with each line containing a path to a data file
      const std::string& sample_map, // see SampleMap struct
      const std::string& mapping_file, // see GenomicMap struct
      bool position_major // sample major (false) or position major (true)
    );
    ~OmicsLoader() {
      tiledb_close_array();
    }
    void import();// import data from callsets
    virtual void create_schema() = 0; //
    void initialize(); // cannot be part of constructor because it invokes create_schema, which is virtual
  protected:
    std::shared_ptr<SampleMap> m_sample_map;
    int tiledb_write_buffers();
    // data that will be given to tiledb
    std::vector<std::vector<uint8_t>> offset_buffers; // encodes offset for variable length attributes, and data for constant length ones
    std::vector<std::vector<uint8_t>> var_buffers; // entries for constant length attributes will be empty
    std::vector<size_t> coords_buffer; // keeps 3d coords
    std::vector<size_t> attribute_offsets; // persists between writes
    size_t buffer_size = 10240; // size of tiledb buffers

    // insert relevant information in offset_buffers
    // and var_buffers (for variable fields)
    void buffer_cell(const OmicsCell& cell, int level = 0);

    virtual void add_reader(const std::string& filename) = 0; // construct a derived class of OmicsFileReader and insert in m_files
    std::string m_file_list;
    std::vector<std::shared_ptr<OmicsFileReader>> m_files;
    typedef std::shared_ptr<OmicsFileReader> omics_fptr;
    std::priority_queue<OmicsCell, std::vector<OmicsCell>, std::function<bool(OmicsCell, OmicsCell)>> m_pq;
    int m_idx;
    static bool comparitor(OmicsCell _l, OmicsCell _r) { // returns which cell is first in schema order
      auto l = _l.coords;
      auto r = _r.coords;
      return (l[0] > r[0]) || (l[0] == r[0] && l[1] > r[1]);
    }

    bool less_than(const std::array<int64_t, 2>& l, const std::array<int64_t, 2>& r) {
      return (l[0] < r[0]) || (l[0] == r[0] && l[1] < r[1]);
    }
    bool less_than(const OmicsCell& l, const OmicsCell& r) {
      return less_than(l.coords, r.coords);
    }
    void push_from_idxs(const std::set<int>& idxs);
    void push_file_from_cell(const OmicsCell& cell);
    void push_from_all_files();
};

// used to ingest SAM files
class ReadCountLoader : public OmicsLoader {
  public:
    ReadCountLoader(
      const std::string& workspace,
      const std::string& array,
      const std::string& file_list,
      const std::string& sample_map,
      const std::string& mapping_file,
      bool position_major
    ) : OmicsLoader(workspace, array, file_list, sample_map, mapping_file, position_major) {
    }
    virtual void create_schema() override;
  protected:
    virtual void add_reader(const std::string& filename) override;
};

// used to ingest Bed and Matrix files
class TranscriptomicsLoader : public OmicsLoader {
  public:
    TranscriptomicsLoader(
      const std::string& workspace,
      const std::string& array,
      const std::string& file_list,
      const std::string& sample_map,
      const std::string& mapping_file,
      const std::string& gene_mapping_file,
      bool position_major
    ) : OmicsLoader(workspace, array, file_list, sample_map, mapping_file, position_major), m_gene_id_map(std::make_shared<GeneIdMap>(gene_mapping_file, m_schema)) {
    }
    virtual void create_schema() override;
  protected:
    virtual void add_reader(const std::string& filename) override;
    std::shared_ptr<GeneIdMap> m_gene_id_map;
};
