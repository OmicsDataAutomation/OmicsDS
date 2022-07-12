#pragma once
//#include <htslib/sam.h>
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
#include "tiledb.h"
#include "tiledb_utils.h"
#include "tiledb_storage.h"

#define CHECK_RC(...)                                      \
do {                                                       \
  int rc = __VA_ARGS__;                                    \
  if (rc) {                                                \
    printf("%s", &tiledb_errmsg[0]);                       \
    printf("[Examples::%s] Runtime Error.\n", __FILE__);   \
    return rc;                                             \
  }                                                        \
} while (false)

void read_sam_file(std::string filename);

// split str into tokens by sep
// similar to java/python split
std::vector<std::string> split(std::string str, std::string sep);

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
    static std::string cigar_to_string(const uint32_t* cigar, size_t n_cigar);

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

// common base class for OmicsLoader and OmicsExporter with various common functionalities
class OmicsModule {
  public:
    OmicsModule(const std::string& workspace, const std::string& array) : m_workspace(workspace), m_array(array) {}
    OmicsModule(const std::string& workspace, const std::string& array, const std::string& mapping_file, bool position_major) : m_workspace(workspace), m_array(array), m_schema(std::make_shared<OmicsSchema>(mapping_file, position_major)) {}
    void serialize_schema(std::string path) { m_schema->serialize(path); }
    void serialize_schema() { serialize_schema(m_workspace + "/" + m_array + "/omics_schema"); }
    void deserialize_schema(std::string path) { m_schema.reset(); m_schema = std::make_shared<OmicsSchema>(); m_schema->create_from_file(path); }
    void deserialize_schema() { deserialize_schema(m_workspace + "/" + m_array + "/omics_schema"); };

  protected:
    std::string m_workspace;
    std::string m_array;
    // tiledb* functions return tiledb return code
    int tiledb_create_array(const std::string& workspace, const std::string& array_name, const OmicsSchema& schema);
    int tiledb_create_array() { return tiledb_create_array(m_workspace, m_array, *m_schema); }
    int tiledb_open_array(const std::string& workspace, const std::string& array_name, bool write = true);
    int tiledb_open_array(bool write = true) { return tiledb_open_array(m_workspace, m_array, write); }
    int tiledb_close_array();
    TileDB_CTX* m_tiledb_ctx;
    TileDB_Array* m_tiledb_array;
    std::shared_ptr<OmicsSchema> m_schema;
    int m_array_descriptor;
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

// used to query from OmicsDS
class OmicsExporter : public OmicsModule {
  public:
    OmicsExporter(const std::string& workspace, const std::string& array) : OmicsModule(workspace, array) {
      deserialize_schema();
      tiledb_open_array(false);
    }
  
    typedef std::function<void (const std::array<uint64_t, 3>& coords, const std::vector<OmicsFieldData>& data)> process_function;
    // used to query given range
    // will use proc as callback if specified, otherwise will default to process
    void query(std::array<int64_t, 2> sample_range = {0, std::numeric_limits<int64_t>::max()}, std::array<int64_t, 2> position_range = {0, std::numeric_limits<int64_t>::max()}, process_function proc = 0);

  protected:
    // coords are in standard order SAMPLE, POSITION, COLLISION INDEX
    virtual void process(const std::array<uint64_t, 3>& coords, const std::vector<OmicsFieldData>& data);
    std::vector<std::vector<uint8_t>> m_buffers_vector;
    std::pair<std::vector<void*>, std::vector<size_t>> prepare_buffers();
    size_t m_buffer_size = 10240;
    void check(const std::string& name, const OmicsFieldInfo& inf); // check that an attribute exists in schema (useful for specific data e.g. ensure that the data is actually from readcounts/sam files before exporting
};

// for exporting data as SAM files
// will create one per row with name sam_output[row idx].sam
// rows with no data in query range will not appear in output
class SamExporter : public OmicsExporter { // for exporting data as SAM files
  public:
    SamExporter(const std::string& workspace, const std::string& array);
    void export_sams(std::array<int64_t, 2> sample_range = {0, std::numeric_limits<int64_t>::max()}, std::array<int64_t, 2> position_range = {0, std::numeric_limits<int64_t>::max()}, const std::string& ouput_prefix = "sam_output");

  protected:
    // callback to write to sam files
    void sam_interface(std::map<int64_t, std::shared_ptr<std::ofstream>>& files, const std::string& output_prefix, const std::array<uint64_t, 3>& coords, const std::vector<OmicsFieldData>& data);
};
