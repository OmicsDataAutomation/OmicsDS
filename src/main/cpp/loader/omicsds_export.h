#pragma once

#include "omicsds_storage.h"

#include <functional>

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
