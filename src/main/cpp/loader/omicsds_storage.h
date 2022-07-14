/**
 * src/main/cpp/loader/omicsds_storage.cc
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
 * Header files for OmicsDS storage.
 */

#pragma once

#include "omicsds_schema.h"

#include <string>

#include "tiledb.h"
#include "tiledb_storage.h"

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



