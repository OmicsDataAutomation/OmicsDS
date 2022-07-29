/**
 * src/test/cpp/test_omicsds_loader.cc
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
 * Test Matrix Loader
 */

#include <catch2/catch.hpp>
#include "test_base.h"

#include "omicsds_loader.h"
#include "omicsds_export.h"

TEST_CASE_METHOD(TempDir, "test MatrixLoader class", "[matrix]") {
  CHECK(TileDBUtils::is_dir(get_temp_dir()));

  std::string gene_map = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/map_ends.gbed";
  FileUtility gene_reader = FileUtility(gene_map);
  std::string retval;
  REQUIRE(gene_reader.generalized_getline(retval));
  std::vector<std::string> tokens = split(retval, "\t");
  std::string valid_gene_name = tokens[0];

  std::string sample_map = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/small_map";
  FileUtility sample_reader = FileUtility(sample_map);
  REQUIRE(sample_reader.generalized_getline(retval));
  tokens = split(retval, "\t");
  std::string valid_sample_name = tokens[0];

  std::string gene_mapping_file = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/human_g1k_v37.fasta.fai";
  SECTION("test invalid mapping files", "[matrix error]") {
    SECTION("missing sample id", "[matrix error sample]") {
      std::string error_workspace_path = append("error-workspace");
      std::string missing_sample_file = append("missing-sample.resort");
      std::string missing_sample_file_list = append("missing-sample-list");
      REQUIRE(FileUtility::write_file(missing_sample_file_list, missing_sample_file) == TILEDB_OK);

      std::string file_content = "SAMPLE\tMISSING_SAMPLE\nMISSING_GENE\t1";
      REQUIRE(FileUtility::write_file(missing_sample_file, file_content) == TILEDB_OK);

      {
        TranscriptomicsLoader l(error_workspace_path, "array", missing_sample_file_list, sample_map, gene_mapping_file, gene_map, true);
        l.initialize();
        l.import();
      }

      std::vector<std::string> fragments = TileDBUtils::get_fragment_names(error_workspace_path);
      std::vector<std::string> files = TileDBUtils::get_files(error_workspace_path + "/array/" + fragments[0]);

      REQUIRE(files.size() == 2); // Only the metadata should exist
    }
  }
  SECTION("test valid import", "[matrix]") {
    std::string workspace_path = append("workspace");
    std::string good_file = append("good.resort");
    std::string good_file_list = append("good-file-list");
    REQUIRE(FileUtility::write_file(good_file_list, good_file) == TILEDB_OK);

    std::string file_content = "SAMPLE\t" + valid_sample_name + "\n" + valid_gene_name + "\t1";
    REQUIRE(FileUtility::write_file(good_file, file_content) == TILEDB_OK);

    {
      TranscriptomicsLoader l(workspace_path, "array", good_file_list, sample_map, gene_mapping_file, gene_map, true);
      l.initialize();
      l.import();
    }

    OmicsExporter r(workspace_path, "array");
    std::array<int64_t, 2> query_range = {0, std::numeric_limits<long>::max()};

    size_t size = 0;
    r.query(query_range, query_range, [&size](auto a, auto b){size++;});
    REQUIRE(size == 2);
  }
}
