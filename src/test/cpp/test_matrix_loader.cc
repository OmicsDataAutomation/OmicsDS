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
#include "omicsds_logger.h"

TEST_CASE_METHOD(TempDir, "test MatrixLoader class", "[matrix]") {
  CHECK(TileDBUtils::is_dir(get_temp_dir()));
  SECTION("test invalid mapping files", "[matrix error]") {
    std::string sample_map = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/small_map";
    std::string gene_map = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/map_ends.gbed";
    std::string gene_mapping_file = std::string(OMICSDS_TEST_INPUTS)+"OmicsDSTests/human_g1k_v37.fasta.fai";
    SECTION("missing sample id", "[matrix error sample]") {
      std::string error_workspace_path = append("error-workspace");
      std::string missing_sample_file = append("missing-sample.resort");
      std::string missing_sample_file_list = append("missing-sample-list");
      REQUIRE(FileUtility::write_file(missing_sample_file_list, missing_sample_file) == TILEDB_OK);

      FileUtility gene_reader = FileUtility(gene_map);
      std::string retval;
      REQUIRE(gene_reader.generalized_getline(retval));
      std::vector<std::string> tokens = split(retval, "\t");
      std::string valid_gene_name = tokens[0];

      std::string file_content = "SAMPLE\tMISSING_SAMPLE\n" + valid_gene_name + "\t1";
      REQUIRE(FileUtility::write_file(missing_sample_file, file_content) == TILEDB_OK);

      TranscriptomicsLoader l(error_workspace_path, "array", missing_sample_file_list, sample_map, gene_mapping_file, gene_map, true);
      l.initialize();
      l.import();
    }
  }
}
