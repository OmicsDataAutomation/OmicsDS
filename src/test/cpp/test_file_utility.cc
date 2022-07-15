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
 * Test generic SAM reader
 */

#include <catch2/catch.hpp>
#include "test_base.h"

#include "omicsds_schema.h"
#include "tiledb_constants.h"

TEST_CASE_METHOD(TempDir, "test FileUtility", "[utility FileUtility]") {
  std::string test_text = "line1\nline2\n";
  SECTION("test writes", "[utility FileUtility write]") {
    std::string tmp_file = append("write-test");
    REQUIRE(FileUtility::write_file(tmp_file, test_text) == TILEDB_OK);
  }
  SECTION("test reads", "[utility FileUtility read]") {
    SECTION("test iterated getline", "[utility FileUtility read getline]") {
      std::string tmp_file = append("write-test");
      REQUIRE(FileUtility::write_file(tmp_file, test_text) == TILEDB_OK);
      FileUtility fu = FileUtility(tmp_file);
      
      std::string retval;
      fu.generalized_getline(retval);
      REQUIRE(retval == "line1");
      REQUIRE(fu.str_buffer == "line2\n");
      REQUIRE(fu.chars_read == test_text.length());
      fu.generalized_getline(retval);
      REQUIRE(retval == "line2");
      REQUIRE(fu.str_buffer == "");
      fu.generalized_getline(retval);
      REQUIRE(retval == "");
    }
    SECTION("test file bigger than buffer", "[utility FileUtility read big-file]") {
      std::string tmp_file = append("big-write-test");
      std::string large_string = "";
      FileUtility tmp_fu = FileUtility(tmp_file);
      for(int i = 0; i < 2*tmp_fu.buffer_size; i++) {
        large_string = large_string + "a\n";
      }
      REQUIRE(FileUtility::write_file(tmp_file, large_string) == TILEDB_OK);
      FileUtility fu = FileUtility(tmp_file); // Need to recreate the FileUtility class to get new write

      std::string retval;
      fu.generalized_getline(retval);
      REQUIRE(retval == "a");
      REQUIRE(fu.chars_read == fu.str_buffer.size() + (retval.length() + 1));
      REQUIRE(fu.chars_read == fu.buffer_size);
    }
  }
}