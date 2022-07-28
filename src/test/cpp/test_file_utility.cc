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
    std::string tmp_file = append("read-test");
    SECTION("test common reads", "[utility FileUtility read]") {
      SECTION("test no newline") {
        REQUIRE(FileUtility::write_file(tmp_file, test_text.substr(0, test_text.size()-1)) == TILEDB_OK);
        FileUtility fu = FileUtility(tmp_file);
        
        std::string retval;
        
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == "line1");
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == "line2");
      }

      REQUIRE(FileUtility::write_file(tmp_file, test_text) == TILEDB_OK);

      SECTION("test iterated getline", "[utility FileUtility read getline]") {
        FileUtility fu = FileUtility(tmp_file);
        
        std::string retval;
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == "line1");
        REQUIRE(fu.str_buffer == "line2\n");
        REQUIRE(fu.chars_read == test_text.length());
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == "line2");
        REQUIRE(fu.str_buffer == "");
        REQUIRE(!fu.generalized_getline(retval));
        REQUIRE(retval == "");
      }
      SECTION("test read file", "[utility FileUtility readfile]") {
        FileUtility fu = FileUtility(tmp_file);

        SECTION("test small read") {
          char buf[256];
          size_t read_size = 5;
          REQUIRE(fu.read_file(buf, read_size) == TILEDB_OK);
          REQUIRE(strncmp(buf, test_text.c_str(), read_size) == 0);
          REQUIRE(fu.chars_read == read_size);
          fu.chars_read = 0; // Reset for next test
        }
        SECTION("test full read") {
          char buf[256];
          REQUIRE(fu.read_file(buf, test_text.size()) == TILEDB_OK);
          REQUIRE(strncmp(buf, test_text.c_str(), test_text.size()) == 0);
          REQUIRE(fu.chars_read == test_text.size());
          fu.chars_read = 0; // Reset for next test
        }
        SECTION("test over-extended read") {
          char buf[256];
          REQUIRE(fu.read_file(buf, 256) == TILEDB_ERR);
          REQUIRE(strncmp(buf, test_text.c_str(), test_text.size()) == 0);
        }
        SECTION("test multiple reads") {
          char buf[256];
          size_t read_size = 4;
          REQUIRE(fu.read_file(buf, read_size) == TILEDB_OK);
          REQUIRE(strncmp(buf, test_text.c_str(), read_size) == 0);
          REQUIRE(fu.read_file(buf, read_size) == TILEDB_OK);
          REQUIRE(strncmp(buf, test_text.c_str() + read_size, read_size) == 0);
        }
      }
      SECTION("test combindation reads") {
        SECTION("test read_file -> getline") {
          FileUtility fu = FileUtility(tmp_file);
          char buf[10];
          fu.read_file(buf, strlen("line1\n"));
          REQUIRE(strncmp(buf, "line1\n", strlen("line1\n")) == 0);

          std::string retval;
          REQUIRE(fu.generalized_getline(retval));
          REQUIRE(retval == "line2");
        }
        SECTION("test getline -> read_file") {
          FileUtility fu = FileUtility(tmp_file);

          std::string retval;
          REQUIRE(fu.generalized_getline(retval));
          REQUIRE(retval == "line1");

          char buf[10];
          REQUIRE(fu.read_file(buf, strlen("line2\n")) == TILEDB_OK);
          REQUIRE(strncmp(buf, "line2\n", strlen("line2\n")) == 0);
          REQUIRE(fu.str_buffer == "");
        }
      }
    }
    SECTION("test big writes", "[utility FileUtility read big-file]") {
      SECTION("test file bigger than buffer", "[utility FileUtility read big-file]") {
        std::string large_string = "";
        FileUtility tmp_fu = FileUtility(tmp_file);
        for(int i = 0; i < 2*tmp_fu.buffer_size; i++) {
          large_string = large_string + "a\n";
        }
        REQUIRE(FileUtility::write_file(tmp_file, large_string) == TILEDB_OK);
        FileUtility fu = FileUtility(tmp_file); // Need to recreate the FileUtility class to pick up the write

        std::string retval;
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == "a");
        REQUIRE(fu.chars_read == fu.str_buffer.size() + (retval.length() + 1));
        REQUIRE(fu.chars_read == fu.buffer_size);
        
        SECTION("test read_line after generalized_getline with big write") {
          char buf[1024];
          size_t read_size = fu.str_buffer.size() + 10;
          size_t chars_read = fu.chars_read;
          REQUIRE(fu.read_file(buf, read_size) == TILEDB_OK);
          REQUIRE(chars_read + 10 == fu.chars_read);
          REQUIRE(strncmp(buf, large_string.substr(2, read_size).c_str(), read_size) == 0);
        }
      }
      SECTION("test file bigger than buffer without newline", "[utility FileUtility read big-file]") {
        std::string large_string = "";
        FileUtility tmp_fu = FileUtility(tmp_file);
        for(int i = 0; i < 2*tmp_fu.buffer_size; i++) {
          large_string = large_string + "a";
        }
        large_string = large_string + "\n";
        REQUIRE(FileUtility::write_file(tmp_file, large_string) == TILEDB_OK);
        FileUtility fu = FileUtility(tmp_file); // Need to recreate the FileUtility class to get new write

        std::string retval;
        REQUIRE(fu.generalized_getline(retval));
        REQUIRE(retval == large_string.substr(0, large_string.size() - 1));
      }
    }
  }
}