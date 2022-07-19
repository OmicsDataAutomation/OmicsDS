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
 * Test OmicsFieldData
 */

#include <catch2/catch.hpp>
#include "test_base.h"

#include "omicsds_schema.h"

typedef struct {
  int x;
  float y;
  size_t z;
} test_struct;

TEST_CASE("test OmicsFieldData class", "[omicsfielddata]") {
  OmicsFieldData ofd;
  SECTION("test push_back", "[omicsfielddata push_back]") {
    ofd.push_back(5);
    ofd.push_back('c');
    test_struct my_struct = { .x = 42, .y = 3.14, .z = 15903};
    ofd.push_back(my_struct);

    std::string string_arr[3] = {"string0", "string1", "string2"};
    ofd.push_pointer_back(string_arr, 3);

    SECTION("test array gets", "[omicsfielddata operator-bracket]") {
      size_t byte_offset = 0;
      REQUIRE((int) ofd[byte_offset] == 5);
      byte_offset += sizeof(int);

      REQUIRE((char) ofd[byte_offset] == 'c');
      byte_offset += sizeof(char);

      test_struct* ret_struct = (test_struct*) (ofd.get_ptr<uint8_t>() + byte_offset);
      REQUIRE(ret_struct->x == 42);
      REQUIRE(ret_struct->y == 3.14f);
      REQUIRE(ret_struct->z == 15903);
      byte_offset += sizeof(test_struct);

      std::string str_buf[3];
      memcpy(str_buf, ofd.get_ptr<uint8_t>() + byte_offset, 3*sizeof(std::string));
      REQUIRE(str_buf[0] == "string0");
      REQUIRE(str_buf[1] == "string1");
      REQUIRE(str_buf[2] == "string2");

      byte_offset += sizeof(str_buf);
      REQUIRE(ofd.size() == byte_offset);
    }

    SECTION("test typed gets bounds checking", "[omicsfielddata bounds get]") {
      REQUIRE_THROWS(ofd.get<uint8_t>(ofd.size()));

      // Test to exercise accessing an element when the start index position is not past the end of the array, but
      // there is not enough data between the start and the end for the given type. The previous pushes should leave us
      // in this state, but the first REQUIRE will confirm this.
      size_t test_struct_size = ofd.typed_size<test_struct>();
      REQUIRE(test_struct_size % sizeof(test_struct) != 0);
      REQUIRE_THROWS(ofd.get<test_struct>(ofd.typed_size<test_struct>()));
    }

    SECTION("test typed gets", "[omicsfielddata get]") {
      // typed gets are hard to test if the data isn't aligned properly so for testing align it
      OmicsFieldData aligned_ofd;
      char char_buf[4] = {'a', 'b', 'c', 'd'};
      int int_buf[3] = {0, 1, 2}; 
      aligned_ofd.push_pointer_back(char_buf, 4);
      aligned_ofd.push_pointer_back(int_buf, 3);

      REQUIRE(aligned_ofd.get<char>(0) == 'a');
      REQUIRE(aligned_ofd.get_ptr<char>()[1] == 'b');

      REQUIRE(aligned_ofd.get<int>(1) == 0);
      REQUIRE(aligned_ofd.get_ptr<int>()[2] == 1);
    }

    SECTION("test typed sizes", "[omicsfielddata sizes]") {
      size_t size = sizeof(int) + sizeof(char) + sizeof(test_struct) + 3*sizeof(std::string);
      REQUIRE(ofd.typed_size<int>() == size / sizeof(int));
      REQUIRE(ofd.typed_size<char>() == size / sizeof(char));
      REQUIRE(ofd.typed_size<test_struct>() == size / sizeof(test_struct));
      REQUIRE(ofd.typed_size<std::string>() == size / sizeof(std::string));
    }
  }
}
