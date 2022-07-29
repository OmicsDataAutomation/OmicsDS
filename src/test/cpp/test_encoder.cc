/**
 * @file src/test/cpp/test_encoder.h
 *
 * @section LICENSE
 *
 * The MIT License
 * 
 * @copyright Copyright (c) 2022 Omics Data Automation, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * 
 * @section DESCRIPTION
 *
 * Test gtf ids encoding and decoding
 */

#include <catch2/catch.hpp>

#include "omicsds_encoder.h"


TEST_CASE("test encoder", "[test_encoder]") {
  std::string gtf_id = "ENST00000456328";
  auto encoded = encode_gtf_id(gtf_id);
  CHECK(encoded.first == 456328);
  CHECK(encoded.second == 0);
  CHECK(decode_gtf_id(encoded) == gtf_id);

  gtf_id = "ENST00000456328.111";
  encoded = encode_gtf_id(gtf_id);
  CHECK(encoded.first == 456328);
  CHECK(encoded.second == 111);
  CHECK(decode_gtf_id(encoded) == gtf_id);

  gtf_id = "ENSG00000456328.111";
  encoded = encode_gtf_id(gtf_id);
  CHECK(encoded.first == ((uint64_t)(/*encode_id_type["ENSE"]*/1) << 48 | 456328));
  CHECK(encoded.second == 111);
  CHECK(decode_gtf_id(encoded) == gtf_id);

  gtf_id = "ENSEMU00000456328.111";
  encoded = encode_gtf_id(gtf_id);
  CHECK(encoded.first == ((uint64_t)(/*encode_id_type["ENSE"]*/2) << 48
                          | (uint64_t)(/*encode_kind_of_organism["MU"]*/1) << 56 | 456328));
  CHECK(encoded.second == 111);
  CHECK(decode_gtf_id(encoded) == gtf_id);

  // Test not valid encodings
  encoded = encode_gtf_id("gibberish");
  CHECK(encoded.first == 0);
  CHECK(encoded.second == 0);

  encoded = encode_gtf_id("xxxx00000004444");
  CHECK(encoded.first == 0);
  CHECK(encoded.second == 0);

  encoded.first = 0xFFFFFFFFFFFFFFFF;
  CHECK(decode_gtf_id(encoded) == "");
}
