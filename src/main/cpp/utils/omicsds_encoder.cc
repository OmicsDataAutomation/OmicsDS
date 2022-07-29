/**
 * @file   omicsds_encoder.cc
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
 * Implementation for encoding and decoding utilities for ids from gtf files
 */

#include "omicsds_encoder.h"
#include "omicsds_logger.h"

#include <cmath>
#include <regex>

/** 
 * The gene/transcript ids in the Homo_sapiens.GRCh37.87.gtf seem to follow a pattern. e.g. ENSG00000223972 for gene_id,
 * ENST00000456328 for transcript id, ENSE00002234944 for exon_id, etc. The assumption here is that the pattern is a 
 * sequence of letters followed by 11 digits and a version for ids in the matrix files 
 */

// For fast lookup
static std::map<std::string, uint8_t> encode_id_type = {{"ENST", 0}, {"ENSG", 1}, {"ENSE", 2}};
static std::map<uint8_t, std::string> decode_id_type = {{0, "ENST"}, {1, "ENSG"}, {2, "ENSE"}};

static std::map<std::string, uint8_t> encode_kind_of_organism = {{"", 0}, {"MU", 1}};
static std::map<uint8_t, std::string> decode_kind_of_organism = {{0, ""}, {1, "MU"}};

static uint64_t eleven_digit_max = pow(10, 12);

gtf_encoding_t encode_gtf_id(std::string gtf_id) {
  std::regex pattern("^([A-Za-z]{4})([A-Za-z]*)(\\d{11})([.]*)(\\d*)$");
  std::smatch match;
  if (std::regex_match(gtf_id, match, pattern) && match.ready() && !match.empty() && match.size() == 6) {
    // The first sub_match is the whole string
    assert(gtf_id == match[0].str());
    try {
      uint8_t type_of_id = encode_id_type.at(match[1].str());
      uint8_t kind_of_organism = encode_kind_of_organism.at(match[2].str());
      assert(match[3].str().length() == 11);
      uint64_t gtf_id = std::stoll(match[3].str());
      assert(gtf_id < eleven_digit_max);
      uint8_t version = 0;
      if (!match[4].str().empty()) {
        if (match[5].str().length() > 3) {
          logger.error("The version in gtf id {} is not parseable with the current algorithm", gtf_id);
          return {0, 0};
        }
        version = std::stol(match[5].str());
      }
      return { (uint64_t)kind_of_organism << 56 | (uint64_t)type_of_id << 48 | gtf_id, version};
    } catch (const std::exception& ex) {
      logger.error("exception thrown {} : gtf id {} is not parseable", ex.what(), gtf_id);
    }
  } else {
    logger.error("The gtf id {} is not parseable with the current algorithm", gtf_id); 
  }
  return {0, 0};
}

std::string decode_gtf_id(const gtf_encoding_t& encoded_id) {
  std::string gtf_id = "";
  try {
    gtf_id = decode_id_type.at((uint64_t)encoded_id.first >> 48 & 0xFF)
                     + decode_kind_of_organism.at((uint64_t)encoded_id.first >> 56 & 0xFF)
                     + logger.format("{:011}", encoded_id.first & 0xFFFFFFFFFFF);
    if (encoded_id.second) {
      gtf_id += logger.format(".{}", encoded_id.second);
    }
  } catch (const std::exception& ex) {
    logger.error("exception thrown {} : encoded gtf id {} could not be decoded", ex.what(), encoded_id.first);
  }
  return gtf_id;
}
