/**
 * @file   omicsds_encoder.h
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
 * Header file for encoding and decoding utilities for ids from gtf files
 */

#pragma once

#include <map>
#include <stdint.h>
#include <string>
#include <utility>

/** 
 * The gene/transcript ids in the Homo_sapiens.GRCh37.87.gtf seem to follow a pattern. e.g. ENSG00000223972 for gene_id,
 * ENST00000456328 for transcript id, ENSE00002234944 for exon_id, etc. The assumption here is that the pattern is a 
 * sequence of letters followed by 11 digits and a version for ids in the matrix files 
 */

/* The first of the pair points to the encoded gene/transcript id and the second to the version with '0' as the default */
typedef std::pair<uint64_t, uint8_t>  gtf_encoding_t;

/** 
 * Encodes gene/transcript ids, by mapping the id prefixes(ENSG, ENST, ENSE) to numbers and then using a composite of this 
 * number with the rest of the id that is assumed to be a 11 digit numeral. The version is stored separately as the second
 * in a pair, see gtf_encoding_t above.
 */
gtf_encoding_t encode_gtf_id(std::string gtf_id);

/**
 * Decodes the given composite encoded id to its original string format.
 */
std::string decode_gtf_id(const gtf_encoding_t& encoded_id);
