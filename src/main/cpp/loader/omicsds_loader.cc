/**
 * src/main/cpp/loader/omicsds_loader.cc
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
 * Implementation for a generic SAM reader
 */

#include "omicsds_loader.h"
#include "omicsds_export.h"

std::vector<std::string> split(std::string str, std::string sep) {
  std::vector<std::string> retval;
  int index;

  if(str.length() >= 2) {
    if(str[0] == '[') {
       str = str.substr(1, str.length() - 2);
    }
  }

  while((index = str.find(sep)) != std::string::npos) {
    retval.push_back(str.substr(0, index));
    str.erase(0, index + 1);
  }
  retval.push_back(str);

  return retval;
}

bool FileUtility::generalized_getline(std::string& retval) {
  retval = "";

  while(chars_read < file_size || str_buffer.size()) {
    int idx = str_buffer.find('\n');
    if(idx != std::string::npos) {
      retval = retval + str_buffer.substr(0, idx); // exclude newline
      str_buffer.erase(0, idx + 1); // erase newline
      return true;
    }

    retval = retval + str_buffer;
    str_buffer.clear();

    int chars_to_read = std::min<ssize_t>(buffer_size, file_size - chars_read);

    if(chars_to_read) {
      TileDBUtils::read_file(filename, chars_read, buffer, chars_to_read);
      chars_read += chars_to_read;
    }

    str_buffer.insert(str_buffer.end(), buffer, buffer + chars_to_read);
  }

  return false;
}

int FileUtility::read_file(void* buffer, size_t chars_to_read) {
  size_t buf_position = 0;
  if(str_buffer.size() > 0) {
    buf_position = read_from_str_buffer(buffer, chars_to_read);
    chars_to_read -= buf_position;
  }
  int rcode = TileDBUtils::read_file(filename, chars_read, (char*)buffer + buf_position, chars_to_read);
  chars_read += chars_to_read;
  CHECK_RC(rcode);
  return rcode;
}

size_t FileUtility::read_from_str_buffer(void* buffer, size_t chars_to_read) {
  size_t readable_chars = std::min<size_t>(str_buffer.size(), chars_to_read);
  memcpy(buffer, str_buffer.c_str(), readable_chars);
  str_buffer.erase(0, readable_chars);
  return readable_chars;
}

void read_sam_file(std::string filename) {
  std::cerr << "SAM file is " << filename << std::endl;

  samFile *fp_in = hts_open(filename.c_str(),"r"); //open bam file
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
  bam1_t *aln = bam_init1(); //initialize an alignment

  if(!bamHdr) {
    std::cerr << "header is null" << std::endl;
  } else {
    std::cerr << "header is NOT null" << std::endl;
  }
  
  // header parse
  // uint32_t *tar = bamHdr->text ;
  // uint32_t *tarlen = bamHdr->target_len ;

  // printf("%d\n",tar);
  
  int rc;
  std::cerr << "before while" << std::endl;
  while(!(rc = sam_read1(fp_in,bamHdr,aln))){
          
    int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordinate (+1)
    char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    uint32_t len = aln->core.l_qseq; //length of the read.
    
    uint8_t *q = bam_get_seq(aln); //quality string
    uint32_t q2 = aln->core.qual ; //mapping quality
    
    char* qname = bam_get_qname(aln);    
    uint16_t flag = aln->core.flag;
    uint32_t* cigar = bam_get_cigar(aln);
    uint32_t n_cigar = aln->core.n_cigar;
    char cigar_codes[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
    //uint8_t* qual = bam_get_qual(aln);
    char* qual = (char*)bam_get_qual(aln);
    uint8_t mapq = aln->core.qual;
    //char* seq = (char*)bam_get_seq(aln);
    int32_t rnext = aln->core.mtid;
    int32_t pnext = aln->core.mpos;
    int32_t tlen = aln->core.isize;

    char *qseq = (char *)malloc(len);
    
    for(int i=0; i<len; i++){
      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
    }
    
    //printf("chr=%s\tpos=%d\tlen=%d\tqseq=%s\tq=%s\tq2=%d\n",chr,pos,len,qseq,q,q2);
    printf("qname=%s\tflag=%d\tchr=%s\tpos=%d\tmapq=%d\tlen=%d\tqseq=%s\tq2=%d\n",qname,flag,chr,pos,mapq,len,qseq,q2);
    std::cout << "cigar=";
    for(uint32_t i = 0; i < n_cigar; i++) {
      auto op_len = cigar[i] >> 4;
      auto op = cigar_codes[cigar[i] & 0b1111];

      std::cout << op_len << op << ", ";
    }
    std::cout << std::endl;
    std::cout << "rnext=" << rnext << "\tpnext=" << pnext << "\ttlen=" << tlen << std::endl;
    std::cerr << "after while rc is " << rc << std::endl;

    std::cout << "qual=" << std::endl;
    for(uint32_t i = 0; i < len; i++) {
      double q = qual[i];
      std::cout << pow(10, (q/-10)) << ", ";
    }
    std::cout << std::endl << std::endl;

    free(qseq);
  }
  bam_destroy1(aln);
  sam_close(fp_in);
}

GenomicMap::GenomicMap(const std::string& mapping_file): GenomicMap(std::make_shared<FileUtility>(mapping_file)) { }

GenomicMap::GenomicMap(std::shared_ptr<FileUtility> mapping_reader): m_mapping_reader(mapping_reader) {
  std::string str;
  int line_num = -1;
  while(m_mapping_reader->generalized_getline(str)) {
    line_num++;

    auto toks = split(str, "\t");
    if(toks.size() < 3) {
      std::cerr << "Warning: line " << line_num << " of mapping file " << m_mapping_reader->filename << " does not contain enough fields (at least 3), make sure the file is tab separated" << std::endl;
      std::cerr << "Note: mapping file may be subset of serialized schema" << std::endl;
      continue;
    }

    std::string contig_name;
    uint64_t length;
    uint64_t starting_index;

    try {
      contig_name = toks[0];
      length = std::stol(toks[1]);
      starting_index = std::stol(toks[2]);
    }
    catch(...) {
      std::cerr << "Warning: line " << line_num << " of mapping file " << m_mapping_reader->filename << " could not be parsed (2nd or 3rd field was not a valid uint64_t)" << std::endl;
      std::cerr << "Note: mapping file may be subset of serialized schema" << std::endl;
      continue;
    }

    contigs.emplace_back(contig_name, length, starting_index);
  }

  idxs_name.resize(contigs.size());
  std::iota(idxs_name.begin(), idxs_name.end(), 0);
  std::sort(idxs_name.begin(), idxs_name.end(), [&](auto l, auto r) { return contigs[l].name < contigs[r].name; });

  idxs_position.resize(contigs.size());
  std::iota(idxs_position.begin(), idxs_position.end(), 0);
  std::sort(idxs_position.begin(), idxs_position.end(), [&](auto l, auto r) { return contigs[l].starting_index < contigs[r].starting_index; });

  std::cout << "REMOVE" << std::endl;
  for(auto& c : contigs) {
    std::cout << "(" << c.name << ", " << c.length << ", " << c.starting_index << ")" << std::endl;
  }
  std::cout << std::endl << std::endl;
}

uint64_t GenomicMap::flatten(std::string contig_name, uint64_t offset) {
  auto it = std::lower_bound(idxs_name.begin(), idxs_name.end(), contig_name, [&](auto l, auto r) { return contigs[l].name < r; });

  //int idx = std::distance(idxs_name.begin(), it);

  if(it != idxs_name.end() && contigs[*it].name == contig_name) {
    if(offset < contigs[*it].length) {
      return contigs[*it].starting_index + offset;
    }
    else {
      std::cerr << "Error, contig " << contig_name << " is only length " << contigs[*it].length << ", " << offset << "is out of bounds" << std::endl;
      exit(1);
    }
  }
  else {
    std::cerr << "Error, contig " << contig_name << " not found in mapping file " << m_mapping_reader->filename << std::endl;
    exit(1);
  }
}

bool equivalent_schema(const OmicsSchema& l, const OmicsSchema& r) {
  if(l.attributes.size() != r.attributes.size()) return false;

  for(auto li = l.attributes.begin(), ri = r.attributes.begin(); li != l.attributes.end(), ri != r.attributes.end(); li++, ri++) {
    if(li->first != ri->first) return false;
    if(li->second.type != ri->second.type) return false;
  }

  return true;
}

void GenomicMap::serialize(std::string path) {
  for(auto& c : contigs) {
    c.serialize(path);
  }
}

void OmicsSchema::serialize(std::string path) {
  if(TileDBUtils::is_file(path)) {
    TileDBUtils::delete_file(path);
  }

  auto write = [&](std::string str) {
    return FileUtility::write_file(path, str);
  };

  write("v1\n"); // version
  std::string order_str = (order == POSITION_MAJOR)? "POSITION_MAJOR\n" : "SAMPLE_MAJOR\n";
  write(order_str); // order
  std::string num_attributes_str = std::to_string(attributes.size()) + "\tattributes\n";
  write(num_attributes_str);
  // attributes
  for(auto& a : attributes) {
    write(a.first); // attribute name
    write("\t");
    write(a.second.type_to_string()); // attribute type as string
    write("\t");
    write(a.second.length_to_string()); // attribute length as string
    write("\n");
  }

  genomic_map.serialize(path);
}

bool OmicsSchema::create_from_file(const std::string& filename) {
  if(!TileDBUtils::is_file(filename)) {
    std::cerr << "Error: cannot deserialize OmicsSchema, " << filename << " does not exist" << std::endl;
    return false;
  }

  std::shared_ptr<FileUtility> reader = std::make_shared<FileUtility>(filename);
  std::string str;

  // version
  if(!reader->generalized_getline(str)) {
    std::cerr << "Error: cannot deserialize OmicsSchema, " << filename << " is empty" << std::endl;
    return false;
  }
  if(str != "v1") {
    std::cerr << "Note: while deserializing OmicsSchema encountered version " << str << ", only v1 is supported" << std::endl;
  }

  // order
  if(!reader->generalized_getline(str)) {
    std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", issue reading order" << std::endl;
    return false;
  }
  if(str == "POSITION_MAJOR") {
    order = POSITION_MAJOR;
  }
  else if(str == "SAMPLE_MAJOR") {
    order = SAMPLE_MAJOR;
  }
  else {
    std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", issue reading order" << std::endl;
    return false;
  }

  if(!reader->generalized_getline(str)) {
    std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", issue reading number of attributes" << std::endl;
    return false;
  }
  int num_attributes = -1;
  try {
    num_attributes = std::stoi(split(str, "\t")[0]);
    if(num_attributes < 0) { throw std::runtime_error("number of attributes is negative"); }
  }
  catch (...) {
    std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", issue reading number of attributes" << std::endl;
    return false;
  }

  attributes.clear();
  for(int i = 0; i < num_attributes; i++) {
    if(!reader->generalized_getline(str)) {
      std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", fewer attributes than reported (" << i << " of " << num_attributes << ")" << std::endl;
      return false;
    }
    auto tokens = split(str, "\t");
    if(tokens.size() < 3) {
      std::cerr << "Error: cannot deserialize OmicsSchema from " << filename << ", issue reading attribute " << i << std::endl;
      return false;
    }
    int length;
    try {
      length = std::stoi(tokens[2]);
    }
    catch(...) {
      length = -1;
    }

    attributes.emplace(tokens[0], OmicsFieldInfo(tokens[1], length));
  }

  genomic_map = GenomicMap(reader);
  return true;
}

// FIXME could potentially index map instead of doing this, only an issue if very many attributes
int OmicsSchema::index_of_attribute(const std::string& name) {
  auto it = attributes.find(name);  
  if(it == attributes.end()) {
    return -1;
  }

  return std::distance(attributes.begin(), it);
}

OmicsCell OmicsCell::create_invalid_cell() {
  OmicsCell rv;
  rv.coords = {-1, -1};
  return rv;
}

bool OmicsCell::is_invalid_cell(const OmicsCell& cell) {
  return cell.coords[0] < 0 || cell.coords[1] < 0;
}

SamReader::SamReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, int file_idx) : OmicsFileReader(filename, schema, sample_map, file_idx) {
  m_fp = hts_open(filename.c_str(),"r"); //open bam file
  m_hdr = sam_hdr_read(m_fp); //read header
  m_align = bam_init1(); //initialize an alignment

  if(!m_hdr) {
    std::cout << "SamReader header is null" << std::endl;
  }

  auto toks = split(filename, "/");
  std::string sample_name = "";
  if(toks.size()) {
    sample_name = toks.back();
  }
  else {
    std::cerr << "Error, could not deduce file name for file at path " << filename << std::endl;
  }
  if(m_sample_map->count(sample_name)) {
    m_row_idx = (*m_sample_map)[sample_name];
  }
  else {
    std::cerr << "Error, no entry " << sample_name << " in sample map" << std::endl;
    exit(1);
  }

  assert((bool)schema);
}

SamReader::~SamReader() {
  //std::cout << "REMOVE SamReader::~SamReader" << std::endl;
  //std::cout << "FIXME uncomment cleanup functions" << std::endl;
  bam_destroy1(m_align);
  sam_close(m_fp);
}

std::vector<OmicsCell> SamReader::get_next_cells() {
  std::vector<OmicsCell> cells;

  int rc;
  if(!(rc = sam_read1(m_fp,m_hdr,m_align))){
    int32_t pos = m_align->core.pos +1; //left most position of alignment in zero based coordinate (+1)
    char *chr = m_hdr->target_name[m_align->core.tid] ; //contig name (chromosome)
    uint32_t len = m_align->core.l_qseq; //length of the read.

    uint8_t *q = bam_get_seq(m_align); //quality string
    uint32_t q2 = m_align->core.qual ; //mapping quality

    char* qname = bam_get_qname(m_align);
    uint16_t flag = m_align->core.flag;
    uint32_t* cigar = bam_get_cigar(m_align);
    uint32_t n_cigar = m_align->core.n_cigar;
    char cigar_codes[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
    //uint8_t* qual = bam_get_qual(m_align);
    
    char* qual = (char*)bam_get_qual(m_align);
    uint8_t mapq = m_align->core.qual;
    //char* seq = (char*)bam_get_seq(m_align);
    int32_t rnext = m_align->core.mtid;
    int32_t pnext = m_align->core.mpos;
    int32_t tlen = m_align->core.isize;
    std::vector<char> qseq(len);

    for(int i=0; i< len ; i++){
      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
    }

    std::string sample = get_filename();

    int64_t position = m_schema->genomic_map.flatten(chr, pos);

    std::cerr << "\t\t\t\tREMOVE rname len " << std::strlen(chr) << std::endl;
    std::cerr << "\t\t\t\tREMOVE cigar len " << n_cigar << std::endl;

    OmicsCell cell({ (int64_t)m_row_idx, position}, m_schema, m_file_idx);
    cell.add_field_ptr("QNAME", qname, std::strlen(qname));
    cell.add_field("FLAG", flag);
    cell.add_field_ptr("RNAME", chr, std::strlen(chr));
    cell.add_field("POS", pos);
    cell.add_field("MAPQ", mapq);
    cell.add_field_ptr("CIGAR", cigar, n_cigar);
    cell.add_field("RNEXT", rnext);
    cell.add_field("PNEXT", pnext);
    cell.add_field("TLEN", tlen);
    cell.add_field_ptr("SEQ", qseq.data(), (int)qseq.size());
    cell.add_field_ptr("QUAL", qual, std::strlen(qual));
  
    // FIXME REMOVE
    cell.add_field_ptr("SAMPLE_NAME", (char*)sample.c_str(), (int)sample.length());


    OmicsCell end_cell = cell;
    end_cell.file_idx = -1;
    int end_offset = std::abs(tlen) - 1; // FIXME figure out negative template length
    end_cell.coords[1] += end_offset;
    if(tlen) {  // if end cell is in same position, only create one cell
      std::cout << "REMOVE SamReader::get_next_cells return " << cell.coords[1] << ", " << end_cell.coords[1] << std::endl;
      return { cell, end_cell };
    }
    std::cout << "REMOVE SamReader::get_next_cells return " << cell.coords[1] << std::endl;
    return { cell };
  }
  else {
    std::cout << "REMOVE sam_read1 was " << rc << std::endl;
    return { OmicsCell::create_invalid_cell() };
  }
}

BedReader::BedReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, int file_idx): OmicsFileReader(filename, schema, sample_map, file_idx) {
  std::string line;
  if(!m_reader_util->generalized_getline(line)) {
    std::cerr << "Error file " << filename << " is empty, skipping" << std::endl;
    return;
  }

  std::smatch m;
  std::regex exp("description\\s*=\\s*\"(.*?)\"");
  std::regex_search(line, m, exp);

  if(m.size() < 2) {
    std::cerr << "Error, could not find sample name for file " << filename << std::endl;
    std::cerr << "First line of Bed file should contain 'description = \"[sample_name]\"'" << std::endl;
    exit(1);
  }

  m_sample_name = m[1];
  std::cout << "************************************** BED sample is " << m_sample_name << std::endl;

  if(!m_sample_map->count(m_sample_name)) {
    std::cerr << "Error, sample " << m_sample_name << " from file " << filename << " not found in sample map" << std::endl;
    exit(1);
  }

  m_row_idx = (*m_sample_map)[m_sample_name];
}

std::vector<OmicsCell> BedReader::get_next_cells() {
  std::string line;
  while(m_reader_util->generalized_getline(line)) {
    std::stringstream ss(line);
    std::string str;
    std::vector<std::string> fields;
    while(ss >> str) {
      fields.push_back(str);
    }
    if(fields.size() < 5) { continue; }

    std::string chrom = fields[0], gene = "N/A", name = fields[3];
    uint64_t start, end, flattened_start, flattened_end;
    float score;
    
    try {
      start = std::stoul(fields[1]);
      flattened_start = m_schema->genomic_map.flatten(chrom, start);
      end = std::stoul(fields[2]);
      flattened_end = m_schema->genomic_map.flatten(chrom, end);
      score = std::stof(fields[4]);
    }
    catch(...) {
      continue;
    }

    OmicsCell cell({ (int64_t)m_row_idx, (int64_t)flattened_start }, m_schema, m_file_idx);
    cell.add_field_ptr("CHROM", chrom.c_str(), chrom.length());
    cell.add_field("START", start);
    cell.add_field("END", end);
    cell.add_field("SCORE", score);
    cell.add_field_ptr("GENE", gene.c_str(), gene.length());
    cell.add_field_ptr("SAMPLE_NAME", m_sample_name.c_str(), m_sample_name.length());
    cell.add_field_ptr("NAME", name.c_str(), name.length());

    OmicsCell end_cell = cell;
    end_cell.file_idx = -1;
    end_cell.coords[1] = flattened_end;

    if(flattened_start == flattened_end) {
      return { cell };
    }
    return { cell, end_cell };
  }
  return { };
}

MatrixReader::MatrixReader(std::string filename, std::shared_ptr<OmicsSchema> schema, std::shared_ptr<SampleMap> sample_map, std::shared_ptr<GeneIdMap> gene_id_map, int file_idx): OmicsFileReader(filename, schema, sample_map, file_idx), m_gene_id_map(gene_id_map) {
  std::string line;

  if(!m_reader_util->generalized_getline(line)) {
    std::cerr << "Error, file " << filename << " is empty" << std::endl;
    exit(1);
  }

  auto toks = split(line, "\t"); // TODO maybe support commas as well
  if(toks[0] == "GENE") {
    m_position_major = false;
  }
  else if(toks[0] == "SAMPLE") {
    m_position_major = true;
  }
  else {
    std::cerr << "Error reading file " << filename << std::endl;
    std::cerr << "Matrix files should be tab separated in format: " << std::endl;
    std::cerr << "GENE\t[gene1]\t[gene2]" <<  std::endl;
    std::cerr << "[sample1]\t[score1]\t[score2]" << std::endl;
    std::cerr << "Or the transpose (keyword SAMPLE instead of GENE, genes across rows)" << std::endl;
  }
  if(m_position_major != m_schema->position_major()) {
    std::cerr << "Error order of matrix file " << filename << " does not match that of schema (must both be position or sample major for now, WIP)" << std::endl;
    exit(1);
  }
  m_columns = std::vector<std::string>(toks.begin() + 1, toks.end());
  m_row_scores = std::vector<float>(m_columns.size(), 0);
  m_column_idx = m_columns.size(); // to force parsing next line
}

bool MatrixReader::parse_next(std::string& sample, std::string& gene, float& score) {
  while(m_column_idx >= m_row_scores.size()) {
    std::string line;
    if(!m_reader_util->generalized_getline(line)) {
      return false;
    }

    auto toks = split(line, "\t");
    if(toks.size() != m_columns.size() + 1) { // sample/gene id followed by scores
      continue;
    }

    m_current_token = toks[0];

    try {
      for(int i = 0; i < m_columns.size(); i++) {
        m_row_scores[i] = std::stof(toks[i + 1]);
      }
    }
    catch (...) {
      continue;
    }
    m_column_idx = 0;
    break;
  }

  if(m_position_major) {
    sample = m_columns[m_column_idx];
    gene = m_current_token;
  }
  else {
    sample = m_current_token;
    gene = m_columns[m_column_idx];
  }
  score = m_row_scores[m_column_idx];
  m_column_idx++;
  return true;
}

std::vector<OmicsCell> MatrixReader::get_next_cells() {
  std::string sample_name, gene_name; // TODO avoid repeated lookups for row token
  float score;
  uint64_t row_idx;

  while(parse_next(sample_name, gene_name, score)) {
    if(!m_sample_map->count(sample_name) || !m_gene_id_map->count(gene_name)) { continue; }
    
    row_idx = (*m_sample_map)[sample_name];
    GeneIdMap::Gene& gene = (*m_gene_id_map)[gene_name];

    std::string name = "N/A";

    OmicsCell cell({ (int64_t)row_idx, (int64_t)gene.flattened_start }, m_schema, m_file_idx);
    cell.add_field_ptr("CHROM", gene.chrom.c_str(), gene.chrom.length());
    cell.add_field("START", gene.start);
    cell.add_field("END", gene.end);
    cell.add_field("SCORE", score);
    cell.add_field_ptr("GENE", gene_name.c_str(), gene_name.length());
    cell.add_field_ptr("SAMPLE_NAME", sample_name.c_str(), sample_name.length());
    cell.add_field_ptr("NAME", name.c_str(), name.length());

    OmicsCell end_cell = cell;
    end_cell.file_idx = -1;
    end_cell.coords[1] = gene.flattened_end;

    if(gene.flattened_start == gene.flattened_end) {
      return { cell };
    }
    return { cell, end_cell };
  }

  return {};
}

int OmicsModule::tiledb_create_array(const std::string& workspace, const std::string& array_name, const OmicsSchema& schema) {
  // initialize with the default configuration parameters
  TileDB_CTX* tiledb_ctx;
  CHECK_RC(tiledb_ctx_init(&tiledb_ctx, NULL));

  // Create a workspace
  CHECK_RC(tiledb_workspace_create(tiledb_ctx, workspace.c_str()));

  std::string full_name = workspace + "/" + array_name;

  // Prepare parameters for array schema
  std::vector<const char*> attributes_vec;
  std::vector<int32_t> cell_val_num_vec;
  std::vector<int32_t> types_vec;
  for(auto&p : schema.attributes) {
    std::cerr << "REMOVE aname is " << p.first << std::endl;
    std::cerr << "REMOVE c_str is " << p.first.c_str() << std::endl;
    std::cerr << "REMOVE address is " << (int64_t)p.first.c_str() << std::endl;
    attributes_vec.push_back(p.first.c_str());
    cell_val_num_vec.push_back(p.second.length);
    types_vec.push_back(p.second.tiledb_type());
    std::cerr << "REMOVE vec 1 " << container_to_string(attributes_vec) << std::endl;
  }
  types_vec.push_back(TILEDB_INT64); // coords

  std::cerr << "REMOVE vec 2 " << container_to_string(attributes_vec) << std::endl;
  std::cerr << "REMOVE first attribute address is " << (int64_t)(attributes_vec[0]) << std::endl;
  std::cerr << "REMOVE first attribute is " << attributes_vec[0] << std::endl;
  std::cerr << "REMOVE len is " << strlen(attributes_vec[0]) << std::endl;

  const char** attributes = attributes_vec.data();
  const int* cell_val_num = cell_val_num_vec.data();
  const int* types = types_vec.data();

  int32_t order = TILEDB_ROW_MAJOR; // different orders are implemented by reordering coordinates

  const char* dimensions[3];
  dimensions[2] = "LEVEL";
  if(schema.position_major()) {
    dimensions[0] = "POSITION";
    dimensions[1] = "SAMPLE";
  }
  else {
    dimensions[0] = "SAMPLE";
    dimensions[1] = "POSITION";
  }

  int64_t domain[] =
  {
      0, std::numeric_limits<int64_t>::max(),            // 1st dimension limits (SAMPLE or POSITION based on order)
      0, std::numeric_limits<int64_t>::max(),            // 2nd dimension limits
      0, std::numeric_limits<int64_t>::max()             // LaVEL limits
  };

  std::vector<int32_t> compression_vec(schema.attributes.size() + 1, TILEDB_NO_COMPRESSION); // plus 1 for coordinates
  const int* compression = compression_vec.data();

  std::vector<int32_t> offsets_compression_vec(schema.attributes.size(), TILEDB_NO_COMPRESSION);
  const int* offsets_compression = offsets_compression_vec.data();

  int64_t tile_extents[] =
  {
      1,                                                 // 1st dimension extents
      1,                                                 // 2nd dimension  extents
      1                                                  // LEVEL extents
  };

  // Set array schema
  TileDB_ArraySchema array_schema;
  tiledb_array_set_schema(
      &array_schema,                                     // Array schema struct
      full_name.c_str(),                                 // Array name
      attributes,                                        // Attributes
      schema.attributes.size(),                          // Number of attributes
      std::numeric_limits<int64_t>::max(),               // Capacity
      //50,                                                // Capacity
      order,                                             // Cell order
      cell_val_num,                                      // Number of cell values per attribute
      compression,                                       // Compression
      NULL,                                              // Compression level - Use defaults
      offsets_compression,                               // Offsets compression
      NULL,                                              // Offsets compression level
      0,                                                 // Sparse array
      dimensions,                                        // Dimensions
      3,                                                 // Number of dimensions
      domain,                                            // Domain
      6*sizeof(int64_t),                                 // Domain length in bytes
      tile_extents,                                      // Tile extents
      4*sizeof(int64_t),                                 // Tile extents length in bytes
      order,                                             // Tile order
      types                                              // Types
  );

  // Create array
  CHECK_RC(tiledb_array_create(tiledb_ctx, &array_schema));

  // Free array schema
  CHECK_RC(tiledb_array_free_schema(&array_schema));

  /* Finalize context. */
  CHECK_RC(tiledb_ctx_finalize(tiledb_ctx));

  return 0;
}

int OmicsModule::tiledb_open_array(const std::string& workspace, const std::string& array_name, bool write) {
  CHECK_RC(tiledb_ctx_init(&m_tiledb_ctx, NULL));

  std::string path = workspace + "/" + array_name;

  // Initialize array
  CHECK_RC(tiledb_array_init(
      m_tiledb_ctx,                                      // Context
      &m_tiledb_array,                                   // Array object
      path.c_str(),                                      // Array name
      write ? TILEDB_ARRAY_WRITE : TILEDB_ARRAY_READ,    // Mode
      NULL,                                              // Entire domain
      NULL,                                              // All attributes
      0));                                               // Number of attributes

  return 0;
}

int OmicsModule::tiledb_close_array() {
  // Finalize array
  if(m_tiledb_array) CHECK_RC(tiledb_array_finalize(m_tiledb_array));
  m_tiledb_array = 0;

  // Finalize context
  if(m_tiledb_ctx) CHECK_RC(tiledb_ctx_finalize(m_tiledb_ctx));
  m_tiledb_ctx = 0;

  return 0;
}

int OmicsLoader::tiledb_write_buffers() {
  std::vector<void*> buffers_vec;
  std::vector<size_t> buffer_sizes_vec;

  int i = 0;
  for(auto it = m_schema->attributes.begin(); it != m_schema->attributes.end(); it++, i++) {
    buffers_vec.push_back(offset_buffers[i].data());
    buffer_sizes_vec.push_back(offset_buffers[i].size());

    if(it->second.length == TILEDB_VAR_NUM) {
      buffers_vec.push_back(var_buffers[i].data());
      buffer_sizes_vec.push_back(var_buffers[i].size());
    }
  }

  buffers_vec.push_back(coords_buffer.data());
  buffer_sizes_vec.push_back(sizeof(int64_t) * coords_buffer.size());  

  // Write to array
  CHECK_RC(tiledb_array_write(m_tiledb_array, const_cast<const void**>(buffers_vec.data()), buffer_sizes_vec.data()));

  return 0;
}

void OmicsLoader::buffer_cell(const OmicsCell& cell, int level) {
  auto fiter = cell.fields.begin();
  auto aiter = m_schema->attributes.begin();
  int i = 0;
  for(; fiter != cell.fields.end() && aiter != m_schema->attributes.end(); fiter++, aiter++, i++) {
    auto& data = fiter->data;

    int length = aiter->second.length;
    int size = aiter->second.element_size();
    if(length == TILEDB_VAR_NUM) { // variable length
      OmicsFieldData::push_back<size_t>(offset_buffers[i], attribute_offsets[i]);
      attribute_offsets[i] += data.size();
      for(auto& c : data) {
        var_buffers[i].push_back(c);
      }
    }
    else {
      assert(data.size() == size * length);
      for(auto& c : data) {
        offset_buffers[i].push_back(c);
      }
    }
  }
  coords_buffer.push_back(cell.coords[0]);
  coords_buffer.push_back(cell.coords[1]);
  coords_buffer.push_back(level);
}

void ReadCountLoader::create_schema() {
  m_schema->attributes.emplace("SAMPLE_NAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("QNAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("FLAG", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint16_t, 1));
  m_schema->attributes.emplace("RNAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("POS", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_int32_t, 1));
  m_schema->attributes.emplace("MAPQ", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint8_t, 1));
  m_schema->attributes.emplace("CIGAR", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint32_t, -1));
  m_schema->attributes.emplace("RNEXT", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_int32_t, 1));
  m_schema->attributes.emplace("PNEXT", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_int32_t, 1));
  m_schema->attributes.emplace("TLEN", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_int32_t, 1));
  m_schema->attributes.emplace("SEQ", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("QUAL", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
}

void ReadCountLoader::add_reader(const std::string& filename) {
  if(std::regex_match(filename, std::regex("(.*)(sam)($)"))) {
    m_files.push_back(std::make_shared<SamReader>(filename, m_schema, m_sample_map, m_files.size()));
  }
}

void TranscriptomicsLoader::create_schema() {
  m_schema->attributes.emplace("CHROM", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("START", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint64_t, 1));
  m_schema->attributes.emplace("END", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint64_t, 1));
  m_schema->attributes.emplace("SCORE", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_float_t, 1));
  m_schema->attributes.emplace("GENE", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1)); // Will be N/A for bed files
  m_schema->attributes.emplace("SAMPLE_NAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1));
  m_schema->attributes.emplace("NAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char, -1)); // Field in bed files, will be N/A for matrix files
}

void TranscriptomicsLoader::add_reader(const std::string& filename) {
  if(std::regex_match(filename, std::regex("(.*)(bed)($)"))) {
    m_files.push_back(std::make_shared<BedReader>(filename, m_schema, m_sample_map, m_files.size()));
  }
  else if(std::regex_match(filename, std::regex("(.*)(resort)($)"))) {
    m_files.push_back(std::make_shared<MatrixReader>(filename, m_schema, m_sample_map, m_gene_id_map, m_files.size()));
  }
}

SampleMap::SampleMap(const std::string& sample_map) {
  FileUtility file(sample_map);

  std::string line;
  while(file.generalized_getline(line)) {
    auto toks = split(line, "\t");
    if(toks.size() < 2) continue;
    try {
      std::string name = toks[0];
      size_t idx = std::stoul(toks[1]);
      if(!map.count(name)) {
        map[name] = idx;
      }
    }
    catch(...) {
      continue;
    }
  }
}

GeneIdMap::GeneIdMap(const std::string& gene_map, std::shared_ptr<OmicsSchema> schema, bool use_transcript, bool drop_version): schema(schema) {
  auto toks = split(gene_map, ".");
  if(!toks.size()) {
    std::cerr << "Error creating gene map, cannot deduce file format from filename " << gene_map << std::endl;
    exit(1);
  }

  std::string format = toks.back();

  if(format == "gtf" || format == "gi") {
    create_from_gtf(gene_map, use_transcript, drop_version);
    return;
  }
  if(format == "gi") {
    create_from_gi(gene_map);
    return;
  }
  if(format == "gbed") {
    create_from_gbed(gene_map);
    return;
  }
  std::cerr << "Error creating gene map, unsupported file format for filename " << gene_map << " (expected gtf/gff or gi)" << std::endl;
  exit(1);
}

void GeneIdMap::create_from_gtf(const std::string& gene_map, bool use_transcript, bool drop_version) {
  std::string str;

  int ind = -1;

  FileUtility file(gene_map);

  while(file.generalized_getline(str)) {
    if(str[0] == '#') {
      continue;
    }

    ++ind;

    std::stringstream ss(str);
    std::string field;
    std::vector<std::string> fields;

    // read the 8 fields before attribute
    for(int i = 0; i < 8; i++) {
      ss >> field;
      fields.push_back(field);
    }

    if(fields[2] != "transcript") {
      continue;
    }

    std::string chrom = fields[0];
    uint64_t start, end, flattened_start, flattened_end;

    try {
      start = std::stoul(fields[3]);
      end = std::stoul(fields[4]);
      flattened_start = schema->genomic_map.flatten(chrom, start);
      flattened_end = schema->genomic_map.flatten(chrom, end);
    }
    catch (...) {
      continue;
    }

    std::string attributes;
    std::getline(ss, attributes);

    std::string pattern = use_transcript ? "transcript_id" : "gene_id"; // FIXME use regex in case strangeness with whitespace

    std::smatch m;
    std::regex exp(pattern + "\\s*\"(.*?)\"");
    std::regex_search(str, m, exp);

    if(m.size() < 2) {
      continue;
    }
 
    std::string tid = m[1];

    if(drop_version) {
      tid = split(tid, ".")[0];
    }

    if(map.count(tid)) {
      std::cout << "gtf reading: " << tid << " is duplicate, old start/end: " << map[tid].flattened_start << ", " << map[tid].flattened_end << ", new start/end: " << flattened_start << ", " << flattened_end << std::endl;
    }

    map[tid] = Gene(chrom, start, end, flattened_start, flattened_end);
  }
}

void GeneIdMap::create_from_gi(const std::string& gene_map) {
  FileUtility file(gene_map);

  char version;
  file.read_file(&version, 1);
  if(version != 1) {
    std::cout << "version " << (int)version << " not supported" << std::endl;
  }

  int64_t size = 0;
  file.read_file(((char*)&size), 5);

  for(int i = 0; i < size; i++) {
    uint16_t len; // gene name
    file.read_file((char*)&len, 2);
    char name[len];
    file.read_file(name, len);
    uint16_t len2; // chrom
    file.read_file((char*)&len2, 2);
    char chrom[len2];
    file.read_file(chrom, len2);
    uint64_t start = 0;
    uint64_t end = 0;
    file.read_file((char*)&start, 5);
    file.read_file((char*)&end, 5);
    map[std::string(name, len)] = Gene(std::string(chrom, len2), start, end, schema->genomic_map.flatten(chrom, start), schema->genomic_map.flatten(chrom, end));
  }
}

void GeneIdMap::create_from_gbed(const std::string& gene_map) {
  FileUtility file(gene_map);

  std::string str;
  while(file.generalized_getline(str)) {
    auto toks = split(str, "\t");
    if(toks.size() < 4) { continue; }
    std::string chrom = toks[1];
    uint64_t start, end;
    try {
      start = std::stoul(toks[2]);
      end = std::stoul(toks[3]);
    }
    catch(...) {
      continue;
    }
    map[toks[0]] = Gene(chrom, start, end, schema->genomic_map.flatten(chrom, start), schema->genomic_map.flatten(chrom, end));
  }
}

void GeneIdMap::export_as_gi(const std::string& filename) {
  auto write = [&] (const void* buffer, size_t length) {
    return FileUtility::write_file(filename, buffer, length);
  };

  char version = 1;
  write(&version, 1);

  // write 5B number of entries
  int64_t size = map.size();
  write((char*)&size, 5);

  // write 2B name length, name, 2B chrom length, chrom, 5B start, 5B end
  for(auto& a : map) {
    std::string name = a.first;
    int16_t len = name.length();
    uint64_t start = a.second.start;
    uint64_t end = a.second.end;
    std::string chrom = a.second.chrom;
    int16_t len2 = chrom.length();

    write((char*)&len, 2);
    write(name.c_str(), len);
    write((char*)&len2, 2);
    write(chrom.c_str(), len2);
    write((char*)&start, 5);
    write((char*)&end, 5);
  }
}

// FIXME add parallelism of some kind
// FIXME 
OmicsLoader::OmicsLoader(
                         const std::string& workspace,
                         const std::string& array,
                         const std::string& file_list,
                         const std::string& sample_map,
                         const std::string& mapping_file,
                         bool position_major
                        ): OmicsModule(workspace, array, mapping_file, position_major), m_file_list(file_list), m_sample_map(std::make_shared<SampleMap>(sample_map)), m_pq(comparitor) {}

void OmicsLoader::initialize() { // FIXME move file reader creation to somewhere virtual
  create_schema();

  // set up buffers
  offset_buffers = std::vector<std::vector<uint8_t>>(m_schema->attributes.size());
  var_buffers = std::vector<std::vector<uint8_t>>(m_schema->attributes.size());
  coords_buffer.clear();
  attribute_offsets = std::vector<size_t>(m_schema->attributes.size(), 0);

  // create workspace/array
  tiledb_create_array(m_workspace, m_array, *m_schema);
  tiledb_open_array(m_workspace, m_array);

  serialize_schema();

  // add file readers
  FileUtility list(m_file_list);
  std::string filename;
  while(list.generalized_getline(filename)) {
    if(TileDBUtils::is_file(filename)) {
      add_reader(filename);
    }
  }

  // push first cells from all files
  push_from_all_files();
}

void OmicsLoader::push_from_all_files() {
  for(auto& f : m_files) {
    if(!f) continue;

    auto cells = f->get_next_cells();
    std::cout << "REMOVE pushing cells from file " << f->get_filename() << std::endl;
    for(auto& c : cells) {
      if(OmicsCell::is_invalid_cell(c)) continue;
      c.coords = m_schema->swap_order(c.coords); // to schema order
      std::cout << "\tREMOVE " << container_to_string(c.coords) << std::endl;
      m_pq.push(c);
    }
  }
}

void OmicsLoader::push_from_idxs(const std::set<int>& idxs) {
  for(auto& idx : idxs) {
    if(idx < 0) continue;

    auto cells = m_files[idx]->get_next_cells();
    std::cout << "REMOVE pushing cells from file " << idx << std::endl;
    for(auto& c : cells) {
      if(OmicsCell::is_invalid_cell(c)) continue;
      c.coords = m_schema->swap_order(c.coords); // to schema order
      std::cout << "\tREMOVE " << container_to_string(c.coords) << std::endl;
      m_pq.push(c);
    }
  }
}

void OmicsLoader::push_file_from_cell(const OmicsCell& cell) {
  push_from_idxs({ cell.file_idx });
}

void OmicsLoader::import() {
  std::cout << "OmicsLoader::import" << std::endl;
  
  std::cout << "REMOVE beginning m_pq.size() is " << m_pq.size() << std::endl;

  std::array<int64_t, 2> last_coords = { -1, -1 };
  int level = 0;

  while(m_pq.size()) {
    auto cell = m_pq.top();
    std::cout << "\t\t\tREMOVE top coords are " << container_to_string(m_pq.top().coords) << std::endl;
    m_pq.pop();
    push_file_from_cell(cell); // make sure that all files are still represented in pq

    if(cell.coords[0] < 0 || cell.coords[1] < 0) { // invalid cell
      continue;
    }

    if(less_than(m_pq.top(), cell)) {
      std::cerr << "Error, top of priority queue is less than previous cell" << std::endl;
      std::cerr << "prev: " << container_to_string(cell.coords) << std::endl;
      std::cerr << "top: " << container_to_string(m_pq.top().coords) << std::endl;
      exit(1);
    }

    if(cell.coords == last_coords) {
      level++;
    }
    else {
      level = 0;
      last_coords = cell.coords;
    }

    // write cell
    std::cerr << "REMOVE import buffering cell " << std::endl;
    std::cerr << cell.to_string() << std::endl << std::endl; // FIXME write to disk
    buffer_cell(cell, level);
  }
  std::cerr << "REMOVE after while tiledb_write_buffers" << std::endl;
  tiledb_write_buffers();

  //if(valid && current_cell.coords[0] >= 0 && current_cell.coords[1] >= 0) { // see if last cell needs to be written to disk
  //  std::cout << current_cell.to_string() << std::endl << std::endl; // FIXME write to disk
  //}
}
