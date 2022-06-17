#include "omicsds_loader.h"

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

  std::cout << "name orders" << std::endl;
  for(auto& i : idxs_name) {
    std::cout << i << std::endl;
  }
  std::cout << std::endl << std::endl;

  std::cout << "position orders" << std::endl;
  for(auto& i : idxs_position) {
    std::cout << i << std::endl;
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

std::vector<uint8_t> OmicsMultiCell::as_cell() {
  return {};
  //construct cell
  /*std::vector<uint8_t> cell(16);
  *(reinterpret_cast<uint64_t*>(cell.data())) = coords[0]; // write row in cell
  *(reinterpret_cast<uint64_t*>(cell.data()) + 1) = coords[1]; // write position in cell

  // reserve space for cell size
  for(int i = 0; i < sizeof(size_t); i++) {
    cell.push_back(0);
  }

  // attributes
  for(int i = 0; i < schema.attribute_num(); i++) {
    std::string attribute_name = schema.attribute_name(i);

    for(auto& sc : subcells) {
      
      cell.push_back();

      if(attribute_name == "NAME") {
        cell.insert(cell.end(), {0, 0, 0, 0});
        *(reinterpret_cast<uint32_t*>(cell.data() + cell.size() - 4)) = name.length();
        for(auto c : name) {
          cell.push_back(c);
        }
      }
    }
  }
  // fill in cell size
  *(reinterpret_cast<size_t*>(cell.data() + 2*sizeof(int64_t))) = cell.size();*/
}

OmicsMultiCell OmicsMultiCell::create_invalid_cell() {
  OmicsMultiCell rv;
  rv.coords = {-1, -1};
  return rv;
}

bool OmicsMultiCell::is_invalid_cell(const OmicsMultiCell& cell) {
  return cell.coords[0] < 0 || cell.coords[1] < 0;
}

SamReader::SamReader(std::string filename, std::shared_ptr<OmicsSchema> schema, int file_idx) : OmicsFileReader(filename, schema, file_idx) {
  m_fp = hts_open(filename.c_str(),"r"); //open bam file
  m_hdr = sam_hdr_read(m_fp); //read header
  m_align = bam_init1(); //initialize an alignment

  if(!m_hdr) {
    std::cout << "SamReader header is null" << std::endl;
  }

  assert((bool)schema);
}

SamReader::~SamReader() {
  //std::cout << "REMOVE SamReader::~SamReader" << std::endl;
  //std::cout << "FIXME uncomment cleanup functions" << std::endl;
  bam_destroy1(m_align);
  sam_close(m_fp);
}

std::vector<OmicsMultiCell> SamReader::get_next_cells() {
  std::vector<OmicsMultiCell> cells;

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

    OmicsMultiCell cell({ (int64_t)m_file_idx, position}, m_schema);
    cell.push_empty_cell(m_file_idx);
    cell.subcells[0].add_field_ptr("QNAME", qname, std::strlen(qname));
    cell.subcells[0].add_field("FLAG", flag);
    cell.subcells[0].add_field_ptr("RNAME", chr, std::strlen(chr));
    cell.subcells[0].add_field("POS", pos);
    cell.subcells[0].add_field("MAPQ", mapq);
    cell.subcells[0].add_field_ptr("CIGAR", cigar, n_cigar);
    cell.subcells[0].add_field("RNEXT", rnext);
    cell.subcells[0].add_field("PNEXT", pnext);
    cell.subcells[0].add_field("TLEN", tlen);
    cell.subcells[0].add_field_ptr("SEQ", qseq.data(), (int)qseq.size());
    cell.subcells[0].add_field_ptr("QUAL", qual, std::strlen(qual));
  

    OmicsMultiCell end_cell = cell;
    end_cell.subcells[0].file_idx = -1;
    int end_offset = std::abs(tlen); // FIXME figure out negative template length
    end_cell.coords[1] += end_offset;
    if(end_offset) {  // if end cell is in same position, only create one cell
      return { cell, end_cell };
    }
    return { cell };
  }
  else {
    std::cout << "REMOVE sam_read1 was " << rc << std::endl;
    return { OmicsMultiCell::create_invalid_cell() };
  }
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

void OmicsLoader::buffer_cell(const OmicsMultiCell& cell, int level) {
  for(auto& sc : cell.subcells) {
    auto fiter = sc.fields.begin();
    auto aiter = m_schema->attributes.begin();
    int i = 0;
    for(; fiter != sc.fields.end() && aiter != m_schema->attributes.end(); fiter++, aiter++, i++) {
      auto& data = fiter->data;

      int length = aiter->second.length;
      int size = aiter->second.element_size();
      if(length == TILEDB_VAR_NUM) { // variable length
        offset_buffers[i].push_back(var_buffers[i].size());
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

// FIXME add parallelism of some kind
// FIXME 
OmicsLoader::OmicsLoader(
                         const std::string& workspace,
                         const std::string& array,
                         const std::string& file_list,
                         const std::string& mapping_file,
                         bool position_major
                        ): OmicsModule(workspace, array, mapping_file, position_major), m_file_list(file_list), m_pq(comparitor) {}

void OmicsLoader::initialize() { // FIXME move file reader creation to somewhere virtual
  create_schema();

  offset_buffers = std::vector<std::vector<size_t>>(m_schema->attributes.size());
  var_buffers = std::vector<std::vector<char>>(m_schema->attributes.size());
  coords_buffer.clear();

  // FIXME don't hardcode
  //tiledb_create_array("/nfs/home/andrei/OmicsDS/build.debug/workspace", "array", *m_schema);
  //tiledb_open_array("/nfs/home/andrei/OmicsDS/build.debug/workspace/array");
  
  tiledb_create_array(m_workspace, m_array, *m_schema);
  tiledb_open_array(m_workspace, m_array);

  std::ifstream f(m_file_list); // initialize OmicsFileReaders from list of filenames
  std::string s;
  while(std::getline(f, s)) {
    if(TileDBUtils::is_file(s)) {
      //m_files.push_back(std::make_shared<OmicsFileReader>(s));
      if(std::regex_match(s, std::regex("(.*)(sam)($)"))) {
        m_files.push_back(std::make_shared<SamReader>(s, m_schema, m_files.size()));
        std::cout << "REMOVE push sam reader " << s << std::endl;
      }
      else {
        std::cout << "REMOVE would push " << s << std::endl;
      }
    }
  }

  push_from_all_files();
}

void OmicsLoader::push_from_all_files() {
  for(auto& f : m_files) {
    if(!f) continue;

    auto cells = f->get_next_cells();
    std::cout << "REMOVE pushing cells from file " << f->get_filename() << std::endl;
    for(auto& c : cells) {
      if(OmicsMultiCell::is_invalid_cell(c)) continue;
      c.coords = m_schema->standard_to_schema_order(c.coords);
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
      if(OmicsMultiCell::is_invalid_cell(c)) continue;
      c.coords = m_schema->standard_to_schema_order(c.coords);
      std::cout << "\tREMOVE " << container_to_string(c.coords) << std::endl;
      m_pq.push(c);
    }
  }
}

void OmicsLoader::push_files_from_cell(const OmicsMultiCell& cell) {
  push_from_idxs(cell.get_file_idxs());
}

void OmicsLoader::import() {
  std::cout << "OmicsLoader::import" << std::endl;
  //OmicsMultiCell current_cell;
  //bool valid = false; // whether value of current_cell is meaningfull
  
  std::cout << "REMOVE beginning m_pq.size() is " << m_pq.size() << std::endl;

  std::array<int64_t, 2> last_coords = { -1, -1 };
  int level = 0;

  while(m_pq.size()) {
    auto cell = m_pq.top();
    std::cout << "\t\t\tREMOVE top coords are " << container_to_string(m_pq.top().coords) << std::endl;
    m_pq.pop();
    push_files_from_cell(cell); // make sure that all files are still represented in pq

    if(cell.coords[0] < 0 || cell.coords[1] < 0) { // invalid cell
      continue;
    }

    // read cells until top of pq is later
    /*while(m_pq.size() && m_pq.top().coords == cell.coords) {
      std::cout << "\t\t\tREMOVE in while, temp coords are " << container_to_string(m_pq.top().coords) << std::endl;
      OmicsMultiCell temp(m_pq.top());
      m_pq.pop();
      push_files_from_cell(temp); // make sure that all files are still represented in pq

      cell.merge(temp);
    }*/

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
    std::cerr << cell.to_string() << std::endl << std::endl; // FIXME write to disk
    buffer_cell(cell, level);
  }
  std::cerr << "REMOVE after while tiledb_write_buffers" << std::endl;
  tiledb_write_buffers();

  //if(valid && current_cell.coords[0] >= 0 && current_cell.coords[1] >= 0) { // see if last cell needs to be written to disk
  //  std::cout << current_cell.to_string() << std::endl << std::endl; // FIXME write to disk
  //}
}

std::pair<std::vector<void*>, std::vector<size_t>> OmicsReader::prepare_buffers() {
  m_buffers_vector.clear();
  std::vector<void*> pointers;
  std::vector<size_t> sizes;

  for(auto&[_, inf] : m_schema->attributes) {
    for(int i = 0; i < 1 + (inf.is_variable()); i++) { // 1 buffer if fixed length, 2 if variable
      m_buffers_vector.emplace_back(m_buffer_size);
      pointers.push_back(m_buffers_vector.back().data());
      sizes.push_back(m_buffers_vector.back().size());
    }
  }
  // coords
  m_buffers_vector.emplace_back(m_buffer_size);
  pointers.push_back(m_buffers_vector.back().data());
  sizes.push_back(m_buffers_vector.back().size());
  
  return { pointers, sizes };
}

void OmicsReader::query() {
  auto[pointers_vec, sizes_vec] = prepare_buffers();
  void** pointers = pointers_vec.data();
  size_t* sizes = sizes_vec.data();

  std::string array_name = m_workspace + "/" + m_array;  

  // int64_t subarray[] = { 3, 4, 2, 4 };

  TileDB_ArrayIterator* tiledb_array_it;
  tiledb_array_iterator_init(
      m_tiledb_ctx,                                  // Context
      &tiledb_array_it,                              // Array iterator
      array_name.c_str(),                            // Array name
      TILEDB_ARRAY_READ,                             // Mode
      0,                                      // Constrain in subarray
      0,                                    // Subset on attributes
      m_schema->attributes.size(),                   // Number of attributes
      pointers,                                       // Buffers used internally
      sizes);                                 // Buffer sizes

 
  const int* a1_v;
  size_t a1_size;
  while(!tiledb_array_iterator_end(tiledb_array_it)) {
    int i = -1;
    std::cout << std::endl << std::endl << "New cell" << std::endl;
    for(auto& a : m_schema->attributes) {
      ++i;
      // Get value
      tiledb_array_iterator_get_value(
          tiledb_array_it,     // Array iterator
          i,                   // Attribute id
          (const void**) &a1_v,// Value
          &a1_size);           // Value size (useful in variable-sized attributes)

      // Print value (if not a deletion)
      if(*a1_v != TILEDB_EMPTY_INT32) { 
        std::cout << "attribute " << a.first << std::endl;
        printf("%3d\n", *a1_v);
        std::cout << a.first << " size " << a1_size << std::endl;
      }
    }


    ++i;
    // FIXME fix i
    i = m_schema->attributes.size();
    uint64_t* coords;
    tiledb_array_iterator_get_value(
        tiledb_array_it,     // Array iterator
        i,                   // Attribute id
        (const void**)&coords,// Value
        &a1_size);           // Value size (useful in variable-sized attributes)

    // Print value (if not a deletion)
    if(*a1_v != TILEDB_EMPTY_INT32) {
      std::cout << "coords " << coords[0] << ", " << coords[1] << ", " << coords[2] << std::endl;
      std::cout << "coords size " << a1_size << std::endl;
    }

    // Advance iterator
    if(tiledb_array_iterator_next(tiledb_array_it)) {
      exit(1);
    }
  }

}
