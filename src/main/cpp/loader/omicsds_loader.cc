#include "omicsds_loader.h"

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
    
    for(int i=0; i< len ; i++){
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

bool equivalent_schema(const OmicsSchema& l, const OmicsSchema& r) {
  if(l.size() != r.size()) return false;

  for(auto li = l.begin(), ri = r.begin(); li != l.end(), ri != r.end(); li++, ri++) {
    if(li->first != ri->first) return false;
    if(li->second.type != ri->second.type) return false;
  }

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

// get line using TileDBUtils api to work with cloud storage as well
// return value indicates if line was read
bool OmicsFileReader::generalized_getline(std::string& retval) {
  retval = "";

  while(m_chars_read < m_file_size || m_str_buffer.size()) {
    int idx = m_str_buffer.find('\n');
    if(idx != std::string::npos) {
      retval = retval + m_str_buffer.substr(0, idx); // exclude newline
      m_str_buffer.erase(0, idx + 1); // erase newline
      return true;
    }

    retval = retval + m_str_buffer;
    m_str_buffer.clear();

    int chars_to_read = std::min<ssize_t>(m_buffer_size, m_file_size - m_chars_read);

    if(chars_to_read) {
      TileDBUtils::read_file(m_filename, m_chars_read, m_buffer, chars_to_read);
       m_chars_read += chars_to_read;
    }

    m_str_buffer.insert(m_str_buffer.end(), m_buffer, m_buffer + chars_to_read);
  }

  return false;
}

/*bool OmicsFileReader::manage_buffer() {
  if(m_cell_buffer.size()) {
    return true;
  }
  auto vec = get_next_cells();
  for(auto& c : vec) {
    m_cell_buffer.push_back(c);
  }
  std::sort(m_cell_buffer.begin(), m_cell_buffer.end());
  return (bool)(m_cell_buffer.size());
}*/

SamReader::SamReader(std::string filename, std::shared_ptr<OmicsSchema> schema, int file_idx) : OmicsFileReader(filename, schema, file_idx) {
  m_fp = hts_open(filename.c_str(),"r"); //open bam file
  m_hdr = sam_hdr_read(m_fp); //read header
  m_align = bam_init1(); //initialize an alignment

  if(!m_hdr) {
    std::cout << "SamReader header is null" << std::endl;
  }
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
    std::vector<uint8_t> qseq(len);

    for(int i=0; i< len ; i++){
      qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
    }

    std::string sample = m_filename;

    int64_t position = std::stol(chr) * 100000000 + pos; // FIXME use acutal mapping

    OmicsMultiCell cell({ (int64_t)m_file_idx, position}, m_schema);
    cell.push_empty_cell(m_file_idx);
    cell.subcells[0].add_field("FLAG", flag);
    cell.subcells[0].add_field_ptr("QNAME", qname, std::strlen(qname));

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

int OmicsLoader::create_array(const std::string& workspace, const std::string& array_name, const OmicsSchema& schema, bool column_major) {
  // Initialize context with home dir if specified in command line, else
  // initialize with the default configuration parameters
  TileDB_CTX* tiledb_ctx;
  /*if (argc > 1) {
    TileDB_Config tiledb_config;
    tiledb_config.home_ = argv[1];
    CHECK_RC(tiledb_ctx_init(&tiledb_ctx, &tiledb_config));
  } else { }*/
  CHECK_RC(tiledb_ctx_init(&tiledb_ctx, NULL));

  // Create a workspace
  CHECK_RC(tiledb_workspace_create(tiledb_ctx, workspace.c_str()));

  std::string group = workspace + "/sparse_arrays";

  // Create two groups in the worskpace
  CHECK_RC(tiledb_group_create(tiledb_ctx, group.c_str()));

  std::string full_name = group + "/" + array_name;

  // Prepare parameters for array schema
  std::vector<const char*> attributes_vec;
  for(auto&p : schema) {
    std::cerr << "REMOVE aname is " << p.first << std::endl;
    std::cerr << "REMOVE c_str is " << p.first.c_str() << std::endl;
    std::cerr << "REMOVE address is " << (int64_t)p.first.c_str() << std::endl;
    attributes_vec.push_back(p.first.c_str());
    std::cerr << "REMOVE vec 1 " << container_to_string(attributes_vec) << std::endl;
  }

  std::cerr << "REMOVE vec 2 " << container_to_string(attributes_vec) << std::endl;
  std::cerr << "REMOVE first attribute address is " << (int64_t)(attributes_vec[0]) << std::endl;
  std::cerr << "REMOVE first attribute is " << attributes_vec[0] << std::endl;
  std::cerr << "REMOVE len is " << strlen(attributes_vec[0]) << std::endl;

  const char** attributes = attributes_vec.data();

  int32_t order = column_major ? TILEDB_COL_MAJOR : TILEDB_ROW_MAJOR;

  const char* dimensions[] = { "SAMPLE" , "POSITION" };  // Two dimensions

  int64_t domain[] =
  {
      0, std::numeric_limits<int64_t>::max(),            // SAMPLE limits
      0, std::numeric_limits<int64_t>::max()             // POSITION limits
  };

  std::vector<int32_t> cell_val_num_vec(schema.size(), TILEDB_VAR_NUM);
  const int* cell_val_num = cell_val_num_vec.data();     // Encode all as variable char to better support superimposing cells

  std::vector<int32_t> compression_vec(schema.size() + 1, TILEDB_NO_COMPRESSION); // plus 1 for coordinates
  const int* compression = compression_vec.data();

  std::vector<int32_t> offsets_compression_vec(schema.size(), TILEDB_NO_COMPRESSION);
  const int* offsets_compression = offsets_compression_vec.data();

  int64_t tile_extents[] =
  {
      1,                                                 // SAMPLE extents
      1                                                  // POSIION extents
  };

  std::vector<int32_t> types_vec(schema.size(), TILEDB_CHAR);
  types_vec.push_back(TILEDB_INT64);
  const int* types = types_vec.data();                   // Encode all as variable char to better support superimposing cells, except for coords

  // Set array schema
  TileDB_ArraySchema array_schema;
  tiledb_array_set_schema(
      &array_schema,                                     // Array schema struct
      full_name.c_str(),                                 // Array name
      attributes,                                        // Attributes
      schema.size(),                                     // Number of attributes
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
      2,                                                 // Number of dimensions
      domain,                                            // Domain
      4*sizeof(int64_t),                                 // Domain length in bytes
      tile_extents,                                      // Tile extents
      2*sizeof(int64_t),                                 // Tile extents length in bytes
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

int OmicsLoader::open_array(const std::string& path) {
  CHECK_RC(tiledb_ctx_init(&m_tiledb_ctx, NULL));

  // Initialize array
  CHECK_RC(tiledb_array_init(
      m_tiledb_ctx,                                      // Context
      &m_tiledb_array,                                   // Array object
      path.c_str(),                                      // Array name
      TILEDB_ARRAY_WRITE,                                // Mode
      NULL,                                              // Entire domain
      NULL,                                              // All attributes
      0));                                               // Number of attributes

  return 0;
}

int OmicsLoader::write_buffers() {
  std::vector<void*> buffers_vec;
  std::vector<size_t> buffer_sizes_vec;

  for(size_t i = 0; i < m_schema->size(); i++) {
    buffers_vec.push_back(offset_buffers[i].data());
    buffers_vec.push_back(var_buffers[i].data());

    buffer_sizes_vec.push_back(offset_buffers[i].size());
    buffer_sizes_vec.push_back(var_buffers[i].size());
  }

  buffers_vec.push_back(coords_buffer.data());

  buffer_sizes_vec.push_back(sizeof(int64_t) * coords_buffer.size());  

  // Write to array
  CHECK_RC(tiledb_array_write(m_tiledb_array, const_cast<const void**>(buffers_vec.data()), buffer_sizes_vec.data()));

  return 0;
}

void OmicsLoader::buffer_cell(const OmicsMultiCell& cell) {
  for(auto& sc : cell.subcells) {
    for(size_t i = 0; i < sc.fields.size(); i++) {
      auto& data = sc.fields[i].data;

      offset_buffers[i].push_back(var_buffers[i].size());
      for(auto& c : data) {
        var_buffers[i].push_back(c);
      }
    }
    coords_buffer.push_back(cell.coords[0]);
    coords_buffer.push_back(cell.coords[1]);
  }
}

// FIXME add parallelism of some kind
// FIXME 
OmicsLoader::OmicsLoader(
                         const std::string& file_list,
                         OmicsStorageOrder order,
                         const bool superimpose
                        ): m_order(order), m_superimpose(superimpose), m_pq(((order == COLUMN_MAJOR) ? m_column_major_comparitor : m_row_major_comparitor)), m_schema(std::make_shared<OmicsSchema>()) {
  // FIXME create schema somewhere virtual
  m_schema->emplace("SAMPLE_NAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char));
  m_schema->emplace("QNAME", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_char));
  m_schema->emplace("FLAG", OmicsFieldInfo(OmicsFieldInfo::OmicsFieldType::omics_uint16_t));

  offset_buffers = std::vector<std::vector<size_t>>(m_schema->size());
  var_buffers = std::vector<std::vector<char>>(m_schema->size());
  coords_buffer.clear();

  create_array("/nfs/home/andrei/OmicsDS/build.debug/workspace", "array", *m_schema, m_order == COLUMN_MAJOR);
  open_array("/nfs/home/andrei/OmicsDS/build.debug/workspace/sparse_arrays/array");

  std::ifstream f(file_list); // initialize OmicsFileReaders from list of filenames
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
  while(m_pq.size()) {
    auto cell = m_pq.top();
    std::cout << "\t\t\tREMOVE top coords are " << container_to_string(m_pq.top().coords) << std::endl;
    m_pq.pop();
    push_files_from_cell(cell); // make sure that all files are still represented in pq

    if(cell.coords[0] < 0 || cell.coords[1] < 0) { // invalid cell
      continue;
    }

    // read cells until top of pq is later
    while(m_pq.size() && m_pq.top().coords == cell.coords) {
      std::cout << "\t\t\tREMOVE in while, temp coords are " << container_to_string(m_pq.top().coords) << std::endl;
      OmicsMultiCell temp(m_pq.top());
      m_pq.pop();
      push_files_from_cell(temp); // make sure that all files are still represented in pq

      cell.merge(temp);
    }

    if(less_than(m_pq.top(), cell)) {
      std::cerr << "Error, top of priority queue is less than previous cell" << std::endl;
      std::cerr << "prev: " << container_to_string(cell.coords) << std::endl;
      std::cerr << "top: " << container_to_string(m_pq.top().coords) << std::endl;
      exit(1);
    }

    // write cell
    std::cerr << cell.to_string() << std::endl << std::endl; // FIXME write to disk
    buffer_cell(cell);

    /*auto reader = m_pq.top();
    m_pq.pop(); // get file reader with next cell in sorted order
    auto ival = reader->peek(); // get interval

    if(ival[0] < 0 || ival[1] < 0) { // end of file, do not replace in m_pq
      continue;
    }

    if(!valid) {
      current_cell = reader->pop(); // get next cell
      valid = true;
      m_pq.push(reader); // replace reader in m_pq
      continue;
    }
    else {
      if(less_than(ival, current_cell.coords)) { // files are all assumed to be in sorted order
        std::cerr << "Error: cell {" << ival[0] << ", " << ival[1] << "} from file " << reader->get_filename() << " preceeds last processed cell {" << current_cell.coords[0] << ", " << current_cell.coords[1] << "}" << std::endl;
        exit(1);
      }
      else if(ival == current_cell.coords) { // next cell overlaps with current, attempt to merge
        if(!current_cell.merge(reader->pop())) {
          std::cerr << "Error: schema of cell from file " << reader->get_filename() << " does not match that of last processed cell" << std::endl;
          exit(1);
        }
        else { // merge was successful, replace reader in m_pq
          m_pq.push(reader);
        }
      }
      else { // next cell is after current, write current cell, replace current cell, and replace reader
        std::cout << current_cell.to_string() << std::endl << std::endl; // FIXME write to disk
        current_cell = reader->pop();
        m_pq.push(reader);
      }
    }*/
    
  }
  std::cerr << "REMOVE after while write_buffers" << std::endl;
  write_buffers();

  //if(valid && current_cell.coords[0] >= 0 && current_cell.coords[1] >= 0) { // see if last cell needs to be written to disk
  //  std::cout << current_cell.to_string() << std::endl << std::endl; // FIXME write to disk
  //}
}
