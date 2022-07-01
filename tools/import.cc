#include <iostream>
#include <getopt.h>
#include "omicsds_loader.h"
// FIXME REMOVE
#include <chrono>
#include <thread>

enum ArgsEnum {
  ARGS_IDX_MAPPER
};

void print_usage() {
  std::cout << "Usage: omicsds_import [options]\n"
            << "where options include:\n"
            << "Generic options\n"
            << "\t \e[1m--workspace\e[0m, \e[1m-w\e[0m Path to workspace\n"
            << "\t \e[1m--array\e[0m, \e[1m-a\e[0m Name of array (should not include path to workspace)\n"
            << "Import options\n"
            << "\t \e[1m--mapping-file\e[0m, \e[1m-m\e[0m Path to file containing information to map from contig/offset pair to flattened coordinates. Currently supports fasta.fai (only needs first 3 columns: contig name, length, and starting index separated by tabs) \n"
            << "\t \e[1m--file-list\e[0m, \e[1m-f\e[0m Path to file containing paths to files to be ingested (one path per line)\n"
            << "\t \e[1m--sample-map\e[0m, \e[1m-s\e[0m Path to file containing information mapping between samples names and row indices in OmicsDS (each line is a sample name and an integer row number separated by a tab)\n"
            << "\t \e[1m--sample-major\e[0m, Indicates that data from input files should be stored in sample major order on disk (cells for the same sample stored contiguously). Default behavior is position major\n"
            << "\t\t Optional\n"
            << "\t \e[1m--read-counts\e[0m, \e[1m-r\e[0m Command to ingest read count related data (file list should contain SAM files) \n"
            << "\t\t Optional\n"
            << "\t \e[1m--transcriptomics\e[0m, \e[1m-t\e[0m Command to transcriptomic data (file list should contain Bed and Matrix files) \n"
            << "\t\t Optional\n"
            << "\t\t \e[1m--gene-mapping-file\e[0m, Path to gtf/gff/gi v1/gbed file for use with transcriptomics \n"
            << "Query options\n"
            << "\t \e[1m--generic-query\e[0m, \e[1m-g\e[0m Command to perform generic query. WIP. Output is probably not useful.\n"
            << "\t\t Optional\n"
            << "\t \e[1m--export-sam\e[0m, Command to export data from query range as sam files, one per sample. Should only be used on data ingested via --read-counts\n"
            << "\t\t Optional\n";
}

int main(int argc, char* argv[]) {
  char* cwd1 = getcwd(0, 0);
  std::cerr << "**************************** first cwd is " << cwd1 << std::endl;

  //read_sam_file("/nfs/home/andrei/benchmarking_requirements/toy.sam");

  enum options_enum { EXPORT_SAM, SAMPLE_MAJOR, GENE_MAPPING_FILE };

  static struct option long_options[] = {
    {"workspace",1,0,'w'},
    {"array",1,0,'a'},
    {"mapping-file",1,0,'m'},
    {"file-list",1,0,'f'},
    {"sample-map",1,0,'s'},
    {"read-counts",0,0,'r'},
    {"generic-query", 0, 0, 'g'},
    {"export-sam", 0, 0, EXPORT_SAM},
    {"transcriptomics", 0, 0, 't'},
    {"sample-major", 0, 0, SAMPLE_MAJOR},
    {"gene-mapping-file", 1, 0, GENE_MAPPING_FILE}
  };

  std::string workspace = "";
  std::string array = "";
  std::string mapping_file = "";
  std::string file_list = "";
  std::string sample_map = "";
  bool read_counts = false;
  bool generic_query = false;
  bool export_sam = false;
  bool transcriptomics = false;
  bool position_major = true;
  std::string gene_mapping_file = "";

  int c;
  while ((c=getopt_long(argc, argv, "w:a:m:f:rs:gt", long_options, NULL)) >= 0) {
    switch (c) {
      case 'w':
        workspace = std::string(optarg);
        break;
      case 'a':
        array = std::string(optarg);
        break;
      case 'm':
        mapping_file = std::string(optarg);
        break;
      case 'f':
        file_list = std::string(optarg);
        break;
      case 's':
        sample_map = std::string(optarg);
        break;
      case 'r':
        read_counts = true;
        break;
      case 'g':
        generic_query = true;
        break;
      case EXPORT_SAM:
        export_sam = true;
        break;
      case 't':
        transcriptomics = true;
        break;
      case SAMPLE_MAJOR:
        position_major = false;
        break;
      case GENE_MAPPING_FILE:
        gene_mapping_file = std::string(optarg);
        break;
      default:
        std::cerr << "Unknown command line argument " << char(c) << "\n";
        print_usage();
        return -1;
    }
  }

  if(workspace == "") {
    std::cerr << "Workspace required\n";
    print_usage();
    return -1;
  }
  if(array == "") {
    std::cerr << "Array required\n";
    print_usage();
    return -1;
  }
  if(read_counts || transcriptomics) {
    if(mapping_file == "") {
      std::cerr << "Mapping file required\n";
      print_usage();
      return -1;
    }
    if(file_list == "") {
      std::cerr << "File list required\n";
      print_usage();
      return -1;
    }
    if(sample_map == "") {
      std::cerr << "Sample map required\n";
      print_usage();
      return -1;
    }
  }

  if(transcriptomics) {
    if(gene_mapping_file == "") {
      std::cerr << "Gene mapping file required\n";
      print_usage();
      return -1;
    }
  }

  std::cout << "Hello there: " << workspace << ", " << array << ", " << mapping_file << std::endl;

  if(read_counts) {
    {
      ReadCountLoader l(workspace, array, file_list, sample_map, mapping_file, position_major);
      l.initialize();
      std::cout << "After ctor in main" << std::endl;
      l.import();
    }
  }

  if(transcriptomics) {
    TranscriptomicsLoader l(workspace, array, file_list, sample_map, mapping_file, gene_mapping_file, position_major);
    l.initialize();
    l.import();
  }

  if(generic_query) {
    std::cout << "===================================== NORMAL QUERY =========================================" << std::endl;
    OmicsExporter r(workspace, array);
    //r.query({0, 1}, {59, 61});
    r.query();
  }

  if(export_sam) {
    std::cout << "===================================== SAM QUERY =========================================" << std::endl;
    SamExporter s(workspace, array);
    s.export_sams();
  }

  /*// FIXME remove
  std::cerr << "FIXME remove end of main reading" << std::endl;

  // ================================== ARRAY READ ======================

  TileDB_CTX* tiledb_ctx;
  TileDB_Array* tiledb_array;

  CHECK_RC(tiledb_ctx_init(&tiledb_ctx, NULL));

  const char array_name[] = "/nfs/home/andrei/OmicsDS/build.debug/workspace/sparse_arrays/array";

  char buffer[1024];

  char* cwd = getcwd(0, 0);
  std::cerr << "**************************** cwd is " << cwd << std::endl;

  // Initialize array
  CHECK_RC(tiledb_array_init(
           tiledb_ctx,                           // Context
           &tiledb_array,                        // Array object
           &array_name[0],                       // Array name
           TILEDB_ARRAY_READ,                    // Mode
           NULL,                                 // Whole domain
           NULL,                                 // All attributes
           0));                                  // Number of attributes

  // Prepare cell buffers
  size_t buffer_sample[50];
  char buffer_sample_var[50];
  size_t buffer_qname[50];
  char buffer_qname_var[50];
  uint16_t buffer_flag[50];
  int64_t buffer_coords[50];
  void* r_buffers[] =
      { buffer_sample, buffer_sample_var, buffer_qname, buffer_qname_var, buffer_flag, buffer_coords };
  size_t r_buffer_sizes[] =
  {
      sizeof(buffer_sample),
      sizeof(buffer_sample_var),
      sizeof(buffer_qname),
      sizeof(buffer_qname_var),
      sizeof(buffer_flag),
      sizeof(buffer_coords)
  };

  // Read from array
  CHECK_RC(tiledb_array_read(tiledb_array, r_buffers, r_buffer_sizes));

  // Print cell values
  int64_t result_num = r_buffer_sizes[0] / sizeof(int);
  printf("%ld results\n", (long)result_num);
  printf("coords\t flag\t   sample\t    qname\n");
  printf("-----------------------\n");
  for(int i=0; i<result_num; ++i) {
    printf("%ld, %ld, %ld", (long)buffer_coords[3*i], (long)buffer_coords[3*i+1], (long)buffer_coords[3*i+2]);

    printf("\t %3d", buffer_flag[i]);

    size_t var_size = (i != result_num-1) ? buffer_sample[i+1] - buffer_sample[i]
                                          : r_buffer_sizes[2] - buffer_sample[i];
    printf("\t %4.*s\n", int(var_size), &buffer_sample_var[buffer_sample[i]]);

    var_size = (i != result_num-1) ? buffer_qname[i+1] - buffer_qname[i]
                                          : r_buffer_sizes[2] - buffer_qname[i];
    printf("\t %4.*s\n", int(var_size), &buffer_qname_var[buffer_qname[i]]);
  }

  std::cerr << "After reading in import" << std::endl;

  // Finalize the array
  CHECK_RC(tiledb_array_finalize(tiledb_array));

  // Finalize context
  CHECK_RC(tiledb_ctx_finalize(tiledb_ctx));

  return 0;*/
}
