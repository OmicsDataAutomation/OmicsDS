#include <iostream>
#include <getopt.h>
#include "omicsds_loader.h"

enum ArgsEnum {
  ARGS_IDX_MAPPER
};

void print_usage() {
  std::cout << "Usage: omicsds_import [options]\n"
            << "where options include:\n"
            << "\t \e[1m--mapping-file\e[0m, \e[1m-m\e[0m File containing information to map from contig/offset pair to flattened coordinates. Currently supports fasta.fai\n";
}

int main(int argc, char* argv[]) {
  //read_sam_file("/nfs/home/andrei/benchmarking_requirements/toy.sam");

  static struct option long_options[] = {
    {"mapping-file",1,0,'m'}
  };

  std::string mapping_file = "";

  int c;
  while ((c=getopt_long(argc, argv, "m:", long_options, NULL)) >= 0) {
    switch (c) {
      case 'm':
        mapping_file = std::string(optarg);
        break;
      default:
        std::cerr << "Unknown command line argument " << char(c) << "\n";
        print_usage();
        return -1;
    }
  }

  if(mapping_file == "") {
    std::cerr << "Mapping file required\n";
    print_usage();
    return -1;
  }

  std::cout << "Hello there" << std::endl;

  {
    ReadCountLoader l("/nfs/home/andrei/benchmarking_requirements/sam_list", mapping_file, true);
    l.initialize();
    std::cout << "After ctor in main" << std::endl;
    l.import();
    l.serialize_schema("/nfs/home/andrei/benchmarking_requirements/schema");
  }

  // FIXME remove
  std::cerr << "FIXME remove end of main reading" << std::endl;

  // ================================== ARRAY READ ======================

  TileDB_CTX* tiledb_ctx;
  TileDB_Array* tiledb_array;

  CHECK_RC(tiledb_ctx_init(&tiledb_ctx, NULL));

  const char array_name[] = "/nfs/home/andrei/OmicsDS/build.debug/workspace/sparse_arrays/array";

  // Initialize array
  CHECK_RC(tiledb_array_init(
           tiledb_ctx,                           // Context
           &tiledb_array,                        // Array object
           array_name,                           // Array name
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

  // Finalize the array
  CHECK_RC(tiledb_array_finalize(tiledb_array));

  // Finalize context
  CHECK_RC(tiledb_ctx_finalize(tiledb_ctx));

  return 0;
}
