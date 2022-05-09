#include <iostream>
#include "omicsds_loader.h"

int main(int argc, char* argv[]) {
  //read_sam_file("/nfs/home/andrei/benchmarking_requirements/toy.sam");

  std::cout << "Hello there" << std::endl;

  {
    ReadCountLoader l("/nfs/home/andrei/benchmarking_requirements/sam_list", true);
    std::cout << "After ctor in main" << std::endl;
    l.import();
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
