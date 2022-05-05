#include <iostream>
#include "omicsds_loader.h"

int main(int argc, char* argv[]) {
  //read_sam_file("/nfs/home/andrei/benchmarking_requirements/toy.sam");

  std::cout << "Hello there" << std::endl;
  //OmicsLoader l("/nfs/home/andrei/benchmarking_requirements/sam_list", OmicsLoader::OmicsStorageOrder::ROW_MAJOR);
  OmicsLoader l("/nfs/home/andrei/benchmarking_requirements/sam_list", OmicsLoader::OmicsStorageOrder::COLUMN_MAJOR);
  std::cout << "After ctor in main" << std::endl;
  l.import();
}
