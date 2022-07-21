#include <catch2/catch.hpp>
#include "test_base.h"

#include "omicsds.h"

TEST_CASE("test version - sanity check", "[version]") {
  CHECK(!OmicsDS::version().empty());
}

TEST_CASE("test basic query", "[basic]") {
  auto handle = OmicsDS::connect(std::string(OMICSDS_TEST_INPUTS)+"interval-level-ws", "array");
  CHECK(handle >= 0);
  std::array<int64_t, 2> sample_range = {0, 1000000};
  std::array<int64_t, 2> position_range = {0, 1000000};
  OmicsDS::query(handle, sample_range, position_range);
  OmicsDS::disconnect(handle);
}
