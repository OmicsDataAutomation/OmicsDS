library(Rcpp)
library(omicsds)

test_that("version is valid", {
    version <- omicsds::version()
    #expect_message(version, "0.0.1")
})

test_that("test that omicsds connects to an existing workspace for queries", {
    print(getwd())
    omicsds_handle <- omicsds::connect(workspace="../../../../../src/test/inputs/interval-level-ws", array="array")
    expect_equal(omicsds_handle, 0)
    omicsds::query(handle=omicsds_handle);
    omicsds::disconnect(handle=omicsds_handle)
})



