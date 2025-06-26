data("simdf")

test_that("sample file can be loaded and normalization performed", {
  path_to_rcc <- system.file("extdata", "GSE146530", package = "rhythmCI")
  meta <- system.file("extdata", "meta_GSE146530.csv", package = "rhythmCI")

  expect_true(file.exists(meta))

  exp_df3 <- loadnormalizercc(dir = path_to_rcc, meta_dir = meta, normalization = "RUVg",
                              hkgenes = c("LOC106606270", "LOC106600007"), k=3)

  expect_type(exp_df3, "list")

})
