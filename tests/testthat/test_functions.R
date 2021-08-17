
library(testthat)
library(vrnmf)


test_that("vol_preprocess() functionality", {
	small_example <- vrnmf:::sim_factors(5, 5, 5)
	vol <- vol_preprocess(t(small_example$X))
	expect_equal(length(vol), 8)
})


test_that("volnmf_main() functionality", {
	small_example <- vrnmf:::sim_factors(5, 5, 5)
	vol <- vol_preprocess(t(small_example$X))
	volres <- volnmf_main(vol)
	expect_equal(length(volres), 10)
})


test_that("AnchorFree() functionality", {
	small_example <- vrnmf:::sim_factors(5, 5, 5)
	vol <- vol_preprocess(t(small_example$X))
	vol.anchor <- AnchorFree(vol)
	expect_equal(length(vol.anchor), 6)
})




