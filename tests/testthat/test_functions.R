
library(testthat)
library(vrnmf)


test_that("vol_preprocess() functionality", {
	small_example <- sim_factors(5, 5, 5)
	vol <- vol_preprocess(t(small_example$X))
	expect_equal(length(vol), 8)
})


test_that("AnchorFree() functionality", {
	small_example <- sim_factors(5, 5, 5)
	vol <- vol_preprocess(t(small_example$X))
	vol.anchor <- AnchorFree(vol)
	expect_equal(length(vol.anchor), 6)
})




