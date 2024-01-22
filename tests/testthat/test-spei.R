test_that("cal_spei works", {
  # x = read.table("data-raw/data.txt")$V1
  r = cal_spei(wb)
  expect_equal(r$coef[[1]], 17.1677346 )
  expect_equal(r$coef[[2]], 9.9856074)
  expect_equal(r$coef[[3]], -0.18829083)

  expect_no_warning({
    r = cal_spei(wb, fit = "max-lik")
    r = cal_spi(wb)
  })
})
