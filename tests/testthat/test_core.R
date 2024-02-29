test_that("core functions work", {
  car_data <- generate_car_simulation()
  all_run <- hyperr8_run(car_data)
  expect_equal(length(unique(all_run$model)), 7)
  expect_equal(min(all_run$deltaAIC), 0)
})
