context("Sample size results")

test_that("Returns same value as Nam 1992", {
  sample_size_test <- function(correct) {
    ceiling(
      power.cmh.test(
        p2 = c(0.75,0.70,0.65,0.60),
        theta = 3,
        power = 0.9,
        t = c(0.10,0.40,0.35,0.15),
        alternative = "greater",
        correct = correct
      )$N
    )
  }
  # Uncorrected
  expect_equal(sample_size_test(FALSE),171)

  # Corrected
  expect_equal(sample_size_test(TRUE),192)
})
