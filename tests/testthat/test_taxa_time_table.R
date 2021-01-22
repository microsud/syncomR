context("test taxa_time_table")

test_that("taxa_time_table sorting works correctly", {
  library(syncomR)
  data(SyncomFiltData)
  ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
  otu.tb <- taxa_time_table(ps1.b5, normalize = TRUE, time.col = "Time_hr", remove.zero = TRUE)

  expect_equal(sum(otu.tb[, 1]), 1)
})
