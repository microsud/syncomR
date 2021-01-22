
context("test stability properties")

test_that("stability_properties works correctly", {
  library(syncomR)
  data(SyncomFiltData)
  ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
  ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 120)
  dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr")

  expect_equal(dat.stab$numberOfGates, 16)
  expect_equal(names(dat.stab$referenceState[1]), "Bacteroides_xylanisolvens")
})
