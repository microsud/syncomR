## code to prepare `syncom` dataset goes here
list.files("data/", pattern = "rda")


usethis::use_data(syncom, overwrite = TRUE)

load("data/SyncomCopyCorrectedCounts.rda")
usethis::use_data(SyncomCopyCorrectedCounts, compress = "xz",overwrite = TRUE)

load("data/SyncomFiltData.rda")
usethis::use_data(SyncomFiltData, compress = "xz",overwrite = TRUE)

load("data/SyncomGMM.rda")
usethis::use_data(SyncomGMM, compress = "xz",overwrite = TRUE)

load("data/SyncomRawCounts.rda")
usethis::use_data(SyncomRawCounts, compress = "xz",overwrite = TRUE)

load("data/SynComRawRNA.rda")
usethis::use_data(SynComRawRNA, compress = "xz",overwrite = TRUE)

load("data/SynComRNAKEGG.rda")
usethis::use_data(SynComRNAKEGG, compress = "xz",overwrite = TRUE)




