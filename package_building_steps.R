## Loding and installing
# https://kbroman.org/pkg_primer/pages/vignettes.html
# https://github.com/ShadeLab/microbiome_trait_succession/blob/master/succession_analysis.R
# https://github.com/fcentler/EcologicalStabilityPropertiesComputation/blob/master/StabilityComputationV1.0.R
# https://github.com/microbiome/microbiome/blob/master/DESCRIPTION
roxygen2::roxygenize() #Documentation
devtools::document() #
usethis::use_tidy_style()
devtools::document() #
devtools::check(vignettes = F)
# "R CMD INSTALL "F:/Post_doc/MM_MDb-MM/syncomR_0.0.0.9000.tar.gz"

# R CMD INSTALL "F:/Post_doc/MM_MDb-MM/syncomR_0.0.9.tar.gz"

# R CMD INSTALL "F:/Post_doc/MM_MDb-MM/syncomR_0.1.0.tar.gz"

# R CMD INSTALL "F:/Post_doc/MM_MDb-MM/syncomR_0.1.1.tar.gz"

devtools::build_manual()
devtools::build_vignettes()

devtools::build(vignettes = T)

devtools::build_vignettes(
  pkg = ".",
  dependencies = "VignetteBuilder",
  clean = TRUE,
  upgrade = "never",
  quiet = TRUE,
  install = TRUE,
  keep_md = TRUE
)

# Run once to configure package to use pkgdown
usethis::use_pkgdown()
# Run to build the website
devtools::build_site()
browseVignettes("syncomR")
vignette(package = "syncomR")

## bug
options(warn = 2)
pkgdown::build_site(new_process = T)
traceback()
