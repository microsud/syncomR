.onLoad <- function(libname = find.package("syncomR"), pkgname = "syncomR") {
  message("loading dependcies packages")
  library(reshape2)
  library(ggrepel)
  library(vegan)
  library(plotrix)
  library(ggpubr)
}
