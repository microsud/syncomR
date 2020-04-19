#' @title Make taxa abundance table
#'
#' @description Converts \code{\link{phyloseq-class}} object to taxa table arranged by time point.
#'
#' @details The output of this function can be used for functions that require taxa tables in matrix form.
#' @param x a otu/ASV/species abundance table
#' @param y list of samples
#' @param tree phylogenetic tree if interested in unifrac.
#' @param method Can be any of the following c('bray','euclidean','unifrac', 'canberra'))
#' @return Normalized or raw counts taxa abundance table (taxaare rows and timepoits columns).
#' @references
#' \itemize{
#' \item{}{Guittar, J., Shade, A., & Litchman, E. (2019). Trait-based community assembly
#' and succession of the infant gut microbiome. Nature communications, 10(1), 512.}
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' otu.tb <- taxa_time_table(ps1.b5, normalize = TRUE, time.col = "Time_hr", remove.zero = TRUE)
#' head(otu.tb)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @keywords Utilities

distance_melt <- function(x, y, tree = NULL,
                          method = c("bray", "euclidean", "manhattan", "unifrac", "canberra")) {
  # meltdist takes a list of samples (x) and an otu table (y) and calculates intersample dissimilarity using (meth)
  require(reshape2)
  # prune the original otu table
  otumat <- y[y$Sample %in% x, !names(y) %in% c("Sample", "time")] %>% as.data.frame()
  row.names(otumat) <- y$Sample[y$Sample %in% x]

  if (method == "unifrac") {

    # create phyloseq object
    j <- phyloseq(
      otu_table(otumat, taxa_are_rows = FALSE),
      drop.tip(tree, tree$tip.label[!tree$tip.label %in% names(otumat)])
    )

    # calculate UniFrac
    j <- UniFrac(j, weighted = TRUE)
  } else {

    # calculate bray-curtis or euclidean dissimilarity
    j <- vegdist(otumat, method = method)
  }

  # melt into a dataframe, append times
  j <- data.frame(melt(as.matrix(j)), stringsAsFactors = FALSE) %>%
    transmute(
      sample1 = Var1,
      sample2 = Var2,
      t1 = y$time[match(sample1, y$Sample)],
      t2 = y$time[match(sample2, y$Sample)],
      dist = value
    ) %>%
    mutate_if(is.factor, as.character)

  return(j)
}
