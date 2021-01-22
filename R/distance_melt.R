#' @title Calculate and Melt Distance Matrix
#'
#' @description This function calculates a user specified distance measure
#'              and returns a long data frame.
#'
#' @param x a OTU/ASV/species abundance table
#'
#' @param y list of samples
#'
#' @param tree phylogenetic tree if interested in unifrac.
#'
#' @param method Can be any of the following c('bray','euclidean','unifrac', 'canberra'))
#'
#' @return Normalized or raw counts taxa abundance table (taxa are rows and timepoints columns).
#'
#' @references
#' \itemize{
#' \item{}{Guittar, J., Shade, A., & Litchman, E. (2019). Trait-based community assembly
#' and succession of the infant gut microbiome. Nature communications, 10(1), 512.}
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @keywords Utilities
#'
#' @importFrom reshape2 melt

distance_melt <- function(x, y, tree = NULL,
                          method = c("bray", "euclidean", "manhattan", "unifrac", "canberra")) {

  # distance_melt takes a list of samples (x) and an otu table (y) and calculates intersample dissimilarity using (meth)

  Var1 <- Var2 <- sample1 <- sample2 <- value <- NULL
  # prune the original otu table
  otumat <- y[y$Sample %in% x, !names(y) %in% c("Sample", "time")] %>% as.data.frame()

  row.names(otumat) <- y$Sample[y$Sample %in% x]

  if (method == "unifrac") {

    # create phyloseq object
    j <- phyloseq::phyloseq(
      otu_table(otumat,
        taxa_are_rows = FALSE
      ),
      ape::drop.tip(
        tree,
        tree$tip.label[!tree$tip.label %in% names(otumat)]
      )
    )

    # calculate UniFrac
    j <- phyloseq::UniFrac(j, weighted = TRUE)
  } else {

    # calculate bray-curtis or euclidean dissimilarity
    j <- vegan::vegdist(otumat, method = method)
  }

  # melt into a dataframe, append times
  j <- data.frame(melt(as.matrix(j)), stringsAsFactors = FALSE) %>%
    dplyr::transmute(
      sample1 = Var1,
      sample2 = Var2,
      t1 = y$time[match(sample1, y$Sample)],
      t2 = y$time[match(sample2, y$Sample)],
      dist = value
    ) %>%
    dplyr::mutate_if(is.factor, as.character)

  return(j)
}
