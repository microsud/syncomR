#' @title Find Taxa That is Highest or Lowest Abundant in Each Sample
#'
#' @description Identifying either taxon with highest or lowest abundance.
#'
#' @details The output of this function can be used for identifying either
#'          taxon with highest or lowest abundance.
#'
#' @param ps a \code{\link{phyloseq-class}} with compositional data
#'
#' @param which.taxa find "most_abund" or "least_abund"
#'
#' @return A data frame with sample_id, taxon highest or lowest abundance,
#'         and abundance value.
#'
#' @references
#' \itemize{
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#'
#' data(SyncomFiltData)
#' ps.rel <- microbiome::transform(SyncomFiltData, "compositional")
#' otu.tb <- find_taxa(ps.rel, which.taxa = "most_abund")
#' head(otu.tb)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @export
#'
#' @keywords Utilities

find_taxa <- function(ps,
                      which.taxa = c("most_abund", "least_abund")) {
  otu.df <- tax <- abund <- sample_id <- taxon <- res <- NULL

  otu.df <- as.data.frame(t(microbiome::abundances(ps)))

  if (which.taxa == "most_abund") {
    tax <- max.col(otu.df, "first")
    abund <- otu.df[cbind(1:nrow(otu.df), tax)]
    sample_id <- cbind(rownames(otu.df), tax)
    taxon <- names(otu.df)[tax]
  } else if (which.taxa == "least_abund") {
    # add positive value to zero to avoid getting zero abundance as lowest
    otu.df[otu.df == 0] <- max(colSums(otu.df)) # trying to figure this out
    tax <- apply(otu.df, 1, which.min)
    abund <- otu.df[cbind(1:nrow(otu.df), tax)]
    sample_id <- cbind(rownames(otu.df), tax)
    taxon <- names(otu.df)[tax]
  }
  res <- data.frame(sample_id, taxon, abund)
  # head(res)
  colnames(res) <- c("sample_id", "tax_index", paste0(which.taxa, "_taxon"), "abundance")
  return(res)
}
