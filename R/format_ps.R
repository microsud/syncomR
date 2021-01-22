#' @title Formatting the ASV-level Phyloseq Object
#'
#' @description Convert ASV-level \code{\link{phyloseq-class}} object to species level data.
#'
#' @details This is a utility function that converts ASV-level \code{\link{phyloseq-class}}
#'          to species level data. Usually the rownames are ASV IDs of Seqs. Here, we convert
#'          it to taxonomic level identities to make sense. It can be used for other taxonomic
#'          levels too. NOTE: Removes unassigned ASVs and returns a \code{\link{phyloseq-class}}
#'          object analysis.
#'
#' @param ps a \code{\link{phyloseq-class}}
#'
#' @param tax.level Taxonomic level to aggregate
#'
#' @return Filtered \code{\link{phyloseq-class}} object.
#'
#' @examples
#' data(SyncomRawCounts)
#' pseq <- format_ps(SyncomRawCounts, tax.level = "Species")
#' print(pseq)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @importFrom phyloseq tax_table
#'
#' @importFrom phyloseq subset_taxa
#'
#' @importFrom microbiome aggregate_taxa
#'
#' @export
#' @keywords Data and filtering

format_ps <- function(ps, tax.level = NULL) {
  if (is.null(tax.level)) {
    stop("Please provide tax.level")
  } else {
    message(paste0(
      "Aggregating input phyloseq object to ",
      tax.level,
      " phyloseq object"
    ))
  }
  # require(microbiome)
  Species <- NULL
  taxic <- as.data.frame(tax_table(ps))
  taxic$OTU <-
    rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
  colnames(taxic) <- c(
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "OTU"
  )
  colnames(taxic)
  taxmat <- as.matrix(taxic) # convert it into a matrix.
  new.tax <-
    phyloseq::tax_table(taxmat) # convert into phyloseq compatible file.
  phyloseq::tax_table(ps) <- new.tax # incroporate into phyloseq Object
  colnames(tax_table(ps))
  ps.sp <- microbiome::aggregate_taxa(ps, tax.level)
  # head(tax_table(ps.sp)[1:5])
  # table(tax_table(ps.sp)[, tax.level],
  #      exclude = NULL)
  ps <- phyloseq::subset_taxa(
    ps.sp,
    Species != "Unknown"
  )

  message("Removed unassigned ASVs")
  return(ps)
}
