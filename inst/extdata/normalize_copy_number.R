#' @title Normalize phyloseq object with user provided 16S rRNA gene copy number(s)
#'
#' @description Provided a copy number table with names matching the taxa in in \code{\link{phyloseq-class}} object to species level data.
#'
#' @details Normalizes the read count of input a \code{\link{phyloseq-class}} by 16S rRNA copy number
#' and returns a \code{\link{phyloseq-class}} object. NOTE: Should be raw counts not relative.
#' @param ps A \code{\link{phyloseq-class}}
#' @param column_with_ids Genus_species i.e column with names to match
#' @param copy_num_tab Table (data.frame) with copynumbers
#' @return Copy number normalized \code{\link{phyloseq-class}} object.
#' @references
#' \itemize{
#' \item{}{"Shetty SA et al (2019-2024)
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomRawCounts)
#' copy_num <- read.table("database/mdbmm_copy_numbers.txt", header = T, row.names = 1, stringsAsFactors = F, sep = "\t")
#' ps_cp_normalized <- normalize_copy_number(SyncomRawCounts,
#'   column_with_ids = "Genus_species",
#'   copy_num_tab = copy_num
#' )
#' head(abundances(ps_cp_normalized))
#' head(abundances(SyncomRawCounts)) # Compare  with normalized counts
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @importFrom microbiome abundances
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq tax_table
#' @export
#' @keywords Data and filtering

normalize_copy_number <-
  function(ps, column_with_ids, copy_num_tab) {
    # Extract strain counts from phyloseq object
    counts <- as.data.frame(abundances(ps))
    copy_num_file <- copy_num_tab
    # column_with_ids = "Genus_species"
    intersect(rownames(counts), copy_num_file[, column_with_ids])

    taxa_names <- rownames(counts)

    copy_num_file <-
      copy_num_file[order(match(copy_num_file$Genus_species, rownames(counts))), ]

    # The copy_number file needs to have a column matching the rownames of the otu_table.
    # The copy_number file needs to have a column matching the rownames of the otu_table.
    suppressWarnings(if (rownames(counts) != copy_num_file[, column_with_ids]) {
      print(
        "The copy_number file needs to have a column with matching values to the rownames of the otu_table."
      )
    })


    rownames(copy_num_file) <- copy_num_file[, column_with_ids]

    copy_num_file <- copy_num_file[taxa_names, ]

    corrected_tab <-
      ceiling((counts) / copy_num_file[, "copy_number"])
    otu_table(ps) <-
      otu_table(as.matrix(corrected_tab), taxa_are_rows = T)
    message("Check normalized and un-normalized counts for sanity")
    # write.csv(corrected_tab, "output/corrected_tab.csv")
    # DT::datatable(otu_table(ps1))
    # message("Check uncorrected table in output/corrected_tab.csv")
    # write.csv(otu_table(ps), "output/not_corrected_tab.csv")
    return(ps)
  }
