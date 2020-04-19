#' @title Make taxa abundance table
#'
#' @description Converts \code{\link{phyloseq-class}} object to taxa table arranged by time point.
#'
#' @details The output of this function can be used for functions that require taxa tables in matrix form.
#' @param ps a \code{\link{phyloseq-class}}
#' @param normalize TRUE or FALSE uses the normalize function from \code{\link{seqtime}} by Faust et al.
#' @param time.col Specifiy column containing time variable
#' @param remove.zero TRUE or FALSE
#' @return Normalized or raw counts taxa abundance table (taxaare rows and timepoits columns).
#' @references
#' \itemize{
#' \item{}{Faust, K., et al. (2018). Signatures of ecological processes in
#' microbial community time series. Microbiome, 6(1), 120. https://doi.org/10.1186/s40168-018-0496-2}
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' remove_T80 <- c("Ferm_1_5_80", "Ferm_1_6_80", "Ferm_1_8_80")
#' ps1.sub <- prune_samples(!(sample_names(ps1.b5) %in% remove_T80), ps1.b5)
#' otu.tb <- taxa_time_table(ps1.sub, normalize = TRUE, time.col = "Time_hr", remove.zero = TRUE)
#' head(otu.tb)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @export
#' @keywords Anlaysis and visualization

taxa_time_table <-
  function(ps,
           normalize = c(TRUE, FALSE),
           time.col = NULL,
           remove.zero = c(TRUE, FALSE)) {
    otu.mat <- otu.mat2 <- otu.mat3 <- ps.ml <- NULL

    otu.mat <- as.data.frame(abundances(ps))

    require(dplyr)
    require(seqtime)
    ps.ml <- psmelt(ps)
    ps.ml$Time <- as.numeric(ps.ml[, time.col])

    # otu.m <-
    # data.table::dcast(ps.ml, Time ~ OTU, value.var = "Abundance")
    otu.m <- ps.ml %>%
      select(OTU, Time, Abundance) %>%
      group_by(OTU, Time) %>%
      tidyr::pivot_wider(id_cols = Time, names_from = OTU, values_from = Abundance) %>%
      arrange(Time)
    otu.mat2 <- as.data.frame(otu.m)

    # otu.mat2 <- otu.m %>% arrange(Time)
    rownames(otu.mat2) <- otu.mat2$Time
    otu.mat3 <- otu.mat2[, -1]

    if (normalize == TRUE) {
      otu.mat4 <-
        normalize_seqtime(t(otu.mat3), removeZero = remove.zero) # compositional for Whisker neutrailty test
    } else {
      otu.mat4 <- otu.mat3
    }

    return(otu.mat4)
  }
