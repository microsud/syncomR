#' @title Plot Coverage of Taxa
#'
#' @description Calculates and plots coverage of taxa from a phyloseq object.
#'
#' @details Calculates and plots coverage of species from a phyloseq object and
#'          returns a ggplot object. Useful for checking co-existence over
#'          different timepoints. NOTE: Should be raw counts not relative. Adapted from
#'          Fukuyama, J et al. 2017 PLoS computational biology see ref below.
#'
#' @param ps A \code{\link{phyloseq-class}}
#'
#' @param coverage Percent coverage limit. i.e. if 0.95 then number of taxa
#'                 accounting for 95 percent relative abundance in each sample
#'                 timepoint is ploted
#'
#' @param time.column Name of column in sample_data i.e. metadata specifying
#'                    time points.
#'
#' @param shape.variable Name of column in sample_data i.e. metadata specifying
#'                       variable to use for shapes
#'
#' @param color.variable Name of column in sample_data specifying variable to
#'                       use for coloring
#'
#' @param color.pal Colors to use for color.variable
#'
#' @param y.breaks Set breaks for y-axis
#'
#' @param y.limits Set limits for y-axis
#'
#' @param plot TRUE or FALSE if TRUE plot is returned else a data frame
#'
#' @return A ggplot object.
#'
#' @references
#' \itemize{
#' \item{}{"Fukuyama, J., Rumker, L., Sankaran, K., Jeganathan, P., Dethlefsen, L.,
#' Relman, D. A., & Holmes, S. P. (2017). Multidomain analyses of a longitudinal
#' human microbiome intestinal cleanout perturbation experiment.
#' PLoS computational biology, 13(8), e1005706.
#' \url{https://doi.org/10.1371/journal.pcbi.1005706}}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @importFrom ggplot2 aes_string
#'
#' @importFrom ggplot2 scale_color_manual
#'
#' @importFrom ggplot2 scale_fill_manual
#'
#' @importFrom ggplot2 scale_y_continuous
#'
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#'
#' @examples
#'
#' data(SyncomFiltData)
#' fasting_cols <- c("#b3de69", "#fb8072", "#80b1d3")
#' pl <- taxa_coverage(SyncomFiltData,
#'   coverage = 0.9999,
#'   time.column = "Time_hr_num",
#'   shape.variable = "Acetate_Feed",
#'   color.variable = "Fasting",
#'   color.pal = fasting_cols,
#'   y.breaks = seq(0, 16, 1),
#'   y.limits = c(9, 16)
#' )
#' pl + facet_wrap(~StudyIdentifier)
#' pl
#' @export
#' @keywords Anlaysis and visualization

taxa_coverage <- function(ps,
                          coverage = 0.95,
                          time.column = "Time_hr_num",
                          shape.variable = NULL,
                          color.variable = NULL,
                          # facet.var = NULL,
                          color.pal = NULL,
                          y.breaks = seq(0, 16, 1),
                          y.limits = c(11, 16),
                          plot = TRUE) {
  abund <- t(phyloseq::otu_table(ps))
  samples <- phyloseq::sample_data(ps)

  calculate_coverage <- function(abund, coverage) {
    relAb <- abund / sum(abund)
    relAb <- sort(relAb, decreasing = TRUE)
    cumprop <- cumsum(relAb)
    min(which(cumprop >= coverage))
  }

  diversity <- apply(abund, 1, calculate_coverage, coverage)

  if (plot == TRUE) {
    p <- ggplot(data.frame(diversity, phyloseq::sample_data(ps))) +
      geom_point(
        aes_string(
          x = time.column,
          color = color.variable,
          y = diversity,
          shape = shape.variable
        ),
        size = 4,
        alpha = 0.7
      )
    p <-
      p + scale_color_manual(values = color.pal)
    p <-
      p + theme(panel.border = element_rect("transparent", size = 0.15))
    p <-
      p + labs("x" = "Time") + ylab("No. of taxa") + scale_y_continuous(
        breaks = y.breaks,
        limits = y.limits
      )
    # geom_smooth(aes(x = Time_hr_num, color = Fasting, y =diversity))
    return(p)
  } else {
    return(data.frame(diversity, sample_data(ps)))
  }
}
