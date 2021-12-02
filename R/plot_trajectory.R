#' @title Plot Temporal Changes in Taxa Abundances
#'
#' @description Temporal taxa abundances are plotted for visual exploration.
#'
#' @param ps A \code{\link{phyloseq-class}}
#'
#' @param taxa.level Taxonomic level
#'
#' @param type Type of plot either all or single. If `type = single` taxa must be specified
#'
#' @param taxa If single taxa name based on taxa.level that's choosen must be specified
#'
#' @param time.col Column specifying time variable
#'
#' @param group.variable Column specifying variable to be colored.
#'                       For e.g. if multiple subjects/bioreactors.
#'
#' @param color.pal List of hex color codes
#'
#' @param transform.counts Any one of \code{\link{transform}}
#'
#' @return ggplot object.
#'
#' @references
#' \itemize{
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' \dontrun{
#' data(SyncomFiltData)
#' ps1.sycom.rel <- microbiome::transform(SyncomFiltData, "compositional")
#' fer_cols <- c(
#'   `Bioreactor A` = "#b2182b",
#'   `Bioreactor B` = "#2166ac",
#'   `Bioreactor C` = "#35978f"
#' )
#'
#' tax.trac <- plot_trajectory(ps1.sycom.rel,
#'   time.col = "Time_hr",
#'   taxa.level = "Species",
#'   type = "single",
#'   taxa = "Akkermansia_muciniphila",
#'   group.variable = "StudyIdentifier",
#'   color.pal = fer_cols,
#'   transform.counts = "compositional"
#' )
#' tax.trac + geom_smooth() + theme_syncom()
#'
#' }
#'
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @importFrom phyloseq psmelt
#'
#' @importFrom ggplot2 geom_line
#'
#' @importFrom ggplot2 geom_smooth
#'
#' @importFrom scales pretty_breaks
#'
#' @importFrom ggplot2 labs
#'
#' @export
#' @keywords Analysis and visualization

plot_trajectory <- function(ps,
                            taxa.level = "Species",
                            type = c("all", "single"),
                            taxa = NULL,
                            time.col = "Time_hr",
                            group.variable = NULL,
                            color.pal = NULL,
                            transform.counts = "compositional") {
  # require(scales)

  time <- Abundance <- taxa.pl <- NULL

  ps <- microbiome::transform(ps, transform.counts)

  ps1.sub.df <- psmelt(ps)

  ps1.sub.df$time <- as.numeric(ps1.sub.df[, time.col])

  varf <- sym(taxa.level)

  if (type == "all") {
    ps1.sub.df$Abundance <- ps1.sub.df$Abundance * 100
    AllPlot <- ggplot(ps1.sub.df, aes(x = time, y = Abundance)) +
      # geom_point(aes(fill = Fasting), size = 1.0) +
      geom_line(aes_string(
        group =
          group.variable, color = group.variable
      )) +
      scale_x_continuous(breaks = pretty_breaks(3)) +
      scale_y_continuous(breaks = pretty_breaks(3)) +
      facet_wrap(varf, scales = "free_y") +
      labs(x = "Time (hr)", y = "Rel. Abundance (%)") +
      # scale_fill_manual(values = fasting_cols) +
      scale_color_manual(values = color.pal) +
      # geom_vline(xintercept = 168, size = 1.2, alpha= 0.6)
      theme(panel.border = element_rect("transparent", size = 0.15))

    return(AllPlot)
  } else {
    ps1.sub.df$Abundance <- ps1.sub.df$Abundance * 100
    colnames(ps1.sub.df)[colnames(ps1.sub.df) == taxa.level] <- "taxa.pl"
    ps1.sub.df.tax <- subset(ps1.sub.df, taxa.pl == taxa)

    IndPlot <- ggplot(ps1.sub.df.tax, aes_string(x = "time", y = "Abundance")) +
      # geom_point(aes(fill = Fasting), size = 1.0) +
      geom_line(aes_string(
        group = group.variable,
        color = group.variable
      )) +
      scale_x_continuous(breaks = pretty_breaks(3)) +
      scale_y_continuous(breaks = pretty_breaks(3)) +
      # facet_wrap(varf, scales = "free_y") +
      labs(x = "Time (hr)", y = "Rel. Abundance (%)") +
      # scale_fill_manual(values = fasting_cols) +
      scale_color_manual(values = color.pal) +
      theme(panel.border = element_rect("transparent", size = 0.15)) +
      labs(subtitle = taxa)

    return(IndPlot)
  }
}
