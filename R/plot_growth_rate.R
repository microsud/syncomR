#' @title Measure "growth rate" of taxa from temporal abundance data
#'
#' @description Calculates a proxy "growth rate" of taxa using temporal abundnace data.
#'
#' @details This approach is a modification of the approach used by Nathan I. Wisnoski and collegues
#' \link{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on the following calculation
#' log(present abundance/past abundance). However, to identify,  "growth rate" the commonly used approach (not bacteria specific)
#' is log(present abundance-past abudnance/past abundance).
#' This approach is useful for identifying long term growth behavious of taxa.
#' @param ps a non-rarefied raw otu/ASV/species abundance table (longitudinal data)
#' @param rarefy TRUE  Kept this option incase users have absolute quantification.
#' @param depth min(sample_sums(ps) If rarefy=TRUE then depth for rarefaction.
#' @return Plot with labelling for max "growth rate" timepoints.
#' @references
#' \itemize{
#' \item{}{Shetty SA. et al. A Minimalist Approach for Deciphering the Ecophysiology of Human Gut Microbes (2020)}
#' \item{}{To cite the package, see citation('syncomR')}
#' \item{}{Nathan I. Wisnoski and collegues https://github.com/nwisnoski/ul-seedbank}
#' }
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' p.g <- plot_growth_rate(ps1.b5, rarefy = TRUE, depth = min(sample_sums(ps1.b5)))
#' print(p.g)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @import tidyr
#' @import ggrepel
#' @import tibble
#' @export
#' @keywords Anlaysis and visualization

plot_growth_rate <- function(ps, rarefy = TRUE, depth = min(sample_sums(ps))) {
  message("This is a proxy growth rate calculation")
  require(tibble)
  getgrs <- function(sbs) {

    # create new matrix to calculate growth rates
    y.s <- matrix(ncol = ncol(sbs), nrow = nrow(sbs) - 1)
    for (i in 1:(nrow(sbs) - 1)) {
      y.s[i, ] <- log(sbs[i + 1, ] - sbs[i, ] / sbs[i, ])
    }
    return(y.s)
  }

  ps.rar <- otu.tb <- grs.all <- grts <- otu.rar <- maxgrs <- p <- NULL
  ps.rar <- rarefy_even_depth(ps)
  # otu.rar <- abundances(ps.rar)
  otu.tb <- taxa_time_table(ps.rar, normalize = F, time.col = "Time_hr", remove.zero = F)
  otu.tb.t <- t(otu.tb)
  grts <- getgrs(t(otu.tb.t) + .1)
  colnames(grts) <- colnames(t(otu.tb.t))
  rownames(grts) <- rownames(t(otu.tb.t))[-1]


  grs.all <- as.data.frame(grts) %>%
    rownames_to_column(var = "sample.name") %>%
    gather(-"sample.name", key = "OTU", value = "growth") %>%
    mutate(Time = sample.name) %>%
    group_by(OTU, Time) %>%
    summarize(
      mean.growth = mean(growth),
      sd.growth = sd(growth)
    )
  maxgrs <- grs.all %>%
    summarize(max.growth = max(mean.growth))
  grs.all <- grs.all %>%
    left_join(maxgrs)
  grs.all <- grs.all %>%
    mutate(ismax = ifelse(mean.growth == max.growth, T, F))

  grs.all$OTU <- gsub("_", " ", grs.all$OTU)
  # head(persistent.grs.all)
  grs.all$Time <- as.numeric(grs.all$Time)
  grs.all$otu.time <- paste0(grs.all$OTU, " ", grs.all$Time, "h")

  p <- ggplot(grs.all, aes(x = Time, group = OTU, col = OTU)) +
    geom_line(aes(y = mean.growth), alpha = 0.6, size = 1) +
    geom_point(
      data = subset(grs.all, ismax == T),
      aes(y = max.growth), alpha = 0.8, size = 3
    ) +
    # coord_flip() +
    geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = -9),
      fill = "#edf3f5", col = "white", alpha = 0.5
    ) +
    theme(legend.position = "bottom", legend.text = element_text(size = 9)) +
    geom_text_repel(
      data = subset(grs.all, ismax == T),
      aes(y = max.growth, label = otu.time),
      nudge_y = 0, box.padding = .5, max.iter = 10000,
      size = 4, color = "black", segment.alpha = .5, fontface = "italic"
    ) +
    geom_hline(yintercept = 0) +
    labs(
      x = "Time (hr)",
      y = "Rate of change in abundance"
    ) #+ scale_color_manual(values = strain.colors)
  #  y = expression(Mean~growth~rate~(hr^-1))

  return(p + scale_x_continuous(expand = c(0, 0)))
}
