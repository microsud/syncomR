#' @title Short-term changes in abundance of taxa from temporal abundance data
#'
#' @description Calculates short-term changes in abundance of taxa using temporal abundnace data.
#'
#' @details This approach is used by Nathan I. Wisnoski and collegues
#' \link{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on the following calculation
#' log(present abundance/past abundance). Also a compositional version using relative abundance similar to
#' Brian W. Ji, Ravi U. Sheth et al., 2020 \link{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#' This approach is useful for identifying short term growth behaviours of taxa.
#' @param ps a non-rarefied raw otu/ASV/species abundance table (longitudinal data).
#' @param rarefy TRUE  Kept this option incase users have absolute quantification.
#' @param depth min(sample_sums(ps) If rarefy=TRUE then depth for rarefaction.
#' @param plot.type If data to be ploted either line or polar else NULL.
#' @param abbreviation For plotting whether to abbreviate names of taxa TRUE or FALSE.
#' @return A data frame with mean growth, max growth or a plot with labelling for max "change in abundance" timepoints.
#' @references
#' \itemize{
#' \item{}{Shetty SA. et al. A Minimalist Approach for Deciphering the Ecophysiology of Human Gut Microbes (2020)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' p.g <- short_term_change(ps1.b5, rarefy = TRUE, depth = min(sample_sums(ps1.b5)), plot.type = "polar")
#' p.g <- p.g + scale_color_manual(values = syncom_colors2("BacterialSpecies")) + ggtitle("Bioreactor A")
#' print(p.g)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @import tidyr
#' @import ggrepel
#' @import tibble
#' @export
#' @keywords Anlaysis and visualization

short_term_change <- function(ps,
                              rarefy = TRUE,
                              depth = min(sample_sums(ps)),
                              plot.type = NULL,
                              abbreviation = T) {
  message("Calculating short term growth behaviours of taxa")

  ps.1 <- otu.tb <- grs.all <- grts <- otu.rar <- maxgrs <- p <- NULL

  if (rarefy == TRUE) {
    ps.1 <- rarefy_even_depth(ps, sample.size = depth)
  } else {
    ps.1 <- microbiome::transform(ps, "compositional")
  }

  otu.tb <- (taxa_time_table(ps.1, normalize = F, time.col = "Time_hr", remove.zero = F) + .1)
  otu.tb$time <- rownames(otu.tb)

  otu.tb.ldf <- reshape2::melt(otu.tb)
  # head(otu.tb)
  otu.tb.ldf$time <- as.numeric(otu.tb.ldf$time)
  # otu.tb.ldf$value <- (otu.tb.ldf$value+0.1)
  # str(otu.tb.ldf)
  grwt <- as_tibble(otu.tb.ldf) %>%
    arrange(time) %>%
    group_by(variable) %>%
    mutate(
      time_lag = time - lag(time), # time lag since sampling is not always equal
      growth_diff = log(value / lag(value)),
      # growth_diff = log(value - lag(value)/lag(value)), # log(present abund-past abund/past abund)
      # growth = log(value - lag(value)/Diff_timepoint)+.1,
      growth_rate = (growth_diff / time_lag) / lag(value),
      # Mortality = log(lead(value) - value/Diff_timepoint)+.1
    )
  grwt$sample.name <- grwt$time
  # colnames(grwt)["variable"]
  colnames(grwt)[colnames(grwt) == "variable"] <- "OTU"

  grs.all <- grwt %>%
    # gather(-"sample.name", key = "OTU", value = "growth_rate") %>%
    # mutate(Time = sample.name) %>%
    group_by(OTU, time) %>%
    # replace_na(growth_diff, "NaN") %>%
    summarize(
      mean.growth = mean(growth_diff),
      sd.growth = sd(growth_diff)
    )
  # head(grs.all)

  maxgrs <- grs.all %>%
    summarize(max.growth = max(mean.growth, na.rm = T))
  grs.all <- grs.all %>%
    left_join(maxgrs)
  grs.all <- grs.all %>%
    mutate(ismax = ifelse(mean.growth == max.growth, T, F))
  # DT::datatable(grs.all)
  grs.all$OTU <- gsub("_", " ", grs.all$OTU)
  if (abbreviation) {
    grs.all$OTUabb <- toupper(abbreviate(grs.all$OTU,
      minlength = 3,
      method = "both.sides"
    ))
    grs.all$otu.time <- paste0(grs.all$OTUabb, " ", grs.all$time, "h")
  } else {
    grs.all$otu.time <- paste0(grs.all$OTU, " ", grs.all$time, "h")
  }
  # head(persistent.grs.all)
  # grs.all$Time <- as.numeric(grs.all$time)


  if (is.null(plot.type) == TRUE) {
    return(grs.all)
  } else if (plot.type == "line") {
    p <- ggplot(grs.all, aes(x = time, group = OTU, col = OTU)) +
      geom_line(aes(y = mean.growth), alpha = 0.6, size = 1) +
      geom_point(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth), alpha = 0.8, size = 3
      ) +
      # coord_flip() +
      geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = min(grs.all$mean.growth)),
        fill = "#edf3f5", col = "white", alpha = 0.5
      ) +
      theme(
        legend.position = "top", legend.text = element_text(size = 9),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank()
      ) +
      geom_text_repel(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth, label = otu.time),
        nudge_y = 0.2, box.padding = .5, max.iter = 10000,
        size = 3, color = "black", segment.alpha = .5,
        fontface = "italic", direction = "both"
      ) +
      geom_hline(yintercept = 0) +
      labs(
        x = "Time (hr)",
        y = expression(paste("Change in abundance ", " µ= ", abund[t + 1] / abund[t]))
      )
    #+ scale_color_manual(values = strain.colors)
    #  y = expression(Mean~growth~rate~(hr^-1))
    # p <- p +scale_x_continuous(expand = c(0, 0))
    return(p + scale_x_continuous(expand = c(0, 0)))
  } else if (plot.type == "polar") {
    p <- ggplot(grs.all, aes(x = time, group = OTU, col = OTU)) +
      geom_line(aes(y = mean.growth), alpha = 0.6, size = 1) +
      geom_point(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth), alpha = 0.8, size = 3
      ) +
      geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = -9),
        fill = "#edf3f5", col = "white", alpha = 0.5
      ) +
      coord_polar(theta = "x") +
      theme(
        legend.position = "right", legend.text = element_text(size = 9),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(
          size = 0.5, linetype = "solid",
          colour = "#CCD1D1"
        ),
        panel.grid.minor = element_line(
          size = 0.5, linetype = "solid",
          colour = "#CCD1D1"
        ),
        legend.key = element_blank()
      ) +
      geom_text_repel(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth, label = otu.time),
        nudge_y = 0.2, box.padding = .5, max.iter = 10000,
        size = 3, color = "black", segment.alpha = .5,
        fontface = "italic", direction = "both"
      ) +
      geom_hline(yintercept = 0) +
      labs(
        x = "Time (hr)",
        y = expression(paste("Change in abundance ", " µ= ", abund[t + 1] / abund[t]))
      )
    #+ scale_color_manual(values = strain.colors)
    #  y = expression(Mean~growth~rate~(hr^-1))
    # p <- p +scale_x_continuous(expand = c(0, 0))
    return(p + scale_x_continuous(expand = c(0, 0)))
  }
}
