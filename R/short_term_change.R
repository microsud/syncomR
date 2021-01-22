#' @title Short term Changes in Abundance
#'
#' @description Calculates short term changes in abundance of taxa
#'              using temporal Abundance data.
#'
#' @details This approach is used by Wisnoski NI and colleagues
#'          \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
#'          the following calculation log(present abundance/past abundance).
#'          Also a compositional version using relative abundance similar to
#'          Brian WJi, Sheth R et al
#'          \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#'          This approach is useful for identifying short term growth behaviors of taxa.
#'
#' @param ps a non rarefied raw OTU or ASV or species abundance table longitudinal data.
#'
#' @param rarefy TRUE  Kept this option in case users have absolute quantification.
#'
#' @param depth If rarefy is TRUE then depth for rarefaction.
#'
#' @param compositional Logical TRUE or FALSE.
#'
#' @param plot.type If data to be plotted either line or polar else NULL.
#'
#' @param abbreviation For plotting whether to abbreviate names of taxa TRUE or FALSE.
#'
#' @return A data frame with mean growth, max growth or a plot with labeling
#'         for max \emph{change in abundance} timepoints.
#'
#' @references
#'
#' \itemize{
#' \item{}{Shetty SA. et al. A Minimalist Approach for Deciphering the
#'         Ecophysiology of Human Gut Microbes 2020}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#'\dontrun{
#' data(SyncomFiltData)
#' short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
#' syncom_ps <- subset_samples(SyncomFiltData, !(Time_label %in% short_time_labels))
#' bioA <- subset_samples(syncom_ps, StudyIdentifier == "Bioreactor A")
#' bioA.lg <- add_time_lag(bioA)
#' bioA.lg <- subset_samples(bioA.lg, time_lag >= 4)
#' p.bioA <- short_term_change(bioA.lg,
#'   rarefy = TRUE,
#'   depth = min(sample_sums(bioA.lg)),
#'   plot.type = "polar"
#' )
#' p.g <- p.bioA +
#'   scale_color_manual(values = syncom_colors2("BacterialSpecies")) +
#'   ggtitle("Bioreactor A")
#' print(p.g)
#' }
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @import tidyr dplyr
#'
#' @importFrom ggrepel geom_label_repel geom_text_repel
#'
#' @export
#' @keywords Analysis and visualization

short_term_change <- function(ps,
                              rarefy = FALSE,
                              compositional = FALSE,
                              depth = NULL,
                              plot.type = NULL,
                              abbreviation = T) {
  # message("Calculating short term growth behaviours of taxa")

  if (is.null(depth) || depth > min(phyloseq::sample_sums(ps))) {
    stop("Depth cannot be NULL or more than
         minimun number of reads in your data")
  }

  time <- variable <- value <- time_lag <- growth_diff <- NULL
  max.growth <- OTU <- ismax <- otu.time <- NULL

  ps.1 <- otu.tb <- grs.all <- grts <- otu.rar <- maxgrs <- p <- NULL

  if (rarefy == TRUE & compositional == FALSE) {
    message("rarefy is set to TRUE, calculating short term change using counts")
    ps.1 <- phyloseq::rarefy_even_depth(ps, sample.size = depth)
  } else if (rarefy == FALSE & compositional == FALSE) {
    message("rarefy is set to FALSE, compositional==FALSE
            calculating short term change using raw counts as provided by user")
    ps.1 <- ps
  } else if (rarefy == FALSE & compositional == TRUE) {
    message("rarefy is set to FALSE, compositional==TRUE,
            using relative abundances")
    ps.1 <- microbiome::transform(ps, "compositional")
  } else if (rarefy == TRUE & compositional == TRUE) {
    stop(message("Both rRarefy and compositional cannot be TRUE"))
  }

  otu.tb <- (taxa_time_table(ps.1,
    normalize = F,
    time.col = "Time_hr",
    remove.zero = F
  ) + .1)

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
      growth_rate = log(value - lag(value) / lag(value)),
      var_abund = value - lag(value) / time_lag
      # growth_diff = log(value - lag(value)/lag(value)), # log(present abund-past abund/past abund)
      # growth = log(value - lag(value)/Diff_timepoint)+.1,
      # Mortality = log(lead(value) - value/Diff_timepoint)+.1
    )

  grwt$sample.name <- grwt$time
  # colnames(grwt)["variable"]
  colnames(grwt)[colnames(grwt) == "variable"] <- "OTU"

  maxgrs <- grwt %>%
    summarize(max.growth = max(growth_diff, na.rm = T))
  colnames(maxgrs)[colnames(maxgrs) == "variable"] <- "OTU"
  grs.all <- grwt %>%
    left_join(maxgrs)
  grs.all <- grs.all %>%
    mutate(ismax = ifelse(growth_diff == max.growth, T, F))
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
      geom_line(aes(y = growth_diff), alpha = 0.6, size = 1) +
      geom_point(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth), alpha = 0.8, size = 3
      ) +
      # coord_flip() +
      geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = min(grs.all$growth_diff)),
        fill = "white", col = "#edf3f5", alpha = 0.7
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
        y = expression(paste("Change in abundance ", " \U00B5 = ", abund[t + 1] / abund[t]))
      )
    #+ scale_color_manual(values = strain.colors)
    #  y = expression(Mean~growth~rate~(hr^-1))
    # p <- p +scale_x_continuous(expand = c(0, 0))
    return(p + scale_x_continuous(expand = c(0, 0)))
  } else if (plot.type == "polar") {
    p <- ggplot(grs.all, aes(x = time, group = OTU, col = OTU)) +
      geom_line(aes(y = growth_diff), alpha = 0.6, size = 1) +
      geom_point(
        data = subset(grs.all, ismax == T),
        aes(y = max.growth), alpha = 0.8, size = 3
      ) +
      geom_ribbon(aes(group = NULL, col = NULL, ymax = 0, ymin = -9),
        fill = "#edf3f5", col = "#edf3f5", alpha = 0.5
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
        y = expression(paste("Change in abundance ", " \U00B5 = ", abund[t + 1] / abund[t]))
      )
    #+ scale_color_manual(values = strain.colors)
    #  y = expression(Mean~growth~rate~(hr^-1))
    # p <- p +scale_x_continuous(expand = c(0, 0))
    return(p + scale_x_continuous(expand = c(0, 0)))
  }
}
