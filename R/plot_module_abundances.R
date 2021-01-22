#' @title Plot Temporal Changes in Module Abundances
#'
#' @description Temporal taxa abundances are plotted for visual exploration.
#'
#' @param x A data.frame with Taxon, Module, Timepoint, Counts labels as value column
#'
#' @param tax.variable Taxa to plot
#'
#' @param mm.variable Module to plot
#'
#' @param color.pal List of hex color codes
#'
#' @param nrow When multiple values for modules and taxa
#'
#' @param ncol When multiple values for modules and taxa
#'
#' @references
#' \itemize{
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' # data(SyncomGMM)
#' # focal.species <- c("Akkermansia_muciniphila",
#' #                    "Bacteroides_xylanisolvens",
#' #                    "Bacteroides_ovatus")
#'
#' # p <- plot_module_abundances(SyncomGMM,
#' #  tax.variable = focal.species,
#' # mm.variable = c("propionate production III", "mucin degradation"),
#' # color.pal = strain.colors,
#' # nrow = 3,
#' # ncol = 1
#' # )
#' # p + theme_syncom() + theme(legend.position = "bottom")
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
#' @importFrom stats sd
#'
#' @export
#'
#' @keywords Analysis and visualization
#'


plot_module_abundances <- function(x,
                                   tax.variable = NULL, # if taxa then taxa
                                   mm.variable = NULL, # if module then taxa
                                   color.pal = NULL,
                                   nrow,
                                   ncol) {
  mod_mean <- value <- Taxon <- Module <- TimePoint <- mod_mean <- mod_ci <- NULL

  df <- filter(x, Taxon %in% tax.variable & Module %in% mm.variable)

  if (length(tax.variable) == 1 && length(mm.variable) == 1) {
    df <- df %>%
      group_by(TimePoint, Module, add = TRUE) %>%
      summarise(
        mod_mean = mean(value),
        mod_ci = 1.96 * sd(value) / sqrt(n())
      )

    mod.plot <-
      ggplot(df,
        aes(
          x = as.factor(TimePoint),
          y = mod_mean,
          group = mm.variable
        ),
        color = Module
      ) +
      geom_line(position = position_dodge(width = 0.2)) +
      geom_errorbar(
        aes(ymin = mod_mean - mod_ci, ymax = mod_mean + mod_ci),
        width = .1,
        position = position_dodge(width = 0.2),
        linetype = 1
      ) +
      geom_point(size = 2, position = position_dodge(width = 0.2))
    labs(
      title = paste(tax.variable, mm.variable, sep = " - "),
      x = "Time (hr)",
      y = "Mean module abundance"
    ) + theme_bw(base_size = 14)

    return(mod.plot)
  } else if (length(tax.variable) > 1 && length(mm.variable) == 1) {

    # Multi taxa and single module
    df <- df %>%
      group_by(TimePoint, Module, Taxon, add = TRUE) %>%
      summarise(
        mod_mean = mean(value),
        mod_ci = 1.96 * sd(value) / sqrt(n())
      )
    # head(SyncomGMM.df)

    mod.plot <- ggplot(df, aes(x = as.factor(TimePoint), y = mod_mean, group = Taxon)) +
      geom_line(position = position_dodge(width = 0.2), aes(color = Taxon)) +
      geom_errorbar(aes(ymin = mod_mean - mod_ci, ymax = mod_mean + mod_ci),
        width = .1, position = position_dodge(width = 0.2), linetype = 1
      ) +
      geom_point(size = 2, position = position_dodge(width = 0.2), aes(color = Taxon)) +
      scale_color_manual(values = color.pal) +
      labs(
        title = mm.variable,
        x = "Time (hr)",
        y = "Mean module abundance"
      ) +
      theme_bw(base_size = 12)

    return(mod.plot)
  } else if (length(tax.variable) > 1 && length(mm.variable) > 1) {
    df <- df %>%
      group_by(TimePoint, Module, Taxon, add = TRUE) %>%
      summarise(
        mod_mean = mean(value),
        mod_ci = 1.96 * sd(value) / sqrt(n())
      )
    # head(df)
    # , color = Taxon
    mod.plot <- ggplot(df, aes(x = as.factor(TimePoint), y = mod_mean, group = Taxon)) +
      geom_line(position = position_dodge(width = 0.2), aes(color = Taxon)) +
      geom_errorbar(aes(ymin = mod_mean - mod_ci, ymax = mod_mean + mod_ci),
        width = .1, position = position_dodge(width = 0.2), linetype = 1
      ) +
      geom_point(size = 2, position = position_dodge(width = 0.2), aes(color = Taxon)) +
      guides(linetype = guide_legend("Species")) +
      scale_color_manual(values = color.pal) +
      facet_wrap(~Module, nrow = nrow, ncol = ncol) +
      labs(
        title = "",
        x = "Time (hr)",
        y = "Mean module abundance"
      ) +
      theme_bw(base_size = 12)
    # Single taxa and multi module
    return(mod.plot)
  } else if (length(tax.variable) == 1 && length(mm.variable) > 1) {
    df <- df %>%
      group_by(TimePoint, Module, Taxon, add = TRUE) %>%
      summarise(
        mod_mean = mean(value),
        mod_ci = 1.96 * sd(value) / sqrt(n())
      )
    # head(df)
    ggplot(df, aes(x = as.factor(TimePoint), y = mod_mean, group = Module)) +
      geom_line(position = position_dodge(width = 0.2), aes(color = Module)) +
      geom_errorbar(aes(ymin = mod_mean - mod_ci, ymax = mod_mean + mod_ci),
        width = .1, position = position_dodge(width = 0.2), linetype = 1
      ) +
      geom_point(size = 2, position = position_dodge(width = 0.2), aes(color = Module)) +
      guides(linetype = guide_legend("Species")) +
      scale_color_manual(values = color.pal) +
      labs(
        title = tax.variable,
        x = "Time (hr)",
        y = "Mean module abundance"
      ) +
      theme_bw(base_size = 12)
  }
}
