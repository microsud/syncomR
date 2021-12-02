#' @title Plot temporal changes in taxa composition
#'
#' @description Temporal community composition.
#'
#' @details Temporal community composition
#' @param ps A \code{\link{phyloseq-class}}
#' @param type  Type of plot either area or bar
#' @param time.col Column specifying time variable
#' @param taxa.level Taxonomic level
#' @param sp.fill.pal List of colors for plotting
#' @param facet.var Column specifying varible to facet by. Here, its by Bioreactors
#' @param ... options to pass for ggplot2 facet_
#' @return ggplot object.
#' @import ggpubr
#' @references
#' \itemize{
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' \dontrun{
#' data(SyncomFiltData)
#' SyncomFiltData.rel <- microbiome::transform(SyncomFiltData, "compositional")
#' strain.colors <- syncom_colors("BacterialSpecies")[1:16]
#' p <- plot_syncom_composition(SyncomFiltData.rel,
#'   type = "bar",
#'   time.col = "Time_hr_num",
#'   taxa.level = "OTU",
#'   sp.fill.pal = strain.colors, facet.var = "StudyIdentifier"
#' )
#' print(p)
#' }
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @importFrom phyloseq psmelt
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_area
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 facet_wrap
#' @export
#' @keywords Anlaysis and visualization

plot_syncom_composition <-
  function(ps,
           type = c("area", "bar"),
           time.col = "Time_hr",
           taxa.level = "OTU",
           sp.fill.pal = NULL,
           facet.var = NULL, ...) {
    ps1.sub.df <- comp.plt <- NULL
    ps1.sub.df <- psmelt(ps)

    ps1.sub.df$Time <- as.numeric(ps1.sub.df[, time.col])

    ps1.sub.df$Time.ch <- as.factor(ps1.sub.df[, time.col])


    comp.plt <- ggplot(ps1.sub.df)
    if (type == "bar") {

      comp.plt <- ggplot(ps1.sub.df)
      comp.plt <-
        comp.plt + geom_bar(aes_string(x = "Time.ch",
                                       y = "Abundance",
                                       fill = taxa.level),
                            position = "stack",
                            stat = "identity") +
        scale_fill_manual("Taxa", values = sp.fill.pal)

    } else if (type == "area") {
      comp.plt <- ggplot(ps1.sub.df)
      comp.plt <- comp.plt + geom_area(aes_string(x = "Time",
                                                  y = "Abundance",
                                                  fill = taxa.level)) +
        scale_fill_manual("Taxa", values = sp.fill.pal)
    }

    guide_italics <- guides(fill = guide_legend(label.theme =
                                                  element_text(size = 15,
                                                               face = "italic",
                                                               colour = "Black",
                                                               angle = 0)))

    comp.plt <- comp.plt +
      theme_minimal() +
      # scale_x_datetime(date_labels="Day %d") +
      ylab("Proportions")

    if (!is.null(facet.var)) {
      varf <- sym(facet.var)
      comp.plt <- comp.plt + facet_wrap(varf, nrow = 3)
    }

    comp.plt <- comp.plt +
      # guides(fill = element_text(size=16))+
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = rel(1), colour = "grey30"),
        axis.title.y = element_text(size = rel(1), colour = "grey30"),
        axis.text.y = element_text(size = rel(1), colour = "grey30"),
        strip.text.x = element_text(size = rel(1), colour = "grey30"),
        legend.text = element_text(size = rel(1), colour = "grey30", face = "italic"),
        legend.title = element_text(size = rel(1), colour = "grey30")
      ) + guide_italics
    return(comp.plt)
  }
