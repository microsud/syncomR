#' @title Line Plot Species
#'
#' @description Temporal community composition using line plot.
#'
#' @details Temporal community composition using line plot
#' @param ps A \code{\link{phyloseq-class}}
#' @param type  Type of plot either area, line or bar
#' @param group.var Grouping variable like bioreactor
#' @param time.col Column specifying time variable
#' @param taxa.level Taxonomic level
#' @param ... For facet_wrap
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
#' }
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @importFrom phyloseq psmelt sample_variables
#' @importFrom ggplot2 geom_line facet_wrap geom_pointrange
#' @importFrom rlang sym
#' @export
#' @keywords Anlaysis and visualization

plot_species_lineplot <- function(ps,
                                  type = c("grouped", "individual"),
                                  group.var = NULL,
                                  time.col = "Time_hr",
                                  taxa.level = "OTU", ...) {

  ps1.sub.df <- comp.plt <- Time <- mean.ab <- Abundance <- NULL

  ps1.sub.df <- phyloseq::psmelt(ps)

  ps1.sub.df$Time <- as.numeric(ps1.sub.df[, time.col])

  ps1.sub.df$Time.ch <- as.factor(ps1.sub.df[, time.col])

  comp.plt <- ggplot(ps1.sub.df)

  if(!is.null(group.var) && group.var %in% sample_variables(ps))

    if(type=="grouped"){
      comp.plt <- ps1.sub.df %>%
        dplyr::group_by(!!sym(taxa.level), Time) %>%
        dplyr::summarise(mean.ab=mean(Abundance, na.rm=TRUE),
                         sd.ab=sd(Abundance, na.rm=TRUE)) %>%
        ggplot(aes(x=Time, y=mean.ab)) +
        geom_pointrange(aes_string(ymin = "mean.ab" - "sd.ab",
                                   ymax = "mean.ab" + "sd.ab",
                                   color=taxa.level)) +
        # color="#69b3a2", size=2, alpha=0.9, linetype=2
        geom_line(aes_string(color=taxa.level)) +
        theme_syncom()
      return(comp.plt)
    }

  if(type=="individual"){
    var.facet <- sym(group.var)

    comp.plt <- ps1.sub.df %>%
      ggplot(aes(x=Time, y=Abundance)) +
      geom_point(aes_string(color=taxa.level)) +
      # color="#69b3a2", size=2, alpha=0.9, linetype=2
      facet_wrap(var.facet, ...) +
      geom_line(aes_string(color= taxa.level)) +
      theme_syncom()
    return(comp.plt)
  }

}
