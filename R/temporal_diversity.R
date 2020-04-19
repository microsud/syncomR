#' @title Measure and plot temporal changes in eveness of community
#'
#' @description Calculate and plot diversity measures over time.
#'
#' @details Can be useful to identify potential trends in temporal community diversities.
#' @param ps a otu/ASV/species abundance table
#' @param time.col Column specifying time variable
#' @param div.measure Currently only Gini as a measure of inequality is used
#' @param plot TRUE or FALSE if TRUE plot is returned else a data frame
#' @return Temporal plot of inquality.
#' @references
#' \itemize{
#' \item{}{To cite the package, see citation('microbiome')}
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' p <- temporal_diversity(SyncomFiltData, time.col = "Time_hr_num", div.measure = "gini")
#' p + facet_grid(~StudyIdentifier) + geom_smooth(fill = "#a6bddb", color = "#d94801")
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @import dplyr
#' @export
#' @keywords Anlaysis and visualization
#'
temporal_diversity <- function(ps, time.col = "Time_hr_num",
                               div.measure = "gini", plot = TRUE) {
  even_df <- as.data.frame(microbiome::inequality(ps))
  colnames(even_df) <- "Inequality"
  meta_df <- meta(ps)
  # head(meta_df)

  metadf <- cbind(meta_df, even_df)
  metadf$time <- metadf[, time.col]
  metadf <- metadf %>% arrange(as.numeric(time))
  # metadf$StudyIdentifier <- gsub("_", " ", metadf$StudyIdentifier)
  if (plot == TRUE) {
    p.ev <- ggplot(metadf, aes(x = time, Inequality))
    p.ev <- p.ev + geom_point() + geom_line() #+ scale_color_manual(values = color.pal)
    p.ev <- p.ev + xlab("Time (hr)")
    return(p.ev)
  } else {
    return(metadf)
  }
  # , color = StudyIdentifier
}
