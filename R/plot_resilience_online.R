#' @title Plot Stability Properties Online Resilience
#'
#' @description Online resilience of the community as described by
#'              Florian Centler and Zishu Liu, UFZ Leipzig, Germany.
#'
#' @details This script was modified to get ggplot object from the paper on
#'          "Ecological Stability Properties.Microbial Communities Assessed by
#'          Flow Cytometry" by Liu et al., 2018
#'          \url{http://msphere.asm.org/content/3/1/e00564-17}. This functions takes
#'          output from the \code{\link{stability_properties}}. The closer it is to zero
#'          the system is deviating continuously away from reference state, while values
#'          closer to 1 indicate the onset of recovery (1 = perfect recovery).
#'
#' @param stab.in The input come from \code{\link{stability_properties}} function.
#'
#' @param col.low color aesthetics.
#'
#' @param col.high color aesthetics.
#'
#' @return ggplot object.
#'
#' @seealso Input for this functions comes from \code{\link{stability_properties}}
#'
#' @references
#' \itemize{
#' \item{}{Liu, Z., et al. (2018). Ecological stability properties of
#'         microbial communities assessed by flow cytometry. \emph{mSphere}, 3(1), e00564-17.
#' \url{http://msphere.asm.org/content/3/1/e00564-17}
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @examples
#'\dontrun{
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 172)
#' dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr")
#' plot_resilience_online(dat.stab, col.low = "#fa9fb5", col.high = "#49006a")
#' }
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @export
#' @keywords Analysis and visualization
#'
plot_resilience_online <- function(stab.in,
                                   col.low = "#fa9fb5",
                                   col.high = "#49006a") {
  Time <- RLcanberraOnline <- NULL
  # titleSize <- 0.95
  tref <- stab.in$tref
  trefId <- stab.in$trefId
  # experimentEnd <- stab.in$experimentEnd
  # maxEuclideanId <- stab.in$maxEuclideanId
  data <- stab.in$data


  # plots over time
  plotData <- data[data[1] >= tref, ]

  # set tref RLcanberraOnline to zero
  plotData$RLcanberraOnline[1] <- 0

  dist.index <- plotData[, c("Time", "RLcanberraOnline")]

  dist.index$tref <- ifelse(dist.index$Time == tref, "tref", "NA")

  textdf <- dist.index[dist.index[1] == tref, ]


  # df.1
  p <- ggplot(dist.index, aes(Time, RLcanberraOnline)) +
    theme_syncom() +
    geom_line(size = 0.6, alpha = 0.8) +
    ylim(c(0, 1)) +
    geom_point(aes(fill = RLcanberraOnline, shape = tref), size = 3) +
    scale_fill_continuous("Resilience", low = col.low, high = col.high) +
    scale_shape_manual(values = c(tref = 24, "NA" = 21), guide = FALSE)
  p <- p + geom_text_repel(data = textdf, label = tref, vjust = 1)
  p <- p + geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(expression(~ Resilience[Online - Canberra]))

  return(p)
}
