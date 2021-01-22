
context("test plot_growth_rate")

test_that("plot_growth_rate works correctly", {
  library(syncomR)
  data(SyncomFiltData)
  ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
  p.g <- plot_growth_rate(ps1.b5, rarefy = TRUE, depth = min(sample_sums(ps1.b5)))
  dat <- p.g$data
  expect_equal(dat$mean.growth[3], 6.27, tolerance = 0.03)
})


#' @title Utility for stability calculations
#' @description Function for plotting community evolution with respect to two taxa
#' @param gates defunct
#' @param title defunct
#' @param xlabel defunct
#' @param ylabel defunct
#' @param ref_radius defunct
#' @seealso Used internally by \code{\link{stability_properties}}
#'
#' @references
#' \itemize{
#' \item{}{"Liu, Z., Cichocki, N., Bonk, F., Günther, S., Schattenberg, F., Harms, H., ... & Müller, S. (2018). Ecological stability properties of
#' microbial communities assessed by flow cytometry. mSphere, 3(1), e00564-17.
#' http://msphere.asm.org/content/3/1/e00564-17
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @keywords Utilities
NULL
plot_2D_projection <- function(gates, title, xlabel, ylabel, ref_radius) {
  # plot states
  plot(data[data[1] >= experimentStart & data[1] <= experimentEnd, gates], col = c("black", "grey")[as.factor(data[data[1] >= experimentStart & data[1] <= experimentEnd, "referencePhase"])], font.main = 1, main = title, xlab = xlabel, ylab = ylabel, cex.main = titleSize)
  # add lines to be able to follow dynamic evolution of the state
  lines(data[data[1] > tref & data[1] <= experimentEnd, gates])
  # ZS: give the direction of evolution from reference state to disturbed state

  # mark tref
  # points(data[data[1]==tref, gates], pch=16, col="lightgreen")
  # mark final state
  points(data[data[1] == experimentEnd, gates], pch = 19, cex = 1.3)
  # mark the reference state
  points(referenceState[gates[1]], y = referenceState[gates[2]], col = "red", pch = 21, bg = "white", cex = 1.3)
  # mark states of maximal deviation from reference state
  points(maxDevStateEuclidean[gates[1]], y = maxDevStateEuclidean[gates[2]], col = "deepskyblue3", pch = 24, cex = 1.2, bg = "white")
  points(maxDevStateCanberra[gates[1]], y = maxDevStateCanberra[gates[2]], col = "brown3", pch = 24, cex = 1.2, bg = "white")
  arrows(referenceState[gates[1]], referenceState[gates[2]], data[, gates[1]][trefId + 1], data[, gates[2]][trefId + 1], length = 0.08, angle = 8, lwd = 1)
  # mark reference space defining
  draw.circle(referenceState[gates[1]], y = referenceState[gates[2]], radius = ref_radius, border = "grey", lty = 2)
}
