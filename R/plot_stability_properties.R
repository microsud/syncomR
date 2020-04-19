#' @title Plot stability properties
#'
#' @description Plot various stability properties decribed by Liu et al., 2018.
#'
#' @details This script was modified to from the paper on "Ecological Stability Properties.
#' Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
#'  \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}.
#'  This functions takes output from the \code{\link{stability_properties}}.
#' \itemize{
#' \item{}{Plot reference states and community evolution}
#' }
#' @param stab.in The input come from \code{\link{stability_properties}} function.
#' @param property c("dist.ot", "resilience.ot.eucl", "resilience.ot.canb", "resilience.oline")
#' for more details check original work by Liu et al., 2018 \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}
#' \itemize{
#' \item{}{dist.ot= Plot dissimilarity over time}
#' \item{}{resilience.ot.eucl= Resilience using elucidean distance}
#' \item{}{resilience.ot.canb= Resilience using Canberra}
#' }
#'
#' @return Plots.
#'
#' @seealso Input for this functions come from \code{\link{stability_properties}}
#' @importFrom graphics matplot
#' @references
#' \itemize{
#' \item{}{Liu, Z., et al. (2018). Ecological stability properties of
#' microbial communities assessed by flow cytometry. mSphere, 3(1), e00564-17.
#' http://msphere.asm.org/content/3/1/e00564-17
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 28)
#' dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr", experimentStart = 52, tref = 152)
#' plot_stability_properties(dat.stab, property = "dist.ot")
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @export
#' @keywords Anlaysis and visualization

plot_stability_properties <- function(stab.in,
                                      property = c(
                                        "dist.ot", "resilience.ot.eucl",
                                        "resilience.ot.canb", "resilience.oline"
                                      )) {
  stab.in <- dat.stab

  titleSize <- 0.95
  numberOfGates <- stab.in$numberOfGates
  tref <- stab.in$tref
  trefId <- stab.in$trefId
  referenceState <- stab.in$referenceState
  experimentStart <- stab.in$experimentStart
  experimentEnd <- stab.in$experimentEnd
  maxEuclideanId <- stab.in$maxEuclideanId
  maxCanberraId <- stab.in$maxCanberraId
  data <- stab.in$data


  # compute radii for reference space
  r_Euclidean <- max(data$euclidean[data$referencePhase == TRUE])
  r_Canberra <- max(data$canberra[data$referencePhase == TRUE])


  maxDevStateEuclidean <- data[maxEuclideanId, 2:(numberOfGates + 1)]
  maxDevStateCanberra <- data[maxCanberraId, 2:(numberOfGates + 1)]

  maxDevStateEuclideanTime <- data[maxEuclideanId, 1]
  maxDevStateCanberraTime <- data[maxCanberraId, 1]
  labelStart <- "Relative abundance in"

  # computing stability properties
  RSeuclidean <- 1.0 - data[maxEuclideanId, "euclidean"]
  RScanberra <- 1.0 - data[maxCanberraId, "canberra"]
  DisSpeedeuclidean <- data[maxEuclideanId, "euclidean"] / (data[maxEuclideanId, 1] - tref)
  DisSpeedcanberra <- data[maxCanberraId, "canberra"] / (data[maxCanberraId, 1] - tref)
  ElasticityEuclidean <- (data[maxEuclideanId, "euclidean"] - data[nrow(data), "euclidean"]) / (data[nrow(data), 1] - data[maxEuclideanId, 1])
  ElasticityCanberra <- (data[maxCanberraId, "canberra"] - data[nrow(data), "canberra"]) / (data[nrow(data), 1] - data[maxCanberraId, 1])

  # plots over time
  plotData <- data[data[1] >= tref, ]
  plotData[1, 2:(numberOfGates + 1)] <- referenceState
  plotData$euclidean[1] <- 0
  plotData$canberra[1] <- 0
  plotData$RLcanberraOnline[1] <- 0

  if (property == "dist.ot") {

    # Distance over time
    matplot(plotData[[1]], cbind(plotData$euclidean, plotData$canberra), ylim = c(0, 1), main = "Resistance and displacement speed\n computing", ylab = expression("Deviation (d) from s"[ref]), xlab = "Time", lty = c(1, 1), pch = c(21, 24), col = c("deepskyblue3", "brown3"), type = "l", font.main = 1, cex.main = titleSize)
    abline(h = c(r_Euclidean, r_Canberra), col = c("deepskyblue3", "brown3"), lty = c(2, 2))

    points(plotData$Time, y = plotData$euclidean, pch = 16, col = "deepskyblue3")
    points(data[maxEuclideanId, 1], y = data[maxEuclideanId, "euclidean"], col = "deepskyblue3", cex = 1.2, pch = 24, bg = "white")
    points(experimentEnd, y = plotData[nrow(plotData), "euclidean"], col = "black", cex = 1.3, pch = 21, bg = "black")
    text(plotData$Time, y = plotData$euclidean, labels = plotData$Time, cex = 0.6, font = 2, pos = 4)
    # spread.labels(plotData$Time, y=plotData$euclidean, labels=plotData$Time,linecol="deepskyblue3",cex=0.)

    points(plotData$Time, y = plotData$canberra, pch = 16, col = "brown3")
    points(data[maxCanberraId, 1], y = data[maxCanberraId, "canberra"], col = "brown3", cex = 1.2, pch = 24, bg = "white")
    points(experimentEnd, y = plotData[nrow(plotData), "canberra"], col = "black", cex = 1.3, pch = 21, bg = "black")
    text(plotData$Time, y = plotData$canberra, labels = plotData$Time, cex = 0.6, font = 2, pos = 4)
    # spread.labels(plotData$Time, y=plotData$canberra, labels=plotData$Time, linecol = "brown3", cex=0.6)

    points(plotData$Time[1], plotData$euclidean[1], col = "red", pch = 21, bg = "white", cex = 1.3)
    points(plotData$Time[1], plotData$canberra[1], col = "red", pch = 21, bg = "white", cex = 1.3)
    legend("topright", c(paste("Euclidean, RS=", toString(round(RSeuclidean, digits = 4)), ", DS=", toString(round(DisSpeedeuclidean, digits = 4)), sep = ""), paste("Canberra, RS=", toString(round(RScanberra, digits = 4)), ", DS=", toString(round(DisSpeedcanberra, digits = 4)), sep = "")), pch = c(19, 19), col = c("deepskyblue3", "brown3"), cex = 0.7)
  } else if (property == "resilience.ot.eucl") {
    plot(data[maxEuclideanId:nrow(data), 1], y = data$RLeuclidean[maxEuclideanId:nrow(data)], ylim = c(0, 1), font.main = 1, main = "Resilience and elasticity computing\n(based on Euclidean distance)", ylab = "Resilience (RL)", xlab = "Time", xlim = c(plotData$Time[1], experimentEnd), type = "l", col = "deepskyblue3", cex.main = titleSize)
    abline(h = 0, col = "black", lty = 2)
    points(data[maxEuclideanId:nrow(data), c(colnames(data)[1], "RLeuclidean")], pch = 19, col = "deepskyblue3")
    points(data[maxEuclideanId, 1], y = data[maxEuclideanId, "RLeuclidean"], col = "deepskyblue3", cex = 1.2, pch = 24, bg = "white")
    points(experimentEnd, y = plotData[nrow(plotData), "RLeuclidean"], col = "black", cex = 1.3, pch = 21, bg = "black")
    legend("topright", paste("RL=", toString(round(plotData$RLeuclidean[nrow(plotData)], digits = 4)), ", E=", toString(round(ElasticityEuclidean, digits = 4))), cex = 0.7)
  } else if (property == "resilience.ot.canb") {
    plot(data[maxCanberraId:nrow(data), 1], y = data$RLcanberra[maxCanberraId:nrow(data)], col = "brown3", ylim = c(0, 1), font.main = 1, main = "Resilience and elasticity computing\n(based on Canberra distance)", ylab = "Resilience (RL)", xlab = "Time", xlim = c(plotData$Time[1], experimentEnd), type = "l", cex.main = titleSize)
    abline(h = 0, col = "black", lty = 2)

    points(data[maxCanberraId:nrow(data), c(colnames(data)[1], "RLcanberra")], pch = 19, col = "brown3")
    points(data[maxCanberraId, 1], y = data[maxCanberraId, "RLcanberra"], col = "brown3", cex = 1.2, pch = 24, bg = "white")
    # text(plotData$Time, y=plotData$RLcanberraOnline, labels=plotData$Time, cex=0.6, font=2, pos=4)
    points(experimentEnd, y = plotData[nrow(plotData), "RLcanberra"], col = "black", cex = 1.3, pch = 21, bg = "black")
    legend("topright", paste("RL=", toString(round(plotData$RLcanberra[nrow(plotData)], digits = 4)), ", E=", toString(round(ElasticityCanberra, digits = 4))), cex = 0.7)
  } else if (property == "resilience.oline") {

    # Online resilience
    plot(plotData[[1]], y = plotData$RLcanberraOnline, ylim = c(0, 1), main = "Online resilience computing\n(based on Canberra distance)", ylab = "Resilience", xlab = "Time", col = "black", pch = 21, bg = "white", font.main = 1, type = "l", cex.main = titleSize)
    abline(h = 0, col = "black", lty = 2)
    text(plotData$Time, y = plotData$RLcanberraOnline, labels = plotData$Time, cex = 0.6, font = 2, pos = 4)
    points(plotData[[1]], y = plotData$RLcanberraOnline, pch = 21, bg = "white", col = "black")
    points(plotData$Time[1], plotData$RLcanberraOnline[1], col = "red", pch = 21, bg = "white", cex = 1.3)
  }
}
