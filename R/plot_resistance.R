#' @title Plot stability properties Resistance
#'
#' @description Resistance of the community overtime as decribed by Florian Centler and Zishu Liu, UFZ Leipzig, Germany.
#'
#' @details This script was modified to get ggplot object from the paper on "Ecological Stability Properties.
#' Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
#'  \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}.
#'  This functions takes output from the \code{\link{stability_properties}}.
#' @param stab.in The input come from \code{\link{stability_properties}} function.
#' @return ggplot object.
#' @seealso Input for this functions comes from \code{\link{stability_properties}}
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
#' ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 172)
#' dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr")
#' plot_resistance(dat.stab)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @export
#' @keywords Anlaysis and visualization
#'
plot_resistance <- function(stab.in) {


  # titleSize <- 0.95
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


  # maxDevStateEuclidean <- data[maxEuclideanId, 2:(numberOfGates + 1)]
  # maxDevStateCanberra <- data[maxCanberraId, 2:(numberOfGates + 1)]

  # maxDevStateEuclideanTime <- data[maxEuclideanId, 1]
  # maxDevStateCanberraTime <- data[maxCanberraId, 1]
  # labelStart <- "Relative abundance in"

  # computing stability properties
  RSeuclidean <- 1.0 - data[maxEuclideanId, "euclidean"]
  RScanberra <- 1.0 - data[maxCanberraId, "canberra"]
  DisSpeedeuclidean <- data[maxEuclideanId, "euclidean"] / (data[maxEuclideanId, 1] - tref)
  DisSpeedcanberra <- data[maxCanberraId, "canberra"] / (data[maxCanberraId, 1] - tref)

  # plots over time
  plotData <- data[data[1] >= tref, ]
  plotData[1, 2:(numberOfGates + 1)] <- referenceState
  plotData$euclidean[1] <- 0
  plotData$canberra[1] <- 0
  plotData$RLcanberraOnline[1] <- 0

  dist.index <- plotData[, c("Time", "euclidean", "canberra")]
  df.1 <- dist.index %>% gather(dist.index, Distance, -Time)

  df.1$max.dis <- ifelse(df.1$Time == data[maxEuclideanId, 1] & df.1$dist.index == "euclidean", "MaxDeviation",
    ifelse(df.1$Time == data[maxCanberraId, 1] & df.1$dist.index == "canberra", "MaxDeviation",
      ifelse(df.1$Time == experimentEnd, "ExperimentEnd", "No")
    )
  )

  colnames(df.1) <- c("Time", "DistanceMethod", "Distance", "Variable")
  df.1$DistanceMethod <- as.factor(df.1$DistanceMethod)
  # df.1
  p <- ggplot(df.1, aes(Time, Distance)) +
    geom_line(aes(color = DistanceMethod), size = 0.6, alpha = 0.8) +
    ylim(c(0, 1)) +
    geom_point(aes(color = DistanceMethod, shape = Variable), size = 2) +
    scale_color_manual(values = c(
      euclidean = "deepskyblue3",
      canberra = "brown3"
    ))
  p <- p + scale_shape_manual(values = c(
    MaxDeviation = 2,
    No = 16,
    ExperimentEnd = 0
  ), guide = FALSE)
  lab.txt <- paste(
    " Euclidean, RS=", round(RSeuclidean, digits = 4), "| DS=", round(DisSpeedeuclidean, digits = 4), "\n",
    "Canberra, RS=", round(RScanberra, digits = 4), "| DS=", round(DisSpeedcanberra, digits = 4)
  )
  p <- p + geom_hline(aes(yintercept = r_Euclidean),
    colour = "deepskyblue3", linetype = "dashed"
  ) + ylab("Resistance (RS)") +
    geom_hline(aes(yintercept = r_Canberra),
      colour = "brown3", linetype = "dashed"
    ) + ggtitle(label = "", subtitle = lab.txt)
  return(p)
}
