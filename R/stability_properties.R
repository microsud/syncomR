#' @title Comunity-level stability properties
#'
#' @description Calculate community stability properties as decribed by Florian Centler and Zishu Liu, UFZ Leipzig, Germany.
#'
#' @details This script was modified from the paper on "Ecological Stability Properties.
#' Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
#'  \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}.
#'  Please cite this work when using this function. A value of -1 will let the script select these time points automatically,
#'  assuming that the file contains a single disturbance experiment (i.e. start and
#'  end of the experiment refer to the first and last entry, respectively). Otherwise, specify a time value which is actually present in the data, i.e.
#'  a value which appears in the first column of the input file.
#'
#' \itemize{
#' \item{}{Identify reference states}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @param ps A \code{\link{phyloseq-class}}. The input \code{\link{phyloseq}} object must contain a single disturbance experiment (i.e. start and end of the
#' experiment refer to the first and last entry, respectively).
#' @param experimentStart Start of the experiment (experimentStart). A value of -1 will let
#' the script select these time points automatically. Otherwise, specify a time value which is actually present in the data, i.e.
#'  a value which appears in the first column of the input file.
#' @param experimentEnd End of the experiment (experimentEnd). A value of -1 will let
#' the script select these time points automatically. Otherwise, specify a time value which is actually present in the data, i.e.
#' a value which appears in the first column of the input file.
#' @param overrideMaxEuclidean Maximal deviation (d_max) using Euclidean distance (overrideMaxEuclidean).
#' A value of -1 will let the script select these time points automatically. Otherwise, specify a time value which is actually present in the data, i.e.
#'  a value which appears in the first column of the input file.
#' @param overrideMaxCanberra Maximal deviation (d_max) using Canberra distance (overrideMaxCanberra).
#' A value of -1 will let the script select these time points automatically. Otherwise, specify a time value which is actually present in the data, i.e.
#'  a value which appears in the first column of the input file.
#' @param tref Disturbance event (tref). Otherwise, specify a time value which is actually present in the data, i.e.
#'  a value which appears in the first column of the input file.
#' @param time.col Name of column in sample_data i.e. metadata specifying
#' time points.
#' @return A data frame with stability proprties.
#'
#' @seealso To be used internally by \code{\link{stability_properties}}
#' @importFrom vegan vegdist
#' @references
#' \itemize{
#' \item{}{Liu, Z., et al., (2018). Ecological stability properties of
#' microbial communities assessed by flow cytometry. mSphere, 3(1), e00564-17.
#' http://msphere.asm.org/content/3/1/e00564-17
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' ps1.sub <- subset_samples(ps1.b5, Time_hr_num >= 120)
#' dat.stab <- stability_properties(ps1.sub, time.col = "Time_hr")
#' kable(dat.stab)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @export
#' @keywords Anlaysis and visualization

stability_properties <- function(ps,
                                 experimentStart = -1,
                                 experimentEnd = -1,
                                 overrideMaxEuclidean = -1,
                                 overrideMaxCanberra = -1,
                                 tref = -1,
                                 time.col = "Time_hr") {
  tax_count <- data <- NULL

  message(paste0("No. of taxa ", ntaxa(ps)))
  message("Normalizing counts to relative abundances")
  tax_count <- taxa_time_table(ps,
    time.col = "Time_hr",
    normalize = TRUE,
    remove.zero = TRUE
  )

  # head(tax_count)
  message(paste0("No. of timepoints ", nrow(tax_count)))

  tax_count <- as.data.frame(t(tax_count))
  tax_count$Time <- as.numeric(rownames(tax_count))
  tax_count <- tax_count %>% arrange(Time)
  # head(otu.rel.b5.t)
  # data <- read.table(filename, header=TRUE, stringsAsFactors = FALSE)

  data <- tax_count[, c(ncol(tax_count), 1:ntaxa(ps))]

  if (experimentStart == -1) {
    experimentStart <- data[1, 1]
  }
  if (experimentEnd == -1) {
    experimentEnd <- data[nrow(data), 1]
  }
  if (tref == -1) {
    tref <- get_reference_phase_end(data, experimentStart, experimentEnd) # time point of the disturbance, estimated from change in community structure
  }

  message("Characterizing the reference space")
  # Characterizing the reference space

  numberOfGates <- ncol(data) - 1
  data$referencePhase <- FALSE
  data[data[1] >= experimentStart & data[1] <= tref, "referencePhase"] <- TRUE

  referenceState <- colMeans(data[data$referencePhase == TRUE, 2:(numberOfGates + 1)])
  gateSDs <- apply(data[data$referencePhase == TRUE, 2:(numberOfGates + 1)], 2, sd)

  # (Radii will be computed later)

  # Sudarshan: Still need to work on this aspect
  # This chunk Start
  # get two most abundant gates
  if (!exists("plotGatesStart")) {
    topGates <- names(sort(referenceState, decreasing = TRUE))[1:2]
  } else {
    topGates <- plotGatesStart
  }

  if (!exists("plotGatesEnd")) {
    finalTopGates <- names(sort(data[data[1] == experimentEnd, 2:(numberOfGates + 1)], decreasing = TRUE))[1:2]
  } else {
    finalTopGates <- plotGatesEnd
  }

  # This chunk end

  # restricting data set to one experiment
  data <- data[data[1] >= experimentStart & data[1] <= experimentEnd, ]

  message("Computing deviation from reference state")
  # compute deviation from reference state
  data$euclidean <- apply(data[2:(numberOfGates + 1)], 1, function(x) dist(rbind(referenceState, x), method = "euclidean"))
  data$euclidean <- data$euclidean / sqrt(2)
  data$canberra <- apply(data[2:(numberOfGates + 1)], 1, function(x) dist(rbind(referenceState, x), method = "canberra"))
  data$canberra <- data$canberra / numberOfGates

  message("Computing radii for reference space")
  # compute radii for reference space
  r_Euclidean <- max(data$euclidean[data$referencePhase == TRUE])
  r_Canberra <- max(data$canberra[data$referencePhase == TRUE])

  if (overrideMaxEuclidean == -1) {
    maxEuclideanId <- order(data$euclidean, decreasing = TRUE)[1]
  } else {
    maxEuclideanId <- match(overrideMaxEuclidean, data[[1]])
  }
  if (overrideMaxCanberra == -1) {
    maxCanberraId <- order(data$canberra, decreasing = TRUE)[1]
  } else {
    maxCanberraId <- match(overrideMaxCanberra, data[[1]])
  }

  maxDevStateEuclidean <- data[maxEuclideanId, 2:(numberOfGates + 1)]
  maxDevStateCanberra <- data[maxCanberraId, 2:(numberOfGates + 1)]

  maxDevStateEuclideanTime <- data[maxEuclideanId, 1]
  maxDevStateCanberraTime <- data[maxCanberraId, 1]

  message("Computing online version of resilience")
  # compute online version of resilience
  data$maxCanberraOnline[1] <- data$canberra[1]
  for (i in 2:nrow(data)) {
    if (data$canberra[i] < data$maxCanberraOnline[i - 1]) {
      data$maxCanberraOnline[i] <- data$maxCanberraOnline[i - 1]
    } else {
      data$maxCanberraOnline[i] <- data$canberra[i]
    }
  }

  message("Computing Resilience RL")
  # compute Resilience RL
  data$RLeuclidean <- apply(data, 1, function(x) computeRL(x["euclidean"], data[maxEuclideanId, "euclidean"]))
  data$RLcanberra <- apply(data, 1, function(x) computeRL(x["canberra"], data[maxCanberraId, "canberra"]))
  data$RLcanberraOnline <- apply(data, 1, function(x) computeRL(x["canberra"], x["maxCanberraOnline"]))
  trefId <- match(tref, data[[1]])
  return(list(
    "data" = data, "numberOfGates" = numberOfGates, "referenceState" = referenceState,
    "tref" = tref, "trefId" = trefId,
    "experimentStart" = experimentStart,
    "experimentEnd" = experimentEnd,
    "maxEuclideanId" = maxEuclideanId, "maxCanberraId" = maxCanberraId
  ))
}
