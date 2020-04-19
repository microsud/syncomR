#' @title Plot stability properties Resilience
#'
#' @description Resilience of the community as decribed by Florian Centler and Zishu Liu, UFZ Leipzig, Germany.
#'
#' @details This script was modified to get ggplot object from the paper on "Ecological Stability Properties.
#' Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
#'  \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}.
#'  This functions takes output from the \code{\link{stability_properties}}.
#' @param stab.in The input come from \code{\link{stability_properties}} function.
#' @param method Either euclidean or canberra
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
#' plot_resilience(dat.stab, method = "euclidean")
#' plot_resilience(dat.stab, method = "canberra")
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @export
#' @keywords Anlaysis and visualization
#'
plot_resilience <- function(stab.in, method = c("euclidean", "canberra")) {
  if (method == "euclidean") {
    p <- plot_resilience_euclidean(stab.in)
    return(p)
  } else if (method == "canberra") {
    p <- plot_resilience_canberra(stab.in)

    return(p)
  }
}

plot_resilience_euclidean <- function(stab.in) {
  experimentEnd <- stab.in$experimentEnd
  maxEuclideanId <- stab.in$maxEuclideanId
  # maxCanberraId <- stab.in$maxCanberraId
  data <- stab.in$data
  ElasticityEuclidean <- (data[maxEuclideanId, "euclidean"] - data[nrow(data), "euclidean"]) / (data[nrow(data), 1] - data[maxEuclideanId, 1])


  df.1 <- subset(data, Time %in% data[maxEuclideanId:nrow(data), 1])

  df.2 <- df.1[, c("Time", "RLeuclidean")]

  # df.2 <- df.resl %>% gather(df.resl, RLDistance, -Time)

  df.2$max.dis <- ifelse(df.2$Time == data[maxEuclideanId, 1], "MaxDeviation",
    ifelse(df.2$Time == experimentEnd, "ExperimentEnd", "No")
  )

  colnames(df.2) <- c("Time", "RLeuclidean", "Variable")
  p <- ggplot(df.2, aes(Time, RLeuclidean, label = Time)) +
    theme_syncom() +
    geom_line(size = 0.6, alpha = 0.8, color = "deepskyblue3") +
    ylim(c(0, 1)) +
    geom_point(aes(shape = Variable), color = "deepskyblue3", size = 2) +
    scale_color_manual(values = c(
      RLeuclidean = "deepskyblue3"
    )) +
    geom_hline(
      yintercept = 0,
      colour = "black", linetype = "dashed"
    )
  p <- p + scale_shape_manual(values = c(
    MaxDeviation = 2,
    No = 16,
    ExperimentEnd = 0
  ), guide = FALSE) + ylab(expression(~ Resilience[Euclidean])) #+ geom_label_repel(subset(df.2, Time== df.2[1,1]), aes(Time,RLeuclidean))

  lab.txt <- paste(
    "RL=", round(data$RLeuclidean[nrow(data)], digits = 4),
    ", E=", round(ElasticityEuclidean, digits = 4)
  )
  p <- p + ggtitle(label = "", subtitle = lab.txt)

  return(p)
}
# expression(paste("Community evolve apart ", s[ref]))
plot_resilience_canberra <- function(stab.in) {
  experimentEnd <- stab.in$experimentEnd
  # maxEuclideanId <- stab.in$maxEuclideanId
  maxCanberraId <- stab.in$maxCanberraId
  data <- stab.in$data
  ElasticityCanberra <- (data[maxCanberraId, "canberra"] - data[nrow(data), "canberra"]) / (data[nrow(data), 1] - data[maxCanberraId, 1])

  df.1 <- subset(data, Time %in% data[maxCanberraId:nrow(data), 1])

  df.2 <- df.1[, c("Time", "RLcanberra")]

  # df.2 <- df.resl %>% gather(df.resl, RLDistance, -Time)

  df.2$max.dis <- ifelse(df.2$Time == data[maxCanberraId, 1], "MaxDeviation",
    ifelse(df.2$Time == experimentEnd, "ExperimentEnd", "No")
  )

  colnames(df.2) <- c("Time", "RLcanberra", "Variable")
  p <- ggplot(df.2, aes(Time, RLcanberra, label = Time)) +
    theme_syncom() +
    geom_line(size = 0.6, alpha = 0.8, color = "brown3") +
    ylim(c(0, 1)) +
    geom_point(aes(shape = Variable), color = "brown3", size = 2) +
    scale_color_manual(values = c(
      RLeuclidean = "brown3"
    )) +
    geom_hline(
      yintercept = 0,
      colour = "black", linetype = "dashed"
    )
  p <- p + scale_shape_manual(values = c(
    MaxDeviation = 2,
    No = 16,
    ExperimentEnd = 0
  ), guide = FALSE) + ylab(expression(~ Resilience[Canberra]))
  lab.txt <- paste(
    "RL=", round(data$RLcanberra[nrow(data)], digits = 4),
    ", E=", round(ElasticityCanberra, digits = 4)
  )
  p <- p + ggtitle(label = "", subtitle = lab.txt)

  return(p)
}
