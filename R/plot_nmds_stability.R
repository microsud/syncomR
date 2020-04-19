#' @title Plot stability properties NMDS
#'
#' @description NMDS plot compared to reference states as decribed by Florian Centler and Zishu Liu, UFZ Leipzig, Germany.
#'
#' @details This script was modified to get ggplot object from the paper on "Ecological Stability Properties.
#' Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
#'  \url{\link{http://msphere.asm.org/content/3/1/e00564-17}}.
#'  This functions takes output from the \code{\link{stability_properties}}.
#' \itemize{
#' \item{}{Plot reference states and community evolusion}
#' }
#' @param stab.in The input come from \code{\link{stability_properties}} function.
#' @param refphase Color for reference phase. DEFAULT = grey60
#' @param refstate Color for reference state DEFAULT = red
#' @param evolution.points Color for to track evolution points DEFAULT = yellow
#' @param max.distant.states.euc Color for sites/timepoint that are max distant, Euclidean DEFAULT = steelblue
#' @param max.distant.states.can Color for sites/timepoint that are max distant, Canberra. DEFAULT = brown3
#' @param final.state Color for last time state DEFAULT = black
#' @note
#' \itemize{
#' \item{}{refphase= Color for reference phase. DEFAULT = grey60}
#' \item{}{refstate= Color for reference state DEFAULT = red}
#' \item{}{evolution.points= Color for to track evolution points DEFAULT = yellow}
#' \item{}{max.distant.states.euc= Color for sites/timepoint that are max distant, Euclidean DEFAULT = steelblue}
#' \item{}{max.distant.states.can = Color for sites/timepoint that are max distant, Canberra. DEFAULT = brown3}
#' \item{}{final.state= Color for last time state DEFAULT = black}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#'
#' @return A data frame and ggplot object.
#'
#' @seealso Input for this functions comes from \code{\link{stability_properties}}
#' @importFrom vegan metaMDS
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
#' nmds.dt <- plot_nmds_stability(dat.stab)
#' print(nmds.dt$plot)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @export
#' @keywords Anlaysis and visualization

plot_nmds_stability <- function(stab.in,
                                refphase = "grey60",
                                refstate = "red",
                                evolution.points = "yellow",
                                max.distant.states.euc = "steelblue",
                                final.state = "black",
                                max.distant.states.can = "brown3") {
  numberOfGates <- stab.in$numberOfGates
  trefId <- stab.in$trefId
  referenceState <- stab.in$referenceState
  experimentStart <- stab.in$experimentStart
  experimentEnd <- stab.in$experimentEnd
  maxEuclideanId <- stab.in$maxEuclideanId
  maxCanberraId <- stab.in$maxCanberraId
  data <- stab.in$data


  mds.out <- metaMDS(rbind(data[, 2:(numberOfGates + 1)], referenceState), distance = "bray", autotransform = FALSE, zerodist = "add")
  nmd.dat <- data.frame(MDS1 = mds.out$point[, 1], MDS2 = mds.out$point[, 2])
  arow <- rbind(nmd.dat[length(nmd.dat$MDS1), 1:2], nmd.dat[(trefId + 1), 1:2])
  p <- ggplot(nmd.dat, aes(MDS1, MDS2)) +
    ylim(c(-0.5, 0.5)) +
    theme_syncom() +
    geom_hline(aes(yintercept = 0), alpha = 0.1) +
    geom_vline(aes(xintercept = 0), alpha = 0.1) +
    geom_path(
      data = nmd.dat[c((trefId + 1):(length(nmd.dat[, 1]) - 1)), ],
      aes(MDS1, MDS2), arrow = arrow(
        angle = 15, length = unit(0.2, "inches"),
        ends = "last", type = "closed"
      ),
      size = 1.0, alpha = 0.6
    ) +
    geom_path(
      data = arow,
      aes(MDS1, MDS2), linejoin = "mitre",
      arrow = arrow(
        angle = 15, length = unit(0.2, "inches"),
        ends = "last", type = "closed"
      ), size = 1.0,
      color = refstate, alpha = 0.6
    ) +
    geom_point(shape = 21, size = 4) +
    xlim(c(-0.5, 0.5))
  # refphase
  p <- p + geom_point(data = nmd.dat[c(1:trefId), ], aes(MDS1, MDS2), fill = refphase, size = 4, shape = 21)

  # refstate
  p <- p + geom_point(
    data = nmd.dat[c(length(nmd.dat[, 1]), length(nmd.dat[, 1])), ],
    aes(MDS1, MDS2), fill = refstate, size = 4, shape = 21
  )

  # evolution mds.out$points
  p <- p + geom_point(
    data = nmd.dat[c((trefId + 1):(length(nmd.dat[, 1]) - 1)), ],
    aes(MDS1, MDS2), fill = evolution.points, shape = 21, size = 4
  )

  # mark most distant states
  p <- p + geom_point(
    data = nmd.dat[c(maxEuclideanId, maxEuclideanId), ],
    aes(MDS1, MDS2), fill = max.distant.states.euc, shape = 24, size = 4
  )

  p <- p + geom_point(
    data = nmd.dat[c(maxCanberraId, maxCanberraId), ],
    aes(MDS1, MDS2), fill = max.distant.states.can, shape = 24, size = 4
  )
  # mark final state
  p <- p + geom_point(
    data = nmd.dat[c(length(nmd.dat[, 1]) - 1, length(nmd.dat[, 2]) - 1), ],
    aes(MDS1, MDS2), fill = final.state, shape = 21, size = 4
  )
  p <- p + labs(title = "", subtitle = expression(paste("Community evolve apart ", s[ref])))
  # p <- p + geom_path(nmd.dat, aes(MDS1,MDS2), colour="black")


  return(list("nmds_data" = mds.out, "plot" = p))
}
