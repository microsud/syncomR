#' @title Niche Overlap Estimation Quantitative Traits from Temporal Abundance Data
#'
#' @description Compute the niche overlap using quantitative traits Mouillot et al. 2005.
#'
#' @details This approach was developed by Mouillot et al. 2005 and The code was
#'          obtained and modified from
#'          \url{https://github.com/umr-marbec/nicheoverlap/blob/master/nicheoverlap.R}.
#'          Computation is done between N species along a trait axis (functional or
#'          environmental) generalization of the Niche Overlap index (NO) based on kernel
#'          density estimation.  The index is non-parametric and assumes no normal distribution.
#'          This is particularly important as -omics data are known for challenges w.r.t normality.
#'          The N densities along the finite range of values considered are corrected according
#'          to their respective integrals, so that integrals of corrected densities along the range
#'          equal 1. If plot is T,each species niche have its own color, their overlap being shaded
#'          by lines.
#'
#' @param trab a matrix (x, N+1) with the abundances of the N species at x values of the
#'             considered trait, NAs are not allowed in traits values.
#'
#' @param graph TRUE or FALSE  indicating whether a graphical illustration of niche overlap is wanted.
#'
#' @param colorsN vector (length=N) for colors of niches in hexadecimal code (if NA default colors).
#'
#' @return A plot or a data frame with the value of the overlap between corrected
#'         densities: 0 (no overlap) to 1(N identic densities).
#' @references
#' \itemize{
#' \item{}{Mouillot D, et al. Niche overlap estimates based on quantitative functional traits: a new
#' family of non-parametric indices. Oecologia. 2005 Sep 1;145(3):345-53.
#' \url{https://link.springer.com/article/10.1007/s00442-005-0151-z}}
#' \item{}{Shetty SA. et al. A Minimalist Approach for Deciphering the Ecophysiology of Human Gut Microbes (2020)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#'
#' @keywords Analysis and visualization
#'
#' @importFrom graphics abline arrows axis legend lines mtext points polygon rect text title
#'
#' @importFrom stats density integrate sd splinefun time
#'
#' @importFrom utils data
#'
#' @export
#'
niche_overlap <- function(trab, graph = T, colorsN = NA) {

  # data extraction and checking

  # traits values
  tr <- trab[, 1]
  if (length(which(is.na(tr))) != 0) stop("error : NA are not allowed in traits values")
  mintr <- min(tr)
  maxtr <- max(tr)

  # abundances
  abN <- trab[, -1]
  if (dim(abN)[2] < 2) stop("error : there must be at least 2 species")

  # number of species
  N <- dim(abN)[2]


  ###############################################################################################################
  #                     function to compute kernel density estimates

  nichedens <- function(trab1) {

    # data extraction
    tr1 <- trab1[, 1]
    ab1 <- trab1[, 2]

    # transformation of abundances into relative abundances and then conversion into whole numbers (precision of 0.001)

    abrel <- ab1 / sum(ab1, na.rm = T)
    abrelW <- round(abrel, 3) * 1000 # chaneg 3 to 6

    # transformation of data: traits values are replicated according to their abundances
    trrep <- vector()
    for (i in 1:length(abrelW)) {
      if (abrelW[i] != 0 & is.na(abrelW[i]) == F) trrep <- c(trrep, rep(tr1[i], abrelW[i]))
    }

    # Bandwidth
    bandwdth <- 1.06 * sd(trrep) * (length(trrep)^-0.2)

    # computation of density based on gaussian kernel along the range of trait
    resdens <- density(trrep,
      bw = bandwdth, kernel = "gaussian",
      from = min(tr1), to = max(tr1), n = 1024,
      na.rm = TRUE
    )
    kx <- resdens$x
    ky <- resdens$y

    dens <- cbind(kx, ky)

    return(dens)
  } # end of function nichedens
  ###############################################################################################################

  # computation of the N densities along the traits axis
  for (i in 1:N)
  {
    densi <- nichedens(cbind(tr, abN[, i]))
    if (i == 1) {
      resdens <- densi[, 1:2]
    } else {
      resdens <- cbind(resdens, densi[, 2])
    }
  } # end of i

  # minimum of the N densities and interpolation
  x <- resdens[, 1]
  densN <- resdens[, 2:(N + 1)]

  # computation of the N integrals and correction of densities
  integ <- vector()
  densNcorr <- densN

  for (i in 1:N)
  {
    fdensk <- splinefun(x, densN[, i])
    integ[i] <- integrate(fdensk, min(x), max(x), subdivisions = 1000)$value
    densNcorr[, i] <- densN[, i] / integ[i]
  } # end of i

  # overlap computation
  minN <- apply(densNcorr, 1, min)
  fminN <- splinefun(x, minN)
  overlap <- integrate(fminN, min(x), max(x), subdivisions = 1000)$value
  NO <- round(overlap, 3)

  #############################################################################################################
  # graphic

  if (graph == T) {
    if (N < 1) stop("error : min 2 densities for graphical representation")

    # colors:  red,       blue,     green,      grey,     violet,    brown
    if (is.na(colorsN[1]) == T) color <- syncom_colors("BacterialSpecies")
    if (is.na(colorsN[1]) == F) color <- colorsN
    if (is.na(colorsN[1]) == F & length(colorsN) != N) stop("'colorsN' must be of same length than number of species")
    transp <- 50


    # axes scale
    labx <- pretty(x, n = 5, min.n = 3)
    laby <- pretty(c(0, as.vector(densNcorr)), n = 5, min.n = 3)

    # plot of axes
    plot(mean(labx), mean(laby), type = "n", xlim = range(x), ylim = range(laby), xlab = "", ylab = "", xaxt = "n", yaxt = "n", axes = F)
    title(sub = paste("Niche overlap=", round(NO, 3), sep = ""), line = 2, cex.sub = 1.2)
    axis(side = 1, labx[labx <= max(x) & labx >= min(x)], tcl = -0.2, pos = 0, labels = F)
    mtext(labx[labx <= max(x) & labx >= min(x)], side = 1, at = labx[labx <= max(x) & labx >= min(x)], las = 1, line = -0.2, cex = 1)

    axis(side = 2, laby, tcl = -0.2, pos = min(x), labels = F)
    mtext(laby, side = 2, at = laby, las = 1, line = -0.4, cex = 1)

    # plot of densities
    for (k in 1:N) {
      polygon(rbind(c(min(x), 0), c(min(x), densNcorr[1, k]), cbind(x, densNcorr[, k]), c(max(x), 0)), border = color[k], col = paste(color[k], transp, sep = ""))
    } # end of k

    # shading lines
    polygon(rbind(c(min(x), 0), c(min(x), minN[1]), cbind(x, minN), c(max(x), 0)), density = 10)

    # border
    rect(min(x), min(laby), max(x), max(laby))
  } # end of if (graph==T)

  #############################################################################################################

  return(NO)
} # end of function nicheoverlap

################################################################################

## TO DO Make it work on long data frame with multiple samples or timepoints i.e taxa*traits*samples
