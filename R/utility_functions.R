#' @title Utility for stability calculations
#'
#' @description List taxa for stability properties analysis as well as identifies reference states.
#'
#' @details These are utility functions that are required for \code{\link{stability_properties}} analysis
#' @param mergeTable input from stability_properties
#' @param idx time id for stability_properties
#' @return Taxa list and reference states as calulated by user provided distance metric from \code{\link{vegan}}.
#'
#' @seealso To be used internally by \code{\link{stability_properties}}
#' @keywords Utilities
#' @references
#' \itemize{
#' \item{}{"Liu, Z., Cichocki, N., Bonk, F., Günther, S., Schattenberg, F., Harms, H., ... & Müller, S. (2018). Ecological stability properties of
#' microbial communities assessed by flow cytometry. mSphere, 3(1), e00564-17.
#' http://msphere.asm.org/content/3/1/e00564-17
#' }
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @export
#'
get_members <- function(mergeTable, idx) {
  if (idx < 0) {
    return(-idx)
  } else {
    aList <- get_members(mergeTable, mergeTable[idx, 1])
    bList <- get_members(mergeTable, mergeTable[idx, 2])

    return(c(aList, bList))
  }
}

#' @title Utility for stability calculations
#'
#' @description List taxa for stability properties analysis as well as identifies reference states.
#'
#' @details Fundamental niche values are calculated on the basis of a multi-variate
#' @param data internal stability_properties
#' @param begin internal stability_properties
#' @param end internal stability_properties
#' @return Taxa list and reference states as calulated by user provided distance metric from \code{\link{vegan}}.
#'
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
#' @importFrom stats dist
#' @importFrom stats hclust
#' @export
#'
get_reference_phase_end <-
  function(data, begin, end) {
    myData <- data[data[1] >= begin & data[1] <= end, ]
    rownames(myData) <- myData[[1]]
    myData <- myData[-1]
    myD <- dist(myData, method = "canberra")
    hc <- hclust(myD)

    # plot(hc, cex=0.6, main="Reference state detection",
    # xlab=NA, sub=NA, font.main=1, cex.main=titleSize)

    startTime <- min(as.numeric(hc$labels))

    index <- match(startTime, hc$labels)

    nMerger <- nrow(hc$merge)
    startIdx <- 0

    for (i in 1:nMerger) {
      if (hc$merge[i, 1] == -index) {
        startIdx <- hc$merge[i, 2]
        break
      } else {
        if (hc$merge[i, 2] == -index) {
          startIdx <- hc$merge[i, 1]
          break
        }
      }
    }

    if (startIdx == 0) {
      print("Error!")
    }

    myList <- get_members(hc$merge, startIdx)

    value <- max(as.numeric(hc$labels[myList]))
    return(value)
  }



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


#' @title Utility for stability calculations resilience
#' @description Function for computing resilience
#' @param dx0xt internal
#' @param dx0x1 internal
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
#' @export computeRL

computeRL <- function(dx0xt, dx0x1) {
  return(2.0 * dx0x1 / (dx0x1 + dx0xt) - 1.0)
}


#' @title Utility for stability calculations
#' @description Function for computing resilience
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
#'

#' @title Default ggplot2 theme
#' @description Default ggplot2 theme for consistency in publications.
#' @param base_size Similiar to ggplot2 themes
#' @param base_family Similiar to ggplot2 themes
#' @param base_line_size Similiar to ggplot2 themes
#' @param base_rect_size Similiar to ggplot2 themes
#' @author Sudarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @references See citation("syncomR")
#' @export
#' @keywords Theme elements


theme_syncom <- function(base_size = 12,
                         base_family = "",
                         base_line_size = base_size / 170,
                         base_rect_size = base_size / 170) {
  theme_bw(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size
  ) %+replace%
    theme(
      panel.background = element_rect(fill = "white"),
      # panel.border = element_rect(size = 1),
      # plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "line"),
      # axis.ticks = element_blank(),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      # legend.key = element_blank(),
      complete = TRUE
    )
}




#' @title Default Color Scheme
#' @description Default colors for different variables.
#' @param x Name of the variable type
#' @param v Optional. Vector of elements to color.
#' @return Named character vector of default colors
#' @author Sudarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @references See citation("syncomR")
#' @export
#' @examples
#' \dontrun{
#' col <- syncom_colors("Phylum")
#' }
#' @keywords Theme elements
syncom_colors <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    syncom_colors <- c(
      Akkermansia_muciniphila = "#e6194b", Bifidobacterium_adolescentis = "#3cb44b",
      Collinsella_aerofaciens = "#1f1f1f", Bacteroides_ovatus = "#4363d8",
      Bacteroides_xylanisolvens = "#f58231", Agathobacter_rectalis = "#911eb4",
      Anaerobutyricum_soehngenii = "#46f0f0", Eubacterium_siraeum = "#f032e6",
      Blautia_hydrogenotrophica = "#bcf60c", Coprococcus_catus = "#fabebe",
      Flavonifractor_plautii = "#008080", Roseburia_intestinalis = "#e6beff",
      Faecalibacterium_prausnitzii = "#9a6324", Blautia_obeum = "#e4cd05",
      Ruminococcus_bromii = "#800000", Subdoligranulum_variabile = "#aaffc3",
      Acetate = "#808000", Butyrate = "#ffd8b1",
      Formate = "#000075", Propionate = "#808080",
      Lactate = "#ffffff", Succinate = "#000000"
    )
  }

  syncom_colors
}

#' @title Default Color Scheme
#' @description Default colors for different variables.
#' @param x Name of the variable type
#' @param v Optional. Vector of elements to color.
#' @return Named character vector of default colors
#' @author Sudarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @references See citation("syncomR")
#' @export
#' @examples
#' \dontrun{
#' col <- syncom_colors2("Species")
#' }
#' @keywords Theme elements
syncom_colors2 <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    syncom_colors2 <- c(
      `Akkermansia muciniphila` = "#e6194b", `Bifidobacterium adolescentis` = "#3cb44b",
      `Collinsella aerofaciens` = "#1f1f1f", `Bacteroides ovatus` = "#4363d8",
      `Bacteroides xylanisolvens` = "#f58231", `Agathobacter rectalis` = "#911eb4",
      `Anaerobutyricum soehngenii` = "#46f0f0", `Eubacterium siraeum` = "#f032e6",
      `Blautia hydrogenotrophica` = "#bcf60c", `Coprococcus catus` = "#fabebe",
      `Flavonifractor plautii` = "#008080", `Roseburia intestinalis` = "#e6beff",
      `Faecalibacterium prausnitzii` = "#9a6324", `Blautia obeum` = "#e4cd05",
      `Ruminococcus bromii` = "#800000", `Subdoligranulum variabile` = "#aaffc3",
      Acetate = "#808000", Butyrate = "#ffd8b1",
      Formate = "#000075", Propionate = "#808080",
      Lactate = "#ffffff", Succinate = "#000000"
    )
  }

  syncom_colors2
}
