#' @title Utility for stability calculations
#'
#' @description List taxa for stability properties analysis as well as identifies reference states.
#'
#' @details These are utility functions that are required for \code{\link{stability_properties}} analysis
#' @param mergeTable input from stability_properties
#' @param idx time id for stability_properties
#' @return Taxa list and reference states as calculated by user provided distance metric from \code{\link{vegan}}.
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
#'
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
#'
#' @keywords Utilities
#' @importFrom stats dist
#' @importFrom stats hclust
#' @export

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
#'
#' @export computeRL
computeRL <- function(dx0xt, dx0x1) {
  return(2.0 * dx0x1 / (dx0x1 + dx0xt) - 1.0)
}


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
#' @import ggplot2

theme_syncom <- function(base_size = 11,
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
      #panel.border = element_rect(fill = NA, colour = "#303030", size = 1),
      # plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "line"),
      # axis.ticks = element_blank(),
      legend.title = element_text(colour = "#303030", size = rel(1.2), face = "bold", hjust=0),
      legend.text = element_text(size = rel(1.2), colour = "#303030"),
      legend.key = element_rect(colour = NA, fill = NA),
      #legend.background = element_rect(colour = NA, fill = NA),
      axis.line = element_blank(),
      axis.text = element_text(colour = "#303030",size = rel(1)),
      axis.title = element_text(colour = "#303030",size = rel(1.2)),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text = element_text(colour = "#303030", size = rel(1.2)),
      plot.title = element_text(colour = "#303030", size = rel(1.2),hjust = 0, margin=margin(0,0,10,0)),
      plot.subtitle = element_text(colour = "#303030", size = rel(1),hjust = 0, margin=margin(0,0,10,0)),
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
#' @examples
#' \dontrun{
#' col <- syncom_colors("OTU")
#' }
#' @keywords Theme elements
#' @rdname colors
#' @name syncom_colors
#' @export
syncom_colors <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    syncom_colors <- c(
      Akkermansia_muciniphila = "#e6194b", Bifidobacterium_adolescentis = "#3cb44b",
      Collinsella_aerofaciens = "#1f1f1f", Bacteroides_ovatus = "#7789c9",
      Bacteroides_xylanisolvens = "#f58231", Agathobacter_rectalis = "#911eb4",
      Anaerobutyricum_soehngenii = "blue4", Eubacterium_siraeum = "#b097ad",
      Blautia_hydrogenotrophica = "#bcf60c", Coprococcus_catus = "#fabebe",
      Flavonifractor_plautii = "#008080", Roseburia_intestinalis = "#e6beff",
      Faecalibacterium_prausnitzii = "#9a6324", Blautia_obeum = "#e4cd05",
      Ruminococcus_bromii = "#800000", Subdoligranulum_variabile = "#88cf9e",
      Acetate = "#808000", Butyrate = "#ffd8b1",
      Formate = "#000075", Propionate = "#808080",
      Lactate = "#ffffff", Succinate = "#000000"
    )
  }

  syncom_colors
}


#' @examples
#' \dontrun{
#' col <- syncom_colors2("Species")
#' }
#' @keywords Theme elements
#' @rdname colors
#' @name syncom_colors2
#' @export
syncom_colors2 <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    syncom_colors2 <- c(
      `Akkermansia muciniphila` = "#e6194b", `Bifidobacterium adolescentis` = "#3cb44b",
      `Collinsella aerofaciens` = "#1f1f1f", `Bacteroides ovatus` = "#7789c9",
      `Bacteroides xylanisolvens` = "#f58231", `Agathobacter rectalis` = "#911eb4",
      `Anaerobutyricum soehngenii` = "blue4", `Eubacterium siraeum` = "#b097ad",
      `Blautia hydrogenotrophica` = "#bcf60c", `Coprococcus catus` = "#fabebe",
      `Flavonifractor plautii` = "#008080", `Roseburia intestinalis` = "#e6beff",
      `Faecalibacterium prausnitzii` = "#9a6324", `Blautia obeum` = "#e4cd05",
      `Ruminococcus bromii` = "#800000", `Subdoligranulum variabile` = "#88cf9e",
      Acetate = "#808000", Butyrate = "#ffd8b1",
      Formate = "#000075", Propionate = "#808080",
      Lactate = "#ffffff", Succinate = "#000000"
    )
  }

  syncom_colors2
}


#' @examples
#' \dontrun{
#' col <- bugs.colors("Species")
#' }
#' @keywords Theme elements
#' @rdname colors
#' @name bugs.colors
#' @export
#'
bugs.colors <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    bugs.colors <- c(
      `Akkermansia muciniphila` = "#e6194b", `Bifidobacterium adolescentis` = "#3cb44b",
      `Collinsella aerofaciens` = "#1f1f1f", `Bacteroides ovatus` = "#7789c9",
      `Bacteroides xylanisolvens` = "#f58231", `Agathobacter rectalis` = "#911eb4",
      `Anaerobutyricum soehngenii` = "blue4", `Eubacterium siraeum` = "#b097ad",
      `Blautia hydrogenotrophica` = "#bcf60c", `Coprococcus catus` = "#fabebe",
      `Flavonifractor plautii` = "#008080", `Roseburia intestinalis` = "#e6beff",
      `Faecalibacterium prausnitzii` = "#9a6324", `Blautia obeum` = "#e4cd05",
      `Ruminococcus bromii` = "#800000", `Subdoligranulum variabile` = "#88cf9e"
    )
  }

  bugs.colors
}

#' @examples
#' \dontrun{
#' col <- bugs.colors2("Species")
#' }
#' @keywords Theme elements
#' @rdname colors
#' @name bugs.colors2
#' @export
#'
bugs.colors2 <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    bugs.colors2 <- c(
      `Akkermansia.muciniphila` = "#e6194b", `Bifidobacterium.adolescentis` = "#3cb44b",
      `Collinsella.aerofaciens` = "#1f1f1f", `Bacteroides.ovatus` = "#7789c9",
      `Bacteroides.xylanisolvens` = "#f58231", `Agathobacter.rectalis` = "#911eb4",
      `Anaerobutyricum.soehngenii` = "blue4", `Eubacterium.siraeum` = "#b097ad",
      `Blautia.hydrogenotrophica` = "#bcf60c", `Coprococcus.catus` = "#fabebe",
      `Flavonifractor.plautii` = "#008080", `Roseburia.intestinalis` = "#e6beff",
      `Faecalibacterium.prausnitzii` = "#9a6324", `Blautia.obeum` = "#e4cd05",
      `Ruminococcus.bromii` = "#800000", `Subdoligranulum.variabile` = "#88cf9e"
    )
  }

  bugs.colors2
}


#' @examples
#' \dontrun{
#' col <- bugs.colors2("Species")
#' }
#' @keywords Theme elements
#' @rdname colors
#' @name bugs.colors3
#' @export
#'
bugs.colors3 <- function(x, v = NULL) {
  if (x == "Species" | x == "Taxon" | x == "Bacteria" | x == "BacterialStrain" | x == "BacterialSpecies" | x == "OTU") {
    # http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    # https://www.r-graph-gallery.com/42-colors-names/
    bugs.colors3 <- c(
      `Akkermansia_muciniphila` = "#e6194b", `Bifidobacterium_adolescentis` = "#3cb44b",
      `Collinsella_aerofaciens` = "#1f1f1f", `Bacteroides_ovatus` = "#7789c9",
      `Bacteroides_xylanisolvens` = "#f58231", `Agathobacter_rectalis` = "#911eb4",
      `Anaerobutyricum_soehngenii` = "blue4", `Eubacterium_siraeum` = "#b097ad",
      `Blautia_hydrogenotrophica` = "#bcf60c", `Coprococcus_catus` = "#fabebe",
      `Flavonifractor_plautii` = "#008080", `Roseburia_intestinalis` = "#e6beff",
      `Faecalibacterium_prausnitzii` = "#9a6324", `Blautia_obeum` = "#e4cd05",
      `Ruminococcus_bromii` = "#800000", `Subdoligranulum_variabile` = "#88cf9e"
    )
  }

  bugs.colors3
}
