#' syncomR: Synthetic Microbial Community Data Exploration, Analysis and Visualization using R.
#'
#' The syncomR package is developed to analyse and visualize temporal microbial community profiles. The development
#' was done for and using the human gut synthetic microbial community designed to investigate ecophysiological properties of
#' core human gut microbes by \href{https://github.com/microsud/syncomR/tree/master}{Shetty SA et al., 2020.
#' Minimalist Approach for Deciphering the Ecophysiology of Human Gut Microbes}.
#' The package requires \code{\link{phyloseq-class}} object of microbial community profile for most of the
#' analysis aimed at the ecological aspects. In addition, package data includes metatranscriptomics data, metabolites data,
#' and gut micorbiota module expression data. The focus is to analyse temporal dynamics of  microbial community,
#' test changes in abundances, overall community variability, resistance and resilience properties of the microbial communities. The ecological stability properties are calculated using the modified R code
#' from \href{http://msphere.asm.org/content/3/1/e00564-17}{Liu et al., 2018. Ecological Stability Properties of
#' Microbial Communities Assessed by Flow Cytometry}, while community level temporal variability and succession
#' is calculated with modified R code from \href{https://www.nature.com/articles/s41467-019-08377-w}{Guittar J et al., 2019.
#' Trait-based community assembly and succession of the infant gut microbiome}.
#' Additional functions are provided for data handling and formatting. Visualization tools are provided to
#' explore the microbial community dynamics, including area, bar and line plot for composition, ordinaion plots
#' for visualizing movement of community in and out of reference phase.
#'
#' \tabular{ll}{
#' Package: \tab syncomR\cr
#' Type: \tab Package\cr
#' Version: \tab See sessionInfo() or DESCRIPTION file\cr
#' Date: \tab 2019-2024\cr
#' License: \tab FreeBSD\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Synthetic Microbial Community Data Exploration, Analysis and Visualization using R
#'
#' @name syncomR-package
#' @aliases syncomR
#' @docType package
#' @title syncomR: Synthetic Microbial Community Data Exploration, Analysis and Visualization using R
#'
#' @author Sudarshan Shetty et al. 2020
#'   \email{sudarshanshetty9@gmail.com}
#' @references
#' \itemize{
#' \item{}{Liu, Z., et al., (2018). Ecological stability properties of microbial communities assessed by flow cytometry. mSphere, 3(1), e00564-17.
#' \href{http://msphere.asm.org/content/3/1/e00564-17}{link to paper}
#' }
#' \item{}{Guittar, J., Shade, A., & Litchman, E. (2019). Trait-based community assembly
#' and succession of the infant gut microbiome. Nature communications, 10(1), 512.
#' \href{https://www.nature.com/articles/s41467-019-08377-w}{link to paper}
#' }
#' \item{}{'Shetty SA et al., (2019-2020) Minimalist Approach for Deciphering the Ecophysiology of Human Gut
#' Microbes XXXX}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' \url{https://github.com/microsud/syncomR/tree/master}
#' @examples
#' citation("syncomR")
#' @keywords package
NULL
