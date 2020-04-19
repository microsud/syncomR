#' @title Normalize a matrix
#'
#' @description Normalize a matrix column-wise by dividing each entry by its corresponding column sum. Transfered
#' the function from seqtime R package \url{https://github.com/hallucigenia-sparsa/seqtime/blob/master/R/normalize.R}
#' becasue of issues with igraph dependency.
#'
#' @details Columns summing to zero are removed by default.
#'
#' @param x A matrix
#' @param removeZero Remove columns summing to zero.
#' @return A normalized matrix
#' @export
#' @keywords Data and filtering

normalize_seqtime <- function(x, removeZero = TRUE) {
  # remove columns with only zeros from matrix, to avoid dividing by a zero
  colsums <- apply(x, 2, sum)
  if (removeZero == TRUE) {
    zero.col.indices <- which(colsums == 0)
    # print(length(zero.col.indices))
    if (length(zero.col.indices) > 0) {
      colsums <- colsums[setdiff(1:ncol(x), zero.col.indices)]
      x <- x[setdiff(1:ncol(x), zero.col.indices), ]
    }
  }
  for (i in 1:ncol(x)) {
    x[, i] <- x[, i] / colsums[i]
  }
  x
}
