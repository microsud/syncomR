#' @title Utility for formatting SynComGMM and SynComRNAKEGG data
#'
#' @description Group by GMM ID or KO ID.
#'
#' @details These are utility functions.
#' @param x A data.frame to formate. Data frames similar to `SynComRNAKEGG`
#' @param by Either KOID or Module
#' @return A data frame.
#'
#' @seealso To be used internally
#' @keywords Utilities
#' @references
#' \itemize{
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @export
#' @keywords Utilities
#'
group_by_variable <- function(x = NULL, by = c("KOID", "Module")) {
  xdf <- xdf2 <- xdf3 <- NULL
  # x and reference should have one matiching column to merge
  # module_df <- merge(x, reference, by.x=column)
  if (by == "KOID") {
    xdf <- x %>%
      group_by(KOID) %>%
      summarize_if(is.numeric, sum, na.rm = TRUE)

    xdf2 <- as.data.frame(xdf)
    rownames(xdf2) <- xdf2$KOID
    # mods_sum.vfa <- as.data.frame(filter(mods_sum.vfa, colnames(mods_sum.vfa) %in% rownames(hplc)))

    xdf3 <- as.data.frame(t(xdf2[, -1]))
  }
  else if (by == "Module") {
    xdf <- x %>%
      group_by(Module) %>%
      summarize_if(is.numeric, sum, na.rm = TRUE)

    xdf2 <- as.data.frame(xdf)
    rownames(xdf2) <- xdf2$Module
    # mods_sum.vfa <- as.data.frame(filter(mods_sum.vfa, colnames(mods_sum.vfa) %in% rownames(hplc)))

    xdf3 <- as.data.frame(t(xdf2[, -1]))
  }
  return(xdf3)
}
