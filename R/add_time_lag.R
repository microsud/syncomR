#' @title Add time lags between sampling
#' @description For time-series data, knowing the time lags can be useful for analysis.
#' @param ps A \code{\link{phyloseq-class}}.
#' @param time.col Name of column in sample_data i.e. metadata specifying
#' time points.
#' @param sam.name.col Name of column in sample_data i.e. metadata specifying
#' sample names. This should match the output of sample_names(ps).
#' @return A \code{\link{phyloseq-class}} with addition of time span between two sampling
#' events to the metadata.
#' @author Sudarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @references See citation("syncomR")
#' @export
#' @examples
#' \dontrun{
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' ps1.b5.lag <- add_time_lag(ps1.b5, time.col = "Time_hr_num", sam.name.col = "FAIR_Labels")
#' head(meta(ps1.b5.lag))
#' }
#' @keywords Utilities

add_time_lag <- function(ps, time.col = "Time_hr_num", sam.name.col = "FAIR_Labels") {
  ps.df <- NULL
  ps.df <- meta(ps)
  ps.df$Time.2 <- ps.df[, time.col]
  ps.df <- ps.df %>%
    arrange(Time.2) %>%
    mutate(time_lag = Time.2 - lag(Time.2))
  ps.df <- select(ps.df, -(Time.2))
  rownames(ps.df) <- ps.df[, sam.name.col]
  ps.df$time_lag <- ps.df$time_lag
  message("column `time_lag` added to `sample_data`")
  sample_data(ps) <- sample_data(ps.df)
  return(ps)
}
