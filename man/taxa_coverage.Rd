% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_coverage.R
\name{taxa_coverage}
\alias{taxa_coverage}
\title{Plot Coverage of Taxa}
\usage{
taxa_coverage(
  ps,
  coverage = 0.95,
  time.column = "Time_hr_num",
  shape.variable = NULL,
  color.variable = NULL,
  color.pal = NULL,
  y.breaks = seq(0, 16, 1),
  y.limits = c(11, 16),
  plot = TRUE
)
}
\arguments{
\item{ps}{A \code{\link{phyloseq-class}}}

\item{coverage}{Percent coverage limit. i.e. if 0.95 then number of taxa
accounting for 95 percent relative abundance in each sample
timepoint is ploted}

\item{time.column}{Name of column in sample_data i.e. metadata specifying
time points.}

\item{shape.variable}{Name of column in sample_data i.e. metadata specifying
variable to use for shapes}

\item{color.variable}{Name of column in sample_data specifying variable to
use for coloring}

\item{color.pal}{Colors to use for color.variable}

\item{y.breaks}{Set breaks for y-axis}

\item{y.limits}{Set limits for y-axis}

\item{plot}{TRUE or FALSE if TRUE plot is returned else a data frame}
}
\value{
A ggplot object.
}
\description{
Calculates and plots coverage of taxa from a phyloseq object.
}
\details{
Calculates and plots coverage of species from a phyloseq object and
         returns a ggplot object. Useful for checking co-existence over
         different timepoints. NOTE: Should be raw counts not relative. Adapted from
         Fukuyama, J et al. 2017 PLoS computational biology see ref below.
}
\examples{

data(SyncomFiltData)
fasting_cols <- c("#b3de69", "#fb8072", "#80b1d3")
pl <- taxa_coverage(SyncomFiltData,
  coverage = 0.9999,
  time.column = "Time_hr_num",
  shape.variable = "Acetate_Feed",
  color.variable = "Fasting",
  color.pal = fasting_cols,
  y.breaks = seq(0, 16, 1),
  y.limits = c(9, 16)
)
pl + facet_wrap(~StudyIdentifier)
pl
}
\references{
\itemize{
\item{}{"Fukuyama, J., Rumker, L., Sankaran, K., Jeganathan, P., Dethlefsen, L.,
Relman, D. A., & Holmes, S. P. (2017). Multidomain analyses of a longitudinal
human microbiome intestinal cleanout perturbation experiment.
PLoS computational biology, 13(8), e1005706.
\url{https://doi.org/10.1371/journal.pcbi.1005706}}
\item{}{To cite the package, see citation('syncomR')}
}
}
\author{
Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
}
\keyword{Anlaysis}
\keyword{and}
\keyword{visualization}
