% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance_melt.R
\name{distance_melt}
\alias{distance_melt}
\title{Calculate and Melt Distance Matrix}
\usage{
distance_melt(
  x,
  y,
  tree = NULL,
  method = c("bray", "euclidean", "manhattan", "unifrac", "canberra")
)
}
\arguments{
\item{x}{a OTU/ASV/species abundance table}

\item{y}{list of samples}

\item{tree}{phylogenetic tree if interested in unifrac.}

\item{method}{Can be any of the following c('bray','euclidean','unifrac', 'canberra'))}
}
\value{
Normalized or raw counts taxa abundance table (taxa are rows and timepoints columns).
}
\description{
This function calculates a user specified distance measure
             and returns a long data frame.
}
\references{
\itemize{
\item{}{Guittar, J., Shade, A., & Litchman, E. (2019). Trait-based community assembly
and succession of the infant gut microbiome. Nature communications, 10(1), 512.}
\item{}{'Shetty SA et al (2019-2024)}
\item{}{To cite the package, see citation('syncomR')}
}
}
\author{
Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
}
\keyword{Utilities}
