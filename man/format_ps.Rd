% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_ps.R
\name{format_ps}
\alias{format_ps}
\title{Formatting the ASV-level Phyloseq Object}
\usage{
format_ps(ps, tax.level = NULL)
}
\arguments{
\item{ps}{a \code{\link{phyloseq-class}}}

\item{tax.level}{Taxonomic level to aggregate}
}
\value{
Filtered \code{\link{phyloseq-class}} object.
}
\description{
Convert ASV-level \code{\link{phyloseq-class}} object to species level data.
}
\details{
This is a utility function that converts ASV-level \code{\link{phyloseq-class}}
         to species level data. Usually the rownames are ASV IDs of Seqs. Here, we convert
         it to taxonomic level identities to make sense. It can be used for other taxonomic
         levels too. NOTE: Removes unassigned ASVs and returns a \code{\link{phyloseq-class}}
         object analysis.
}
\examples{
data(SyncomRawCounts)
pseq <- format_ps(SyncomRawCounts, tax.level = "Species")
print(pseq)
}
\author{
Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
}
\keyword{Data}
\keyword{and}
\keyword{filtering}
