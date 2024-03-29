% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/short_term_change.R
\name{short_term_change}
\alias{short_term_change}
\title{Short term Changes in Abundance}
\usage{
short_term_change(
  ps,
  rarefy = FALSE,
  compositional = FALSE,
  depth = NULL,
  plot.type = NULL,
  abbreviation = T
)
}
\arguments{
\item{ps}{a non rarefied raw OTU or ASV or species abundance table longitudinal data.}

\item{rarefy}{TRUE  Kept this option in case users have absolute quantification.}

\item{compositional}{Logical TRUE or FALSE.}

\item{depth}{If rarefy is TRUE then depth for rarefaction.}

\item{plot.type}{If data to be plotted either line or polar else NULL.}

\item{abbreviation}{For plotting whether to abbreviate names of taxa TRUE or FALSE.}
}
\value{
A data frame with mean growth, max growth or a plot with labeling
        for max \emph{change in abundance} timepoints.
}
\description{
Calculates short term changes in abundance of taxa
             using temporal Abundance data.
}
\details{
This approach is used by Wisnoski NI and colleagues
         \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
         the following calculation log(present abundance/past abundance).
         Also a compositional version using relative abundance similar to
         Brian WJi, Sheth R et al
         \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
         This approach is useful for identifying short term growth behaviors of taxa.
}
\examples{
\dontrun{
data(SyncomFiltData)
short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
syncom_ps <- subset_samples(SyncomFiltData, !(Time_label \%in\% short_time_labels))
bioA <- subset_samples(syncom_ps, StudyIdentifier == "Bioreactor A")
bioA.lg <- add_time_lag(bioA)
bioA.lg <- subset_samples(bioA.lg, time_lag >= 4)
p.bioA <- short_term_change(bioA.lg,
  rarefy = TRUE,
  depth = min(sample_sums(bioA.lg)),
  plot.type = "polar"
)
p.g <- p.bioA +
  scale_color_manual(values = syncom_colors2("BacterialSpecies")) +
  ggtitle("Bioreactor A")
print(p.g)
}
}
\references{
\itemize{
\item{}{Shetty SA. et al. A Minimalist Approach for Deciphering the
        Ecophysiology of Human Gut Microbes 2020}
\item{}{To cite the package, see citation('syncomR')}
}
}
\author{
Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
}
\keyword{Analysis}
\keyword{and}
\keyword{visualization}
