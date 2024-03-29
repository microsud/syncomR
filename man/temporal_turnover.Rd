% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporal_turnover.R
\name{temporal_turnover}
\alias{temporal_turnover}
\title{Measure and Plot Temporal Turnover of Community}
\usage{
temporal_turnover(
  ps,
  time.col = NULL,
  tree = NULL,
  method = c("bray", "euclidean", "unifrac", "canberra"),
  compositional = c(TRUE, FALSE),
  compared.to = "start",
  plot = TRUE
)
}
\arguments{
\item{ps}{a OTU or ASV or Species abundance table}

\item{time.col}{Column specifying time variable. Must be numeric.}

\item{tree}{phylogenetic tree if interested in unifrac.}

\item{method}{Can be any of the following c('bray','euclidean','unifrac', 'canberra')}

\item{compositional}{TRUE or FALSE If relative abundance to be used as input for dist calculations}

\item{compared.to}{final, subsequent or start}

\item{plot}{TRUE or FALSE if TRUE plot is returned else a data frame}
}
\value{
Plot of community dissimilarity over time.
}
\description{
Calculates dissimilarity overtime compared to either starting,
             subsequent or final timepoint.
}
\details{
Dissimilarity between subsequent samples can help in identifying succession
         of the community. Dissimilarity between longitudinal samples and start or final
         timepoint can help in identifying long term succession of the community.
}
\examples{
data(SyncomFiltData)
ps1.b5 <- subset_samples(
  SyncomFiltData,
  StudyIdentifier == "Bioreactor A"
)
p <- temporal_turnover(ps1.b5,
  tree = NULL, time.col = "Time_hr_num",
  method = "canberra", compositional = TRUE, compared.to = "start"
)
p + theme_syncom() #+ geom_label_repel()
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
\keyword{Analysis}
\keyword{and}
\keyword{visualization}
