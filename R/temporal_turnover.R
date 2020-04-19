#' @title Measure and plot long term temporal turnover of community
#'
#' @description Calculates dissimilarity overtime compared to either strating, subsequent or final timepoint.
#'
#' @details Dissimilarity between subsequent samples can help in identifying succession of the community.
#' Dissimilarity between longitudinal samples and start or final timepoint can help in identifying long term
#' succession of the community.
#' @param ps a otu/ASV/species abundance table
#' @param time.col Column specifying time variable. Must be numeric.
#' @param tree phylogenetic tree if interested in unifrac.
#' @param method Can be any of the following c('bray','euclidean','unifrac', 'canberra')
#' @param compositional TRUE or FALSE If relative abundance to be used as input for dist calculations
#' @param compared.to final, subsequent or start
#' @param plot TRUE or FALSE if TRUE plot is returned else a data frame
#' @return Plot of community dissimilalrity over time.
#' @references
#' \itemize{
#' \item{}{Guittar, J., Shade, A., & Litchman, E. (2019). Trait-based community assembly
#' and succession of the infant gut microbiome. Nature communications, 10(1), 512.}
#' \item{}{'Shetty SA et al (2019-2024)}
#' \item{}{To cite the package, see citation('syncomR')}
#' }
#' @examples
#' data(SyncomFiltData)
#' ps1.b5 <- subset_samples(SyncomFiltData, StudyIdentifier == "Bioreactor A")
#' p <- temporal_turnover(ps1.b5,
#'   tree = NULL, time.col = "Time_hr_num",
#'   method = "canberra", compositional = TRUE, compared.to = "start"
#' )
#' p + geom_label_repel() + theme_syncom()
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@gmail.com}
#' @import tidyr
#' @importFrom vegan vegdist
#' @export
#' @keywords Anlaysis and visualization
#'
temporal_turnover <- function(ps, time.col = NULL, tree = NULL,
                              method = c("bray", "euclidean", "unifrac", "canberra"),
                              compositional = c(TRUE, FALSE),
                              compared.to = "start",
                              plot = TRUE) {
  require(tidyr)
  require(vegan)
  require(ggrepel)
  if (is.null(time.col)) {
    stop("Please provide column name for time variable")
  }
  if (is.null(method)) {
    stop("Please provide method either 'bray','euclidean','unifrac', 'canberra'")
  }


  if (compositional == TRUE) {
    ps <- microbiome::transform(ps, "compositional")
  } else {
    ps <- ps
  }
  otus <- t(abundances(ps))
  metadf_b5 <- meta(ps)
  ps.df <- psmelt(ps)

  # head(otus)
  otus <- as.data.frame(otus)
  otus$Sample <- rownames(otus)
  otus$time <- metadf_b5[, time.col]

  j5 <- ps.df %>%
    spread(OTU, Abundance, fill = 0) %>%
    do(distance_melt(
      x = .$Sample,
      y = otus,
      method = method
    )) %>%
    mutate(method = "Taxa-based dissimilarity", data = "obs")


  if (compared.to == "subsequent") {
    j5a <- j5 %>%
      group_by(method, t1) %>%
      filter(t2 > t1) %>%
      filter(t2 == min(t2)) %>%
      mutate(group = "Dissimilarity to next sample")
    j5a$com <- paste0(j5a$t1, "-", j5a$t2)
    j5a$time.lag <- j5a$t2 - j5a$t1
  } else if (compared.to == "final") {
    j5a <- j5 %>%
      group_by(method, t1) %>%
      filter(t2 == max(t2) & t1 != max(t2)) %>%
      mutate(group = "Dissimilarity to final sample")
  } else if (compared.to == "start") {
    j5a <- j5 %>%
      group_by(method, t1) %>%
      filter(t2 == min(t2) & t1 != min(t2)) %>%
      mutate(group = "Dissimilarity to starting sample")
  }

  # j5b$StudyIdentifier <- "Bioreactor A"
  j5b <- j5a %>% arrange(t1)
  if (plot == TRUE) {
    p1 <-
      ggplot(j5b, aes(t1, dist)) +
      geom_line() +
      geom_point()
    p1 <-
      p1 + geom_smooth() + ylab(paste0("Community dissimilarity\nto ", compared.to, " sample(s)")) + xlab("Time (hr)")
    #+ geom_label_repel() +

    p1 <- p1 + scale_linetype_discrete(name = "")

    return(p1)
  } else {
    return(j5b)
  }
}
