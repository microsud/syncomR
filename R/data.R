#' @title SynComm raw ASV data
#' @description The SyncomRawCounts data set contains high-through taxonomic profiling data from
#' 210 samples taken from the syncom experiment and 5269 ASV-level counts. The data set includes control samples
#' samples from feed, reagent and DNA extraction procedure.
#' For details, see \url{XXXX}.
#' @name SyncomRawCounts
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SyncomRawCounts)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#' \itemize{
#' \item HPLC: Quantification of metabolites majorly short chain fatty acids
#' \item microbes: 16S rRNA gene counts of 5269 ASV-level taxa across 210 samples
#' \item meta: Sample metadata including time point, StudyIdentifier, Perturbations, etc
#' \item phyloseq The syn comm data set converted into a
#' \code{\link{phyloseq-class}} object.
#' }
#'
#' @keywords Data
NULL


#' @title SynComm Copy Number Corrected Species-level data
#' @description The SyncomCopyCorrectedCounts data set contains high-through taxonomic profiling data from
#' 210 samples taken from the syncom experiment and 16 species-level counts. The data set includes control samples
#' samples from feed, reagent and DNA extraction procedure. For raw ASV-level counts check SyncomRawCounts
#' For details, see \url{XXXX}.
#' @name SyncomCopyCorrectedCounts
#' @return Loads the data set in R.
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SyncomCopyCorrectedCounts)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#' \itemize{
#' \item HPLC: Quantification of metabolites majorly short chain fatty acids
#' \item microbes: 16S rRNA gene counts corrected for copy numbers of 16 species-level counts across 210 samples
#' \item meta: Sample metadata including time point, StudyIdentifier, Perturbations, etc
#' \item phyloseq The syn comm data set converted into a
#' \code{\link{phyloseq-class}} object.
#' }
#'
#' @keywords data
NULL


#' @title SynComm Filtered Copy Number Corrected Species-level data
#' @description The SyncomFiltData data set contains high-through taxonomic profiling data from
#' 183 samples taken from the syncom experiment and 16 species-level counts. The data set contains subset of data with
#' only time points sampled for each bioreacotr. No control samples
#' samples from feed, reagent and DNA extraction procedure, for this check either SyncomRawCounts or SyncomCopyCorrectedCounts.
#' For details, see \url{XXXX}.
#' @name SyncomFiltData
#' @return Loads the data set in R.
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SyncomFiltData)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#' \itemize{
#' \item HPLC: Quantification of metabolites majorly short chain fatty acids
#' \item microbes: 16S rRNA gene counts corrected for copy numbers of 16 species-level counts across 207 samples
#' \item meta: Sample metadata including time point, StudyIdentifier, Perturbations, etc
#' \item phyloseq The syn comm data set converted into a
#' \code{\link{phyloseq-class}} object.
#' }
#'
#' @keywords data
NULL

#' @title SynComm raw locus tag count from raw read processing
#' @description The `SynComRawRNA` consist of raw counts from RNAseq for 36 samples (51138 locus tags)
#' stored as data.frame object and made available in `syncomR` package. For details, see \url{XXXX}.
#' @name SynComRawRNA
#' @return Loads the data set in R.
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SynComRawRNA)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#'
#' @keywords data
NULL

#' @title SynComm KEGG linking of SynComRawRNA
#' @description The `SynComRNAKEGG` consists 36 samples and locus tags linked to KEGG IDs with following columns
#' KO, LocusTag, , BacterialStrain, GeneName, Level_1, Level_2, Level_3, KOID, KOGeneName
#' stored as data.frame object and made available in `syncomR` package. For details, see \url{XXXX}.
#' @name SynComRNAKEGG
#' @return Loads the data set in R.
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SynComRNAKEGG)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#'
#' @keywords data
NULL



#' @title SyncomGMM Gut Metabolic Modules (GMMs)
#' @description The `SyncomGMM` consists 36 samples and abundances of GMMs linked with taxa. This was obtained by running
#' the omixerRpm on KEGG KO abundances and stored as data.frame object and made available in `syncomR` package. For details, see \url{XXXX}.
#' @name SyncomGMM
#' @return Loads the data set in R.
#' @docType data
#' @author Sudarshan Shetty \email{sudarshanshetty9@gmail.com}
#' @references Shetty SA et al. (20XX) XXX 1:e32
#' \url{http://XXX/XXX/XXX}
#' @usage data(SyncomGMM)
#' @return Loads the data set in R.
#' @format List of the following data matrices as described in detail in
#' Shetty SA et al. (20XX):
#'
#' @keywords data
NULL
