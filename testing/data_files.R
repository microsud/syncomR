# Hardrive G:\Post_doc\MM_MDb-MM\symcommR-master\03_amplicon\output\rds
SyncomRawCounts <- readRDS("ps_fr.rds")
sample_data(SyncomRawCounts)$StudyIdentifier <- gsub("Fermentor_5", "Bioreactor A",sample_data(SyncomRawCounts)$StudyIdentifier)
sample_data(SyncomRawCounts)$StudyIdentifier <- gsub("Fermentor_6", "Bioreactor B",sample_data(SyncomRawCounts)$StudyIdentifier)
sample_data(SyncomRawCounts)$StudyIdentifier <- gsub("Fermentor_8", "Bioreactor C",sample_data(SyncomRawCounts)$StudyIdentifier)

save(SyncomRawCounts, file = "SyncomRawCounts.rda")

SyncomCopyCorrectedCounts <- readRDS("CP_corrected_ps1.rds")
sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier <- gsub("Fermentor_5", "Bioreactor A",sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier)
sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier <- gsub("Fermentor_6", "Bioreactor B",sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier)
sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier <- gsub("Fermentor_8", "Bioreactor C",sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier)

save(SyncomCopyCorrectedCounts, file = "SyncomCopyCorrectedCounts.rda")

controls <- c("FeedBottle", "PCRcontrol", "PositiveSeqControl", "ReagentControlDNAkit")
SyncomFiltData <- prune_samples( !(sample_data(SyncomCopyCorrectedCounts)$StudyIdentifier %in% controls), SyncomCopyCorrectedCounts)
#ps1.sub <- readRDS("inst/extdata/ps1.sub.rds")
remove_T80 <- c("Ferm_1_5_80", "Ferm_1_6_80", "Ferm_1_8_80")
SyncomFiltData <- prune_samples( !(sample_names(SyncomFiltData) %in% remove_T80), SyncomFiltData)
sample_data(SyncomFiltData)$Time_hr_num <- as.numeric(sample_data(SyncomFiltData)$Time_hr)
save(SyncomFiltData, file = "SyncomFiltData.rda")
