
# Author: Kevin Blighe
# Contact: kevin.blighe@ranchobiosciences.com
# Company: Rancho Biosciences

  setwd('/Kev/CollegeWork/12ClinBio2015_/Projects/Rancho/BMS/CLL17p_Chris_Work')

  merge <- readRDS('qc_kb/WGS_output/merge.AnyImpact.Rds')
  merge$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge$HetSamples))
  merge.tumoronly <- readRDS('qc_kb/WGS_output/merge.TumorOnly.AnyImpact.Rds')
  merge.tumoronly$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge.tumoronly$HetSamples))
  # Before just binding the datasets, filter out any tumor-only variants already found in T-N paireds
  # any already present, will have to be intelligently merged for any samples where the variant was not previously identified.
    idx1 <- which(merge.tumoronly$Variant %in% merge$Variant)
    idx2 <- which(merge$Variant %in% merge.tumoronly$Variant)
    df.tumoronly <- merge.tumoronly[idx1,]
    df <- merge[idx2,]
    merge.tumoronly <- merge.tumoronly[-idx1,]
    merge <- merge[-idx2,]
    df$Variant == df.tumoronly$Variant
    hetsamples <- c()
    nhet <- c()
    for (i in 1:nrow(df)) {
      vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
        unlist(strsplit(df.tumoronly$HetSamples[i], ',')))))
      nhet <- c(nhet, length(vec))
      hetsamples <- c(hetsamples, paste(vec, collapse = ','))
    }
    df$nHet <- nhet
    df$HetSamples <- hetsamples
  merge <- rbind(df, merge, merge.tumoronly)

  write.table(subset(merge, Gene == 'TP53'), 'TP53mutations.AnyImpact.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

