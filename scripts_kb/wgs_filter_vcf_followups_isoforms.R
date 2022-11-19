
  setwd('/shared/CLL')

  # tumor-only
    var <- data.table::fread('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.IsoformLevel.tsv', data.table = FALSE)
    var <- subset(var, Impact1 == 'HIGH|MODERATE')

    for (i in 1:nrow(var)) {
      x <- unlist(strsplit(var[i,'Impact2'], ';'))

      idx_high <- which(x %in% 'HIGH')
      idx_moderate <- which(x %in% 'MODERATE')
      idx_modifier <- which(x %in% 'MODIFIER')
      idx_low <- which(x %in% 'LOW')
      idx <- c(idx_high, idx_moderate, idx_modifier, idx_low)

      var[i,'Impact2,'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Function'], ';'))
      var[i,'Function'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Symbol'], ';'))
      var[i,'Symbol'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Ensembl Gene'], ';'))
      var[i,'Ensembl Gene'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Ensembl Transcript'], ';'))
      var[i,'Ensembl Transcript'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Biotype'], ';'))
      var[i,'Biotype'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'HGVS'], ';'))
      var[i,'HGVS'] <- paste(x[idx], collapse = ';')

      if (i %% 1000 == 0) {
        message('--', i, ' of ', nrow(var))
      }
     }

    write.table(var,
      'qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.IsoformLevel.HighModerateOrdered.tsv',
      sep = '\t', quote = FALSE, row.names = FALSE)

  # paired
    var <- data.table::fread('qc_kb/WGS_Merged.Filtered_Variant_Stats.IsoformLevel.tsv', data.table = FALSE)
    var <- subset(var, Impact1 == 'HIGH|MODERATE')

    for (i in 1:nrow(var)) {
      x <- unlist(strsplit(var[i,'Impact2'], ';'))

      idx_high <- which(x %in% 'HIGH')
      idx_moderate <- which(x %in% 'MODERATE')
      idx_modifier <- which(x %in% 'MODIFIER')
      idx_low <- which(x %in% 'LOW')
      idx <- c(idx_high, idx_moderate, idx_modifier, idx_low)

      var[i,'Impact2,'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Function'], ';'))
      var[i,'Function'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Symbol'], ';'))
      var[i,'Symbol'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Ensembl Gene'], ';'))
      var[i,'Ensembl Gene'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Ensembl Transcript'], ';'))
      var[i,'Ensembl Transcript'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'Biotype'], ';'))
      var[i,'Biotype'] <- paste(x[idx], collapse = ';')

      x <- unlist(strsplit(var[i,'HGVS'], ';'))
      var[i,'HGVS'] <- paste(x[idx], collapse = ';')

      if (i %% 1000 == 0) {
        message('--', i, ' of ', nrow(var))
      }
     }

    write.table(var,
      'qc_kb/WGS_Merged.Filtered_Variant_Stats.IsoformLevel.HighModerateOrdered.tsv',
      sep = '\t', quote = FALSE, row.names = FALSE)