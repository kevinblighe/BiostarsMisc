# Author: Kevin Blighe
# Contact: kevin.blighe@ranchobiosciences.com
# Company: Rancho Biosciences

  setwd('/Kev/CollegeWork/12ClinBio2015_/Projects/Rancho/BMS/CLL17p_Chris_Work')

  tn_unfilt <- data.table::fread('qc_kb/WGS_Merged.UnFiltered_Variant_Stats.tsv', header = TRUE, data.table = FALSE)
  t_unfilt <- data.table::fread('qc_kb/WGS_Merged.UnFiltered.TumorOnly_Variant_Stats.tsv', header = TRUE, data.table = FALSE)

  est <- data.frame(readxl::read_xlsx('Metadata/Esteban/mutations found in targeted sequencing KB.xlsx'))
  paste0(est$CHROM_hg38, ':', est$POS_hg38, ',', est$REF, '>', est$ALT)

  # filter for Esteban's mutation list
    tn_unfilt <- tn_unfilt[which(tn_unfilt$Variant %in% paste0(est$CHROM_hg38, ':', est$POS_hg38, ',', est$REF, '>', est$ALT)),]
    t_unfilt <- t_unfilt[which(t_unfilt$Variant %in% paste0(est$CHROM_hg38, ':', est$POS_hg38, ',', est$REF, '>', est$ALT)),]
    gc(); gc();

  merge <- rbind(tn_unfilt, t_unfilt)
  merge <- merge[!duplicated(merge$Variant),]
  write.table(
    merge,
    'WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv',
    sep = '\t', row.names = FALSE)
  saveRDS(merge, 'WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.Rds')

