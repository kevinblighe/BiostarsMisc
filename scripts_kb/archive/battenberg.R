
  setwd('/shared/CLL/WGS/AnalyzedData/')

  cn <- rbind(
    readRDS('20200715/DA882_20200715_battenberg_hg19liftFromhg38_tidy_overlaps_20220525.rds'),
    readRDS('20210330/DA882_20210330_battenberg_hg19liftFromhg38_tidy_overlaps_20220525.rds'),
    readRDS('20210616/DA882_20210616_battenberg_hg19liftFromhg38_tidy_overlaps_20220525.rds'))
