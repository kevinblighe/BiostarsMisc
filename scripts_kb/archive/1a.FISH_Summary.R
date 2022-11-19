
samorder <- data.frame(readxl::read_xlsx('DA0000882/FISH_Summary.xlsx'))[-1,1]
df <- data.frame(readxl::read_xlsx('CLL-Clinical/project2final58_20211109.xlsx'))
df <- df[,which(colnames(df) %in% c('patient','fishcat','c17p','z_fish11q','z_fish12tri','z_fish13q','z_fish6','z_fish14'))]
df <- df[match(gsub('\\-T|_CD19', '', samorder), df$patient),]
df$patient == gsub('\\-T|_CD19', '', samorder)
write.table(df, 'DA0000882/FISH_Summary_Align.tsv', row.names = FALSE, quote = FALSE)
