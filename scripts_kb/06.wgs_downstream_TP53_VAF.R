
  require(ggplot2)
  require(ggrepel)

  # TP53, tumour
    vars1 <- data.table::fread('qc_kb/WGS_Merged.Filtered_Variant_Stats.VAF.tsv', header = TRUE, data.table = FALSE)
    vars1 <- subset(vars1, Impact == 'HIGH|MODERATE')

    vars2 <- data.table::fread('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.VAF.tsv', header = TRUE, data.table = FALSE)
    vars2 <- vars2[-grep('\\-N', vars2[,'HetSamples']),]
    vars2 <- subset(vars2, Impact == 'HIGH|MODERATE')

    ggdata <- data.frame(
      AF = sort(c(
        as.numeric(do.call(c, strsplit(subset(vars1, Gene == 'TP53')[,'AF'], ','))),
        as.numeric(do.call(c, strsplit(subset(vars2, Gene == 'TP53')[,'AF'], ','))))) * 100)
    unique(sort(c(
      sub('\\-T', '', do.call(c, strsplit(subset(vars1, Gene == 'TP53')[,'HetSamples'], ','))),
      sub('\\-T', '', do.call(c, strsplit(subset(vars2, Gene == 'TP53')[,'HetSamples'], ','))))))

    ggdata$label <- make.unique(as.character(ggdata$AF))
    ggdata$label <- factor(ggdata$label, levels = rev(ggdata$label))
    ggplot(aes(x = label, y = AF, label = AF), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'royalblue') +
      theme_minimal() +
      #geom_label_repel(max.overlaps = 100) +
      theme(
        plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 16, vjust = 0.5),
        axis.title = element_text(size = 18, face = 'bold')) +
      xlab(NULL) + ylab('Allele Frequency (%)') +
      labs(
        title = 'TP53 somatic mutation allele frequencies',
        subtitle = paste0('Date generated: ', Sys.Date()),
        caption = paste0(
          'min, ', min(ggdata$AF), '%\n',
          'max, ', max(ggdata$AF), '%\n',
          'mean, ', round(mean(ggdata$AF), digits = 2), '%')) +
      coord_cartesian(ylim = c(0, 100)) +
      scale_y_continuous(breaks = seq(0, 100, 5))

  # TP53, normal
    vars2 <- data.table::fread('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.VAF.tsv', header = TRUE, data.table = FALSE)
    vars2 <- vars2[-grep('\\-T|_CD19', vars2[,'HetSamples']),]
    vars2 <- subset(vars2, Impact == 'HIGH|MODERATE')

    write.table(subset(vars2, Gene == 'TP53'), 'WGS TumorOnly NormalVariants TP53.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

    ggdata <- data.frame(
      AF = sort(c(
        as.numeric(do.call(c, strsplit(subset(vars2, Gene == 'TP53')[,'AF'], ','))))) * 100)
    unique(sort(c(
      sub('\\-N', '', do.call(c, strsplit(subset(vars2, Gene == 'TP53')[,'HetSamples'], ','))))))

    ggdata$label <- make.unique(as.character(ggdata$AF))
    ggdata$label <- factor(ggdata$label, levels = rev(ggdata$label))
    ggplot(aes(x = label, y = AF, label = AF), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'royalblue') +
      theme_minimal() +
      #geom_label_repel(max.overlaps = 100) +
      theme(
        plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 16, vjust = 0.5),
        axis.title = element_text(size = 18, face = 'bold')) +
      xlab(NULL) + ylab('Allele Frequency (%)') +
      labs(
        title = 'TP53 variants in normal samples, allele frequencies',
        subtitle = paste0('Date generated: ', Sys.Date()),
        caption = paste0(
          'min, ', min(ggdata$AF), '%\n',
          'max, ', max(ggdata$AF), '%\n',
          'mean, ', round(mean(ggdata$AF), digits = 2), '%')) +
      coord_cartesian(ylim = c(0, 100)) +
      scale_y_continuous(breaks = seq(0, 100, 5))

