
  require(ggplot2)

  data <- data.frame(readxl::read_xlsx('WGS/FISH_Summary.xlsx', skip = 2))
  colnames(data) <- gsub('\\.\\.', '', sub('\\.\\.\\.[0-9]*$', '', colnames(data)))
  colnames(data)[12:14] <- paste0('chr17p_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data)[15:17] <- paste0('chr11q_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data)[18:20] <- paste0('chr12_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data)[21:23] <- paste0('chr13q_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data)[24:26] <- paste0('chr6_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data)[27:29] <- paste0('chr14_', c('Locus', 'RawCN', 'CNploidy'))
  colnames(data) <- make.unique(colnames(data))

  data[,'z_fish11q'] <- factor(ifelse(data[,'z_fish11q'] == 1, 'chr11q', 'WT'), levels = c('WT','chr11q'))
  data[,'z_fish12tri'] <- factor(ifelse(data[,'z_fish12tri'] == 1, 'Trisomy chr12', 'WT'), levels = c('WT','Trisomy chr12'))
  data[,'z_fish13q'] <- factor(ifelse(data[,'z_fish13q'] == 1, 'chr13q', 'WT'), levels = c('WT','chr13q'))
  data[,'z_fish6'] <- factor(ifelse(data[,'z_fish6'] == 1, 'chr6', 'WT'), levels = c('WT','chr6'))
  data[,'z_fish14'] <- factor(ifelse(data[,'z_fish14'] == 1, 'chr14', 'WT'), levels = c('WT','chr14'))
  data$fishcat <- factor(sub('5\\=', '', data$fishcat))

  custom_theme <- theme_bw(base_size=24) + theme(
    legend.position = 'top',
    legend.background=element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle=element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption=element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(angle = 0, size = 16, vjust = 0.5),
    axis.title=element_text(size = 16),
    legend.key=element_blank(),
    legend.key.size=unit(0.5, 'cm'),
    legend.text=element_text(size = 14),
    title=element_blank())

  data$chr6_CNploidy <- ifelse(data$chr6_CNploidy < 0, 0, data$chr6_CNploidy)
  p1 <- ggplot(data = data, aes(x = z_fish6, y = chr6_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = z_fish6)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,6) +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  data$chr11q_CNploidy <- ifelse(data$chr11q_CNploidy < 0, 0, data$chr11q_CNploidy)
  p2 <- ggplot(data = data, aes(x = z_fish11q, y = chr11q_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = z_fish11q)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,6) +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  data$chr12_CNploidy <- ifelse(data$chr12_CNploidy < 0, 0, data$chr12_CNploidy)
  p3 <- ggplot(data = data, aes(x = z_fish12tri, y = chr12_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = z_fish12tri)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,6) +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  data$chr13q_CNploidy <- ifelse(data$chr13q_CNploidy < 0, 0, data$chr13q_CNploidy)
  p4 <- ggplot(data = data, aes(x = z_fish13q, y = chr13q_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = z_fish13q)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,6) +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  data$chr14_CNploidy <- ifelse(data$chr14_CNploidy < 0, 0, data$chr14_CNploidy)
  p5 <- ggplot(data = data, aes(x = z_fish14, y = chr14_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = z_fish14)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,6) +
    labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  data$chr17p_CNploidy <- ifelse(data$chr17p_CNploidy < 0, 0, data$chr17p_CNploidy)
  p6 <- ggplot(data = data, aes(x = fishcat, y = chr17p_RawCN)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      outlier.shape = 17,
      outlier.colour = 'red',
      outlier.size = 0.1,
      aes(fill = fishcat)) +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.6, colour="black") +
    custom_theme +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    xlab('') + ylab('Copy number') +
    ylim(0,85) +
    labs(title = NULL, subtitle = NULL, caption = 'Copy number estimates not adjusted for ploidy') +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank())

  cairo_pdf('WGS/Battenberg_x_FISH_Analysis.pdf', width = 10, height = 7)
    cowplot::plot_grid(
      p1, p2, p3, p4, p5, p6,
      ncol = 3)
  dev.off()

