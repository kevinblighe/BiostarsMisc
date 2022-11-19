# Author: Kevin Blighe
# Contact: kevin.blighe@ranchobiosciences.com
# Company: Rancho Biosciences

  setwd('/Kev/CollegeWork/12ClinBio2015_/Projects/Rancho/BMS/CLL17p_Chris_Work')
  dir.create('WGS_output/', showWarnings = FALSE)

  require(ggplot2)
  require(ggrepel)
  require(ComplexHeatmap)
  library(trackViewer)
  library(RColorBrewer)

  # generate histogram of coverage, and other coverage metrics
    qc <- read.table('qc_kb/WGS_picard.qcstats', header = TRUE, sep = '\t')
    qc <- qc[order(qc$MEAN_COVERAGE, decreasing = TRUE),]
    qc$SAMPLE <- factor(qc$SAMPLE, levels = qc$SAMPLE)
    ggplot(qc, aes(x = MEAN_COVERAGE)) + 
      geom_histogram(aes(y =..density..),
        breaks = seq(30, 60, by = 3), 
        colour = 'black', 
        fill = 'white') +
      stat_function(fun = dnorm,
        args = list(
          mean = mean(qc$MEAN_COVERAGE),
          sd = sd(qc$MEAN_COVERAGE))) +
      theme_minimal() +
      theme(
        plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5, hjust = 0.0),#element_blank(),
        axis.text.y = element_blank(),#element_text(angle = 0, size = 16, vjust = 0.5),
        axis.title = element_text(size = 18, face = 'bold')) +
      xlab('Mean coverage') + ylab('Density') +
      labs(title = NULL, subtitle = NULL,
        caption = paste0('Date generated: ', Sys.Date())) +
      coord_cartesian(xlim = c(30, 60)) +
      scale_x_continuous(breaks = seq(30, 60, 5))

    qc <- read.table('qc_kb/WGS_picard.qcstats', header = TRUE, sep = '\t')
    mean(qc$PCT_1X); sd(qc$PCT_1X)
    mean(qc$PCT_5X); sd(qc$PCT_5X)
    mean(qc$PCT_10X); sd(qc$PCT_10X)
    mean(qc$PCT_15X); sd(qc$PCT_15X)
    mean(qc$PCT_20X); sd(qc$PCT_20X)
    mean(qc$PCT_25X); sd(qc$PCT_25X)
    mean(qc$PCT_30X); sd(qc$PCT_30X)
    mean(qc$PCT_40X); sd(qc$PCT_40X)
    mean(qc$PCT_50X); sd(qc$PCT_50X)
    mean(qc$PCT_60X); sd(qc$PCT_60X)
    mean(qc$PCT_70X); sd(qc$PCT_70X)
    mean(qc$PCT_80X); sd(qc$PCT_80X)
    mean(qc$PCT_90X); sd(qc$PCT_90X)
    mean(qc$PCT_100X); sd(qc$PCT_100X)



  # stats for different classes of variants
    vars <- read.table('qc_kb/WGS_Filtered_Variant_Stats.tsv', header = TRUE, sep = '\t')
    vars$High.Moderate.Impact <- vars$High.Impact + vars$Medium.Impact
    par(mar = c(4,4,4,4), mfrow = c(2,2), cex = 1.0)
    #p1 <- hist(vars$SNVs, breaks = 25, plot = FALSE)
    p2 <- hist(vars$Protein.Coding.Region, breaks = 25, plot = FALSE)
    p3 <- hist(vars$Missense.Variants, breaks = 25, plot = FALSE)
    p4 <- hist(vars$High.Moderate.Impact, breaks = 25, plot = FALSE)

    #plot(p1, col = 'grey', xlab = NULL, main = NULL)
    plot(p2, col = 'royalblue', xlab = NULL, main = NULL)
    plot(p3, col = 'forestgreen', xlab = NULL, main = NULL, add = FALSE)
    legend('topright', cex = 1.0, c('Protein Coding', 'Missense', 'High | Moderate Impact'),
      fill = c('royalblue', 'forestgreen', 'red2'), bty = 'n')
    plot(p4, col = 'red2', xlab = NULL, main = paste0('Date generated: ', Sys.Date()), add = FALSE)

    #median(vars$SNVs); sd(vars$SNVs)
    median(vars$Protein.Coding.Region); sd(vars$Protein.Coding.Region)
    median(vars$Missense.Variants); sd(vars$Missense.Variants)
    median(vars$High.Moderate.Impact); sd(vars$High.Moderate.Impact)

    ggdata <- vars[order(vars$SNVs, decreasing = TRUE),]
    ggdata$Sample <- unlist(lapply(strsplit(ggdata$Sample, '-|_'), function(x) x[1]))
    ggdata$Sample <- factor(ggdata$Sample, levels = ggdata$Sample)
    p1 <- ggplot(aes(x = Sample, y = SNVs, label = Sample), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'grey') +
      theme_minimal() +
      theme(
        plot.title = element_text(angle = 0, size = 22, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_text(angle = 90, size = 10, face = 'plain', vjust = 0.5, hjust = 0.0),#element_blank(),
        axis.text.y = element_text(angle = 0, size = 13, vjust = 0.5),
        axis.title = element_text(size = 20, face = 'bold')) +
      xlab('') + ylab('Count') +
      labs(title = 'SNVs', subtitle = NULL,
        caption = '') +
      coord_cartesian(ylim = c(0, 10000)) +
      scale_y_continuous(breaks = seq(0, 10000, 1000))

    ggdata <- vars[order(vars$Protein.Coding.Region, decreasing = TRUE),]
    ggdata$Sample <- unlist(lapply(strsplit(ggdata$Sample, '-|_'), function(x) x[1]))
    ggdata$Sample <- factor(ggdata$Sample, levels = ggdata$Sample)
    p2 <- ggplot(aes(x = Sample, y = Protein.Coding.Region, label = Sample), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'royalblue') +
      theme_minimal() +
      theme(
        plot.title = element_text(angle = 0, size = 22, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_text(angle = 90, size = 10, face = 'plain', vjust = 0.5, hjust = 0.0),#element_blank(),
        axis.text.y = element_text(angle = 0, size = 13, vjust = 0.5),
        axis.title = element_text(size = 20, face = 'bold')) +
      xlab('') + ylab('Count') +
      labs(title = 'Protein coding', subtitle = NULL,
        caption = '') +
      coord_cartesian(ylim = c(0, 4000)) +
      scale_y_continuous(breaks = seq(0, 4000, 500))

    ggdata <- vars[order(vars$Missense.Variants, decreasing = TRUE),]
    ggdata$Sample <- unlist(lapply(strsplit(ggdata$Sample, '-|_'), function(x) x[1]))
    ggdata$Sample <- factor(ggdata$Sample, levels = ggdata$Sample)
    p3 <- ggplot(aes(x = Sample, y = Missense.Variants, label = Sample), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'forestgreen') +
      theme_minimal() +
      theme(
        plot.title = element_text(angle = 0, size = 22, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_text(angle = 90, size = 10, face = 'plain', vjust = 0.5, hjust = 0.0),#element_blank(),
        axis.text.y = element_text(angle = 0, size = 13, vjust = 0.5),
        axis.title = element_text(size = 20, face = 'bold')) +
      xlab('') + ylab('Count') +
      labs(title = 'Missense', subtitle = NULL,
        caption = '') +
      coord_cartesian(ylim = c(0, 60)) +
      scale_y_continuous(breaks = seq(0, 60, 5))

    ggdata <- vars[order(vars$High.Moderate.Impact, decreasing = TRUE),]
    ggdata$Sample <- unlist(lapply(strsplit(ggdata$Sample, '-|_'), function(x) x[1]))
    ggdata$Sample <- factor(ggdata$Sample, levels = ggdata$Sample)
    p4 <- ggplot(aes(x = Sample, y = High.Moderate.Impact, label = Sample), data = ggdata) +
      geom_bar(stat = 'identity', fill = 'red2') +
      theme_minimal() +
      theme(
        plot.title = element_text(angle = 0, size = 22, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
        plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
        axis.text.x = element_text(angle = 90, size = 10, face = 'plain', vjust = 0.5, hjust = 0.0),#element_blank(),
        axis.text.y = element_text(angle = 0, size = 13, vjust = 0.5),
        axis.title = element_text(size = 20, face = 'bold')) +
      xlab('') + ylab('Count') +
      labs(title = 'High | Moderate impact', subtitle = NULL,
        caption = paste0('Date generated: ', Sys.Date())) +
      coord_cartesian(ylim = c(0, 100)) +
      scale_y_continuous(breaks = seq(0, 100, 10))
    cowplot::plot_grid(p1, p2, p3, p4)



  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  for (i in 1:length(dp)) {
    dp[[i]]$Sample <- rep(names(dp)[i], nrow(dp[[i]]))
    dp[[i]] <- dp[[i]][,c('Read.Depth','Sample')]
  }
  ggdata <- do.call(rbind, dp)
  p1 <- ggplot(data = subset(ggdata, Sample %in% names(dp)[1:27]), aes(x = Sample, y = Read.Depth)) +
    geom_violin(
      stat = 'ydensity',
      position = 'dodge',
      draw_quantiles = NULL,
      trim = FALSE,
      scale = 'area',
      na.rm = TRUE,
      orientation = NA,
      aes(fill = Sample)) +
    #facet_grid(Sample ~ ., scales = 'free_y') +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.05, colour = 'black') +
    theme_bw(base_size=24) + theme(
      legend.position = 'none',
      legend.background=element_rect(),
      plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      axis.text.x = element_text(angle = 45, size = 14, face = 'bold', vjust = 1.4, hjust = 1.3),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 0, size = 14, vjust = 0.5),
      axis.title = element_text(size = 18),
      legend.key = element_blank(), #removes the border
      legend.key.size = unit(0.75, 'cm'), #Sets overall area/size of the legend
      legend.text = element_text(size = 18), #Text size
      title = element_blank(),#element_text(size = 12), #Title text size
      strip.text.x = element_text(size = 18, face = 'bold'),
      strip.text.y = element_text(size = 18, face = 'bold', margin = margin(), angle = 90),
      strip.background = element_rect(fill = 'white', colour = 'white'), 
      strip.text = element_text(size = 18, face = 'bold', colour = 'black', margin = margin()),
      strip.switch.pad.grid = unit(0, 'cm')) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    xlab('') + ylab('Position Read Depth') +
    #ylim(2, 10.5) +
    labs(
      title = NULL,
      subtitle = NULL,
      caption = NULL) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    coord_cartesian(ylim = c(0, 300)) +
    scale_y_continuous(breaks = seq(0, 300, 50))
  p2 <- ggplot(data = subset(ggdata, Sample %in% names(dp)[28:54]), aes(x = Sample, y = Read.Depth)) +
    geom_violin(
      stat = 'ydensity',
      position = 'dodge',
      draw_quantiles = NULL,
      trim = FALSE,
      scale = 'area',
      na.rm = TRUE,
      orientation = NA,
      aes(fill = Sample)) +
    #facet_grid(Sample ~ ., scales = 'free_y') +
    stat_summary(
      geom = 'crossbar',
      width = 0.8,
      fatten = 2,
      color = 'black',
      fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
    geom_jitter(position=position_jitter(width = 0.4), size = 0.05, colour = 'black') +
    theme_bw(base_size=24) + theme(
      legend.position = 'none',
      legend.background=element_rect(),
      plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      axis.text.x = element_text(angle = 45, size = 14, face = 'bold', vjust = 1.4, hjust = 1.3),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 0, size = 14, vjust = 0.5),
      axis.title = element_text(size = 18),
      legend.key = element_blank(), #removes the border
      legend.key.size = unit(0.75, 'cm'), #Sets overall area/size of the legend
      legend.text = element_text(size = 18), #Text size
      title = element_blank(),#element_text(size = 12), #Title text size
      strip.text.x = element_text(size = 18, face = 'bold'),
      strip.text.y = element_text(size = 18, face = 'bold', margin = margin(), angle = 90),
      strip.background = element_rect(fill = 'white', colour = 'white'), 
      strip.text = element_text(size = 18, face = 'bold', colour = 'black', margin = margin()),
      strip.switch.pad.grid = unit(0, 'cm')) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    xlab('') + ylab('Position Read Depth') +
    #ylim(2, 10.5) +
    labs(
      title = NULL,
      subtitle = NULL,
      caption = paste0('Date generated: ', Sys.Date())) +
    theme(axis.line = element_line(
      size = 1.0, colour = 'black'),
      panel.border = element_blank(),
      panel.background = element_blank()) +
    coord_cartesian(ylim = c(0, 300)) +
    scale_y_continuous(breaks = seq(0, 300, 50))
  #cowplot::plot_grid(p1, p2, nrow = 2)
  p1
  p2



  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  merge <- read.table('qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv', sep = '\t', header = TRUE)
  saveRDS(merge, 'WGS_output/merge.AnyImpact.Rds')
  merge <- subset(merge, Impact == 'HIGH|MODERATE')
  saveRDS(merge, 'WGS_output/merge.Rds')

  ggdata <- reshape2::melt(sort(table(merge$Gene)))
  colnames(ggdata) <- c('Gene', 'Count')
  ggdata$Gene <- factor(ggdata$Gene, levels = ggdata$Gene)
  ggplot(aes(x = Gene, y = Count, label = Gene), data = subset(ggdata, Count > 2)) +
    geom_bar(stat = 'identity', fill = 'royalblue') +
    theme_minimal() +
    #geom_label_repel(max.overlaps = 10) +
    theme(
      plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 14, face = 'plain', vjust = 1),
      axis.text.x = element_text(angle = 0, size = 14, face = 'bold', vjust = 0.5, hjust = 0.0),#element_blank(),
      axis.text.y = element_text(angle = 0, size = 16, vjust = 0.5),
      axis.title = element_text(size = 18, face = 'bold')) +
    xlab(NULL) + ylab('Number of unique mutations') +
    labs(title = 'Number of unique High|Moderate impact mutations per gene', subtitle = NULL,
      caption = paste0('Date generated: ', Sys.Date())) +
    coord_cartesian(ylim = c(0, 15)) +
    scale_y_continuous(breaks = seq(0, 15, 1)) +
    coord_flip()

  genes <- sort(unique(merge$Gene))
  merge <- do.call(rbind,
    lapply(genes, function(x) {
      tab <- subset(merge, Gene == x)
      data.frame(
        Gene = tab$Gene[1],
        Impact = tab$Impact[1],
        nHet = sum(tab$nHet),
        HetSamples = paste(tab$HetSamples, collapse = ','))
    }))
  topgenes <- subset(merge, nHet > 2)$Gene
  saveRDS(topgenes, 'WGS_output/TN-PairedOncoPlot.Rds')
  onco <- data.frame(row.names = merge$Gene)
  names <- names(dp)
  for (i in 1:length(names)) {
    vec <- c()
    for (j in 1:nrow(merge)) {
      if (grepl(names[i], merge$HetSamples[j])) {
        vec <- c(vec, 'Somatic')
      } else {
        vec <- c(vec, '')
      }
    }
    onco <- cbind(onco, vec)
    colnames(onco)[i] <- names[i]
  }
  onco <- as.matrix(onco)
  rownames(onco) <- merge$Gene
  onco <- onco[topgenes,]
  alter_fun <- list(
    background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
     #grid.rect(x, y, w-unit(0.15, 'mm'), h-unit(0.15, 'mm'), gp = gpar(fill = 'white', col = NA))
    },
    Somatic = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, 'mm'), h-unit(1.5, 'mm'), gp = gpar(fill = 'red2', col = 'red2'))
    },
    TP53 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2.5, 'mm'), h-unit(2.5, 'mm'), gp = gpar(fill = 'black', col = 'black'))
    })
  cols <- c('Somatic' = 'red2', 'TP53' = 'black')
  ht_opt$message = TRUE
  geneticprint <- oncoPrint(
    onco,
    get_type = function(x) strsplit(x, ';')[[1]],
    name = 'Oncoprint',
    alter_fun = alter_fun,
    col = cols,
    #row_order = NULL,
    #column_order = ,
    remove_empty_columns = TRUE,
    row_title = 'Gene',
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16, fontface = 'plain'),
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10, fontface = 'plain'),
    row_names_max_width = unit(3, 'cm'),
    column_title = 'Samples',
    column_title_side = 'top',
    column_title_gp = gpar(fontsize = 16, fontface = 'plain'),
    column_title_rot = 0,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 10),
    pct_gp = gpar(fontsize = 10, fontface = 'plain', fill = 'white', col = 'white'),
    #bottom_annotation = annBlack,
    heatmap_legend_param = list(
      title = 'Mutation',
      at = c('Somatic', 'TP53'),
      labels = c('Somatic', 'TP53'),
      nrow = 1,
      title_position = 'topcenter'))
  draw(geneticprint, heatmap_legend_side = 'top', annotation_legend_side = 'bottom', newpage = TRUE)



  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  merge <- read.table('qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv', sep = '\t', header = TRUE)
  merge <- subset(merge, Gene == 'TP53' & Impact == 'HIGH|MODERATE')
  onco <- data.frame(row.names = merge$Variant)
  names <- names(dp)
  for (i in 1:length(names)) {
    vec <- c()
    for (j in 1:nrow(merge)) {
      if (grepl(names[i], merge$HetSamples[j])) {
        vec <- c(vec, 'TP53')
      } else {
        vec <- c(vec, '')
      }
    }
    onco <- cbind(onco, vec)
    colnames(onco)[i] <- names[i]
  }
  onco <- as.matrix(onco)

  # collapse the data to be just a single row
    onco <- t(data.frame(apply(onco, 2, function(x) ifelse(any(x == 'TP53'), 'TP53', ''))))
    rownames(onco) <- 'TP53'

  meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
  meta <- meta[match(colnames(onco), meta$patient),]
  all(colnames(onco) == meta$patient)

  idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
  meta <- meta[idx,]
  onco <- onco[,idx]

  # TP53 mutation
    tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)
    names(tp53) <- meta$patient

  ann <- HeatmapAnnotation(
    'p53 result\n(clinical)\n\n' = factor(tp53, levels = c('Abnormal', 'Normal', 'Unknown')),
    #c17p = anno_points(as.numeric(meta$c17p), gp = gpar(col = 'black'), ylim = c(0, 100), axis = TRUE, pch = '.', size = unit(5.0, 'mm')),
    col = list(
      'p53 result\n(clinical)\n\n' = c('Abnormal' = 'red', 'Normal' = 'royalblue', 'Unknown' = 'grey')),
    na_col = 'white',
    #annotation_height = c(2, 2, 2, 2, 2, 2, 2, 10),
    gap = unit(1.5, 'mm'),
    annotation_legend_param = list(
      'p53 result\n(clinical)\n\n' = list(title = 'p53 result')))

  alter_fun <- list(
    background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
     #grid.rect(x, y, w-unit(0.15, 'mm'), h-unit(0.15, 'mm'), gp = gpar(fill = 'white', col = NA))
    },
    Somatic = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, 'mm'), h-unit(1.5, 'mm'), gp = gpar(fill = 'red2', col = 'red2'))
    },
    TP53 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(2.5, 'mm'), h-unit(2.5, 'mm'), gp = gpar(fill = 'black', col = 'black'))
    })
  cols <- c('Somatic' = 'red2', 'TP53' = 'black')
  ht_opt$message = TRUE
  geneticprint <- oncoPrint(
    t(data.matrix(onco)),
    get_type = function(x) strsplit(x, ';')[[1]],
    name = 'Oncoprint',
    alter_fun = alter_fun,
    col = cols,
    row_order = NULL,
    column_order = match(names(onco), names(tp53)),
    remove_empty_columns = FALSE,
    row_title = expression(italic(TP53)~mutation),
    row_title_side = 'left',
    row_title_gp = gpar(fontsize = 16, fontface = 'plain'),
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 14, fontface = 'plain'),
    row_names_max_width = unit(3, 'cm'),
    column_title = 'Sample',
    column_title_side = 'bottom',
    column_title_gp = gpar(fontsize = 16, fontface = 'plain'),
    column_title_rot = 0,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 14),
    pct_gp = gpar(fontsize = 14, fontface = 'plain', fill = 'white', col = 'white'),
    top_annotation = ann,
    heatmap_legend_param = list(
      title = 'Mutation',
      at = c('Somatic', 'TP53'),
      labels = c('Somatic', 'TP53'),
      nrow = 1,
      title_position = 'topcenter'))
  draw(geneticprint, heatmap_legend_side = 'top', annotation_legend_side = 'right', newpage = TRUE)



  # retrieve co-ordinates from tiny 'Download Track Data' button at: https://www.ncbi.nlm.nih.gov/gene/7157
  merge <- read.table('qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv', sep = '\t', header = TRUE)
  merge <- subset(merge, Gene == 'TP53' & Impact == 'HIGH|MODERATE')
  features <- GRanges('chr17',
    IRanges(c(7668421, 7670609, 7673535, 7673701, 7674181, 7674859, 7675053, 7675994, 7676382, 7676521, 7687377),
      width = c(1269, 106, 73, 136, 109, 112, 183, 278, 21, 101, 113),
      names = paste0('exon ', 1:11)),
    fill = brewer.pal(11, "Spectral"),
    height = rep(0.1, 11))
  SNP <- as.numeric(sub(',[A-Z]*>[A-Z]*$', '', sub('chr17:', '',  merge$Variant)))
  TP53 <- GRanges('chr17', IRanges(SNP, width = 1, names = merge$Variant),
    color = sample.int(length(SNP), length(SNP), replace = TRUE),
    score = merge$nHet)
  lolliplot(TP53, features)

