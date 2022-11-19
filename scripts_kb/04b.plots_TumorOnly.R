# Author: Kevin Blighe
# Contact: kevin.blighe@ranchobiosciences.com
# Company: Rancho Biosciences

  setwd('/Kev/CollegeWork/12ClinBio2015_/Projects/Rancho/BMS/CLL17p_Chris_Work')

  require(ggplot2)
  require(ggrepel)
  require(ComplexHeatmap)
  library(trackViewer)
  library(RColorBrewer)

  vars <- data.table::fread('qc_kb/WGS_Filtered.TumorOnly_Variant_Stats.tsv', header = TRUE, data.table = FALSE)
  colnames(vars) <- make.names(colnames(vars))
  vars$High.Moderate.Impact <- vars$High.Impact + vars$Medium.Impact
  vars[,1] <- sub('WGS\\/FilteredData\\.TumorOnly\\/', '', sub('\\-null\\.filt\\.vcf\\.gz', '', vars[,1]))
  vars.N <- vars[grep('\\-N$', vars[,1]),]
  vars.T <- vars[grep('\\-T$', vars[,1]),]

  par(mar = c(4,4,4,4), mfrow = c(2,4), cex = 1.0)
  p1n <- hist(vars.N$SNVs, breaks = 10, plot = FALSE)
  p2n <- hist(vars.N$Protein.Coding.Region, breaks = 10, plot = FALSE)
  p3n <- hist(vars.N$Missense.Variants, breaks = 10, plot = FALSE)
  p4n <- hist(vars.N$High.Moderate.Impact, breaks = 10, plot = FALSE)
  p1t <- hist(vars.T$SNVs, breaks = 10, plot = FALSE)
  p2t <- hist(vars.T$Protein.Coding.Region, breaks = 10, plot = FALSE)
  p3t <- hist(vars.T$Missense.Variants, breaks = 10, plot = FALSE)
  p4t <- hist(vars.T$High.Moderate.Impact, breaks = 10, plot = FALSE)

  #plot(p1n, col = 'grey', xlab = NULL, main = NULL)
  plot(p2n, col = 'royalblue', xlab = NULL, main = NULL, add = FALSE)
  plot(p3n, col = 'forestgreen', xlab = NULL, main = NULL, add = FALSE)
  legend('topright', cex = 1.0, c('SNVs', 'Protein Coding', 'Missense', 'High | Moderate Impact'),
    fill = c('grey', 'royalblue', 'forestgreen', 'red2'), bty = 'n')
  plot(p4n, col = 'red2', xlab = NULL, main = paste0('Date generated: ', Sys.Date()), add = FALSE)
  #plot(p1t, col = 'grey', xlab = NULL, main = NULL)
  plot(p2t, col = 'royalblue', xlab = NULL, main = NULL, add = FALSE)
  plot(p3t, col = 'forestgreen', xlab = NULL, main = NULL, add = FALSE)
  legend('topright', cex = 1.0, c('SNVs', 'Protein Coding', 'Missense', 'High | Moderate Impact'),
    fill = c('grey', 'royalblue', 'forestgreen', 'red2'), bty = 'n')
  plot(p4t, col = 'red2', xlab = NULL, main = paste0('Date generated: ', Sys.Date()), add = FALSE)

  #median(vars.N$SNVs); sd(vars.N$SNVs)
  median(vars.N$Protein.Coding.Region); sd(vars.N$Protein.Coding.Region)
  median(vars.N$Missense.Variants); sd(vars.N$Missense.Variants)
  median(vars.N$High.Moderate.Impact); sd(vars.N$High.Moderate.Impact)
  #median(vars.T$SNVs); sd(vars.T$SNVs)
  median(vars.T$Protein.Coding.Region); sd(vars.T$Protein.Coding.Region)
  median(vars.T$Missense.Variants); sd(vars.T$Missense.Variants)
  median(vars.T$High.Moderate.Impact); sd(vars.T$High.Moderate.Impact)



  files <- list.files('qc_kb/', pattern = '*.TumorOnly.DP.txt', full.names = TRUE)
  files <- files[grep('\\-T', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  for (i in 1:length(dp)) {
    dp[[i]]$Sample <- rep(names(dp)[i], nrow(dp[[i]]))
    dp[[i]] <- dp[[i]][,c('Read.Depth','Sample')]
  }
  ggdata <- do.call(rbind, dp)
  p1 <- ggplot(data = subset(ggdata, Sample %in% names(dp)[1:20]), aes(x = Sample, y = Read.Depth)) +
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
  p2 <- ggplot(data = subset(ggdata, Sample %in% names(dp)[21:39]), aes(x = Sample, y = Read.Depth)) +
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



  files <- list.files('qc_kb/', pattern = '*.TumorOnly.DP.txt', full.names = TRUE)
  files <- files[grep('\\-T', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  merge <- read.table('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.NoNormal.tsv', sep = '\t', header = TRUE)

  # filter for 50 driver genes from Landau et al. (https://www.nature.com/articles/nature15395), any impact type
    landau <- c('SF3B1','ATM','TP53','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
      'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
      'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
      'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')
    # add in drivers from Knisbacher et al. (2022
      knisbacher <- data.frame(readxl::read_xlsx('Doc/41588_2022_1140_MOESM1_ESM.xlsx', sheet = 'Driver_Genes'))[,6]
      knisbacher <- knisbacher[!is.na(knisbacher)]
    #drivers <- sort(unique(c(landau, knisbacher)))
    drivers <- sort(unique(c(landau)))
    write.table(
      merge[which(merge$Gene %in% drivers),],
      'WGS_output/WGS_Merged.Filtered.TumorOnly_Variant_Stats.NoNormal.Drivers.AnyImpact.tsv',
      sep = '\t', row.names = FALSE)

  # filter for drivers, high|moderate impact
    saveRDS(merge, 'WGS_output/merge.TumorOnly.AnyImpact.Rds')
    merge <- subset(merge, Impact == 'HIGH|MODERATE')
    merge <- merge[which(merge$Gene %in% drivers),]
    write.table(merge, 'WGS_output/WGS_Merged.Filtered.TumorOnly_Variant_Stats.NoNormal.Drivers.High_ModImpact.tsv',
      sep = '\t', row.names = FALSE)

  saveRDS(merge, 'WGS_output/merge.TumorOnly.Rds')

  ggdata <- reshape2::melt(sort(table(merge$Gene)))
  colnames(ggdata) <- c('Gene', 'Count')
  ggdata$Gene <- factor(ggdata$Gene, levels = ggdata$Gene)
  ggplot(aes(x = Gene, y = Count, label = Gene), data = subset(ggdata, Count > 0)) +
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
    coord_cartesian(ylim = c(0, 40)) +
    scale_y_continuous(breaks = seq(0, 40, 2)) +
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
  topgenes <- subset(merge, nHet > 0)$Gene
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
    remove_empty_columns = FALSE,
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
      title = 'Mutation type',
      at = c('Somatic', 'TP53'),
      labels = c('Somatic', 'TP53'),
      nrow = 1,
      title_position = 'topcenter'))
  draw(geneticprint, heatmap_legend_side = 'top', annotation_legend_side = 'bottom', newpage = TRUE)



  files <- list.files('qc_kb/', pattern = '*.TumorOnly.DP.txt', full.names = TRUE)
  files <- files[grep('\\-T', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  merge <- read.table('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.NoNormal.tsv', sep = '\t', header = TRUE)
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
    #row_order = NULL,
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
      title = 'Mutation type',
      at = c('Somatic', 'TP53'),
      labels = c('Somatic', 'TP53'),
      nrow = 1,
      title_position = 'topcenter'))
  draw(geneticprint, heatmap_legend_side = 'top', annotation_legend_side = 'right', newpage = TRUE)



  # retrieve co-ordinates from tiny 'Download Track Data' button at: https://www.ncbi.nlm.nih.gov/gene/7157
  merge <- read.table('qc_kb/WGS_Merged.Filtered.TumorOnly_Variant_Stats.NoNormal.tsv', sep = '\t', header = TRUE)
  merge <- subset(merge, Gene == 'TP53' & Impact == 'HIGH|MODERATE')

  features <- GRanges('chr17',
    IRanges(c(7668421, 7670609, 7673535, 7673701, 7674181, 7674859, 7675053, 7675994, 7676382, 7676521, 7687377),
      width = c(1269, 106, 73, 136, 109, 112, 183, 278, 21, 101, 113),
      names = paste0('exon ', 1:11)),
    fill = brewer.pal(11, "Spectral"),
    height = rep(0.1, 11))
  SNP <- as.numeric(sub(',[A-Z,]*>[A-Z,]*$', '', sub('chr17:', '',  merge$Variant)))
  TP53 <- GRanges('chr17', IRanges(SNP, width = 1, names = merge$Variant),
    color = sample.int(length(SNP), length(SNP), replace = TRUE),
    score = merge$nHet)
  lolliplot(TP53, features)

 
