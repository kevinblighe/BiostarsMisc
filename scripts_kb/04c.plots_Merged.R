
# Author: Kevin Blighe
# Contact: kevin.blighe@ranchobiosciences.com
# Company: Rancho Biosciences

  setwd('/Kev/CollegeWork/12ClinBio2015_/Projects/Rancho/BMS/CLL17p_Chris_Work')
  require(ggplot2)
  library(survminer)
  require(ComplexHeatmap)
  require(survival)

  merge <- readRDS('WGS_output/merge.Rds')
  merge$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge$HetSamples))
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

  # save IGH-V mutated samples for later
    #ighv <- merge[grep('IGHV', merge$Gene),]$HetSamples
    ighv <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))[,c(1,6)]
    ighv <- ighv[ifelse(is.na(ighv$z_mut), 'unknown', ighv$z_mut) == "0=Mutated",1]

  sum(merge$nHet)
  nrow(merge)

  ggdata <- reshape2::melt(sort(table(merge$Gene)))
  colnames(ggdata) <- c('Gene', 'Count')
  ggdata$Gene <- factor(ggdata$Gene, levels = ggdata$Gene)
  p1 <- ggplot(aes(x = Gene, y = Count, label = Gene), data = subset(ggdata, Count > 2)) +
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
    scale_y_continuous(breaks = seq(0, 52, 2)) +
    coord_flip()

  genes <- sort(unique(merge$Gene))
  merge <- do.call(rbind,
    lapply(genes, function(x) {
      tab <- subset(merge, Gene == x)
      data.frame(
        Gene = tab$Gene[1],
        Impact = tab$Impact[1],
        Function = paste(sort(unique(do.call(c, (strsplit(paste(tab$Function, collapse = ';'), ';|&'))))), collapse = ','),
        nHet = sum(tab$nHet),
        HetSamples = paste(tab$HetSamples, collapse = ','))
    }))
  #topgenes <- subset(merge, nHet > 2)$Gene
  onco <- data.frame(row.names = merge$Gene)

  merge$Function <- gsub('synonymous_variant|intron_variant|disruptive_inframe_deletion|disruptive_inframe_insertion|downstream_gene_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|protein_protein_contact|sequence_feature|conservative_inframe_deletion|structural_interaction_variant|upstream_gene_variant', '', merge$Function)
  merge$Function <- gsub('3_prime_UTR_variant|5_prime_UTR_variant', 'UTR', merge$Function)
  merge$Function <- gsub('frameshift_variant', 'frameshift', merge$Function)
  #merge$Function <- gsub('intron_variant', 'intronic', merge$Function)
  merge$Function <- gsub('missense_variant', 'missense', merge$Function)
  merge$Function <- gsub('splice_acceptor_variant|splice_donor_variant|splice_region_variant', 'splicing', merge$Function)
  merge$Function <- gsub('start_lost|stop_gained|stop_lost', 'nonsense', merge$Function)
  #merge$Function <- gsub('synonymous_variant', 'synonymous', merge$Function)
  merge$Function <- sub('^,', '', sub(',$', '', gsub(',+', ',', merge$Function )))

  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  names <- names(dp)
  for (i in 1:length(names)) {
    vec <- c()
    for (j in 1:nrow(merge)) {
      if (grepl(names[i], merge$HetSamples[j])) {
        vec <- c(vec, paste(sort(unique(do.call(c, strsplit(merge[which(merge$Gene == rownames(onco)[j]),'Function'], ',')))), collapse = ';'))
        #vec <- c(vec, 'Somatic')
      } else {
        vec <- c(vec, '')
      }
    }
    onco <- cbind(onco, vec)
    colnames(onco)[i] <- names[i]
  }
  onco <- as.matrix(onco)
  rownames(onco) <- merge$Gene
  #onco <- onco[topgenes,]

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
    onco <- onco[which(rownames(onco) %in% drivers),]

  # annotation
    # Richter's
      richters <- data.frame(readxl::read_xlsx('Metadata/project2final58_richter6_20211201.xlsx'))$patient

    # IGH-V
      ighv

    meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
    meta <- meta[match(colnames(onco), meta$patient),]
    all(colnames(onco) == meta$patient)
    idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
    meta <- meta[idx,]
    onco <- onco[,idx]
    
    # TP53 mutation
      tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)

    ann <- data.frame(
      '6q-deleted' = meta$z_fish6,
      '11q-deleted' = meta$z_fish11q,
      '13q-deleted' = meta$z_fish13q,
      #'17p-deleted' = rep(1, nrow(meta)),
      #'14q translocation' = meta$z_fish14,
      'Trisomy 12' = meta$z_fish12tri,
      #'p53 result (Sanger)' = factor(tp53, levels = c('Abnormal', 'Normal', 'Unknown')),
      'IGHV-mutated' = ifelse(meta$patient %in% ighv, 1, 0),
      'Richter\'s' = ifelse(meta$patient %in% richters, 1, 0))
    #colnames(ann) <- c('6q-deleted','11q-deleted','13q-deleted','17p-deleted','14q translocation',
    #  'Trisomy 12','p53 result (Sanger)','IGHV-mutated','Richter\'s')
    colnames(ann) <- c('6q-deleted','11q-deleted','13q-deleted','Trisomy 12','IGHV-mutated','Richter\'s')
    colAnn <- HeatmapAnnotation(
      df = ann,
      col = list(
        '6q-deleted' = c('1' = 'black', '0' = 'white'),
        '11q-deleted' = c('1' = 'black', '0' = 'white'),
        '13q-deleted' = c('1' = 'black', '0' = 'white'),
        #'17p-deleted' = c('1' = 'black', '0' = 'white'),
        #'14q translocation' = c('1' = 'grey', '0' = 'white'),
        'Trisomy 12' = c('1' = 'blue', '0' = 'white'),
        'IGHV-mutated' = c('1' = 'forestgreen', '0' = 'white'),
        'Richter\'s' = c('1' = 'gold', '0' = 'white')),
        #'p53 result (Sanger)' = c('Abnormal' = 'red', 'Normal' = 'royalblue', 'Unknown' = 'grey')),
      na_col = 'white',
      #annotation_height = c(2, 2, 2, 2, 2, 2, 2, 10),
      gap = unit(1.5, 'mm'))
      #annotation_legend_param = list(
      #  'p53 result (Sanger)' = list(title = 'p53 result (Sanger)')))

  alter_fun <- list(
    background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
     #grid.rect(x, y, w-unit(0.15, 'mm'), h-unit(0.15, 'mm'), gp = gpar(fill = 'white', col = NA))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, 'mm'), h-unit(0.5, 'mm'), gp = gpar(fill = 'red2', col = 'red2'))
    },
    frameshift = function(x, y, w, h) {
      grid.points(x, y, pch = 4)
    },
    nonsense = function(x, y, w, h) {
      grid.rect(x, y, w, h, gp = gpar(lwd = 2))
    },
    splicing = function(x, y, w, h) {
      grid.points(x, y, pch = '*')
      #grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5, gp = gpar(lwd = 2))
    },
    UTR = function(x, y, w, h) {
      grid.points(x, y, pch = 3)
    })
  test_alter_fun(alter_fun)
  cols <- c('missense' = 'black', 'frameshift' = 'black',
    'nonsense' = 'black', 'splicing' = 'black', 'UTR' = 'black')
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
    top_annotation = colAnn,
    heatmap_legend_param = list(
      title = 'Mutation type',
      at = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      labels = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      nrow = 1,
      title_position = 'topcenter'))
  pdf('WGS_output/OncoPrint.pdf', width = 9, height = 7)
    draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
  dev.off()

  # freq of genes associated with del 17p
    (table(onco['TP53',] != '')['TRUE'] / 54) * 100
    (table(onco['SF3B1',] != '')['TRUE'] / 54) * 100
    (table(onco['RPS15',] != '')['TRUE'] / 54) * 100
    (table(onco['NRAS',] != '')['TRUE'] / 54) * 100
    (table(onco['NOTCH1',] != '')['TRUE'] / 54) * 100
    (table(onco['MYD88',] != '')['TRUE'] / 54) * 100
    (table(onco['KRAS',] != '')['TRUE'] / 54) * 100
    (table(onco['IGLL5',] != '')['TRUE'] / 54) * 100
    (table(onco['GPS2',] != '')['TRUE'] / 54) * 100
    (table(onco['DDX3X',] != '')['TRUE'] / 54) * 100
    (table(onco['ATM',] != '')['TRUE'] / 54) * 100






  # without knisbacher variants added in (just 50 CLL drivers)
  merge <- readRDS('WGS_output/merge.Rds')
  merge$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge$HetSamples))
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

  # save IGH-V mutated samples for later
    ighv <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))[,c(1,6)]
    ighv <- ighv[ifelse(is.na(ighv$z_mut), 'unknown', ighv$z_mut) == "0=Mutated",1]

  genes <- sort(unique(merge$Gene))
  merge$Function <- gsub('\\;[A-Za-z0-9_;&]*$|\\&[A-Za-z0-9_;&]*$', '', merge$Function)
  onco <- data.frame(row.names = genes)

  merge$Function <- gsub('synonymous_variant|intron_variant|disruptive_inframe_deletion|disruptive_inframe_insertion|downstream_gene_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|protein_protein_contact|sequence_feature|conservative_inframe_deletion|structural_interaction_variant|upstream_gene_variant|3_prime_UTR_variant|5_prime_UTR_variant', 'other', merge$Function)
  merge$Function <- gsub('frameshift_variant', 'frameshift', merge$Function)
  merge$Function <- gsub('missense_variant', 'missense', merge$Function)
  merge$Function <- gsub('splice_acceptor_variant|splice_donor_variant|splice_region_variant', 'splicing', merge$Function)
  merge$Function <- gsub('start_lost|stop_gained|stop_lost', 'nonsense', merge$Function)
  merge$Function <- gsub('^,', '', gsub(',$', '', gsub(',+', ',', merge$Function)))

  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  names <- names(dp)
  for (i in 1:length(names)) {
    for (j in 1:nrow(onco)) {
      merge.sub <- subset(merge, Gene == rownames(onco)[j])
      if (any(grepl(names[i], merge.sub$HetSamples))) {
        conseq <- sub('^,', '', paste(unique(merge.sub[grep(names[i], merge.sub$HetSamples),'Function']), collapse = ','))
      } else {
        conseq <- ''
      }
      onco[j,i] <- conseq
    }
    colnames(onco)[i] <- names[i]
  }
  onco <- as.matrix(onco)

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
    onco <- onco[which(rownames(onco) %in% drivers),]

    # export data for oncoplot
      merge.wout <- merge[which(merge$Gene %in% drivers),]
      merge.wout <- data.frame(
        chr = unlist(lapply(strsplit(merge.wout$Variant, ':'), function(x) x[1])),
        pos = unlist(lapply(strsplit(merge.wout$Variant, ':|,'), function(x) x[2])),
        ref = unlist(lapply(strsplit(merge.wout$Variant, ':|,|>'), function(x) x[3])),
        alt = unlist(lapply(strsplit(merge.wout$Variant, ':|,|>'), function(x) x[3])),
        merge.wout[,-c(1,7,8,10,11)])
      write.table(merge.wout, 'WGS_output/OncoPrintV2.tsv',
        sep ='\t', row.names = FALSE, quote = FALSE)

  # remove empty rows
    onco <- onco[apply(onco, 1, function(x) !all(x == '')),]

  # annotation
    # Richter's
      richters <- data.frame(readxl::read_xlsx('Metadata/project2final58_richter6_20211201.xlsx'))$patient

    # IGH-V
      ighv

    meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
    meta <- meta[match(colnames(onco), meta$patient),]
    all(colnames(onco) == meta$patient)
    idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
    meta <- meta[idx,]
    onco <- onco[,idx]
    
    # TP53 mutation
      tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)

    ann <- data.frame(
      '6q-deleted' = meta$z_fish6,
      '11q-deleted' = meta$z_fish11q,
      '13q-deleted' = meta$z_fish13q,
      'Trisomy 12' = meta$z_fish12tri,
      'IGHV-mutated' = ifelse(meta$patient %in% ighv, 1, 0),
      'Richter\'s' = ifelse(meta$patient %in% richters, 1, 0))
    colnames(ann) <- c('6q-deleted','11q-deleted','13q-deleted','Trisomy 12','IGHV-mutated','Richter\'s')
    colAnn <- HeatmapAnnotation(
      df = ann,
      col = list(
        '6q-deleted' = c('1' = 'black', '0' = 'white'),
        '11q-deleted' = c('1' = 'black', '0' = 'white'),
        '13q-deleted' = c('1' = 'black', '0' = 'white'),
        'Trisomy 12' = c('1' = 'blue', '0' = 'white'),
        'IGHV-mutated' = c('1' = 'forestgreen', '0' = 'white'),
        'Richter\'s' = c('1' = 'gold', '0' = 'white')),
      na_col = 'white',
      gap = unit(1.5, 'mm'))

  for (i in 1:nrow(onco)) {
    rownames(onco)[i] <- paste0(rownames(onco)[i], ' (',
      round((table(onco[i,]!='')/54) * 100, digits = 2)['TRUE'], '%)')
  }

  alter_fun <- list(
    background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'red2', col = NA))
    },
    frameshift = function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = 'blue2', col = NA))
    },
    nonsense = function(x, y, w, h) {
      grid.points(x, y, pch = 20)
    },
    splicing = function(x, y, w, h) {
      grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(lwd = 2))
    },
    other = function(x, y, w, h) {
      grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(lwd = 2))
    })
  test_alter_fun(alter_fun)
  cols <- c('missense' = 'red2', 'frameshift' = 'blue2',
    'nonsense' = 'black', 'splicing' = 'black', 'other' = 'black')
  ht_opt$message = TRUE
  geneticprint <- oncoPrint(
    onco,
    get_type = function(x) strsplit(x, ',')[[1]],
    name = 'Oncoprint',
    alter_fun = alter_fun,
    col = cols,
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
    top_annotation = colAnn,
    heatmap_legend_param = list(
      title = 'Mutation type',
      at = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      labels = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      nrow = 1,
      title_position = 'topcenter'))
  pdf('WGS_output/OncoPrintV2.pdf', width = 9, height = 7)
    draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
  dev.off()

  # as above, but separate oncoplot by IGHV mutation status
    geneticprint <- oncoPrint(
      onco,
      get_type = function(x) strsplit(x, ',')[[1]],
      column_split = ann[,'IGHV-mutated'],
      name = 'Oncoprint',
      alter_fun = alter_fun,
      col = cols,
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
      top_annotation = colAnn,
      heatmap_legend_param = list(
        title = 'Mutation type',
        at = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        labels = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        nrow = 1,
        title_position = 'topcenter'))
    pdf('WGS_output/OncoPrintV2_IGHV.pdf', width = 9, height = 7)
      draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
    dev.off()

    # check mutation differences between IGHV mutated versus not
      t1 <- table(onco[,which(ann[,'IGHV-mutated'] == 1)])
      t2 <- table(onco[,which(ann[,'IGHV-mutated'] == 0)])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(ann[,'IGHV-mutated'] == 1))
      t2 <- sum(t2[names(t2) != '']) / length(which(ann[,'IGHV-mutated'] == 0))

      t1 <- table(onco[grep('TP53', rownames(onco)),which(ann[,'IGHV-mutated'] == 1)])
      t2 <- table(onco[grep('TP53', rownames(onco)),which(ann[,'IGHV-mutated'] == 0)])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(ann[,'IGHV-mutated'] == 1))
      t2 <- sum(t2[names(t2) != '']) / length(which(ann[,'IGHV-mutated'] == 0))

      tab <- table(ann[,'IGHV-mutated'], onco[grep('NOTCH1', rownames(onco)),])
      tab[1,2] <- sum(tab[1,2:5])
      tab[2,2] <- sum(tab[2,2:5])
      tab <- tab[,1:2]
      colnames(tab) <- c('NOTCH1-unmutated', 'NOTCH1-mutated')
      chisq.test(tab)

  # Oncoplot for bi-allelic vs mono-allelic P53 (P53 mutation profile as the as the top track with all mono-allelics on one side etc)
    geneticprint <- oncoPrint(
      onco,
      get_type = function(x) strsplit(x, ',')[[1]],
      column_split = ifelse(onco[grep('TP53', rownames(onco)),] == '', 'Monoallelic', 'Biallelic'),
      name = 'Oncoprint',
      alter_fun = alter_fun,
      col = cols,
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
      top_annotation = colAnn,
      heatmap_legend_param = list(
        title = 'Mutation type',
        at = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        labels = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        nrow = 1,
        title_position = 'topcenter'))
    pdf('WGS_output/OncoPrintV3_MonoBiAllelic.pdf', width = 9, height = 7)
      draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
    dev.off()

    # check mutation differences between IGHV mutated versus not
      t1 <- table(onco[,which(ifelse(onco[grep('TP53', rownames(onco)),] == '', 'Monoallelic', 'Biallelic') == 'Biallelic')])
      t2 <- table(onco[,which(ifelse(onco[grep('TP53', rownames(onco)),] == '', 'Monoallelic', 'Biallelic') == 'Monoallelic')])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(ifelse(onco[grep('TP53', rownames(onco)),] == '', 'Monoallelic', 'Biallelic') == 'Biallelic'))
      t2 <- sum(t2[names(t2) != '']) / length(which(ifelse(onco[grep('TP53', rownames(onco)),] == '', 'Monoallelic', 'Biallelic') == 'Monoallelic'))

  # Omcoplot to compare CK ≥3  to CK<3 ( this will be limited to 35 patients that have CK data available – have CK as the top track and segregate based on patients with Ck≥ 3 on one side and lower CK on the other)
    ck <- data.frame(readxl::read_xlsx(
      'Metadata/BMS CpG research cases for Preeti.xlsx',
      sheet = 'PROJECT 1 & 2 KB'))
    ck <- ck[,c(1,3,4,5)]                                                   
    colnames(ck) <- c('SampleID', 'Aberrations', 'ISCN', 'Description')
    ck$SampleID <- unlist(lapply(strsplit(
      unlist(lapply(strsplit(ck$SampleID, '-'), function(x) x[1])), ' '), function(x) rev(rev(x)[1])))
    ck$Description <- gsub('\\;\\ \\;\\ ', '; ', gsub('\\r|\\n', '; ', ck$Description))
    ck$Aberrations_Three <- factor(
      ifelse(ck$Aberrations >= 3, '≥3','<3'))

    idx <- which(colnames(onco) %in% ck$SampleID)
    onco.sub <- onco[,idx]
    ann.sub <- ann[idx,]
    ck <- ck[match(colnames(onco.sub), ck$SampleID),]
    ck$SampleID == colnames(onco.sub)
    ann.sub$ComplexKaryotype <- ck$Aberrations_Three
    colnames(ann.sub)[7] <- 'CK aberrations'
    colAnn.sub <- HeatmapAnnotation(
      df = ann.sub,
      col = list(
        'CK aberrations' = c('≥3' = 'red2', '<3' = 'white'),
        '6q-deleted' = c('1' = 'black', '0' = 'white'),
        '11q-deleted' = c('1' = 'black', '0' = 'white'),
        '13q-deleted' = c('1' = 'black', '0' = 'white'),
        'Trisomy 12' = c('1' = 'blue', '0' = 'white'),
        'IGHV-mutated' = c('1' = 'forestgreen', '0' = 'white'),
        'Richter\'s' = c('1' = 'gold', '0' = 'white')),
      na_col = 'white',
      gap = unit(1.5, 'mm'))

    geneticprint <- oncoPrint(
      onco.sub,
      get_type = function(x) strsplit(x, ',')[[1]],
      column_split = ann.sub[,'CK aberrations'],
      name = 'Oncoprint',
      alter_fun = alter_fun,
      col = cols,
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
      top_annotation = colAnn.sub,
      heatmap_legend_param = list(
        title = 'Mutation type',
        at = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        labels = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        nrow = 1,
        title_position = 'topcenter'))
    cairo_pdf('WGS_output/OncoPrintV4_CK.pdf', width = 7, height = 8)
      draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
    dev.off()

    chisq.test(table(ann.sub[,'CK aberrations'], ann.sub[,'IGHV-mutated']))

    # check mutation differences between IGHV mutated versus not
      t1 <- table(onco.sub[,which(ann.sub[,'CK aberrations'] == '≥3')])
      t2 <- table(onco.sub[,which(ann.sub[,'CK aberrations'] == '<3')])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(ann.sub[,'CK aberrations'] == '≥3'))
      t2 <- sum(t2[names(t2) != '']) / length(which(ann.sub[,'CK aberrations'] == '<3'))


  # Oncoplot for patients never treated vs treated at some point ( there are 15 untreated patients and 39 treated)
    geneticprint <- oncoPrint(
      onco,
      get_type = function(x) strsplit(x, ',')[[1]],
      column_split = ifelse(meta$z_treatdesc == '0=Not\ntreated', 'Untreated', ifelse(meta$z_treatdesc == '1=Treated', 'Treated', NA)),
      name = 'Oncoprint',
      alter_fun = alter_fun,
      col = cols,
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
      top_annotation = colAnn,
      heatmap_legend_param = list(
        title = 'Mutation type',
        at = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        labels = c('missense', 'frameshift',
          'nonsense', 'splicing', 'UTR'),
        nrow = 1,
        title_position = 'topcenter'))
    pdf('WGS_output/OncoPrintV6_Treatment.pdf', width = 9, height = 7)
      draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
    dev.off()

    # check mutation differences between IGHV mutated versus not
      t1 <- table(onco[,which(meta$z_treatdesc == '1=Treated')])
      t2 <- table(onco[,which(meta$z_treatdesc == '0=Not\ntreated')])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(meta$z_treatdesc == '1=Treated'))
      t2 <- sum(t2[names(t2) != '']) / length(which(meta$z_treatdesc == '0=Not\ntreated'))

      t1 <- table(onco[grep('TP53', rownames(onco)),which(meta$z_treatdesc == '1=Treated')])
      t2 <- table(onco[grep('TP53', rownames(onco)),which(meta$z_treatdesc == '0=Not\ntreated')])
      chisq.test(rbind(c(sum(t1[names(t1) != '']), sum(t1[names(t1) == ''])), c(sum(t2[names(t2) != '']), sum(t2[names(t2) == '']))))
      t1 <- sum(t1[names(t1) != '']) / length(which(meta$z_treatdesc == '1=Treated'))
      t2 <- sum(t2[names(t2) != '']) / length(which(meta$z_treatdesc == '0=Not\ntreated'))

      tab <- table(meta$z_treatdesc, onco[grep('NOTCH1', rownames(onco)),])
      tab[1,2] <- sum(tab[1,2:5])
      tab[2,2] <- sum(tab[2,2:5])
      tab <- tab[,1:2]
      colnames(tab) <- c('NOTCH1-unmutated', 'NOTCH1-mutated')
      chisq.test(tab)






  # with knisbacher variants added in
  merge <- readRDS('WGS_output/merge.Rds')
  merge$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge$HetSamples))
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

  # save IGH-V mutated samples for later
    ighv <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))[,c(1,6)]
    ighv <- ighv[ifelse(is.na(ighv$z_mut), 'unknown', ighv$z_mut) == "0=Mutated",1]

  genes <- sort(unique(merge$Gene))
  merge$Function <- gsub('\\;[A-Za-z0-9_;&]*$|\\&[A-Za-z0-9_;&]*$', '', merge$Function)
  onco <- data.frame(row.names = genes)

  merge$Function <- gsub('synonymous_variant|intron_variant|disruptive_inframe_deletion|disruptive_inframe_insertion|downstream_gene_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|protein_protein_contact|sequence_feature|conservative_inframe_deletion|structural_interaction_variant|upstream_gene_variant|3_prime_UTR_variant|5_prime_UTR_variant', 'other', merge$Function)
  merge$Function <- gsub('frameshift_variant', 'frameshift', merge$Function)
  merge$Function <- gsub('missense_variant', 'missense', merge$Function)
  merge$Function <- gsub('splice_acceptor_variant|splice_donor_variant|splice_region_variant', 'splicing', merge$Function)
  merge$Function <- gsub('start_lost|stop_gained|stop_lost', 'nonsense', merge$Function)
  merge$Function <- gsub('^,', '', gsub(',$', '', gsub(',+', ',', merge$Function)))

  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  names <- names(dp)
  for (i in 1:length(names)) {
    for (j in 1:nrow(onco)) {
      merge.sub <- subset(merge, Gene == rownames(onco)[j])
      if (any(grepl(names[i], merge.sub$HetSamples))) {
        conseq <- sub('^,', '', paste(unique(merge.sub[grep(names[i], merge.sub$HetSamples),'Function']), collapse = ','))
      } else {
        conseq <- ''
      }
      onco[j,i] <- conseq
    }
    colnames(onco)[i] <- names[i]
  }
  onco <- as.matrix(onco)

  # filter for 50 driver genes from Landau et al. (https://www.nature.com/articles/nature15395), any impact type
    landau <- c('SF3B1','ATM','TP53','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
      'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
      'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
      'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')
    # add in drivers from Knisbacher et al. (2022
      knisbacher <- data.frame(readxl::read_xlsx('Doc/41588_2022_1140_MOESM1_ESM.xlsx', sheet = 'Driver_Genes'))[,6]
      knisbacher <- knisbacher[!is.na(knisbacher)]
    drivers <- sort(unique(c(landau, knisbacher)))
    #drivers <- sort(unique(c(landau)))
    onco <- onco[which(rownames(onco) %in% drivers),]

  # remove empty rows
    onco <- onco[apply(onco, 1, function(x) !all(x == '')),]

  # annotation
    # Richter's
      richters <- data.frame(readxl::read_xlsx('Metadata/project2final58_richter6_20211201.xlsx'))$patient

    # IGH-V
      ighv

    meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
    meta <- meta[match(colnames(onco), meta$patient),]
    all(colnames(onco) == meta$patient)
    idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
    meta <- meta[idx,]
    onco <- onco[,idx]
    
    # TP53 mutation
      tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)

    ann <- data.frame(
      '6q-deleted' = meta$z_fish6,
      '11q-deleted' = meta$z_fish11q,
      '13q-deleted' = meta$z_fish13q,
      'Trisomy 12' = meta$z_fish12tri,
      'IGHV-mutated' = ifelse(meta$patient %in% ighv, 1, 0),
      'Richter\'s' = ifelse(meta$patient %in% richters, 1, 0))
    colnames(ann) <- c('6q-deleted','11q-deleted','13q-deleted','Trisomy 12','IGHV-mutated','Richter\'s')
    colAnn <- HeatmapAnnotation(
      df = ann,
      col = list(
        '6q-deleted' = c('1' = 'black', '0' = 'white'),
        '11q-deleted' = c('1' = 'black', '0' = 'white'),
        '13q-deleted' = c('1' = 'black', '0' = 'white'),
        'Trisomy 12' = c('1' = 'blue', '0' = 'white'),
        'IGHV-mutated' = c('1' = 'forestgreen', '0' = 'white'),
        'Richter\'s' = c('1' = 'gold', '0' = 'white')),
      na_col = 'white',
      gap = unit(1.5, 'mm'))

  for (i in 1:nrow(onco)) {
    rownames(onco)[i] <- paste0(rownames(onco)[i], ' (',
      round((table(onco[i,]!='')/54) * 100, digits = 2)['TRUE'], '%)')
  }

  alter_fun <- list(
    background = function(x, y, w, h) {
     grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'red2', col = NA))
    },
    frameshift = function(x, y, w, h) {
      grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = 'blue2', col = NA))
    },
    nonsense = function(x, y, w, h) {
      grid.points(x, y, pch = 20)
    },
    splicing = function(x, y, w, h) {
      grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(lwd = 2))
    },
    other = function(x, y, w, h) {
      grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(lwd = 2))
    })
  test_alter_fun(alter_fun)
  cols <- c('missense' = 'red2', 'frameshift' = 'blue2',
    'nonsense' = 'black', 'splicing' = 'black', 'other' = 'black')
  ht_opt$message = TRUE
  geneticprint <- oncoPrint(
    onco,
    get_type = function(x) strsplit(x, ',')[[1]],
    name = 'Oncoprint',
    alter_fun = alter_fun,
    col = cols,
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
    top_annotation = colAnn,
    heatmap_legend_param = list(
      title = 'Mutation type',
      at = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      labels = c('missense', 'frameshift',
        'nonsense', 'splicing', 'UTR'),
      nrow = 1,
      title_position = 'topcenter'))
  pdf('WGS_output/OncoPrint_knisbacher.pdf', width = 9, height = 10)
    draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
  dev.off()






  merge <- readRDS('WGS_output/merge.Rds')
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  #merge <- read.table('qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv', sep = '\t', header = TRUE)
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

  meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
  meta <- meta[match(colnames(onco), meta$patient),]
  all(colnames(onco) == meta$patient)

  idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
  meta <- meta[idx,]
  onco <- onco[,idx]
    
  # TP53 mutation
    tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)

  ann <- HeatmapAnnotation(
    'p53 result (clinical)' = factor(tp53, levels = c('Abnormal', 'Normal', 'Unknown')),
    #c17p = anno_points(as.numeric(meta$c17p), gp = gpar(col = 'black'), ylim = c(0, 100), axis = TRUE, pch = '.', size = unit(5.0, 'mm')),
    col = list(
      'p53 result (clinical)' = c('Abnormal' = 'red', 'Normal' = 'royalblue', 'Unknown' = 'grey')),
    na_col = 'white',
    #annotation_height = c(2, 2, 2, 2, 2, 2, 2, 10),
    gap = unit(1.5, 'mm'),
    annotation_legend_param = list(
      'p53 result (clinical)' = list(title = 'p53 result')))

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






  modeling <- data.frame(
    NGS = ifelse(apply(onco, 2, function(x) any(x == 'TP53')) == TRUE, 'TP53', 'WildType'),
    p53 = tp53,
    OS = meta$z_os,
    Death = ifelse(meta$z_vital == '1=Dead', 1, ifelse(meta$z_vital == '0=Alive', 0, NA)),
    Tx1 = ifelse(grepl('Ibrutinib', meta$txtype), 'Ibrutinib', ifelse(is.na(meta$txtype), 'Unknown', 'Non-Ibrutinib')),
    Tx2 = meta$z_saptxcat,
    Tx3 = meta$txcat_main)                 

  pdf('WGS_output/Survival_p53.Unknown.Removed.pdf', width = 15, height = 9)
    arrange_ggsurvplots(list(
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown')),
        data = subset(modeling, p53 != 'Unknown'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Overall',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx1 == 'Ibrutinib')),
        data = subset(modeling, p53 != 'Unknown' & Tx1 == 'Ibrutinib'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Ibrutinib',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx1 == 'Non-Ibrutinib')),
        data = subset(modeling, p53 != 'Unknown' & Tx1 == 'Non-Ibrutinib'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Non-Ibrutinib',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'BTKi-based')),
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'BTKi-based'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'BTKi-based',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'Monoclonal\nantibody-based')),
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'Monoclonal\nantibody-based'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Monoclonal antibody-based',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'Untreated')),
        data = subset(modeling, p53 != 'Unknown' & Tx2 == 'Untreated'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Untreated',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx3 == 'Antibody only')),
        data = subset(modeling, p53 != 'Unknown' & Tx3 == 'Antibody only'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Antibody only',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, p53 != 'Unknown' & Tx3 == 'RTK (Receptor Tyrosine\nKinase) inhibitor or PI3K\ninhibitor')),
        data = subset(modeling, p53 != 'Unknown' & Tx3 == 'RTK (Receptor Tyrosine\nKinase) inhibitor or PI3K\ninhibitor'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'RTK or PI3K inhibitor',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE)),
      ncol = 4,
      nrow = 2)
  dev.off()

  pdf('WGS_output/Survival.pdf', width = 15, height = 9)
    arrange_ggsurvplots(list(
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = modeling),
        data = modeling,
        risk.table = TRUE,
        pval = TRUE,
        title = 'Overall',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx1 == 'Ibrutinib')),
        data = subset(modeling, Tx1 == 'Ibrutinib'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Ibrutinib',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx1 == 'Non-Ibrutinib')),
        data = subset(modeling, Tx1 == 'Non-Ibrutinib'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Non-Ibrutinib',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx2 == 'BTKi-based')),
        data = subset(modeling, Tx2 == 'BTKi-based'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'BTKi-based',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx2 == 'Monoclonal\nantibody-based')),
        data = subset(modeling, Tx2 == 'Monoclonal\nantibody-based'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Monoclonal antibody-based',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx2 == 'Untreated')),
        data = subset(modeling, Tx2 == 'Untreated'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Untreated',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),

      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx3 == 'Antibody only')),
        data = subset(modeling, Tx3 == 'Antibody only'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'Antibody only',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE),
      ggsurvplot(survfit(Surv(OS, Death) ~ NGS,
        data = subset(modeling, Tx3 == 'RTK (Receptor Tyrosine\nKinase) inhibitor or PI3K\ninhibitor')),
        data = subset(modeling, Tx3 == 'RTK (Receptor Tyrosine\nKinase) inhibitor or PI3K\ninhibitor'),
        risk.table = TRUE,
        pval = TRUE,
        title = 'RTK or PI3K inhibitor',
        break.time.by = 2,
        ggtheme = theme_minimal(),
        risk.table.y.text.col = TRUE,
        risk.table.y.text = FALSE)),
      ncol = 4,
      nrow = 2)
  dev.off()






  merge <- readRDS('WGS_output/merge.Rds')
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

  files <- list.files('qc_kb/', pattern = '*.txt', full.names = TRUE)
  files <- files[-grep('TumorOnly', files)]
  dp <- lapply(sort(files),
    function(x) read.table(x, sep = '\t', header = TRUE))
  names(dp) <- sub('kb\\/\\/', '', unlist(lapply(strsplit(files, '-|_'), function(x) x[2])))
  #merge <- read.table('qc_kb/WGS_Merged.Filtered_Variant_Stats.tsv', sep = '\t', header = TRUE)
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

  # esteban
    est <- data.frame(readxl::read_xlsx('Metadata/Esteban/mutations found in targeted sequencing KB.xlsx'))
    est <- subset(est, Gene == 'TP53')
    est <- est[match(colnames(onco), est$SAMPLE),]
    all(colnames(onco) == est$SAMPLE)

  idx <- c(which(meta$p53_result == 'Abnormal'), which(meta$p53_result == 'Normal'), which(is.na(meta$p53_result)))
  meta <- meta[idx,]
  onco <- onco[,idx]
  est <- est[idx,]

  # TP53 mutation
    tp53 <- ifelse(is.na(meta$p53_result), 'Unknown', meta$p53_result)
    names(tp53) <- meta$patient

  # Esteban
    est <- ifelse(is.na(est$WGS), 'Unknown', est$WGS)
    names(est) <- meta$patient

  # Venn diagram
    VennDiagram::draw.triple.venn(
      length(which(tp53 == 'Abnormal')),
      length(which(onco == 'TP53')),
      length(which(est == 'Yes' | est == 'No')),
        n12 = length(which(tp53 == 'Abnormal' & onco == 'TP53')),
        n23 = length(which(onco == 'TP53' & est == 'Yes')),
        n13 = length(which(tp53 == 'Abnormal' & est == 'Yes')),
        n123 = length(which(tp53 == 'Abnormal' & onco == 'TP53' & est == 'Yes')),
      #category = NULL,#c('p53 Sanger\n(abnormal)', 'NGS (TP53\nmutation)', 'Esteban (TP53 mutation)'),
      col = c('red2', 'black', 'purple'))

  ann <- HeatmapAnnotation(
    'p53 result (clinical)\n' = factor(tp53, levels = c('Abnormal', 'Normal', 'Unknown')),
    'Esteban\'s Data\n' = est,
    #c17p = anno_points(as.numeric(meta$c17p), gp = gpar(col = 'black'), ylim = c(0, 100), axis = TRUE, pch = '.', size = unit(5.0, 'mm')),
    col = list(
      'p53 result (clinical)\n' = c('Abnormal' = 'red', 'Normal' = 'royalblue', 'Unknown' = 'grey'),
      'Esteban\'s Data\n' = c('Unknown' = 'grey', 'Yes' = 'purple', 'No' = 'pink')),
    na_col = 'white',
    #annotation_height = c(2, 2, 2, 2, 2, 2, 2, 10),
    gap = unit(1.5, 'mm'),
    annotation_legend_param = list(
      'p53 result (clinical)\n' = list(title = 'p53 result')))

  alter_fun <- list(
    #background = function(x, y, w, h) {
    # grid.rect(x, y, w-unit(0.15, 'mm'), height = unit(0.15, 'mm'), gp = gpar(fill = 'white', col = 'white'))
    # grid.rect(x, y, w-unit(0.15, 'mm'), h-unit(0.15, 'mm'), gp = gpar(fill = 'white', col = NA))
    background = function(x, y, w, h) {
      grid.points(x, y, pch = 0)
    },
    Somatic = function(x, y, w, h) {
      grid.rect(x, y, w-unit(1.5, 'mm'), h-unit(1.5, 'mm'), gp = gpar(fill = 'red2', col = 'red2'))
    },
    TP53 = function(x, y, w, h) {
      grid.points(x, y, pch = 15)
      #grid.rect(x, y, w-unit(2.5, 'mm'), h-unit(2.5, 'mm'), gp = gpar(fill = 'black', col = 'black'))
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
    #column_order = match(names(onco), names(tp53)),
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
  library(trackViewer)
  library(RColorBrewer)
  merge <- readRDS('WGS_output/merge.Rds')
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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

  # pull in variants called by Esteban
    est <- read.table('WGS_output/WGS_Merged.UnFiltered.Esteban_Variant_Stats.tsv', sep = '\t', header = TRUE)
    est <- est[,-which(colnames(est) == 'AF')]

    # Before just binding the datasets, filter out any variants already found in the merged T-N paireds + tumor-only dataset
      idx1 <- which(est$Variant %in% merge$Variant)
      idx2 <- which(merge$Variant %in% est$Variant)
      df.est <- est[idx1,]
      df <- merge[idx2,]
      est <- est[-idx1,]
      merge <- merge[-idx2,]
      df.est <- df.est[match(df$Variant, df.est$Variant),]
      df$Variant == df.est$Variant
      hetsamples <- c()
      nhet <- c()
      for (i in 1:nrow(df)) {
        vec <- sort(unique(c(unlist(strsplit(df$HetSamples[i], ',')),
          unlist(strsplit(df.est$HetSamples[i], ',')))))
        nhet <- c(nhet, length(vec))
        hetsamples <- c(hetsamples, paste(vec, collapse = ','))
      }
      df$nHet <- nhet
      df$HetSamples <- hetsamples
    merge <- rbind(df, merge, est)

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






  # retrieve co-ordinates from tiny 'Download Track Data' button at: https://www.ncbi.nlm.nih.gov/gene/7157
  library(trackViewer)
  library(RColorBrewer)
  merge <- readRDS('WGS_output/merge.Rds')
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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
  merge <- subset(merge, Gene == 'TP53' & Impact == 'HIGH|MODERATE')
  features <- GRanges('chr17',
    IRanges(c(7668421, 7670609, 7673535, 7673701, 7674181, 7674859, 7675053, 7675994, 7676382, 7676521, 7687377),
      width = c(1269, 106, 73, 136, 109, 112, 183, 278, 21, 101, 113),
      names = paste0('exon ', 1:11)),
    fill = brewer.pal(11, "Spectral"),
    height = rep(0.1, 11))
  cbio <- read.csv('WGS_output/cBioportal_TP53.tsv', header = TRUE, sep = '\t')[,c('Protein.Change','Chromosome','Start.Pos','Ref','Var')]
  cbio$variant <- paste0(cbio$Chromosome, ':', cbio$Start.Pos, ',', cbio$Ref, '>', cbio$Var)
  SNP <- as.numeric(sub(',[A-Z,]*>[A-Z,]*$', '', sub('chr17:', '',  merge$Variant)))
  merge$Variant_original <- merge$Variant
  merge$Variant <- cbio[match(merge$Variant, cbio$variant),'Protein.Change']
  merge$Variant <- ifelse(is.na(merge$Variant), 'indel', merge$Variant)
  merge[,c('Variant','Variant_original')]
  # manually change indel annotation
    merge$Variant[merge$Variant_original == 'chr17:7674184,GA>G'] <- 'c.778_779delinsC'
    merge$Variant[merge$Variant_original == 'chr17:7674903,TTC>T'] <- 'c.626_628delinsA'
    merge$Variant[merge$Variant_original == 'chr17:7675166,GAATC>G'] <- 'c.442_446delinsC'
    merge$Variant[merge$Variant_original == 'chr17:7676250,ATTGCT>A'] <- 'c.114_119delinsT'
    merge$Variant[merge$Variant_original == 'chr17:7670684,CG>C'] <- 'c.1024_1025delinsG'
    merge$Variant[merge$Variant_original == 'chr17:7673588,A>AG'] <- 'c.940delinsCT'
    merge$Variant[merge$Variant_original == 'chr17:7674184,G>GT'] <- 'c.779delinsAC'
    merge$Variant[merge$Variant_original == 'chr17:7674212,TG>T'] <- 'c.750_751delinsA'
    merge$Variant[merge$Variant_original == 'chr17:7674261,G>GT'] <- 'c.702delinsAC'
    merge$Variant[merge$Variant_original == 'chr17:7674289,ACCTAGGAG>A'] <- 'c.673_674delinsT'
    merge$Variant[merge$Variant_original == 'chr17:7674910,ATCCAAATAC>A'] <- 'c.612_621delinsT'
    merge$Variant[merge$Variant_original == 'chr17:7675216,CTTG>C'] <- 'c.393_396delinsG'
    merge$Variant[merge$Variant_original == 'chr17;7675997,G>GC'] <- 'c.372delinsGC'
    merge$Variant[merge$Variant_original == 'chr17:7675997,G>GC'] <- 'c.372delinsGC'
    merge$Variant[merge$Variant_original == 'chr17:7673193,CATTTTCAACTTACAAT>C'] <- 'r.spl'
    merge$Variant[merge$Variant_original == 'chr17:7667260,TAA>TAAAA,TAAA,TA,T'] <- 'g.7667260TAA>TAAAA'
    merge$Variant[merge$Variant_original == 'chr17:7667260,T>TA'] <- 'g.7667260T>TA'
    merge$Variant[merge$Variant_original == 'chr17:7667261,A>T'] <- 'g.7667261A>T'
    merge[,c('Variant','Variant_original')]
    subset(merge[,c('Variant','Variant_original')], Variant %in% c('indel',''))
  TP53 <- GRanges('chr17', IRanges(SNP, width = 1, names = merge$Variant),
    color = sample.int(length(SNP), length(SNP), replace = TRUE),
    score = merge$nHet)
  lolliplot(TP53, features)






  # retrieve co-ordinates from tiny 'Download Track Data' button at: https://www.ncbi.nlm.nih.gov/gene/7157
  library(trackViewer)
  library(RColorBrewer)
  merge <- readRDS('WGS_output/merge.Rds')
  merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
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
  merge <- subset(merge, Gene == 'MED12' & Impact == 'HIGH|MODERATE')
  features <- GRanges('chrX',
    IRanges(c(71118596,71119373,71119686,71120014,71120971,71121327,71121562,71122200,71122508,
      71122738,71123095,71123594,71124159,71124764,71124976,71125351,71125663,71126036,71126341,71126969,
      71127336,71127893,71128296,71128598,71129114,71129316,71129680,71130035,71131550,71132073,71132377,
      71132845,71133123,71134357,71134713,71135092,71136281,71136879,71137187,71137558,71137726,71140635,
      71141230,71141883,71142175),
      width = c(257,104,191,156,181,110,254,146,99,136,131,126,229,80,170,144,50,118,143,163,131,227,144,120,101,113,
        175,179,71,133,161,111,89,109,135,161,374,150,196,77,217,222,140,81,275),
      names = paste0('exon ', 1:45)),
    fill = brewer.pal(11, "Spectral"),
    height = rep(0.1, 11))
  cbio <- read.csv('WGS_output/cBioportal_MED12.tsv', header = TRUE, sep = '\t')[,c('Protein.Change','Chromosome','Start.Pos','Ref','Var')]
  cbio$Chromosome <- paste0('chr', cbio$Chromosome)
  cbio$variant <- paste0(cbio$Chromosome, ':', cbio$Start.Pos, ',', cbio$Ref, '>', cbio$Var)
  SNP <- as.numeric(sub(',[A-Z,]*>[A-Z,]*$', '', sub('chrX:', '',  merge$Variant)))
  merge$Variant_original <- merge$Variant
  merge$Variant <- cbio[match(merge$Variant, cbio$variant),'Protein.Change']
  merge$Variant <- ifelse(is.na(merge$Variant), 'indel', merge$Variant)
  merge[,c('Variant','Variant_original')]
  # manually change indel annotation
    merge$Variant[merge$Variant_original == 'chrX:71138754,C>CA'] <- 'c.6044+811C>CA'
    merge$Variant[merge$Variant_original == 'chrX:71139898,C>CA'] <- 'c.6045-737C>CA'
    merge$Variant[merge$Variant_original == 'chrX:71140080,CTTT>C'] <- 'r.spl'
    merge$Variant[merge$Variant_original == 'chrX:71141301,A>ACAGCAACACCAG'] <- 'g.71141301A>ACAGCAACACCAG'
    merge[,c('Variant','Variant_original')]
    subset(merge[,c('Variant','Variant_original')], Variant %in% c('indel',''))
  MED12 <- GRanges('chrX', IRanges(SNP, width = 1, names = merge$Variant),
    color = sample.int(length(SNP), length(SNP), replace = TRUE),
    score = merge$nHet)
  lolliplot(MED12, features)






  # are all key variants that are identifed in T-N paired workflwo also identified in Tumor-only?
    landau <- c('SF3B1','ATM','TP53','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
      'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
      'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
      'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')

    knisbacher <- data.frame(readxl::read_xlsx('Doc/41588_2022_1140_MOESM1_ESM.xlsx', sheet = 'Driver_Genes'))[,6]
    knisbacher <- knisbacher[!is.na(knisbacher)]

    merge <- readRDS('WGS_output/merge.Rds')
    merge$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge$HetSamples))
    merge.tumoronly <- readRDS('WGS_output/merge.TumorOnly.Rds')
    merge.tumoronly$HetSamples <- gsub('_CD19', '', gsub('-T', '', merge.tumoronly$HetSamples))
    merge.subset <- merge[-which(merge$Variant %in% merge.tumoronly$Variant),]
    subset(merge.subset, Gene == 'TP53')
    subset(merge.subset, Gene %in% landau)
    subset(merge.subset, Gene %in% knisbacher)

