
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
    drivers <- sort(unique(c(landau, knisbacher)))
    #drivers <- sort(unique(c(landau)))
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
  pdf('WGS_output/OncoPrint.NewDrivers.pdf', width = 9, height = 9)
    draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
  dev.off()
  png('WGS_output/OncoPrint.NewDrivers.png', units = 'in', res = 300, width = 9, height = 9)
    draw(geneticprint, heatmap_legend_side = 'bottom', annotation_legend_side = 'top', newpage = TRUE)
  dev.off()

