
  # as per 05a.wgs_Battenberg.R, but excluding samples with erroneous TP53 CN

  require(ggplot2)
  require(ggrepel)
  require(RegParallel)
  require(survminer)
  require(survival)

  mytheme <- theme(
    legend.position = 'top',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'plain', vjust = 1),
    axis.line = element_line(size = 1.0, colour = 'black'),
    axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5),
    axis.title = element_text(size = 16, face = 'bold'),
    legend.key = element_blank(),
    legend.key.size = unit(0.75, 'cm'),
    legend.text = element_text(size = 16),
    title = element_text(size = 16),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_text(size = 14, face = 'bold'),
    strip.text.y = element_text(size = 14, face = 'bold', margin = margin(), angle = 90),
    strip.background = element_rect(fill = 'white', colour = 'white'), 
    strip.text = element_text(size = 14, face = 'bold', colour = 'black', margin = margin()),
    strip.switch.pad.grid = unit(0, 'cm'))



  # import data
    cn <- list(
      readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
      readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
      readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal)
    genes <- unique(c(rownames(cn[[1]]), rownames(cn[[2]]), rownames(cn[[3]])))
    for (i in 1:length(cn)) {
      cn[[i]] <- cn[[i]][match(genes, rownames(cn[[i]])),]
    }
    all(rownames(cn[[1]]) == rownames(cn[[2]]))
    all(rownames(cn[[1]]) == rownames(cn[[3]]))
    all(rownames(cn[[2]]) == rownames(cn[[3]]))
    cn <- do.call(cbind, cn)
    colnames(cn) <- sub('_CD19|\\-T', '', colnames(cn))
    wout <- data.frame(SampleID = rownames(cn), cn)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Gene.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    # determine missing sample in CN data compared to WGS
      muts <- readRDS('WGS_output/merge.Rds')
      muts_ids <- unique(do.call(c, strsplit(gsub('\\-T|_CD19', '', muts$HetSamples), ',')))
      muts_ids[-which(muts_ids %in% colnames(cn))]

    ploidy <- rbind(
      readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_purityPloidy_20220525.RDS'),
      readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_purityPloidy_20220525.RDS'),
      readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_purityPloidy_20220525.RDS'))
    ploidy$celgene_id <- sub('_CD19|\\-T', '', ploidy$celgene_id)
    ploidy <- ploidy[match(colnames(cn), ploidy$celgene_id),]
    ploidy$celgene_id == colnames(cn)
    write.table(ploidy, 'CN_output/p53Exclusions.PurityPloidy.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    cyto <- list(
      readRDS('AnalyzedData/20200715/DA882_20200715.cyto.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
      readRDS('AnalyzedData/20210330/DA882_20210330.cyto.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
      readRDS('AnalyzedData/20210616/DA882_20210616.cyto.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal)
    cytobands <- unique(c(rownames(cyto[[1]]), rownames(cyto[[2]]), rownames(cyto[[3]])))
    for (i in 1:length(cyto)) {
      cyto[[i]] <- cyto[[i]][match(cytobands, rownames(cyto[[i]])),]
    }
    all(rownames(cyto[[1]]) == rownames(cyto[[2]]))
    all(rownames(cyto[[1]]) == rownames(cyto[[3]]))
    all(rownames(cyto[[2]]) == rownames(cyto[[3]]))
    cyto <- do.call(cbind, cyto)
    colnames(cyto) <- sub('_CD19|\\-T', '', colnames(cyto))
    cyto <- cyto[,match(colnames(cn), colnames(cyto))]
    colnames(cyto) == colnames(cn)
    wout <- data.frame(SampleID = rownames(cyto), cyto)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Cytobands.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

  # adjust all CN data by ploidy
    cn <- round(t(apply(cn, 1, function(x) x / (ploidy$ploidy/2))), digits = 1)
    cyto <- round(t(apply(cyto, 1, function(x) x / (ploidy$ploidy/2))), digits = 1)
    ploidy$ploidy <- 2
    wout <- data.frame(SampleID = rownames(cn), cn)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Gene.PloidyAdjusted.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
    wout <- data.frame(SampleID = rownames(cyto), cyto)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Cytobands.PloidyAdjusted.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

  # exclude samples with erroneous CN values for TP53 (should all be deleted)
    idx <- which(cn['TP53', ] / (round(ploidy$ploidy, 0)/2) >= 2)
    all(colnames(cn)[idx] == ploidy[idx,]$celgene_id)
    all(colnames(cyto)[idx] == ploidy[idx,]$celgene_id)

    cn <- cn[,-idx]
    cyto <- cyto[,-idx]
    ploidy <- ploidy[-idx,]

    wout <- data.frame(SampleID = rownames(cn), cn)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Gene.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    wout <- data.frame(SampleID = rownames(cyto), cyto)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Cytobands.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

  # compare to FISH and other metadata
    meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
    #meta <- meta[,c('patient','z_fish11q','z_fish12tri','z_fish13q','z_fish6','z_fish14')]

    cn_meta <- data.frame(t(cn), meta[match(colnames(cn), meta$patient),])
    all(rownames(cn_meta) == cn_meta$patient)
    wout <- data.frame(SampleID = rownames(cn_meta), cn_meta)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Gene.PloidyAdjusted.Metadata.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    cyto_meta <- data.frame(t(cyto), meta[match(colnames(cyto), meta$patient),])
    all(rownames(cyto_meta) == cyto_meta$patient)
    wout <- data.frame(SampleID = rownames(cyto_meta), cyto_meta)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/p53Exclusions.CN.Cytobands.PloidyAdjusted.Metadata.tsv', sep = '\t', quote = FALSE, row.names = FALSE)



  # ploidy versus ploidy-adjusted copy number
    ggdata <- data.frame(SampleID = colnames(cn),
      TP53 = unlist(cn['TP53',]), Ploidy = ploidy$ploidy)
    p1a <- ggplot(data = ggdata, aes(x = Ploidy, y = TP53, label = SampleID)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Adjusted ploidy') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
      labs(title = NULL, subtitle = NULL,
        caption = NULL) +
      geom_label_repel(size = 3) +
      coord_cartesian(ylim = c(0, max(round(ggdata$TP53, 0) + 1))) +
      scale_y_continuous(breaks = seq(0, max(round(ggdata$TP53, 0) + 1), 1)) +
      scale_x_continuous(breaks = seq(1, 3, 1))
    p1b <- ggplot(data = ggdata, aes(x = Ploidy, y = TP53, label = SampleID)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Adjusted ploidy') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
      labs(title = NULL, subtitle = 'zoomed 0-5 CN range',
        caption = NULL) +
      #geom_label_repel(size = 3) +
      coord_cartesian(ylim = c(0,5)) +
      scale_x_continuous(breaks = seq(1, 3, 1)) +
      scale_y_continuous(breaks = seq(0, 5, 1))

    ggdata <- data.frame(SampleID = colnames(cn),
      TP53 = unlist(cn['TP53',]), Cellularity = ploidy$cellularity)
    p2a <- ggplot(data = ggdata, aes(x = Cellularity, y = TP53, label = SampleID)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Cellularity') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
      labs(title = NULL, subtitle = NULL,
        caption = NULL) +
      geom_label_repel(size = 3) +
      coord_cartesian(ylim = c(0, max(round(ggdata$TP53, 0) + 1))) +
      scale_y_continuous(
        breaks = seq(0, max(round(ggdata$TP53, 0) + 1), 1))
    p2b <- ggplot(data = ggdata, aes(x = Cellularity, y = TP53, label = SampleID)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Cellularity') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
      labs(title = NULL, subtitle = 'zoomed 0-5 CN range',
        caption = NULL) +
      #geom_label_repel(size = 3) +
      coord_cartesian(ylim = c(0, 5)) +
      scale_y_continuous(
        breaks = seq(0, 5, 1))

    ggdata <- data.frame(SampleID = colnames(cn),
      Ploidy = ploidy$ploidy, Cellularity = ploidy$cellularity)
    p3 <- ggplot(data = ggdata, aes(x = Ploidy, y = Cellularity, label = SampleID)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Adjusted ploidy') + ylab('Cellularity') +
      labs(title = NULL, subtitle = NULL,
        caption = NULL) +
      geom_label_repel(size = 3) +
      coord_cartesian(ylim = c(0,1)) +
      scale_x_continuous(breaks = seq(1, 3, 1)) +
      scale_y_continuous(breaks = seq(0, 1, 0.1))

    pdf('CN_output/p53Exclusions.PurityPloidy.TP53.pdf', width = 8.5, height = 9)
      cowplot::plot_grid(p1a, p2a, p3,
        p1b, p2b, ncol = 3)
    dev.off()



  # ploidy TiN contamination
    tin <- data.frame(readxl::read_xlsx(
      'ChrisHartl_WorkingDir/CLL-Clinical/Preeti_Requested_Del17p tumor contamination in germline fraction_11262021.xlsx',
      skip = 3))[,c(1,3)]
    tin$Sample.ID <- gsub('\\-[0-9]*\\/[0-9]*\\/[0-9A-Z\\-]*', '', tin$Sample.ID)
    tin <- tin[match(ploidy$celgene_id, tin$Sample.ID),]
    all(ploidy$celgene_id == tin$Sample.ID, na.rm = TRUE)
    ggdata <- data.frame(ploidy, tin)
    colnames(ggdata)[6] <- 'TiN'
    p1 <- ggplot(data = ggdata, aes(x = ploidy, y = TiN, label = celgene_id)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Adjusted ploidy') + ylab('Tumor-in-Normal\ncontamination fraction') +
      labs(title = NULL, subtitle = NULL,
        caption = NULL) +
      geom_label_repel(size = 3) +
      scale_x_continuous(breaks = seq(1, 3, 1))
    p2 <- ggplot(data = ggdata, aes(x = cellularity, y = TiN, label = celgene_id)) +
      geom_point(colour = 'black', size = 2) +
      theme_bw(base_size = 24) + mytheme +
      guides(colour = guide_legend(override.aes = list(size = 1.5))) +
      xlab('Cellularity') + ylab('Tumor-in-Normal\ncontamination fraction') +
      labs(title = NULL, subtitle = NULL,
        caption = NULL) +
      geom_label_repel(size = 3)
    pdf('CN_output/p53Exclusions.TiN.pdf', width = 7, height = 5)
      cowplot::plot_grid(p1, p2, ncol = 2)
    dev.off()


  # check against clinical info (FISH): regression
    # FISH prob target info is in Table 1 from 'Standardization of fluorescence in situ hybridization studies on chronic lymphocytic leukemia (CLL) blood and marrow cells by the CLL Research Consortium'
    # z_fish11q (11q22) / ATM
      cyto_meta.11q <- data.frame(cyto_meta[,grepl('11q22', colnames(cyto_meta))], z_fish11q = cyto_meta$z_fish11q)
      cyto_meta.11q$z_fish11q <- factor(ifelse(cyto_meta.11q$z_fish11q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cyto_meta.11q[cyto_meta.11q > 2] <- 'AMP'
      cyto_meta.11q[cyto_meta.11q < 2] <- 'DEL'
      cyto_meta.11q[cyto_meta.11q == 2] <- 'NORMAL'
      for (i in 1:(ncol(cyto_meta.11q) -1)) {
        cyto_meta.11q[,i] <- factor(cyto_meta.11q[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      res.11q <- RegParallel(
        data = cyto_meta.11q, formula = 'z_fish11q ~ [*]',
        FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'),
        FUNtype = 'glm',
        variables = colnames(cyto_meta.11q)[-ncol(cyto_meta.11q)],
        blocksize = 3,
        cores = 2,
        nestedParallel = FALSE,
        conflevel = 95,
        excludeTerms = NULL,
        excludeIntercept = TRUE)
      res.11q <- res.11q[order(res.11q$P),]
      res.11q$Variable <- sub('^X', '', res.11q$Variable)
      res.11q$Term <- sub('^X', '', res.11q$Term)

      cn.11q <- data.frame(ATM = cn[grep('^ATM$', rownames(cn)),], z_fish11q = cyto_meta$z_fish11q)
      cn.11q$z_fish11q <- factor(ifelse(cn.11q$z_fish11q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cn.11q[cn.11q > 2] <- 'AMP'
      cn.11q[cn.11q < 2] <- 'DEL'
      cn.11q[cn.11q == 2] <- 'NORMAL'
      for (i in 1:(ncol(cn.11q) -1)) {
        cn.11q[,i] <- factor(cn.11q[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      summary(glm('z_fish11q ~ ATM', data = cn.11q, family = binomial(link = 'logit'), method = 'glm.fit'))

    # z_fish12tri (12q15)
      cyto_meta.12tri <- data.frame(X12q15 = cyto_meta[,grepl('12q15', colnames(cyto_meta))], z_fish12tri = cyto_meta$z_fish12tri)
      cyto_meta.12tri$z_fish12tri <- factor(ifelse(cyto_meta.12tri$z_fish12tri == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cyto_meta.12tri[cyto_meta.12tri > 2] <- 'AMP'
      cyto_meta.12tri[cyto_meta.12tri < 2] <- 'DEL'
      cyto_meta.12tri[cyto_meta.12tri == 2] <- 'NORMAL'
      for (i in 1:(ncol(cyto_meta.12tri) -1)) {
        cyto_meta.12tri[,i] <- factor(cyto_meta.12tri[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      summary(glm(z_fish12tri ~ X12q15, cyto_meta.12tri, family = binomial(link = 'logit'), method = 'glm.fit'))

      cn.12tri <- data.frame(MDM2 = cn[grep('^MDM2$', rownames(cn)),], z_fish12tri = cyto_meta$z_fish12tri)
      cn.12tri$z_fish12tri <- factor(ifelse(cn.12tri$z_fish12tri == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cn.12tri[cn.12tri > 2] <- 'AMP'
      cn.12tri[cn.12tri < 2] <- 'DEL'
      cn.12tri[cn.12tri == 2] <- 'NORMAL'
      for (i in 1:(ncol(cn.12tri) -1)) {
        cn.12tri[,i] <- factor(cn.12tri[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      summary(glm('z_fish12tri ~ MDM2', data = cn.12tri, family = binomial(link = 'logit'), method = 'glm.fit'))

    # z_fish13q (13q14.3 13q34)
      cyto_meta.13q <- data.frame(cyto_meta[,grepl('13q14.3|13q34', colnames(cyto_meta))], z_fish13q = cyto_meta$z_fish13q)
      cyto_meta.13q$z_fish13q <- factor(ifelse(cyto_meta.13q$z_fish13q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cyto_meta.13q[cyto_meta.13q > 2] <- 'AMP'
      cyto_meta.13q[cyto_meta.13q < 2] <- 'DEL'
      cyto_meta.13q[cyto_meta.13q == 2] <- 'NORMAL'
      for (i in 1:(ncol(cyto_meta.13q) -1)) {
        cyto_meta.13q[,i] <- factor(cyto_meta.13q[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      res.13q <- RegParallel(
        data = cyto_meta.13q, formula = 'z_fish13q ~ [*]',
        FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'),
        FUNtype = 'glm',
        variables = colnames(cyto_meta.13q)[-ncol(cyto_meta.13q)],
        blocksize = 1,
        cores = 2,
        nestedParallel = FALSE,
        conflevel = 95,
        excludeTerms = NULL,
        excludeIntercept = TRUE)
      res.13q <- res.13q[order(res.13q$P),]
      res.13q$Variable <- sub('^X', '', res.13q$Variable)
      res.13q$Term <- sub('^X', '', res.13q$Term)

    # z_fish6 (6q23)
      cyto_meta.6 <- data.frame(cyto_meta[,grepl('6q23', colnames(cyto_meta))], z_fish6 = cyto_meta$z_fish6)
      cyto_meta.6 <- cyto_meta.6[,!grepl('16', colnames(cyto_meta.6))]
      cyto_meta.6$z_fish6 <- factor(ifelse(cyto_meta.6$z_fish6 == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cyto_meta.6[cyto_meta.6 > 2] <- 'AMP'
      cyto_meta.6[cyto_meta.6 < 2] <- 'DEL'
      cyto_meta.6[cyto_meta.6 == 2] <- 'NORMAL'
      for (i in 1:(ncol(cyto_meta.6) -1)) {
        cyto_meta.6[,i] <- factor(cyto_meta.6[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      res.6 <- RegParallel(
        data = cyto_meta.6, formula = 'z_fish6 ~ [*]',
        FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'),
        FUNtype = 'glm',
        variables = colnames(cyto_meta.6)[-ncol(cyto_meta.6)],
        blocksize = 3,
        cores = 2,
        nestedParallel = FALSE,
        conflevel = 95,
        excludeTerms = NULL,
        excludeIntercept = TRUE)
      res.6 <- res.6[order(res.6$P),]
      res.6$Variable <- sub('^X', '', res.6$Variable)
      res.6$Term <- sub('^X', '', res.6$Term)

      cn.6 <- data.frame(MYB = cn[grep('^MYB$', rownames(cn)),], z_fish6 = cyto_meta$z_fish6)
      cn.6$z_fish6 <- factor(ifelse(cn.6$z_fish6 == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cn.6[cn.6 > 2] <- 'AMP'
      cn.6[cn.6 < 2] <- 'DEL'
      cn.6[cn.6 == 2] <- 'NORMAL'
      for (i in 1:(ncol(cn.6) -1)) {
        cn.6[,i] <- factor(cn.6[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      summary(glm('z_fish6 ~ MYB', data = cn.6, family = binomial(link = 'logit'), method = 'glm.fit'))

    # z_fish14 (14q32 and 11q13)
      cyto_meta.14 <- data.frame(cyto_meta[,grepl('14q32|11q13', colnames(cyto_meta))], z_fish14 = cyto_meta$z_fish14)
      cyto_meta.14$z_fish14 <- factor(ifelse(cyto_meta.14$z_fish14 == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      cyto_meta.14[cyto_meta.14 > 2] <- 'AMP'
      cyto_meta.14[cyto_meta.14 < 2] <- 'DEL'
      cyto_meta.14[cyto_meta.14 == 2] <- 'NORMAL'
      for (i in 1:(ncol(cyto_meta.14) -1)) {
        cyto_meta.14[,i] <- factor(cyto_meta.14[,i],
          levels = c('NORMAL','DEL','AMP'))
      }
      res.14 <- RegParallel(
        data = cyto_meta.14, formula = 'z_fish14 ~ [*]',
        FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'),
        FUNtype = 'glm',
        variables = colnames(cyto_meta.14)[-ncol(cyto_meta.14)],
        blocksize = 3,
        cores = 2,
        nestedParallel = FALSE,
        conflevel = 95,
        excludeTerms = NULL,
        excludeIntercept = TRUE)
      res.14 <- res.14[order(res.14$P),]
      res.14$Variable <- sub('^X', '', res.14$Variable)
      res.14$Term <- sub('^X', '', res.14$Term)

    # 17p- P53 (17p13)
      cyto_meta.del17p <- data.frame(cyto_meta[,grepl('17p13', colnames(cyto_meta))])
      cyto_meta.del17p[cyto_meta.del17p > 2] <- 'AMP'
      cyto_meta.del17p[cyto_meta.del17p < 2] <- 'DEL'
      cyto_meta.del17p[cyto_meta.del17p == 2] <- 'NORMAL'
      apply(cyto_meta.del17p, 2, table)
      table(ifelse(cn['TP53',] < 2, 'DEL', ifelse(cn['TP53',] > 2, 'AMP', 'NORMAL')))



  # check against clinical info (FISH): box-and-whiskers
    # FISH prob target info is in Table 1 from 'Standardization of fluorescence in situ hybridization studies on chronic lymphocytic leukemia (CLL) blood and marrow cells by the CLL Research Consortium'
    # z_fish11q (11q22)
      cyto_meta.11q <- data.frame(cyto_meta[,grepl('11q22', colnames(cyto_meta))], z_fish11q = cyto_meta$z_fish11q)
      cyto_meta.11q$z_fish11q <- factor(ifelse(cyto_meta.11q$z_fish11q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cyto_meta.11q)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish11q <- ifelse(ggdata$z_fish11q == 'YES', '11q\nDeletion',
        ifelse(ggdata$z_fish11q == 'NO', 'Cyto\nNormal', NA))
      p1 <- ggplot(data = ggdata, aes(x = z_fish11q, y = value)) +
        geom_boxplot(aes(fill = z_fish11q), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish11q)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish11q (11q22)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))
      table(subset(ggdata, variable == '11q22.3')[,'z_fish11q'])

    # z_fish12tri (12q15)
      cyto_meta.12tri <- data.frame(X12q15 = cyto_meta[,grepl('12q15', colnames(cyto_meta))], z_fish12tri = cyto_meta$z_fish12tri)
      cyto_meta.12tri$z_fish12tri <- factor(ifelse(cyto_meta.12tri$z_fish12tri == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cyto_meta.12tri)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish12tri <- ifelse(ggdata$z_fish12tri == 'YES', 'Trisomy\n12',
        ifelse(ggdata$z_fish12tri == 'NO', 'Cyto\nNormal', NA))
      ggdata$z_fish12tri <- factor(ggdata$z_fish12tri, levels = c('Trisomy\n12','Cyto\nNormal'))
      p2 <- ggplot(data = ggdata, aes(x = z_fish12tri, y = value)) +
        geom_boxplot(aes(fill = z_fish12tri), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish12tri)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish12tri (12q15)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))
      table(subset(ggdata, variable == '12q15')[,'z_fish12tri'])

    # z_fish13q (13q14.3 & 13q34)
      cyto_meta.13q <- data.frame(cyto_meta[,grepl('13q14.3|13q34', colnames(cyto_meta))], z_fish13q = cyto_meta$z_fish13q)
      cyto_meta.13q$z_fish13q <- factor(ifelse(cyto_meta.13q$z_fish13q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cyto_meta.13q)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish13q <- ifelse(ggdata$z_fish13q == 'YES', '13q\nDeletion',
        ifelse(ggdata$z_fish13q == 'NO', 'Cyto\nNormal', NA))
      p3 <- ggplot(data = ggdata, aes(x = z_fish13q, y = value)) +
        geom_boxplot(aes(fill = z_fish13q), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish13q)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish13q (13q14.3 & 13q34)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 4))
      table(subset(ggdata, variable == '13q14.3')[,'z_fish13q'])
      table(subset(ggdata, variable == '13q34')[,'z_fish13q'])

    # z_fish6 (6q23)
      cyto_meta.6 <- data.frame(cyto_meta[,grepl('6q23', colnames(cyto_meta))], z_fish6 = cyto_meta$z_fish6)
      cyto_meta.6 <- cyto_meta.6[,!grepl('16', colnames(cyto_meta.6))]
      cyto_meta.6$z_fish6 <- factor(ifelse(cyto_meta.6$z_fish6 == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cyto_meta.6)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish6 <- ifelse(ggdata$z_fish6 == 'YES', '6q\nDeletion',
        ifelse(ggdata$z_fish6 == 'NO', 'Cyto\nNormal', NA))
      p4 <- ggplot(data = ggdata, aes(x = z_fish6, y = value)) +
        geom_boxplot(aes(fill = z_fish6), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish6)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish6 (6q23)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))

    pdf('CN_output/p53Exclusions.FISH.z_fish11q.pdf', width = 10, height = 6)
      cowplot::plot_grid(p1, ncol = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.FISH.z_fish12tri.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p2, ncol = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.FISH.z_fish13q.pdf', width = 7, height = 6)
      cowplot::plot_grid(p3, ncol = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.FISH.z_fish6.pdf', width = 10, height = 6)
      cowplot::plot_grid(p4, ncol = 1)
    dev.off()



  # check gene copy number against clinical info (FISH): box-and-whiskers
    # z_fish11q (11q22) (ATM)
      cn_meta.11q <- data.frame(ATM = cn_meta[,grepl('^ATM$', colnames(cn_meta))], z_fish11q = cn_meta$z_fish11q)
      cn_meta.11q$z_fish11q <- factor(ifelse(cn_meta.11q$z_fish11q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cn_meta.11q)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish11q <- ifelse(ggdata$z_fish11q == 'YES', '11q (ATM)\nDeletion',
        ifelse(ggdata$z_fish11q == 'NO', 'Cyto\nNormal', NA))
      p1 <- ggplot(data = ggdata, aes(x = z_fish11q, y = value)) +
        geom_boxplot(aes(fill = z_fish11q), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish11q)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish11q (11q22)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))
      table(ggdata$z_fish11q)

    # z_fish12tri (12q15) (MDM2)
      cn_meta.12tri <- data.frame(MDM2 = cn_meta[,grepl('^MDM2$', colnames(cn_meta))], z_fish12tri = cn_meta$z_fish12tri)
      cn_meta.12tri$z_fish12tri <- factor(ifelse(cn_meta.12tri$z_fish12tri == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cn_meta.12tri)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish12tri <- ifelse(ggdata$z_fish12tri == 'YES', 'Trisomy\n12 (MDM2)',
        ifelse(ggdata$z_fish12tri == 'NO', 'Cyto\nNormal', NA))
      ggdata$z_fish12tri <- factor(ggdata$z_fish12tri, levels = c('Trisomy\n12 (MDM2)','Cyto\nNormal'))
      p2 <- ggplot(data = ggdata, aes(x = z_fish12tri, y = value)) +
        geom_boxplot(aes(fill = z_fish12tri), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish12tri)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish12tri (12q15)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))
      table(ggdata$z_fish12tri)

    # z_fish6 (6q23) (MYB)
      cn_meta.6 <- data.frame(MYB = cn_meta[,grepl('^MYB$', colnames(cn_meta))], z_fish6 = cn_meta$z_fish6)
      cn_meta.6 <- cn_meta.6[,!grepl('16', colnames(cn_meta.6))]
      cn_meta.6$z_fish6 <- factor(ifelse(cn_meta.6$z_fish6 == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cn_meta.6)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish6 <- ifelse(ggdata$z_fish6 == 'YES', '6q (MYB)\nDeletion',
        ifelse(ggdata$z_fish6 == 'NO', 'Cyto\nNormal', NA))
      p4 <- ggplot(data = ggdata, aes(x = z_fish6, y = value)) +
        geom_boxplot(aes(fill = z_fish6), outlier.shape = NA) +
        #geom_violin(
        #  stat = 'ydensity',
        #  position = 'dodge',
        #  draw_quantiles = NULL,
        #  trim = FALSE,
        #  scale = 'area',
        #  na.rm = TRUE,
        #  orientation = NA,
        #  aes(fill = z_fish6)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        facet_wrap(. ~ variable, scales = 'free_y', ncol = 4) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('') + ylab('Battenberg ploidy-adjusted\ncopy number') +
        labs(title = 'z_fish6 (6q23)', subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, max(ggdata$value) + 2)) +
        scale_y_continuous(
          breaks = seq(0, max(ggdata$value) + 2, 1))

    pdf('CN_output/p53Exclusions.FISH.z_fish11q.ATM.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p1, ncol = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.FISH.z_fish12tri.MDM2.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p2, ncol = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.FISH.z_fish6.MYB.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p4, ncol = 1)
    dev.off()






  # show higher genomic loss in mono- versus bi-allelic
    tp53 <- data.frame(readxl::read_xlsx('WGS_output/TP53-mutated-samples.xlsx'))
    modeling <- cn
    modeling[modeling > 2] <- 'AMP'
    modeling[modeling < 2] <- 'DEL'
    modeling[modeling == 2] <- 'NORMAL'
    modeling <- do.call(rbind, lapply(apply(modeling, 2, table), function(x) x[c('NORMAL', 'AMP', 'DEL')]))
    modeling_meta <- data.frame(modeling, meta[match(rownames(modeling), meta$patient),])
    all(rownames(modeling_meta) == modeling_meta$patient)

    modeling_meta$OS <- modeling_meta$z_os
    modeling_meta$Death <- ifelse(modeling_meta$z_vital == '1=Dead', 1, ifelse(modeling_meta$z_vital == '0=Alive', 0, NA))
    modeling_meta$Tx1 <- ifelse(grepl('Ibrutinib', modeling_meta$txtype), 'Ibrutinib', ifelse(is.na(modeling_meta$txtype), 'Unknown', 'Non-Ibrutinib'))
    modeling_meta$Tx2 <- modeling_meta$z_saptxcat
    modeling_meta$Tx3 <- modeling_meta$txcat_main                

    all(tp53[match(modeling_meta$patient, tp53$Sample),'Sample'] == modeling_meta$patient)
    tp53 <- tp53[match(modeling_meta$patient, tp53$Sample),]
    tp53$Mutation[tp53$Mutation == 'TP53'] <- 'biallelic'
    tp53$Mutation[is.na(tp53$Mutation)] <- 'monoallelic'
    modeling_meta$tp53 <- factor(tp53$Mutation, levels = c('monoallelic', 'biallelic'))

    quant <- quantile(modeling_meta$AMP, c(0.25,0.5,0.75), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP <= quant[1], 'LOWER',
      ifelse(modeling_meta$AMP >= quant[3], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.25,0.5,0.75), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL <= quant[1], 'LOWER',
      ifelse(modeling_meta$DEL >= quant[3], 'UPPER',
        'MIDDLE'))

    modeling_meta$TP53_AMP_CAT <- paste0(tp53$Mutation, '_', modeling_meta$AMP_CAT)
    modeling_meta$TP53_DEL_CAT <- paste0(tp53$Mutation, '_', modeling_meta$DEL_CAT)

    summary(lm(AMP ~ tp53, data = modeling_meta))
    summary(lm(DEL ~ tp53, data = modeling_meta))
    # deletions compared between mono- and bi-allelics
    # remove outlier sample
      summary(lm(log(DEL) ~ tp53, data = subset(modeling_meta, DEL < 5000),
        family = binomial(link = 'logit')))

      ggplot(data = subset(modeling_meta, DEL < 5000), aes(x = tp53, y = DEL)) +
        #geom_boxplot(aes(fill = z_fish6)) +
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = tp53)) +
        stat_summary(
          geom = 'crossbar',
          width = 0.8,
          fatten = 1.5,
          color = 'black',
          fun.data = function(x){return(c(y = median(x), ymin = median(x), ymax = median(x)))}) +
        geom_jitter(position = position_jitter(width = 0.25),
          size = 1.0, colour = 'black') +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('TP53 status') + ylab('Number of deleted genes') +
        labs(title = NULL, subtitle = NULL,
          caption = paste0('Date: ', Sys.Date())) +
        coord_cartesian(ylim = c(0, 5000)) +
        scale_y_continuous(
          breaks = seq(0, 5000, 250)) +
        annotate(geom = 'text', x = 1.5, y = 4900, size = 6,
          label = 'p=0.0269')

    # pull in CCF
      ccf_major <- unlist(c(
        readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_hg19liftFromhg38_20220525.RDS')$fracMajor['TP53',],
        readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_hg19liftFromhg38_20220525.RDS')$fracMajor['TP53',],
        readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_hg19liftFromhg38_20220525.RDS')$fracMajor['TP53',]))
      ccf_minor <- unlist(c(
        readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_hg19liftFromhg38_20220525.RDS')$fracMinor['TP53',],
        readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_hg19liftFromhg38_20220525.RDS')$fracMinor['TP53',],
        readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_hg19liftFromhg38_20220525.RDS')$fracMinor['TP53',]))






  # survival associations
    modeling <- cn
    modeling[modeling > 2] <- 'AMP'
    modeling[modeling < 2] <- 'DEL'
    modeling[modeling == 2] <- 'NORMAL'
    modeling <- do.call(rbind, lapply(apply(modeling, 2, table), function(x) x[c('NORMAL', 'AMP', 'DEL')]))
    modeling_meta <- data.frame(modeling, meta[match(rownames(modeling), meta$patient),])
    all(rownames(modeling_meta) == modeling_meta$patient)

    modeling_meta$OS <- modeling_meta$z_os
    modeling_meta$Death <- ifelse(modeling_meta$z_vital == '1=Dead', 1, ifelse(modeling_meta$z_vital == '0=Alive', 0, NA))
    modeling_meta$Tx1 <- ifelse(grepl('Ibrutinib', modeling_meta$txtype), 'Ibrutinib', ifelse(is.na(modeling_meta$txtype), 'Unknown', 'Non-Ibrutinib'))
    modeling_meta$Tx2 <- modeling_meta$z_saptxcat
    modeling_meta$Tx3 <- modeling_meta$txcat_main                

    quant <- quantile(modeling_meta$AMP, c(0.25,0.5,0.75), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP <= quant[1], 'LOWER',
      ifelse(modeling_meta$AMP >= quant[3], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.25,0.5,0.75), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL <= quant[1], 'LOWER',
      ifelse(modeling_meta$DEL >= quant[3], 'UPPER',
        'MIDDLE'))

    modeling_meta$Tx4 <- ifelse(modeling_meta$Tx2 == 'BTKi-based', 'BTKi-based',
      ifelse(modeling_meta$Tx2 == 'Untreated', 'Untreated', 'non-BTKi-based'))

    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Quartiles.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ paste0(AMP_CAT, Tx4),
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Quartiles.BTKi.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Quartiles.nonBTKi.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - quartile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()

    med <- median(modeling_meta$AMP, na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < med, 'LOWER', 'UPPER')
    med <- median(modeling_meta$DEL, na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < med, 'LOWER', 'UPPER')
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Median.pdf', width = 12, height = 6)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Median.BTKi.pdf', width = 12, height = 6)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 1)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Median.nonBTKi.pdf', width = 12, height = 6)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - median-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 1)
    dev.off()

    quant <- quantile(modeling_meta$AMP, c(0.333,0.666), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < quant[1], 'LOWER',
      ifelse(modeling_meta$AMP > quant[2], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.333,0.666), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < quant[1], 'LOWER',
      ifelse(modeling_meta$DEL > quant[2], 'UPPER',
        'MIDDLE'))
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Tertiles.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Tertiles.BTKi.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'BTKi-based' & DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Tertiles.nonBTKi.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, Tx4 == 'non-BTKi-based' & DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - tertile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()

    quant <- quantile(modeling_meta$AMP, c(0.1,0.9), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < quant[1], 'LOWER',
      ifelse(modeling_meta$AMP > quant[2], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.1,0.9), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < quant[1], 'LOWER',
      ifelse(modeling_meta$DEL > quant[2], 'UPPER',
        'MIDDLE'))
    pdf('CN_output/p53Exclusions.Survival.CN.Burden.Deciles.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - decile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = modeling_meta),
          data = modeling_meta,
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - decile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE')),
          data = subset(modeling_meta, AMP_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - decile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ DEL_CAT,
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE')),
          data = subset(modeling_meta, DEL_CAT != 'MIDDLE'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Copy Number burden - decile-stratified',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 2,
        nrow = 2)
    dev.off()



  # survival associations at the gene level
    modeling <- cn
    modeling[modeling > 2] <- 'AMP'
    modeling[modeling < 2] <- 'DEL'
    modeling[modeling == 2] <- 'aaNORMAL'
    modeling <- t(modeling)
    modeling_meta <- data.frame(modeling, meta[match(rownames(modeling), meta$patient),])
    all(rownames(modeling_meta) == modeling_meta$patient)
    modeling_meta$OS <- modeling_meta$z_os
    modeling_meta$Death <- ifelse(modeling_meta$z_vital == '1=Dead', 1, ifelse(modeling_meta$z_vital == '0=Alive', 0, NA))
    #vars <- colnames(modeling_meta)[1:19429]
    #keep <- apply(apply(modeling_meta[,vars], 2, function(x) levels(factor(x))), 2, length) > 1
    #vars <- vars[keep]
    #drivers <- c('SF3B1','ATM','TP53','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
    #  'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
    #  'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
    #  'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')
    drivers <- c('SF3B1','ATM','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
      'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
      'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
      'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')
    drivers <- drivers[which(drivers %in% colnames(modeling_meta))]
    res <- RegParallel(
      data = modeling_meta,
      formula = 'Surv(OS, Death) ~ factor([*])',
      FUN = function(formula, data)
        coxph(formula = formula,
          data = modeling_meta,
          ties = 'breslow',
          singular.ok = TRUE),
      FUNtype = 'coxph',
      variables = drivers,
      blocksize = 1,
      cores = 1,
      conflevel = 95,
      excludeTerms = NULL,
      excludeIntercept = TRUE)
    res <- res[order(res$P),]
    res$Term <- sub('\\)', ', ', sub('factor\\(', '', res$Term))
    res$Term <- unlist(lapply(strsplit(sub('\\)', ' ', sub('factor\\(', '', res$Term)), ' '), function(x) x[2]))
    #write.table(res, 'CN_output/p53Exclusions.Survival.CN.Gene.tsv', row.names = FALSE, sep = '\t', quote = FALSE)
    # filter for 50 driver genes from Landau et al. (https://www.nature.com/articles/nature15395), any impact type
    write.table(res[which(res$Variable %in% drivers),],
      'CN_output/p53Exclusions.Survival.CN.Gene.CLL.Drivers.tsv', row.names = FALSE, sep = '\t', quote = FALSE)

    modeling_meta <- data.frame(
      apply(modeling_meta[,which(colnames(modeling_meta) %in% drivers)], 2, function(x) sub('aaNORMAL', 'NORMAL', x)),
      OS = modeling_meta$OS,
      Death = modeling_meta$Death)

    pdf('CN_output/p53Exclusions.Survival.CN.Gene.CLL.Drivers.DEL.pdf', width = 15, height = 11)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ NOTCH1,
          data = subset(modeling_meta, NOTCH1 != 'AMP')),
          data = subset(modeling_meta, NOTCH1 != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(NOTCH1)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ TRAF2,
          data = subset(modeling_meta, TRAF2 != 'AMP')),
          data = subset(modeling_meta, TRAF2 != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(TRAF2)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ RPS15,
          data = subset(modeling_meta, RPS15 != 'AMP')),
          data = subset(modeling_meta, RPS15 != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(RPS15)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ MGA,
          data = subset(modeling_meta, MGA != 'AMP')),
          data = subset(modeling_meta, MGA != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(MGA)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ IGLL5,
          data = subset(modeling_meta, IGLL5 != 'DEL')),
          data = subset(modeling_meta, IGLL5 != 'DEL'),
          palette = c('forestgreen', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(IGLL5)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 3,
        nrow = 2)
    dev.off()



  # survival associations at the gene level
    modeling <- cn
    modeling[modeling > 2] <- 'AMP'
    modeling[modeling < 2] <- 'DEL'
    modeling[modeling == 2] <- 'aaNORMAL'
    modeling <- t(modeling)
    modeling_meta <- data.frame(modeling, meta[match(rownames(modeling), meta$patient),])
    all(rownames(modeling_meta) == modeling_meta$patient)
    modeling_meta$OS <- modeling_meta$z_os
    modeling_meta$Death <- ifelse(modeling_meta$z_vital == '1=Dead', 1, ifelse(modeling_meta$z_vital == '0=Alive', 0, NA))
    vars <- colnames(modeling_meta)[1:19429]
    keep <- unlist(lapply(apply(modeling_meta[,vars], 2, function(x) levels(factor(x))), function(x) length(x) > 1))
    vars <- vars[keep]
    modeling_meta <- data.frame(
      MYC = sub('aaNORMAL', 'NORMAL', modeling_meta[,'MYC']),
      OS = modeling_meta$OS,
      Death = modeling_meta$Death)
    pdf('CN_output/p53Exclusions.Survival.CN.Gene.CLL.Drivers.MYC.pdf', width = 14, height = 6)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ MYC,
          data = modeling_meta),
          data = modeling_meta,
          palette = c('forestgreen', 'red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(MYC)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ MYC,
          data = subset(modeling_meta, MYC != 'DEL')),
          data = subset(modeling_meta, MYC != 'DEL'),
          palette = c('forestgreen', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(MYC)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
        ggsurvplot(survfit(Surv(OS, Death) ~ MYC,
          data = subset(modeling_meta, MYC != 'AMP')),
          data = subset(modeling_meta, MYC != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(MYC)),
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE)),
        ncol = 3,
        nrow = 1)
    dev.off()



  # compare to Complex karyotype data
    ck <- data.frame(readxl::read_xlsx(
      'Metadata/BMS CpG research cases for Preeti.xlsx',
      sheet = 'PROJECT 1 & 2 KB'))
    ck <- ck[,c(1,3,4,5)]                                                   
    colnames(ck) <- c('SampleID', 'Aberrations', 'ISCN', 'Description')
    ck$SampleID <- unlist(lapply(strsplit(
      unlist(lapply(strsplit(ck$SampleID, '-'), function(x) x[1])), ' '), function(x) rev(rev(x)[1])))
    ck$Description <- gsub('\\;\\ \\;\\ ', '; ', gsub('\\r|\\n', '; ', ck$Description))

    ck_meta <- data.frame(ck, meta[match(ck$SampleID, meta$patient),])
    all(ck$SampleID == ck_meta$patient)
    ck_meta$OS <- ck_meta$z_os
    ck_meta$Death <- ifelse(ck_meta$z_vital == '1=Dead', 1, ifelse(ck_meta$z_vital == '0=Alive', 0, NA))

    ck_meta$Aberrations_Zero <- factor(
      ifelse(ck_meta$Aberrations > 0, 'Complex karyotype','Normal karyotype'),
      levels = c('Normal karyotype','Complex karyotype'))
    ck_meta$Aberrations_One <- factor(
      ifelse(ck_meta$Aberrations >= 1, '>=1','<1'))
    ck_meta$Aberrations_Two <- factor(
      ifelse(ck_meta$Aberrations >= 2, '>=2','<2'))
    ck_meta$Aberrations_Three <- factor(
      ifelse(ck_meta$Aberrations >= 3, '>=3','<3'))

    write.table(apply(ck_meta, 2, function(x) gsub('\\\n', ' ', x)),
      'CN_output/p53Exclusions.ComplexKaryotype_Metadata_Merged.tsv',
      sep = '\t', quote = FALSE, row.names = FALSE)

    nick1 <- readRDS('CN_output/DA882.CNAsize20221028.RDS')
    nick2 <- readRDS('CN_output/DA882.CNAsize.noBuffer.20221028.RDS')
    rownames(nick1) <- sub('\\-T|_CD19', '', rownames(nick1))
    rownames(nick2) <- sub('\\-T|_CD19', '', rownames(nick2))
    nick1 <- nick1[match(ck_meta$patient, rownames(nick1)),]
    nick2 <- nick2[match(ck_meta$patient, rownames(nick2)),]
    ck_meta <- cbind(ck_meta, nick1)

    pdf('CN_output/p53Exclusions.ComplexKaryotype.Survival.pdf', width = 10, height = 10)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Zero,
          data = ck_meta),
          data = ck_meta,
          palette = c('royalblue','red'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Complex vs Normal karyotypes',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_One,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n1 aberration',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Two,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n2 aberrations',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Three,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n3 aberrations',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE)),
          ncol = 2,
          nrow = 2)
      dev.off()

    # karyoptype complexity X deletion burden
      modeling <- cn
      modeling[modeling > 2] <- 'AMP'
      modeling[modeling < 2] <- 'DEL'
      modeling[modeling == 2] <- 'NORMAL'
      modeling <- data.frame(do.call(rbind, lapply(apply(modeling, 2, table), function(x) x[c('NORMAL', 'AMP', 'DEL')])))
      modeling$SampleID <- rownames(modeling)

      modeling_meta <- data.frame(ck, modeling[match(ck$SampleID, modeling$SampleID),])
      all(modeling_meta$SampleID == modeling_meta$SampleID.1)
      modeling_meta$AMP_DEL <- modeling_meta$AMP + modeling_meta$DEL

      nick1 <- readRDS('CN_output/DA882.CNAsize20221028.RDS')
      nick2 <- readRDS('CN_output/DA882.CNAsize.noBuffer.20221028.RDS')
      rownames(nick1) <- sub('\\-T|_CD19', '', rownames(nick1))
      rownames(nick2) <- sub('\\-T|_CD19', '', rownames(nick2))
      nick1 <- nick1[match(modeling_meta$SampleID, rownames(nick1)),]
      nick2 <- nick2[match(modeling_meta$SampleID, rownames(nick2)),]
      modeling_meta <- cbind(modeling_meta, nick1)

      p1 <- ggplot(data = subset(modeling_meta, AMP < 5000), aes(x = Aberrations, y = AMP, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg amplifications') +
        labs(title = NULL, subtitle = NULL,
          caption = '\n') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
      p2 <- ggplot(data = subset(modeling_meta, DEL < 5000), aes(x = Aberrations, y = DEL, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg deletions') +
        labs(title = NULL, subtitle = NULL,
          caption = '\n') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
        #scale_y_continuous(breaks = seq(0, 5, 1))
      p3 <- ggplot(data = subset(modeling_meta, AMP_DEL < 5000), aes(x = Aberrations, y = AMP_DEL, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg\namplifications + deletions') +
        labs(title = NULL, subtitle = NULL,
          caption = 'grey shade, 95% CI;\n# Battenberg events limited to 5000 (samples above this excluded)') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
        #scale_y_continuous(breaks = seq(0, 1, 0.1))

      pdf('CN_output/p53Exclusions.ComplexKaryotype.CNBurden.pdf', width = 13, height = 5.5)
        cowplot::plot_grid(p1, p2, p3, ncol = 3)
      dev.off()
      summary(lm(AMP ~ Aberrations, data = subset(modeling_meta, AMP < 5000)))
      summary(lm(DEL ~ Aberrations, data = subset(modeling_meta, DEL < 5000)))
      summary(lm(AMP_DEL ~ Aberrations, data = subset(modeling_meta, AMP_DEL < 5000)))

      # DA882.CNAsize20221028.RDS')
      p1 <- ggplot(data = subset(modeling_meta, TotalCNA < 2000), aes(x = Aberrations, y = Gain, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nGain') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        annotate(geom = 'text', x = 2, y = 590, label = 'p=0.6856', size = 7) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1)) +
        scale_y_continuous(breaks = seq(0, 600, 50))
      p2 <- ggplot(data = subset(modeling_meta, TotalCNA < 2000), aes(x = Aberrations, y = Loss, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nLoss') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        annotate(geom = 'text', x = 2, y = 590, label = 'p=0.00145', size = 7) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1)) +
        scale_y_continuous(breaks = seq(0, 600, 50))
      p3 <- ggplot(data = subset(modeling_meta, TotalCNA < 2000), aes(x = Aberrations, y = TotalCNA, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nTota CNA') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        annotate(geom = 'text', x = 2, y = 590, label = 'p=0.0259', size = 7) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1)) +
        scale_y_continuous(breaks = seq(0, 600, 50))
      cowplot::plot_grid(p1, p2, p3, ncol = 3)
      summary(lm(Gain ~ Aberrations, data = subset(modeling_meta, TotalCNA < 2000)))
      summary(lm(Loss ~ Aberrations, data = subset(modeling_meta, TotalCNA < 2000)))
      summary(lm(TotalCNA ~ Aberrations, data = subset(modeling_meta, TotalCNA < 2000)))

      # DA882.CNAsize.noBuffer.20221028.RDS
      modeling_meta$Gain <- nick2$Gain
      modeling_meta$Loss <- nick2$Loss
      modeling_meta$TotalCNA <- nick2$TotalCNA
      p1 <- ggplot(data = modeling_meta, aes(x = Aberrations, y = Gain, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nGain') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
      p2 <- ggplot(data = modeling_meta, aes(x = Aberrations, y = Loss, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nLoss') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
      p3 <- ggplot(data = modeling_meta, aes(x = Aberrations, y = TotalCNA, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('DA882.CNAsize20221028.RDS\nTota CNA') +
        labs(title = NULL, subtitle = NULL,
          caption = NULL) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations), 1))
      cowplot::plot_grid(p1, p2, p3, ncol = 3)

    mytheme <- theme(
      legend.position = 'top',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 16, face = 'plain', vjust = 1),
      axis.line = element_line(size = 1.0, colour = 'black'),
      axis.text.x = element_text(angle = 45, size = 12, hjust = 1.0, vjust = 1.0),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(angle = 0, size = 12, hjust = 1.0, vjust = 1.0),
      axis.title = element_text(size = 18, face = 'bold'),
      legend.key = element_blank(),
      legend.key.size = unit(0.75, 'cm'),
      legend.text = element_text(size = 16),
      title = element_text(size = 16),
      panel.border = element_blank(),
      panel.background = element_blank(),
      strip.text.x = element_text(size = 14, face = 'bold'),
      strip.text.y = element_text(size = 14, face = 'bold', margin = margin(), angle = 90),
      strip.background = element_rect(fill = 'white', colour = 'white'), 
      strip.text = element_text(size = 14, face = 'bold', colour = 'black', margin = margin()),
      strip.switch.pad.grid = unit(0, 'cm'))   






  # compare to Complex karyotype data (deletions only)
    ck <- data.frame(readxl::read_xlsx(
      'Metadata/BMS CpG research cases for Preeti_del added.xlsx',
      sheet = 'PROJECT 1 & 2'))
    ck <- ck[,c(1,8,7,9)]                                                   
    colnames(ck) <- c('SampleID', 'Aberrations', 'ISCN', 'Description')
    ck <- ck[(ck[,1] != 'PROJECT 2') & (ck[,1] != 'MULTIPLE ATTEMPTS'),]
    ck <- ck[!is.na(ck[,1]),]
    ck[is.na(ck[,2]),2] <- 0
    ck <- ck[-10,]
    ck$SampleID <- unlist(lapply(strsplit(
      unlist(lapply(strsplit(ck$SampleID, '-'), function(x) x[1])), ' '), function(x) rev(rev(x)[1])))
    ck$Description <- gsub('\\;\\ \\;\\ ', '; ', gsub('\\r|\\n', '; ', ck$Description))

    ck <- ck[!is.na(ck$SampleID),]

    ck_meta <- data.frame(ck, meta[match(ck$SampleID, meta$patient),])
    all(ck$SampleID == ck_meta$patient, na.rm = TRUE)
    ck_meta$OS <- ck_meta$z_os
    ck_meta$Death <- ifelse(ck_meta$z_vital == '1=Dead', 1, ifelse(ck_meta$z_vital == '0=Alive', 0, NA))

    ck_meta$Aberrations_Zero <- factor(
      ifelse(ck_meta$Aberrations > 0, 'Complex karyotype','Normal karyotype'),
      levels = c('Normal karyotype','Complex karyotype'))
    ck_meta$Aberrations_One <- factor(
      ifelse(ck_meta$Aberrations >= 1, '>=1','<1'))
    ck_meta$Aberrations_Two <- factor(
      ifelse(ck_meta$Aberrations >= 2, '>=2','<2'))
    ck_meta$Aberrations_Three <- factor(
      ifelse(ck_meta$Aberrations >= 3, '>=3','<3'))

    pdf('CN_output/p53Exclusions.ComplexKaryotype.Survival.DeletionsOnly.pdf', width = 10, height = 10)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Zero,
          data = ck_meta),
          data = ck_meta,
          palette = c('royalblue','red'),
          risk.table = TRUE,
          pval = TRUE,
          title = 'Complex vs Normal karyotypes',
          break.time.by = 2,
          ggtheme = theme_minimal(),
          risk.table.y.text.col = TRUE,
          risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_One,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n1 aberration',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Two,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n2 aberrations',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE),
          ggsurvplot(survfit(Surv(OS, Death) ~ Aberrations_Three,
            data = ck_meta),
            data = ck_meta,
            palette = c('royalblue','red'),
            risk.table = TRUE,
            pval = TRUE,
            title = 'Complex karyotypes\n3 aberrations',
            break.time.by = 2,
            ggtheme = theme_minimal(),
            risk.table.y.text.col = TRUE,
            risk.table.y.text = FALSE)),
          ncol = 2,
          nrow = 2)
      dev.off()

    # karyoptype complexity X deletion burden
      modeling <- cn
      modeling[modeling > 2] <- 'AMP'
      modeling[modeling < 2] <- 'DEL'
      modeling[modeling == 2] <- 'NORMAL'
      modeling <- data.frame(do.call(rbind, lapply(apply(modeling, 2, table), function(x) x[c('NORMAL', 'AMP', 'DEL')])))
      modeling$SampleID <- rownames(modeling)

      modeling_meta <- data.frame(ck, modeling[match(ck$SampleID, modeling$SampleID),])
      all(modeling_meta$SampleID == modeling_meta$SampleID.1, na.rm = TRUE)
      modeling_meta$AMP_DEL <- modeling_meta$AMP + modeling_meta$DEL

      p1 <- ggplot(data = subset(modeling_meta, AMP < 5000), aes(x = Aberrations, y = AMP, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg amplifications') +
        labs(title = NULL, subtitle = NULL,
          caption = '\n') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations, na.rm = TRUE), 1))
      p2 <- ggplot(data = subset(modeling_meta, DEL < 5000), aes(x = Aberrations, y = DEL, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg deletions') +
        labs(title = NULL, subtitle = NULL,
          caption = '\n') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations, na.rm = TRUE), 1))
        #scale_y_continuous(breaks = seq(0, 5, 1))
      p3 <- ggplot(data = subset(modeling_meta, AMP_DEL < 5000), aes(x = Aberrations, y = AMP_DEL, label = SampleID)) +
        geom_point(colour = 'black', size = 2) +
        geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('# complex karyoptying\naberrations') + ylab('Battenberg\namplifications + deletions') +
        labs(title = NULL, subtitle = NULL,
          caption = 'grey shade, 95% CI;\n# Battenberg events limited to 5000 (samples above this excluded)') +
        #geom_label_repel(size = 2.5) +
        #coord_cartesian(ylim = c(0,5)) +
        scale_x_continuous(breaks = seq(0, max(modeling_meta$Aberrations, na.rm = TRUE), 1))
        #scale_y_continuous(breaks = seq(0, 1, 0.1))

      pdf('CN_output/p53Exclusions.ComplexKaryotype.CNBurden.DeletionsOnly.pdf', width = 13, height = 5.5)
        cowplot::plot_grid(p1, p2, p3, ncol = 3)
      dev.off()

      summary(lm(AMP ~ Aberrations, data = subset(modeling_meta, AMP < 5000)))
      summary(lm(DEL ~ Aberrations, data = subset(modeling_meta, DEL < 5000)))
      summary(lm(AMP_DEL ~ Aberrations, data = subset(modeling_meta, AMP_DEL < 5000)))





  # Another late thought.  Can you plot out corrected TP53 copy number vs. the normal TP53 VAF I shared?  I wonder if the 9 samples we are excluding have high contamination by that measure.
    # re-import data but don't remove samples with erroneous TP53 CN
      cn <- list(
        readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
        readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal,
        readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_hg19liftFromhg38_20220525.RDS')$cnvTotal)
      genes <- unique(c(rownames(cn[[1]]), rownames(cn[[2]]), rownames(cn[[3]])))
      for (i in 1:length(cn)) {
        cn[[i]] <- cn[[i]][match(genes, rownames(cn[[i]])),]
      }
      all(rownames(cn[[1]]) == rownames(cn[[2]]))
      all(rownames(cn[[1]]) == rownames(cn[[3]]))
      all(rownames(cn[[2]]) == rownames(cn[[3]]))
      cn <- do.call(cbind, cn)
      colnames(cn) <- sub('_CD19|\\-T', '', colnames(cn))

      ploidy <- rbind(
        readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_purityPloidy_20220525.RDS'),
        readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_purityPloidy_20220525.RDS'),
        readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_purityPloidy_20220525.RDS'))
      ploidy$celgene_id <- sub('_CD19|\\-T', '', ploidy$celgene_id)
      ploidy <- ploidy[match(colnames(cn), ploidy$celgene_id),]
      ploidy$celgene_id == colnames(cn)

      # adjust all CN data by ploidy
        cn <- round(t(apply(cn, 1, function(x) x / (ploidy$ploidy/2))), digits = 1)

      # record samples that were excluded due to erroneous CN
        exclusions <- names(which(cn['TP53', ] / (round(ploidy$ploidy, 0)/2) >= 2))

    # import Nick's VAF data
      vaf <- read.table('WGS_output/TP53_mutations_af_TiN.csv', sep = ',', header = TRUE)
      vaf <- vaf[match(colnames(cn), vaf$Samples),]
      vaf$Samples == colnames(cn)

    # plot TP53 corrected CN vs. normal AF
      ggdata <- data.frame(Sample = colnames(cn), TP53 = cn['TP53',],
        tAF = vaf$tumor_AF, nAF = vaf$normal_AF)
      cowplot::plot_grid(
        ggplot(data = ggdata, aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = '',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,40)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 5)),
        ggplot(data = subset(ggdata, TP53 < 4), aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ggplot(data = subset(ggdata, TP53 < 4 & nAF > 0), aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4 & normal AF > 0.0',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ncol = 3)
      summary(lm(TP53 ~ nAF, data = ggdata))
      summary(lm(TP53 ~ nAF, data = subset(ggdata, TP53 < 4)))
      summary(lm(TP53 ~ nAF, data = subset(ggdata, TP53 < 4 & nAF > 0)))

    # plot TP53 corrected CN vs. tumor AF
      ggdata <- data.frame(Sample = colnames(cn), TP53 = cn['TP53',],
        tAF = vaf$tumor_AF, nAF = vaf$normal_AF)
      cowplot::plot_grid(
        ggplot(data = ggdata, aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = '',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,40)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 5)),
        ggplot(data = subset(ggdata, TP53 < 4), aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ggplot(data = subset(ggdata, TP53 < 4 & tAF > 0), aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4 & tumor AF > 0.0',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ncol = 3)
      summary(lm(TP53 ~ tAF, data = ggdata))
      summary(lm(TP53 ~ tAF, data = subset(ggdata, TP53 < 4)))
      summary(lm(TP53 ~ tAF, data = subset(ggdata, TP53 < 4 & tAF > 0)))

    # plot TP53 corrected CN vs. normal AF (exclusions only)
      ggdata <- data.frame(Sample = colnames(cn), TP53 = cn['TP53',],
        tAF = vaf$tumor_AF, nAF = vaf$normal_AF)
      ggdata <- ggdata[which(ggdata$Sample %in% exclusions),]
      cowplot::plot_grid(
        ggplot(data = ggdata, aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = '',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,40)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 5)),
        ggplot(data = subset(ggdata, TP53 < 4), aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ggplot(data = subset(ggdata, TP53 < 4 & nAF > 0), aes(x = nAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Normal AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4 & normal AF > 0.0',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ncol = 3)
      summary(lm(TP53 ~ nAF, data = ggdata))
      summary(lm(TP53 ~ nAF, data = subset(ggdata, TP53 < 4)))
      summary(lm(TP53 ~ nAF, data = subset(ggdata, TP53 < 4 & nAF > 0)))

    # plot TP53 corrected CN vs. tumor AF (exclusions only)
      ggdata <- data.frame(Sample = colnames(cn), TP53 = cn['TP53',],
        tAF = vaf$tumor_AF, nAF = vaf$normal_AF)
      ggdata <- ggdata[which(ggdata$Sample %in% exclusions),]
      cowplot::plot_grid(
        ggplot(data = ggdata, aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = '',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,40)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 5)),
        ggplot(data = subset(ggdata, TP53 < 4), aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ggplot(data = subset(ggdata, TP53 < 4 & tAF > 0), aes(x = tAF, y = TP53, label = Sample)) +
          geom_point(colour = 'black', size = 2) +
          geom_smooth(method = 'lm', formula= y~x, level = 0.95) +
          theme_bw(base_size = 24) + mytheme +
          guides(colour = guide_legend(override.aes = list(size = 1.5))) +
          xlab('Tumor AF') + ylab(expression(paste(italic('TP53'), ' ploidy-adjusted copy number'))) +
          labs(title = NULL, subtitle = 'TP53 CN < 4 & tumor AF > 0.0',
            caption = '\n') +
          geom_label_repel(size = 2.5) +
          coord_cartesian(xlim = c(0, 0.8), ylim = c(0,10)) +
          scale_x_continuous(breaks = seq(0, 0.8, 0.1)) +
          scale_y_continuous(breaks = seq(0, max(ggdata$TP53, na.rm = TRUE), 2)),
        ncol = 3)
      summary(lm(TP53 ~ tAF, data = ggdata))
      summary(lm(TP53 ~ tAF, data = subset(ggdata, TP53 < 4)))
      summary(lm(TP53 ~ tAF, data = subset(ggdata, TP53 < 4 & tAF > 0)))

