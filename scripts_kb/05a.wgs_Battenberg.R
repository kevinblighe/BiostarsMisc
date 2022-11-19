
  require(ggplot2)
  require(ggrepel)
  require(RegParallel)
  require(survminer)
  require(survival)

  dir.create('CN_output/', showWarnings = FALSE)

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
      'CN_output/CN.Gene.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    ploidy <- rbind(
      readRDS('AnalyzedData/20200715/DA882_20200715.battenberg_purityPloidy_20220525.RDS'),
      readRDS('AnalyzedData/20210330/DA882_20210330.battenberg_purityPloidy_20220525.RDS'),
      readRDS('AnalyzedData/20210616/DA882_20210616.battenberg_purityPloidy_20220525.RDS'))
    ploidy$celgene_id <- sub('_CD19|\\-T', '', ploidy$celgene_id)
    ploidy <- ploidy[match(colnames(cn), ploidy$celgene_id),]
    ploidy$celgene_id == colnames(cn)
    write.table(ploidy, 'CN_output/PurityPloidy.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

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
      'CN_output/CN.Cytobands.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

  # adjust all CN data by ploidy
    cn <- round(t(apply(cn, 1, function(x) x / (ploidy$ploidy/2))), digits = 1)
    cyto <- round(t(apply(cyto, 1, function(x) x / (ploidy$ploidy/2))), digits = 1)
    ploidy$ploidy <- 2
    wout <- data.frame(SampleID = rownames(cn), cn)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/CN.Gene.PloidyAdjusted.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
    write.table(wout,
      'CN_output/CN.Cytobands.PloidyAdjusted.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

  # compare to FISH and other metadata
    meta <- data.frame(readxl::read_xlsx('ChrisHartl_WorkingDir/CLL-Clinical/project2final58_20211109.xlsx'))
    #meta <- meta[,c('patient','z_fish11q','z_fish12tri','z_fish13q','z_fish6','z_fish14')]

    cn_meta <- data.frame(t(cn), meta[match(colnames(cn), meta$patient),])
    all(rownames(cn_meta) == cn_meta$patient)
    wout <- data.frame(SampleID = rownames(cn_meta), cn_meta)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/CN.Gene.PloidyAdjusted.Metadata.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

    cyto_meta <- data.frame(t(cyto), meta[match(colnames(cyto), meta$patient),])
    all(rownames(cyto_meta) == cyto_meta$patient)
    wout <- data.frame(SampleID = rownames(cyto_meta), cyto_meta)
    colnames(wout) <- sub('^X', '', colnames(wout))
    write.table(wout,
      'CN_output/CN.Cytobands.PloidyAdjusted.Metadata.tsv', sep = '\t', quote = FALSE, row.names = FALSE)



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
      coord_cartesian(ylim = c(0, max(round(ggdata$TP53, 0) + 5))) +
      scale_y_continuous(breaks = seq(0, max(round(ggdata$TP53, 0) + 5), 10)) +
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
      coord_cartesian(ylim = c(0, max(round(ggdata$TP53, 0) + 5))) +
      scale_y_continuous(
        breaks = seq(0, max(round(ggdata$TP53, 0) + 5), 10))
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

    pdf('CN_output/PurityPloidy.TP53.pdf', width = 8.5, height = 9)
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
    pdf('CN_output/TiN.pdf', width = 7, height = 5)
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish11q)) +
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish12tri)) +
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

    # z_fish13q (13q14.3 & 13q34)
      cyto_meta.13q <- data.frame(cyto_meta[,grepl('13q14.3|13q34', colnames(cyto_meta))], z_fish13q = cyto_meta$z_fish13q)
      cyto_meta.13q$z_fish13q <- factor(ifelse(cyto_meta.13q$z_fish13q == 1, 'YES', 'NO'),
        levels = c('NO','YES'))
      ggdata <- reshape2::melt(cyto_meta.13q)
      ggdata$variable <- sub('^X', '', ggdata$variable)
      ggdata$z_fish13q <- ifelse(ggdata$z_fish13q == 'YES', '13q\nDeletion',
        ifelse(ggdata$z_fish13q == 'NO', 'Cyto\nNormal', NA))
      p3 <- ggplot(data = ggdata, aes(x = z_fish13q, y = value)) +
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish13q)) +
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish6)) +
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

    pdf('CN_output/FISH.z_fish11q.pdf', width = 10, height = 6)
      cowplot::plot_grid(p1, ncol = 1)
    dev.off()
    pdf('CN_output/FISH.z_fish12tri.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p2, ncol = 1)
    dev.off()
    pdf('CN_output/FISH.z_fish13q.pdf', width = 7, height = 6)
      cowplot::plot_grid(p3, ncol = 1)
    dev.off()
    pdf('CN_output/FISH.z_fish6.pdf', width = 10, height = 6)
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish11q)) +
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish12tri)) +
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
        geom_violin(
          stat = 'ydensity',
          position = 'dodge',
          draw_quantiles = NULL,
          trim = FALSE,
          scale = 'area',
          na.rm = TRUE,
          orientation = NA,
          aes(fill = z_fish6)) +
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

    pdf('CN_output/FISH.z_fish11q.ATM.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p1, ncol = 1)
    dev.off()
    pdf('CN_output/FISH.z_fish12tri.MDM2.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p2, ncol = 1)
    dev.off()
    pdf('CN_output/FISH.z_fish6.MYB.pdf', width = 4.5, height = 6)
      cowplot::plot_grid(p4, ncol = 1)
    dev.off()



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
    pdf('CN_output/Survival.CN.Burden.Quartiles.pdf', width = 12, height = 12)
      arrange_ggsurvplots(list(
        ggsurvplot(survfit(Surv(OS, Death) ~ AMP_CAT,
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

    med <- median(modeling_meta$AMP, na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < med, 'LOWER', 'UPPER')
    med <- median(modeling_meta$DEL, na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < med, 'LOWER', 'UPPER')
    pdf('CN_output/Survival.CN.Burden.Median.pdf', width = 12, height = 6)
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

    quant <- quantile(modeling_meta$AMP, c(0.333,0.666), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < quant[1], 'LOWER',
      ifelse(modeling_meta$AMP > quant[2], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.333,0.666), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < quant[1], 'LOWER',
      ifelse(modeling_meta$DEL > quant[2], 'UPPER',
        'MIDDLE'))
    pdf('CN_output/Survival.CN.Burden.Tertiles.pdf', width = 12, height = 12)
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

    quant <- quantile(modeling_meta$AMP, c(0.1,0.9), na.rm = TRUE)
    modeling_meta$AMP_CAT <- ifelse(modeling_meta$AMP < quant[1], 'LOWER',
      ifelse(modeling_meta$AMP > quant[2], 'UPPER',
        'MIDDLE'))
    quant <- quantile(modeling_meta$DEL, c(0.1,0.9), na.rm = TRUE)
    modeling_meta$DEL_CAT <- ifelse(modeling_meta$DEL < quant[1], 'LOWER',
      ifelse(modeling_meta$DEL > quant[2], 'UPPER',
        'MIDDLE'))
    pdf('CN_output/Survival.CN.Burden.Deciles.pdf', width = 12, height = 12)
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
    vars <- colnames(modeling_meta)[1:19429]
    keep <- apply(apply(modeling_meta[,vars], 2, function(x) levels(factor(x))), 2, length) > 1
    vars <- vars[keep]
    res <- RegParallel(
      data = modeling_meta,
      formula = 'Surv(OS, Death) ~ factor([*])',
      FUN = function(formula, data)
        coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
      FUNtype = 'coxph',
      variables = vars,
      blocksize = 1000,
      cores = 3,
      conflevel = 95,
      excludeTerms = NULL,
      excludeIntercept = TRUE)
    res <- res[order(res$P),]
    res$Term <- sub('\\)', ', ', sub('factor\\(', '', res$Term))
    res$Term <- unlist(lapply(strsplit(sub('\\)', ' ', sub('factor\\(', '', res$Term)), ' '), function(x) x[2]))
    write.table(res, 'CN_output/Survival.CN.Gene.tsv', row.names = FALSE, sep = '\t', quote = FALSE)
    # filter for 50 driver genes from Landau et al. (https://www.nature.com/articles/nature15395), any impact type
      drivers <- c('SF3B1','ATM','TP53','POT1','NOTCH1','XPO1','BIRC3','RPS15','BRAF','EGR2','MYD88','KRAS',
        'MAP2K1','DDX3X','NRAS','CHD2','SAMHD1','FUBP1','FBXW7','DYRK1A','BCOR','HIST1H1E','NXF1','IRF4','PTPN11',
        'MGA','EWSR1','ZMYM3','FAM50A','IKZF3','MED12','TRAF3','IGLL5','BAZ2A','GNB1','ELF4','TRAF2','CARD11','BRCC3',
        'CHEK2','HIST1H1B','XPO4','ASXL1','PIM1')
    write.table(res[which(res$Variable %in% drivers),],
      'CN_output/Survival.CN.Gene.CLL.Drivers.tsv', row.names = FALSE, sep = '\t', quote = FALSE)

    modeling_meta <- data.frame(
      apply(modeling_meta[,which(colnames(modeling_meta) %in% drivers)], 2, function(x) sub('aaNORMAL', 'NORMAL', x)),
      OS = modeling_meta$OS,
      Death = modeling_meta$Death)

    pdf('CN_output/Survival.CN.Gene.CLL.Drivers.DEL.pdf', width = 15, height = 11)
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
        ggsurvplot(survfit(Surv(OS, Death) ~ TRAF3,
          data = subset(modeling_meta, TRAF3 != 'AMP')),
          data = subset(modeling_meta, TRAF3 != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(TRAF3)),
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
        ggsurvplot(survfit(Surv(OS, Death) ~ MYD88,
          data = subset(modeling_meta, MYD88 != 'AMP')),
          data = subset(modeling_meta, MYD88 != 'AMP'),
          palette = c('red', 'royalblue'),
          risk.table = TRUE,
          pval = TRUE,
          title = bquote(Copy~Number~-~~italic(MYD88)),
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
    pdf('CN_output/Survival.CN.Gene.CLL.Drivers.MYC.pdf', width = 14, height = 6)
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

    pdf('CN_output/ComplexKaryotype.Survival.pdf', width = 10, height = 10)
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

      pdf('CN_output/ComplexKaryotype.CNBurden.pdf', width = 13, height = 5.5)
        cowplot::plot_grid(p1, p2, p3, ncol = 3)
      dev.off()

      summary(lm(AMP ~ Aberrations, data = subset(modeling_meta, AMP < 5000)))
      summary(lm(DEL ~ Aberrations, data = subset(modeling_meta, DEL < 5000)))
      summary(lm(AMP_DEL ~ Aberrations, data = subset(modeling_meta, AMP_DEL < 5000)))

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

    # tabulate corroborative results to Battenberg data
      res <- data.frame()

      # ck[1,c('Aberrations', 'Description')]
        idx <- 1

        # psu dic(17;3)(p13;p13) = loss of 17p13pter and loss of 3p13pter
        #target <- '^17p13'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'psu dic(17;3)(p13;p13) = loss of 17p13pter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^3p13'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'psu dic(17;3)(p13;p13) = loss of 3p13pter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # der(13;18)(q10;q10)del(13)(q12q22) = loss of 13q12q22 and loss of 18p
        #target <- '^13q12'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'der(13;18)(q10;q10)del(13)(q12q22) = loss of 13q12',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^13q22'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'der(13;18)(q10;q10)del(13)(q12q22) = loss of 13q22',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^18p'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'der(13;18)(q10;q10)del(13)(q12q22) = loss of 18p',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(7)(q22)  = loss of 7q22qter
        #target <- '^7q22'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(7)(q22)  = loss of 7q22qter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[2,c('Aberrations', 'Description')]
        idx <- 2

        # add(10)(q22) = loss of 10q22qter
        target <- '^10q22'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(10)(q22) = loss of 10q22qter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # +12 = trisomy 12
        target <- '^12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+12 = trisomy 12',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(13)(q32) = loss of 13q32qter
        target <- '^13q32'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(13)(q32) = loss of 13q32qter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -14,der(17)t(14;17)(q13;p11.2) = loss of 14pterq13 and loss of 17p11.2pter
        target <- '^14q13'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-14,der(17)t(14;17)(q13;p11.2) = loss of 14pterq13',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-14,der(17)t(14;17)(q13;p11.2) = loss of 17p11.2pter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # +20 = trisomy 20
        target <- '^20'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+20 = trisomy 20',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # +21 = trisomy 21
        target <- '^21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+21 = trisomy 21',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[3,c('Aberrations', 'Description')]
        idx <- 3

        # add(10)(q22) = loss of 10q22qter
        #target <- '^10q22'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(10)(q22) = loss of 10q22qter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(14)(q22q32) = loss of 14q22q32
        #target <- '^14q22'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'del(14)(q22q32) = loss of 14q22',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^14q32'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'del(14)(q22q32) = loss of 14q32',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -17 = loss of chromosome 17
        #target <- '^17'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], '-17 = loss of chromosome 17',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[4,c('Aberrations', 'Description')]
        idx <- 4

        # der(4;17)(q10;q10) = loss of chromosome 4 short arm and loss of chromosome 17 short arm
        target <- '^4q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(4;17)(q10;q10) = loss of chromosome 4 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(4;17)(q10;q10) = loss of chromosome 17 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[5,c('Aberrations', 'Description')]
        idx <- 5

        # der(3;4)(q10;q10) = loss of chromosome 3 short arm and loss of chromosome 4 short arm
        target <- '^3q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(3;4)(q10;q10) = loss of chromosome 3 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^4q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(3;4)(q10;q10) = loss of chromosome 4 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # +del(11)(q21q23) = gain of 11pterq21 and gain of 11q23qter
        target <- '^11q21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+del(11)(q21q23) = gain of 11pterq21',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^11q23'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+del(11)(q21q23) = gain of 11q23qter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(17)(p11.2) = loss of 17p11.2pter
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of 17p11.2pter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # there is loss of 17p on both chromosome 17's thus homozygous loss of 17p
        target <- '^17p'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'there is loss of 17p on both chromosome 17\'s thus homozygous loss of 17p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[6,c('Aberrations', 'Description')]
        idx <- 6

        # -4,der(17)t(4;17)(q21;p11.2) = loss of 4pter4q21 and loss of 17p11.2-pter
        target <- '^4q21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-4,der(17)t(4;17)(q21;p11.2) = loss of 4pter4q21',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-4,der(17)t(4;17)(q21;p11.2) = loss of 17p11.2-pter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # dic(11;20)(q23;p13) = loss of 11q23qter and loss of 20p13pter
        target <- '^11q23'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'dic(11;20)(q23;p13) = loss of 11q23qter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^20p13'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'dic(11;20)(q23;p13) = loss of 20p13pter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(11)(q21) or del(11)(q21) = loss of 11q21qter
        target <- '^11q21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(11)(q21) or del(11)(q21) = loss of 11q21qter',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[7,c('Aberrations', 'Description')]
        idx <- 7

      # ck[8,c('Aberrations', 'Description')]
        idx <- 8

        # +12 = trisomy 12
        target <- '^12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+12 = trisomy 12',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # i(17)(q10) = loss of 17p (TP53) and gain of 17q.  = isochromosome 17q
        target <- '^17p'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'i(17)(q10) = loss of 17p (TP53) and gain of 17q.  = isochromosome 17q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17q'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'i(17)(q10) = loss of 17p (TP53) and gain of 17q.  = isochromosome 17q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))


      # ck[9,c('Aberrations', 'Description')]
        idx <- 9

        # add(10)(q22) = loss of 10q22qter
        #target <- '^10q22'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(10)(q22) = loss of 10q22qter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # der(17;22)(q10;q10) = loss of 17p and loss of 22p
        #target <- '^17p'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'der(17;22)(q10;q10) = loss of 17p',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^22p'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'der(17;22)(q10;q10) = loss of 22p',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # i(8)(q10) = gain of 8q and loss of  8p
        #target <- '^8q'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'i(8)(q10) = gain of 8q and loss of  8p',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^8p'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'i(8)(q10) = gain of 8q and loss of  8p',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -9, der(3)t(3;9)(p25;q13) = loss of 3p25pter and loss of 9pter9q13
        #target <- '^3p25'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], '-9, der(3)t(3;9)(p25;q13) = loss of 3p25pter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        #target <- '^9q13'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], '-9, der(3)t(3;9)(p25;q13) = loss of 9pter9q13',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -6 = loss of chromosome 6
        #target <- '^6'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], '-6 = loss of chromosome 6',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(8)(p11.2) = loss of 8p11.2pter
        #target <- '^8p11.2'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(8)(p11.2) = loss of 8p11.2pter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(11)(q23) = loss of 11q23qter
        #target <- '^11q23'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(11)(q23) = loss of 11q23qter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(15)(q24) = loss of 15q24qter
        #target <- '^15q24'
        #ggdata <- data.frame(
        #  Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
        #  CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        #res <- rbind(res,
        #  c(ck$SampleID[idx], 'add(15)(q24) = loss of 15q24qter',
        #    mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[10,c('Aberrations', 'Description')]
        idx <- 10

        # NORMAL
        target <- NULL
        ggdata <- data.frame(
          Cytoband = rownames(cyto),
          CN = cyto)
        res <- rbind(res,
          c(ck$SampleID[idx], 'NORMAL',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[11,c('Aberrations', 'Description')]
        idx <- 11

        # add(15)(q24) = additional unidentified material attached to distal 15q
        target <- '^15q24'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(15)(q24) = additional unidentified material attached to distal 15q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -15 = monosomy 15, loss of whle chromosome 15
        target <- '^15'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-15 = monosomy 15, loss of whle chromosome 15',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[12,c('Aberrations', 'Description')]
        idx <- 12

        # NORMAL
        target <- NULL
        ggdata <- data.frame(
          Cytoband = rownames(cyto),
          CN = cyto)
        res <- rbind(res,
          c(ck$SampleID[idx], 'NORMAL',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[13,c('Aberrations', 'Description')]
        idx <- 13

        # add(18)(p11.2) = loss of the chromosome 18 short arm
        target <- '^18p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(18)(p11.2) = loss of the chromosome 18 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[14,c('Aberrations', 'Description')]
        idx <- 14

        # -8 = loss of chromosome 8
        target <- '^8'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-8 = loss of chromosome 8',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(17)(p11.2) = loss of the chromosome 17 short arm (loss of TP53)
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of the chromosome 17 short arm (loss of TP53)',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -18 = loss of chromosome 18
        target <- '^18'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-18 = loss of chromosome 18',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -20 = loss of chromosome 20
        target <- '^20'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-20 = loss of chromosome 20',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(21)(p11.2) = addition of unknown material on short arm of chromosome 21, possibly part of 20q
        target <- '^21p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(21)(p11.2) = addition of unknown material on short arm of chromosome 21, possibly part of 20q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))
        
      # ck[15,c('Aberrations', 'Description')]
        idx <- 15

        # add(3)q25) = loss of distal chromosome 3 long arm and gain of unidentified material
        target <- '^3q25'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(3)q25) = loss of distal chromosome 3 long arm and gain of unidentified material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -9 = loss of chromosome 9 would include loss of CDKN2A/B
        target <- '^9'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-9 = loss of chromosome 9 would include loss of CDKN2A/B',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(13)(q12q14) = deletion in chromosome 13 long arm including miR region and likely RB1
        target <- '^13q12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q12q14) = deletion in chromosome 13 long arm including miR region and likely RB1',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^13q14'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q12q14) = deletion in chromosome 13 long arm including miR region and likely RB1',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -14 = loss of chromosome 14
        target <- '^14'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-14 = loss of chromosome 14',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -15 = loss of chromosome 15
        target <- '^15'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-15 = loss of chromosome 15',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(17)(p11.2) = loss of 17p including TP53, and gain of unidentified material
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of 17p including TP53, and gain of unidentified material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[16,c('Aberrations', 'Description')]
        idx <- 16

        # add(7)(q11.2) = loss of chromosome 7 long arm and addition of unknown material
        target <- '^7q11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(7)(q11.2) = loss of chromosome 7 long arm and addition of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # der(10)t(10;13)(p11.2;q14 = loss of chromosome 10 short arm and likely loss of proximal 13q including miRNA region
        target <- '^10p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(10)t(10;13)(p11.2;q14 = loss of chromosome 10 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^13q14'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(10)t(10;13)(p11.2;q14 = likely loss of proximal 13q including miRNA region',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(13)(q12q32 = loss of part of 13q including miRNA region and likely RB1 region
        target <- '^13q12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q12q32 = loss of part of 13q including miRNA region and likely RB1 region',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^13q32'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q12q32 = loss of part of 13q including miRNA region and likely RB1 region',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(17)(p11.2 = loss of 17p including TP53
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(17)(p11.2 = loss of 17p including TP53',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(1)(p36.1) = loss of very distal chromosome 1 short arm
        target <- '^1p36.1'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(1)(p36.1) = loss of very distal chromosome 1 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[17,c('Aberrations', 'Description')]
        idx <- 17

        # i(17)(q10 = loss of 17p including TP53 and gain of 17q
        target <- '^17p'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'i(17)(q10 = loss of 17p including TP53 and gain of 17q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17q'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'i(17)(q10 = loss of 17p including TP53 and gain of 17q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(17)(p11.2) = loss of 17p
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(17)(p11.2) = loss of 17p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[18,c('Aberrations', 'Description')]
        idx <- 18

        # -20 = loss of chromosome 20
        target <- '^20'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-20 = loss of chromosome 20',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(14)(q22) = deletion of the chromosome 14 long arm
        target <- '^14q22'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(14)(q22) = deletion of the chromosome 14 long arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[19,c('Aberrations', 'Description')]
        idx <- 19

        # add(17)(p11.2) = loss of the short arm including TP53, and gain of unidentified material that could be derived from 8q, 9q or 14q
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of the short arm including TP53, and gain of unidentified material that could be derived from 8q, 9q or 14q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[20,c('Aberrations', 'Description')]
        idx <- 20

        # add(1)(q21) = a large unidentified block of chromosome material
        target <- '^1q21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(1)(q21) = a large unidentified block of chromosome material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # der(8;17)(q10;q10) = loss of the chromosome 8 short arm and the chromosome 17 short arm  (including TP53)
        target <- '^8q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(8;17)(q10;q10) = loss of the chromosome 8 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(8;17)(q10;q10) = loss of the chromosome 17 short arm  (including TP53)',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # +12 = trisomy 12
        target <- '^12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+12 = trisomy 12',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[21,c('Aberrations', 'Description')]
        idx <- 21

        # add(1)(p36.1) = unidentified material attached to distal 1p
        target <- '^1p36.1'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(1)(p36.1) = unidentified material attached to distal 1p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -15 = loss of chromosome 15
        target <- '^15'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-15 = loss of chromosome 15',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(18)(p11.2) = loss of distal 18p & unidentified material attached to distal 18p
        target <- '^18p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(18)(p11.2) = loss of distal 18p & unidentified material attached to distal 18p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[22,c('Aberrations', 'Description')]
        idx <- 22

        # der(8) and der(13) = 13q deletion & no obvious loss of 8p
        target <- '^13q'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(8) and der(13) = 13q deletion & no obvious loss of 8p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[23,c('Aberrations', 'Description')]
        idx <- 23

        # del(6)(q13q23) = deletion of part of 6q
        target <- '^6q13'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(6)(q13q23) = deletion of part of 6q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^6q23'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(6)(q13q23) = deletion of part of 6q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(7)(p11.2) = loss of chromosome 7 short arm
        target <- '^7p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(7)(p11.2) = loss of chromosome 7 short arm',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(14)(q22q24) = 14q interstitial deletion
        target <- '^14q22'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(14)(q22q24) = 14q interstitial deletion',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^14q24'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(14)(q22q24) = 14q interstitial deletion',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(15)(q24) = loss of distal 15q
        target <- '^15q24'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(15)(q24) = loss of distal 15q',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(17)(p11.2) = 17p (TP53) deletion
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(17)(p11.2) = 17p (TP53) deletion',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[24,c('Aberrations', 'Description')]
        idx <- 24

        # +12 = trisomy 12
        target <- '^12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+12 = trisomy 12',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # der(8;17)(q10;q10) = deletion of chromosome 8 & 17 short arms
        target <- '^8q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(8;17)(q10;q10) = deletion of chromosome 8 & 17 short arms',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^17q10'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'der(8;17)(q10;q10) = deletion of chromosome 8 & 17 short arms',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[25,c('Aberrations', 'Description')]
        idx <- 25

        # -8 = monosomy 8
        target <- '^8'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-8 = monosomy 8',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(17)(p11.2) = loss of 17p (TP53) and gain of unidentified material (?11p)
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of 17p (TP53) and gain of unidentified material (?11p)',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(19)(p13.3) = additional material attached to distal 19p
        target <- '^19p13.3'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(19)(p13.3) = additional material attached to distal 19p',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[26,c('Aberrations', 'Description')]
        idx <- 26

        # -17 = loss of chromosome 17 (TP53)
        target <- '^17'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-17 = loss of chromosome 17 (TP53)',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # del(13)(q14q22) - 13q deletion
        target <- '^13q14'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q14q22) - 13q deletion',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        target <- '^13q22'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'del(13)(q14q22) - 13q deletion',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[27,c('Aberrations', 'Description')]
        idx <- 27

      # ck[28,c('Aberrations', 'Description')]
        idx <- 28

        # +12 = trisomy 12
        target <- '^12'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '+12 = trisomy 12',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[29,c('Aberrations', 'Description')]
        idx <- 29

      # ck[30,c('Aberrations', 'Description')]
        idx <- 30

        # -17 = loss of chromosome 17
        target <- '^17'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-17 = loss of chromosome 17',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[31,c('Aberrations', 'Description')]
        idx <- 31

        # add(11)(q21) = loss of distal 11q and gain of unknown material
        target <- '^11q21'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(11)(q21) = loss of distal 11q and gain of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -13 = monosomy (loss) of chromosome 13
        target <- '^13'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-13 = monosomy (loss) of chromosome 13',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(17)(p11.2) = loss of 17p including TP53 and gain of unknown material
        target <- '^17p11.2'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(17)(p11.2) = loss of 17p including TP53 and gain of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(4)(q31.1) = loss of distal 4q and gain of unknown material
        target <- '^4q31.1'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(4)(q31.1) = loss of distal 4q and gain of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # -8 = loss of chromosome 8
        target <- '^8'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], '-8 = loss of chromosome 8',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

        # add(9)(q34) = loss of very distal 9q and gain of unknown material
        target <- '^9q34'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(9)(q34) = loss of very distal 9q and gain of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[32,c('Aberrations', 'Description')]
        idx <- 32

        # add(5)(p15.1) = loss of very distal 5p and gain of unknown material
        target <- '^5p15.1'
        ggdata <- data.frame(
          Cytoband = rownames(cyto)[grep(target, rownames(cyto))],
          CN = cyto[grep(target, rownames(cyto)),ck$SampleID[idx]])
        res <- rbind(res,
          c(ck$SampleID[idx], 'add(5)(p15.1) = loss of very distal 5p and gain of unknown material',
            mean(ggdata[,2]), min(ggdata[,2]), max(ggdata[,2])))

      # ck[33,c('Aberrations', 'Description')]
        idx <- 33

      colnames(res) <- c('SampleID', 'Description', 'Mean CN', 'Min CN', 'Max CN')
      write.table(res, 'CN_output/ComplexKaryotype.tsv', sep = '\t', row.names = FALSE, quote = FALSE)



    # 2530
      # "NORMAL"
      ggdata <- data.frame(Cytoband = rownames(cyto), CN = cyto[,'2530'])

      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 2) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        theme(axis.text.x = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband (genome-wide)') + ylab('Copy number') +
        labs(title = '2530, \"NORMAL\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.2530.pdf', width = 16, height = 4)
        p1
      dev.off()

    # 3833
      # add(17)(p11.2) = loss of the chromosome 17 short arm (loss of TP53) + posibly part of 8p or 20q;
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^17', rownames(cyto))], CN = cyto[grep('^17', rownames(cyto)),'3833'])
      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '3833, \"add(17)(p11.2) = loss of the chromosome 17 short arm (loss of TP53)\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.3833.pdf', width = 11, height = 3.5)
        cowplot::plot_grid(p1, ncol = 1)
      dev.off()

    # 2623
      # del(13)(q12q14) = deletion in chromosome 13 long arm including miR region and likely RB1;
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^13', rownames(cyto))], CN = cyto[grep('^13', rownames(cyto)),'2623'])
      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '2623, \"del(13)(q12q14) = deletion in chromosome 13 long arm\nincluding miR region and likely RB1\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.2623.pdf', width = 12, height = 4.5)
        cowplot::plot_grid(p1, ncol = 1)
      dev.off()

   # 3498
      # del(14)(q22) = deletion of the chromosome 14 long arm
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^14', rownames(cyto))], CN = cyto[grep('^14', rownames(cyto)),'3498'])
      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '3498, \"del(14)(q22) = deletion of the chromosome 14 long arm\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.3498.pdf', width = 10, height = 4.5)
        cowplot::plot_grid(p1, ncol = 1)
      dev.off()

    # 3764
      # der(8;17)(q10;q10) = loss of the chromosome 8 short arm and the chromosome 17 short arm  (including TP53);
      # +12 = trisomy 12;
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^8', rownames(cyto))], CN = cyto[grep('^8', rownames(cyto)),'3764'])
      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '3764, \"der(8;17)(q10;q10) = loss of the chromosome 8 short arm\nand the chromosome 17 short arm  (including TP53)\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^17', rownames(cyto))], CN = cyto[grep('^17', rownames(cyto)),'3764'])
      p2 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '3764, \"der(8;17)(q10;q10) = loss of the chromosome 8 short arm\nand the chromosome 17 short arm  (including TP53)\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^12', rownames(cyto))], CN = cyto[grep('^12', rownames(cyto)),'3764'])
      p3 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '3764, \"+12 = trisomy 12\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.3764.pdf', width = 12, height = 12)
        cowplot::plot_grid(p1, p2, p3, ncol = 1)
      dev.off()

    # 4074
      # 13 del(6)(q13q23) = deletion of part of 6q;
      ggdata <- data.frame(Cytoband = rownames(cyto)[grep('^6', rownames(cyto))], CN = cyto[grep('^6', rownames(cyto)),'4074'])
      p1 <- ggplot(data = ggdata, aes(x = Cytoband, y = CN, label = Cytoband)) +
        geom_point(colour = 'black', size = 3) +
        #geom_line() +
        theme_bw(base_size = 24) + mytheme +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        xlab('Cytoband') + ylab('Copy number') +
        labs(title = '4074, \"13 del(6)(q13q23) = deletion of part of 6q\"', subtitle = NULL,
          caption = NULL) +
        #geom_label_repel(size = 2) +
        coord_cartesian(ylim = c(0, 4)) +
        scale_y_continuous(
          breaks = seq(0, 4, 1))
      pdf('CN_output/ComplexKaryotype.4074.pdf', width = 10, height = 3.5)
        cowplot::plot_grid(p1, ncol = 1)
      dev.off()

