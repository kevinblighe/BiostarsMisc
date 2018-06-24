TCGAanalyze_DEA_Fix_aluesley1 <- function (mat1, mat2, metadata = TRUE, Cond1type, Cond2type, 
    pipeline = "edgeR", method = "exactTest", fdr.cut = 1, logFC.cut = 0, 
    elementsRatio = 30000, batch.factors = NULL, ClinicalDF = data.frame(), 
    paired = FALSE, log.trans = FALSE, voom = FALSE, trend = FALSE, 
    MAT = data.frame(), contrast.formula = "", Condtypes = c()) 
{
    table.code <- c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", 
        "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM", "CELLC", 
        "TRB", "CELL", "XP", "XCL")
    names(table.code) <- c("01", "02", "03", "04", "05", "06", 
        "07", "08", "09", "10", "11", "12", "13", "14", "20", 
        "40", "50", "60", "61")
    if (nrow(MAT) == 0) {
        TOC <- cbind(mat1, mat2)
        Cond1num <- ncol(mat1)
        Cond2num <- ncol(mat2)
    }
    else {
        TOC <- MAT
    }
    if (metadata == TRUE) {
        my_IDs <- get_IDs(TOC)
        Plate <- factor(my_IDs$plate)
        Condition <- factor(my_IDs$condition)
        TSS <- factor(my_IDs$tss)
        Portion <- factor(my_IDs$portion)
        Center <- factor(my_IDs$center)
        Patients <- factor(my_IDs$patient)
    }
    if (paired == TRUE) {
        matched.query <- TCGAquery_MatchedCoupledSampleTypes(my_IDs$barcode, c("TP","NT"))
        my_IDs <- subset(my_IDs, barcode == matched.query)
        TOC <- TOC[, which(colnames(TOC) %in% matched.query)]
    }
    if (nrow(ClinicalDF) > 0) {
        names(ClinicalDF)[names(ClinicalDF) == "bcr_patient_barcode"] <- "patient"
        ClinicalDF$age_at_diag_year <- floor(clinical$age_at_diagnosis/365)
        ClinicalDF$diag_year <- ClinicalDF$age_at_diag_year + 
            clinical$year_of_birth
        diag_yearDF <- ClinicalDF[, c("patient", "diag_year")]
        my_IDs <- merge(my_IDs, ClinicalDF, by = "patient")
        Year <- as.factor(my_IDs$diag_year)
    }
    options <- c("Plate", "TSS", "Year", "Portion", "Center", 
        "Patients")
    if (length(batch.factors) == 0) {
        message("Batch correction skipped since no factors provided")
    }
    else for (o in batch.factors) {
        if (o %in% options == FALSE) 
            stop(paste0(o, " is not a valid batch correction factor"))
        if (o == "Year" & nrow(ClinicalDF) == 0) 
            stop("batch correction using diagnosis year needs clinical info. Provide Clinical Data in arguments")
    }
    additiveformula <- paste(batch.factors, collapse = "+")
    message("----------------------- DEA -------------------------------")
    if (nrow(MAT) == 0) {
        message(message1 <- paste("there are Cond1 type", Cond1type, 
            "in ", Cond1num, "samples"))
        message(message2 <- paste("there are Cond2 type", Cond2type, 
            "in ", Cond2num, "samples"))
        message(message3 <- paste("there are ", nrow(TOC), "features as miRNA or genes "))
    }
    else {
        message(message3 <- paste("there are ", nrow(TOC), "features as miRNA or genes "))
    }
    timeEstimated <- format(ncol(TOC) * nrow(TOC)/elementsRatio, 
        digits = 2)
    message(messageEstimation <- paste("I Need about ", timeEstimated, 
        "seconds for this DEA. [Processing 30k elements /s]  "))
    colnames(TOC) <- paste0("s", 1:ncol(TOC))
    if (length(Condtypes) > 0) {
        tumorType <- factor(x = Condtypes, levels = unique(Condtypes))
    }
    else {
        tumorType <- factor(x = rep(c(Cond1type, Cond2type), 
            c(Cond1num, Cond2num)), levels = c(Cond1type, Cond2type))
    }
    if (length(batch.factors) == 0 & length(Condtypes) > 0) {
        if (pipeline == "edgeR") 
            design <- model.matrix(~tumorType)
        else design <- model.matrix(~0 + tumorType)
    }
    else if (length(batch.factors) == 0 & length(Condtypes) == 
        0) {
        if (pipeline == "edgeR") 
            design <- model.matrix(~tumorType)
        else design <- model.matrix(~0 + tumorType)
    }
    else if (length(batch.factors) > 0 & length(Condtypes) == 
        0) {
        if (pipeline == "edgeR") 
            formula <- paste0("~tumorType+", additiveformula)
        else formula <- paste0("~0+tumorType+", additiveformula)
        design <- model.matrix(eval(parse(text = formula)))
    }
    else if (length(batch.factors) > 0 & length(Condtypes) > 
        0) {
        if (pipeline == "edgeR") 
            formula <- paste0("~tumorType+", additiveformula)
        else formula <- paste0("~0+tumorType+", additiveformula)
        design <- model.matrix(eval(parse(text = formula)))
    }
    if (pipeline == "edgeR") {
        if (method == "exactTest") {
            DGE <- edgeR::DGEList(TOC, group = rep(c(Cond1type, 
                Cond2type), c(Cond1num, Cond2num)))
            disp <- edgeR::estimateCommonDisp(DGE)
            tested <- edgeR::exactTest(disp, pair = c(Cond1type, 
                Cond2type))
            logFC_table <- tested$table
            tableDEA <- edgeR::topTags(tested, n = nrow(tested$table))$table
            tableDEA <- tableDEA[tableDEA$FDR <= fdr.cut, ]
            tableDEA <- tableDEA[abs(tableDEA$logFC) >= logFC.cut, 
                ]
        }
        else if (method == "glmLRT") {
            if (length(unique(tumorType)) == 2) {
                aDGEList <- edgeR::DGEList(counts = TOC, group = tumorType)
                aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, 
                  design)
                aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, 
                  design)
                aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion, 
                  prior.count.total = 0)
                aGlmLRT <- edgeR::glmLRT(aGlmFit, coef = 2)
                tableDEA <- cbind(aGlmLRT$table, FDR = p.adjust(aGlmLRT$table$PValue, 
                  "fdr"))
                tableDEA <- tableDEA[tableDEA$FDR < fdr.cut, 
                  ]
                tableDEA <- tableDEA[abs(tableDEA$logFC) > logFC.cut, 
                  ]
                if (all(grepl("ENSG", rownames(tableDEA)))) 
                  tableDEA <- cbind(tableDEA, map.ensg(genes = rownames(tableDEA))[, 
                    2:3])
            }
            else if (length(unique(tumorType)) > 2) {
                aDGEList <- edgeR::DGEList(counts = TOC, group = tumorType)
                colnames(design)[1:length(levels(tumorType))] <- levels(tumorType)
                prestr = "makeContrasts("
                poststr = ",levels=colnames(design))"
                commandstr = paste(prestr, contrast.formula, 
                  poststr, sep = "")
                commandstr = paste0("limma::", commandstr)
                cont.matrix <- eval(parse(text = commandstr))
                aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, 
                  design)
                aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, 
                  design)
                aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion, 
                  prior.count.total = 0)
                print(cont.matrix)
                tableDEA <- list()
                for (mycoef in colnames(cont.matrix)) {
                  message(paste0("DEA for", " :", mycoef))
                  aGlmLRT <- edgeR::glmLRT(aGlmFit, contrast = cont.matrix[, 
                    mycoef])
                  print("---toptags---")
                  print(topTags(aGlmLRT, adjust.method = "fdr", 
                    sort.by = "PValue"))
                  tt <- aGlmLRT$table
                  tt <- cbind(tt, FDR = p.adjust(aGlmLRT$table$PValue, 
                    "fdr"))
                  tt <- tt[(tt$FDR < fdr.cut & abs(as.numeric(tt$logFC)) > 
                    logFC.cut), ]
                  tableDEA[[as.character(mycoef)]] <- tt
                  if (all(grepl("ENSG", rownames(tableDEA[[as.character(mycoef)]])))) 
                    tableDEA[[as.character(mycoef)]] <- cbind(tableDEA[[as.character(mycoef)]], 
                      map.ensg(genes = rownames(tableDEA[[as.character(mycoef)]]))[, 
                        2:3])
                }
            }
        }
        else stop(paste0(method, " is not a valid DEA method option. Choose 'exactTest' or 'glmLRT' "))
    }
    else if (pipeline == "limma") {
        if (log.trans == TRUE) 
            logCPM <- edgeR::cpm(TOC, log = TRUE, prior.count = 3)
        else logCPM <- TOC
        if (voom == TRUE) {
            message("Voom Transformation...")
            logCPM <- limma::voom(logCPM, design)
        }
        if (length(unique(tumorType)) == 2) {
            colnames(design)[1:2] <- c(Cond1type, Cond2type)
            contr <- paste0(Cond2type, "-", Cond1type)
            cont.matrix <- limma::makeContrasts(contrasts = contr, 
                levels = design)
            fit <- limma::lmFit(logCPM, design)
            fit <- contrasts.fit(fit, cont.matrix)
            if (trend == TRUE) {
                fit <- limma::eBayes(fit, trend = TRUE)
            }
            else {
                fit <- limma::eBayes(fit, trend = FALSE)
            }
            tableDEA <- limma::topTable(fit, coef = 1, adjust.method = "fdr", 
                number = nrow(TOC))
            limma::volcanoplot(fit, highlight = 10)
            index <- which(tableDEA[, 4] < fdr.cut)
            tableDEA <- tableDEA[index, ]
            neg_logFC.cut <- -1 * logFC.cut
            index <- which(abs(as.numeric(tableDEA$logFC)) > 
                logFC.cut)
            tableDEA <- tableDEA[index, ]
        }
        else if (length(unique(tumorType)) > 2) {
            DGE <- edgeR::DGEList(TOC, group = tumorType)
            colnames(design)[1:length(levels(tumorType))] <- levels(tumorType)
            prestr = "makeContrasts("
            poststr = ",levels=colnames(design))"
            commandstr = paste(prestr, contrast.formula, poststr, 
                sep = "")
            commandstr = paste0("limma::", commandstr)
            cont.matrix <- eval(parse(text = commandstr))
            fit <- limma::lmFit(logCPM, design)
            fit <- limma::contrasts.fit(fit, cont.matrix)
            if (trend == TRUE) 
                fit <- limma::eBayes(fit, trend = TRUE)
            else fit <- limma::eBayes(fit, trend = FALSE)
            tableDEA <- list()
            for (mycoef in colnames(cont.matrix)) {
                tableDEA[[as.character(mycoef)]] <- limma::topTable(fit, 
                  coef = mycoef, adjust.method = "fdr", number = nrow(MAT))
                message(paste0("DEA for", " :", mycoef))
                tempDEA <- tableDEA[[as.character(mycoef)]]
                index.up <- which(tempDEA$adj.P.Val < fdr.cut & 
                  abs(as.numeric(tempDEA$logFC)) > logFC.cut)
                tableDEA[[as.character(mycoef)]] <- tempDEA[index.up, 
                  ]
                if (all(grepl("ENSG", rownames(tableDEA[[as.character(mycoef)]])))) 
                  tableDEA[[as.character(mycoef)]] <- cbind(tableDEA[[as.character(mycoef)]], 
                    map.ensg(genes = rownames(tableDEA[[as.character(mycoef)]]))[, 
                      2:3])
            }
        }
    }
    else stop(paste0(pipeline, " is not a valid pipeline option. Choose 'edgeR' or 'limma'"))
    message("----------------------- END DEA -------------------------------")
    return(tableDEA)
}

