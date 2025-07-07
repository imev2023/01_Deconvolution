###### Adapted from 14 Years of Deconv Review 
###### CTVT specific script --> change name of column from which we get the cell type annotations. 
# Helper function to compute MAE
calc_mae <- function(estimate, truth) {
  
  # Join on sampleName and manual.log.anot.fine
  joined <- estimate %>%
    inner_join(truth, by = c("sampleNameEdit", "manual.log.anot.fine"), suffix = c(".est", ".truth"))
  
  # Compute absolute differences of cell types per sampleNameEdit. Take mean MAE per sample. 
  mae <- joined %>%
    mutate(abs_error = abs(cellProportion.est - cellProportion.truth)) %>%
    summarise(mae = mean(abs_error)) %>%
    pull(mae)
  
  return(mae)
}

# Helper function to compute SCorr
calc_scorr <- function(estimate, truth) {
  
  # Join on sampleNameEdit and manual.log.anot.fine
  joined <- estimate %>%
    inner_join(truth, by = c("sampleNameEdit", "manual.log.anot.fine"), suffix = c(".est", ".truth"))
  
  # Compute Spearman correlation for each sample
  scorr_df <- joined %>%
    group_by(sampleNameEdit) %>%
    summarise(
      scorr = cor(cellProportion.est, cellProportion.truth, method = "spearman"),
      .groups = "drop"
    )
  
  # Return average SCorr
  mean(scorr_df$scorr, na.rm = TRUE)
}

calc_ccorr <- function(estimate, truth) {
  # Join on sampleNameEdit and manual.log.anot.fine
  joined <- estimate %>%
    inner_join(truth, by = c("sampleNameEdit", "manual.log.anot.fine"), suffix = c(".est", ".truth"))
  
  # Compute Spearman correlation for each cell type
  ccorr_df <- joined %>%
    group_by(manual.log.anot.fine) %>%
    summarise(
      ccorr = cor(cellProportion.est, cellProportion.truth, method = "spearman"),
      .groups = "drop"
    )
  
  # Return average CCorr
  mean(ccorr_df$ccorr, na.rm = TRUE)
}

# Helper to calculate MAECorr
calc_maecorr <- function(estimate, truth) {
  est_mat <- reshape_to_matrix(estimate)
  truth_mat <- reshape_to_matrix(truth)
  
  # Ensure matching samples
  common_samples <- intersect(rownames(est_mat), rownames(truth_mat))
  est_mat <- est_mat[common_samples, , drop = FALSE]
  truth_mat <- truth_mat[common_samples, , drop = FALSE]
  
  # Compute Pearson correlation matrices (sample-wise)
  cor_est <- cor(t(est_mat), method = "pearson")
  cor_truth <- cor(t(truth_mat), method = "pearson")
  
  # Compute mean absolute error between correlation matrices
  mae_corr <- mean(abs(cor_est - cor_truth), na.rm = TRUE)
  return(mae_corr)
}

# Helper to calculate MAECorr
calc_rmse <- function(estimate, truth) {
  
  # Join on sampleNameEdit and manual.log.anot.fine
  joined <- estimate %>%
    inner_join(truth, by = c("sampleNameEdit", "manual.log.anot.fine"), suffix = c(".est", ".truth"))
  
  # Compute Spearman correlation for each cell type
  rmse_df <- joined %>%
    group_by(manual.log.anot.fine) %>%
    summarise(
      rmse = sqrt(mean(((cellProportion.est - cellProportion.truth) * 100)^2)),
      .groups = "drop"
    )
  
  # Return average CCorr
  mean(rmse_df$rmse, na.rm = TRUE)
}

# Helper to reshape to matrix: rows = samples, cols = cell types
reshape_to_matrix <- function(df, 
                              rownames = "sampleNameEdit",
                              rOrder = NULL, 
                              cOrder = NULL) {
  matrix <- df %>%
    select(all_of(rownames), manual.log.anot.fine, cellProportion) %>%
    pivot_wider(names_from = manual.log.anot.fine, values_from = cellProportion) %>%
    column_to_rownames(var = rownames) %>%
    as.matrix()
  
  # Reorder rows if rOrder is provided
  if (!is.null(rOrder)) {
    matrix <- matrix[rOrder, , drop = FALSE]
  }
  
  # Reorder columns if cOrder is provided
  if (!is.null(cOrder)) {
    matrix <- matrix[, cOrder, drop = FALSE]
  }
  
  return(matrix)
}

unshape_matrix <- function(mat, 
                           rowname_col, 
                           colname_col, 
                           value_col) {
  df <- as.data.frame(mat) %>%
    rownames_to_column(var = rowname_col) %>%
    pivot_longer(
      cols = -all_of(rowname_col),
      names_to = colname_col,
      values_to = value_col
    )
  return(df)
}


absErrorMatrices <- function(plotDf, 
                             keyColumn, 
                             sOrder, 
                             ctypeOrder){
  
  # plotDf is a dataframe derived from masterPropDf outputs of different deconvolution methods. 
  # Importantly, it must contain a column (keyColumn argument - typeName in my visualization scripts with the type of sample (bc, bn, or sn) "_" and the sample ID. 
  # sOrder is the order of the keyColumn (i.e. how you want to organize the rows.)
  # cOrder is the order of the cell types (i.e. how you want to organize the columns. )
  # Output: Three matrices for plotting (sample x cell types.)
  
  # Pre-populating matrices. 
  # Importantly rows will be matched to the rownames so really important to ensure columns are ordered in the same way as cOrder in the per-sample loop. 
  
  # rmseMatrix <- matrix(NA, 
  #                      nrow = length(sOrder), 
  #                      ncol = length(ctypeOrder), 
  #                      dimnames = list(sOrder, 
  #                                      ctypeOrder))
  
  # rpeMatrix <- matrix(NA, 
  #                     nrow = length(sOrder), 
  #                     ncol = length(ctypeOrder), 
  #                     dimnames = list(sOrder, 
  #                                     ctypeOrder))
  
  abseMatrix <- matrix(NA, 
                       nrow = length(sOrder), 
                       ncol = length(ctypeOrder), 
                       dimnames = list(sOrder, 
                                       ctypeOrder))
  
  errMatrix <- matrix(NA, 
                      nrow = length(sOrder), 
                      ncol = length(ctypeOrder), 
                      dimnames = list(sOrder, 
                                      ctypeOrder))
  
  # For each sample (bc or bn), for each cell type, calculate three error metrics.
  # Samples defined by edited name, because changes with bulk/nuc. 
  
  samples <- unique(plotDf$sampleNameEdit) # i.e. 3176T3
  
  for (samp in samples) {
    
    # Take three rows (sn_, bc_, and bn_ for a given sample). Of note, some samples do not have a matching reference so we cannot include them in this matrix. 
    
    tempSampleDf <- plotDf %>%
      filter(sampleNameEdit %in% samp)
    
    if(nrow(tempSampleDf) != ((length(unique(ctypeOrder)))*3) ){ 
      print(paste0("The number of results for sampleNameEdit: ", 
                   samp, 
                   " is not length(unique(ctypeOrder)*3. This sample was likely missing a matched reference in the sn data. Skipping.")
      )
      next
    }
    
    # Create sn, bc, and bn matrices. 
    snSampleRef <- tempSampleDf %>%
      filter(sampleType == "sn")
    
    sampTruth <- reshape_to_matrix(snSampleRef, 
                                   rownames = keyColumn, 
                                   rOrder = NULL, # Because we are only returning one sample at a time.  
                                   cOrder = ctypeOrder)
    
    bcSampleRef <- tempSampleDf %>%
      filter(sampleType == "bc")
    
    bcName <- unique(bcSampleRef[, keyColumn])
    
    sampBc <- reshape_to_matrix(bcSampleRef, 
                                rownames = keyColumn, 
                                rOrder = NULL, 
                                cOrder = ctypeOrder)
    
    bnSampleRef <- tempSampleDf %>%
      filter(sampleType == "bn")
    
    bnName <- unique(bnSampleRef[, keyColumn])
    
    sampBn <- reshape_to_matrix(bnSampleRef, 
                                rownames = keyColumn, 
                                rOrder = NULL, 
                                cOrder = ctypeOrder)
    
    #### RMSE ####
    # rmseBc <- rmse(sampTruth, sampBc)
    # rmseMatrix[bcName, ] <- rmseBc[, ctypeOrder]
    # rm(rmseBc)
    
    # rmseBn <- rmse(sampTruth, sampBn)
    # rmseMatrix[bnName, ] <- rmseBn[, ctypeOrder]
    # rm(rmseBn)
    
    
    #### Absolute error ####
    abseBc <- abs(sampTruth - sampBc)
    abseMatrix[bcName, ] <- abseBc[, ctypeOrder]
    rm(abseBc)
    
    abseBn <- abs(sampTruth - sampBn)
    abseMatrix[bnName, ] <- abseBn[, ctypeOrder]
    rm(abseBn)
    
    #### Error ####
    errBc <- sampTruth - sampBc
    errMatrix[bcName, ] <- errBc[, ctypeOrder]
    rm(errBc)
    
    errBn <- abs(sampTruth - sampBn)
    errMatrix[bnName, ] <- errBn[, ctypeOrder]
    rm(errBn)
    
    # #### RPE (relative prediction error) ####
    # safeTruth <- ifelse(sampTruth == 0, 0.00001, sampTruth)
    # rpeBc <- abs(sampTruth - sampBc) / safeTruth
    # rpeMatrix[bcName, ] <- rpeBc[, ctypeOrder]
    # rm(rpeBc)
    # 
    # rpeBn <- abs(sampTruth - sampBn) / safeTruth
    # rpeMatrix[bnName, ] <- rpeBn[, ctypeOrder]
    # rm(rpeBn)
  }
  
  
  return(list("absErrMatrix" = abseMatrix[sOrder, ctypeOrder], 
              "errMatrix" = errMatrix[sOrder, ctypeOrder]))
  
}

celltypeSampleMetricsMatrices <- function(plotDf, 
                                          keyColumn, 
                                          sOrder, 
                                          ctypeOrder){
  
  # plotDf is a dataframe derived from masterPropDf outputs of different deconvolution methods. 
  # Importantly, it must contain a column (keyColumn argument - typeName in my visualization scripts with the type of sample (bc, bn, or sn) "_" and the sample ID. 
  # sOrder is the order of the keyColumn (i.e. how you want to organize the rows.)
  # cOrder is the order of the cell types (i.e. how you want to organize the columns. )
  # Output: Three matrices for plotting (sample x cell types.)
  
  # Pre-populating matrices. 
  # Importantly rows will be matched to the rownames so really important to ensure columns are ordered in the same way as cOrder in the per-sample loop. 
  # Sample matrices have 1 column and sOrder rows. Cell type matrices have 1 row and cType columns. 
  rmseMatrixS <- matrix(NA,
                        nrow = length(sOrder),
                        ncol = 1,
                        dimnames = list(sOrder,
                                        "RMSE"))
  
  rpeMatrixS <- matrix(NA,
                       nrow = length(sOrder),
                       ncol = 1,
                       dimnames = list(sOrder,
                                       "RPE"))
  
  maeMatrixS <- matrix(NA,
                       nrow = length(sOrder),
                       ncol = 1,
                       dimnames = list(sOrder,
                                       "Mean Abs. Error"))
  
  medMatrixS <- matrix(NA,
                       nrow = length(sOrder),
                       ncol = 1,
                       dimnames = list(sOrder,
                                       "Median Abs. Error"))
  
  rmseMatrixC <- matrix(NA,
                        nrow = 1,
                        ncol = length(ctypeOrder),
                        dimnames = list("RMSE",
                                        ctypeOrder))
  
  rpeMatrixC <- matrix(NA,
                       nrow = 1,
                       ncol = length(ctypeOrder),
                       dimnames = list("RPE",
                                       ctypeOrder))
  
  maeMatrixC <- matrix(NA,
                       nrow = 1,
                       ncol = length(ctypeOrder),
                       dimnames = list("Mean Abs. Error",
                                       ctypeOrder))
  
  medMatrixC <- matrix(NA,
                       nrow = 1,
                       ncol = length(ctypeOrder),
                       dimnames = list("Median Abs. Error",
                                       ctypeOrder))
  
  
  
  
  # For each sample (bc or bn) or for each cell type, calculate four error metrics.
  # Samples defined by edited name, because changes with bulk/nuc. 
  
  samples <- unique(plotDf$sampleNameEdit) # i.e. 3176T3
  
  for (samp in samples) {
    
    # Take three rows (sn_, bc_, and bn_ for a given sample). Of note, some samples do not have a matching reference so we cannot include them in this matrix. 
    
    tempSampleDf <- plotDf %>%
      filter(sampleNameEdit %in% samp)
    
    if(nrow(tempSampleDf) != ((length(unique(ctypeOrder)))*3) ){ 
      print(paste0("The number of results for sampleNameEdit: ", 
                   samp, 
                   " is not length(unique(ctypeOrder)*3. This sample was likely missing a matched reference in the sn data. Skipping.")
      )
      next
    }
    
    ###############################
    ######## Sample-wise ##########
    ###############################
    
    # Create sn, bc, and bn matrices. 
    snSampleRef <- tempSampleDf %>%
      filter(sampleType == "sn")
    
    sampTruth <- reshape_to_matrix(snSampleRef, 
                                   rownames = keyColumn, 
                                   rOrder = NULL, # Because we are only returning one sample at a time.  
                                   cOrder = ctypeOrder)
    
    bcSampleRef <- tempSampleDf %>%
      filter(sampleType == "bc")
    
    bcName <- unique(bcSampleRef[, keyColumn])
    
    sampBc <- reshape_to_matrix(bcSampleRef, 
                                rownames = keyColumn, 
                                rOrder = NULL, 
                                cOrder = ctypeOrder)
    
    bnSampleRef <- tempSampleDf %>%
      filter(sampleType == "bn")
    
    bnName <- unique(bnSampleRef[, keyColumn])
    
    sampBn <- reshape_to_matrix(bnSampleRef, 
                                rownames = keyColumn, 
                                rOrder = NULL, 
                                cOrder = ctypeOrder)
    
    #### RMSE ####
    # rmseBc <- rmse(sampTruth, sampBc)
    # rmseMatrix[bcName, ] <- rmseBc[, ctypeOrder]
    # rm(rmseBc)
    
    # rmseBn <- rmse(sampTruth, sampBn)
    # rmseMatrix[bnName, ] <- rmseBn[, ctypeOrder]
    # rm(rmseBn)
    
    
    #### Absolute error ####
    abseBc <- sampTruth - sampBc
    abseMatrix[bcName, ] <- abseBc[, ctypeOrder]
    rm(abseBc)
    
    abseBn <- abs(sampTruth - sampBn)
    abseMatrix[bnName, ] <- abseBn[, ctypeOrder]
    rm(abseBn)
    
    # #### RPE (relative prediction error) ####
    # safeTruth <- ifelse(sampTruth == 0, 0.00001, sampTruth)
    # rpeBc <- abs(sampTruth - sampBc) / safeTruth
    # rpeMatrix[bcName, ] <- rpeBc[, ctypeOrder]
    # rm(rpeBc)
    # 
    # rpeBn <- abs(sampTruth - sampBn) / safeTruth
    # rpeMatrix[bnName, ] <- rpeBn[, ctypeOrder]
    # rm(rpeBn)
    
    ###############################
    ######## Cell Type-wise ##########
    ###############################
    
  }
  
  
  return(list("rmseMatrix" = rmseMatrix[sOrder, ctypeOrder], 
              "absErrMatrix" = abseMatrix[sOrder, ctypeOrder], 
              "rpeMatrix" = rpeMatrix[sOrder, ctypeOrder]))
  
}



###### Adapted from https://github.com/MedicalGenomicsLab/deconvolution_benchmarking/tree/master
library(Metrics)
library(dplyr)

calculateCtypeMetrics <- function(truth, estimate) {
  
  # Initialize list to store metrics for each cell type
  metrics_list <- list()
  
  # cTypes always defined by ground truth. Make dynamic!!
  cTypes <- unique(truth$manual.log.anot.fine)
  
  for (c_type in cTypes) {
    
    # Extract the columns for the specific cell type
    ctype_truth <- reshape_to_matrix(truth)[,c_type]
    ctype_preds <- reshape_to_matrix(estimate)[,c_type]
    
    # RMSE (scaled by 100)
    rmse_val <- rmse(ctype_truth, ctype_preds)
    
    # MAE (median absolute error)
    mae_val <- median(abs(ctype_truth - ctype_preds))
    
    # RPE (relative prediction error)
    safe_truth <- ifelse(ctype_truth == 0, 0.00001, ctype_truth)
    rpe_val <- median(abs(ctype_truth - ctype_preds) / safe_truth)
    
    # Create named vector for the cell type metrics
    metrics_vec <- c(RMSE = rmse_val, MAE = mae_val, RPE = rpe_val)
    metrics_list[[c_type]] <- metrics_vec
  }
  
  # Convert list of metrics to a data frame (wide format)
  ctype_method_metrics_df <- as.data.frame(do.call(cbind, metrics_list)) %>%
    rownames_to_column(var = "Metric") %>%
    pivot_longer(-Metric, names_to = "Subject", values_to = "Value") %>%
    mutate(Comparison = "perCellType") %>%
    select(Comparison, Subject, Metric, Value)
  
  return(ctype_method_metrics_df)
}


calculateSampleMetrics <- function(truth, estimate) {
  
  # Initialize list to store metrics for each cell type
  metrics_list <- list()
  
  # Samples defined by edited name, because changes with bulk/nuc. 
  samples <- intersect(unique(truth$sampleNameEdit), unique(estimate$sampleNameEdit))
  
  for (samp in samples) {
    
    # Extract the rows of a given sample
    samp_truth <- reshape_to_matrix(truth)[samp,] 
    
    samp_preds <- reshape_to_matrix(estimate)[samp,]
    
    # RMSE (scaled by 100)
    rmse_val <- rmse(samp_truth, samp_preds)
    # rmse_scaled_val <- rmse(samp_truth * 100, samp_preds * 100)
    
    # MAE (median absolute error scaled by 100)
    mae_val <- median(abs(samp_truth - samp_preds)) #* 100
    
    # RPE (relative prediction error)
    safe_truth <- ifelse(samp_truth == 0, 0.00001, samp_truth)
    rpe_val <- median(abs(samp_truth - samp_preds) / safe_truth)
    
    # Create named vector for the cell type metrics
    metrics_vec <- c(RMSE = rmse_val, MAE = mae_val, RPE = rpe_val)
    metrics_list[[samp]] <- metrics_vec
  }
  
  # Convert list of metrics to a data frame
  samp_method_metrics_df <- as.data.frame(do.call(cbind, metrics_list)) %>%
    rownames_to_column(var = "Metric") %>%
    pivot_longer(-Metric, names_to = "Subject", values_to = "Value") %>%
    mutate(Comparison = "perSample") %>%
    select(Comparison, Subject, Metric, Value)
  
  return(samp_method_metrics_df)
}

# Predicted - ground truth
calculateRawPE <- function(truth, estimate){
  samples <- intersect(unique(truth$sampleNameEdit), unique(estimate$sampleNameEdit))
  
  truthReshape <- reshape_to_matrix(truth)[samples,]
  estimateReshape <- reshape_to_matrix(estimate)[samples,]
  rawPE <- estimateReshape - truthReshape
  
}

###### Compute and save all metrics
deconError <- function(data, 
                       parameters, 
                       method, 
                       reference, 
                       seed){
  
  # Filter data into ground truth and estimated sets
  ground_truth <- data %>% 
    filter(sampleType == "sn") 
  
  # Get reference sample names
  groundTruthSampleName <- unique(ground_truth$sampleName)
  groundTruthSampleNameEdit <- unique(ground_truth$sampleNameEdit)
  
  # Estimated bulk cell data (bc)
  est_bc <- data %>% 
    filter(sampleType == "bc")  %>%
    filter(sampleNameEdit %in% groundTruthSampleNameEdit)
  
  # Estimated bulk nucleus data (bn)
  est_bn <- data %>%
    filter(sampleType == "bn") %>%
    filter(sampleName %in% groundTruthSampleName)
  
  mae_bc <- calc_mae(est_bc, ground_truth)
  mae_bn <- calc_mae(est_bn, ground_truth)
  
  mae <- data.frame(bulkCellEstimate = mae_bc, 
                    bulkNucEstimate = mae_bn, 
                    Method = method, 
                    Parameters = parameters, 
                    Reference = reference, 
                    Seed = seed,
                    Metric = "MAE")
  
  scorr_bc <- calc_scorr(est_bc, ground_truth)
  scorr_bn <- calc_scorr(est_bn, ground_truth)
  
  scorr <- data.frame(bulkCellEstimate = scorr_bc, 
                      bulkNucEstimate = scorr_bn, 
                      Method = method, 
                      Parameters = parameters, 
                      Reference = reference, 
                      Seed = seed,
                      Metric = "SCorr")
  
  ccorr_bc <- calc_ccorr(est_bc, ground_truth)
  ccorr_bn <- calc_ccorr(est_bn, ground_truth)
  
  ccorr <- data.frame(bulkCellEstimate = ccorr_bc, 
                      bulkNucEstimate = ccorr_bn, 
                      Method = method, 
                      Parameters = parameters, 
                      Reference = reference, 
                      Seed = seed,
                      Metric = "CCorr")
  
  
  maecorr_bc <- calc_maecorr(est_bc, ground_truth)
  maecorr_bn <- calc_maecorr(est_bn, ground_truth)
  
  maecorr <- data.frame(bulkCellEstimate = maecorr_bc, 
                        bulkNucEstimate = maecorr_bn, 
                        Method = method, 
                        Parameters = parameters, 
                        Reference = reference, 
                        Seed = seed,
                        Metric = "MAECorr")
  
  rmse_bc <- calc_rmse(est_bc, ground_truth)
  rmse_bn <- calc_rmse(est_bn, ground_truth)
  
  rmse <- data.frame(bulkCellEstimate = rmse_bc, 
                     bulkNucEstimate = rmse_bn, 
                     Method = method, 
                     Parameters = parameters, 
                     Reference = reference, 
                     Seed = seed,
                     Metric = "RMSE")
  
  ## New metrics - Tran et al. Deconvolution of TME results. 
  ctype_bc <- calculateCtypeMetrics(ground_truth, est_bc) %>%
    mutate(Estimate = "bulkCell",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  ctype_bn <- calculateCtypeMetrics(ground_truth, est_bn) %>%
    mutate(Estimate = "bulkNuc",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  samp_bc <- calculateSampleMetrics(ground_truth, est_bc) %>%
    mutate(Estimate = "bulkCell",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  samp_bn <- calculateSampleMetrics(ground_truth, est_bn) %>%
    mutate(Estimate = "bulkNuc",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  
  # Calling this rawPE because RPE = relative proportion error. 
  # Predicted - ground truth.
  rawPE_bn <- unshape_matrix(calculateRawPE(ground_truth, est_bn), 
                             value_col = "rawPE", 
                             colname_col = "manual.log.anot.fine",
                             rowname_col = "sampleNameEdit") %>%
    mutate(Estimate = "bulkNuc",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  rawPE_bc <- unshape_matrix(calculateRawPE(ground_truth, est_bc),
                             value_col = "rawPE", 
                             colname_col = "manual.log.anot.fine",
                             rowname_col = "sampleNameEdit") %>%
    mutate(Estimate = "bulkCell",
           Method = method, 
           Parameters = parameters, 
           Reference = reference, 
           Seed = seed)
  
  
  
  resList <- list("overviewMetrics" = rbind(mae, scorr, ccorr, maecorr, rmse), 
                  "granularMetrics" = rbind(ctype_bc, ctype_bn, samp_bc, samp_bn), 
                  "rawPredictionError" = rbind(rawPE_bc, rawPE_bn))
  
  
  return(resList)
  
}
