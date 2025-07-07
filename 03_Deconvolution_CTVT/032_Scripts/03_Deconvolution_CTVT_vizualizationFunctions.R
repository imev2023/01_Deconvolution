rmseByCellTypeMean <- function(df, pred_type) {
  
  dfPred <- df %>%
    filter(sampleType %in% c(pred_type)) %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion, 
           Method_Param_Reference)
  
  dfTruth <- df %>%
    filter(sampleType == "sn") %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion) %>%
    dplyr::rename(truth = cellProportion)
  
  # 7 single-cell samples, 10 bulk cell and 10 bulk nuc samples. 
  intersect <- intersect(dfPred$sampleNameEdit, dfTruth$sampleNameEdit)
  
  dfPredIntersect <- dfPred %>%
    filter(sampleNameEdit %in% c(intersect))
  
  dfTruthIntersect <- dfTruth %>%
    filter(sampleNameEdit %in% c(intersect))
  
  # Join and compute RMSE
  rmseDf <- dfPredIntersect %>%
    left_join(dfTruthIntersect, 
              by = c("sampleNameEdit", "manual.log.anot.fine")) %>%
    group_by(manual.log.anot.fine, 
             Method_Param_Reference) %>%
    summarise(RMSE = sqrt(mean(((cellProportion - truth) * 100)^2)), 
              .groups = "drop")
  
  return(rmseDf)
}

rmseByCellTypeMedian <- function(df, pred_type, col) {
  
  col <- sym(col)
  
  dfPred <- df %>%
    filter(sampleType %in% c(pred_type)) %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion, 
           !!col)
  
  dfTruth <- df %>%
    filter(sampleType == "sn") %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion) %>%
    dplyr::rename(truth = cellProportion)
  
  # 7 single-cell samples, 10 bulk cell and 10 bulk nuc samples. 
  intersect <- intersect(dfPred$sampleNameEdit, dfTruth$sampleNameEdit)
  
  dfPredIntersect <- dfPred %>%
    filter(sampleNameEdit %in% c(intersect))
  
  dfTruthIntersect <- dfTruth %>%
    filter(sampleNameEdit %in% c(intersect))
  
  # Join and compute RMSE
  rmseDf <- dfPredIntersect %>%
    left_join(dfTruthIntersect, 
              by = c("sampleNameEdit", "manual.log.anot.fine")) %>%
    group_by(manual.log.anot.fine, 
             !!col) %>%
    summarise(RMSE = sqrt(median(((cellProportion - truth) * 100)^2)), 
              .groups = "drop") 
  
  return(rmseDf)
}

rmseBySampleMean <- function(df, pred_type) {
  
  dfPred <- df %>%
    filter(sampleType %in% c(pred_type)) %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion, 
           Method_Param_Reference)
  
  dfTruth <- df %>%
    filter(sampleType == "sn") %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion) %>%
    dplyr::rename(truth = cellProportion)
  
  # 7 single-cell samples, 10 bulk cell and 10 bulk nuc samples. 
  intersect <- intersect(dfPred$sampleNameEdit, dfTruth$sampleNameEdit)
  
  dfPredIntersect <- dfPred %>%
    filter(sampleNameEdit %in% c(intersect))
  
  dfTruthIntersect <- dfTruth %>%
    filter(sampleNameEdit %in% c(intersect))
  
  # Join and compute RMSE
  rmseDf <- dfPredIntersect %>%
    left_join(dfTruthIntersect, 
              by = c("sampleNameEdit", "manual.log.anot.fine")) %>%
    group_by(sampleNameEdit, 
             Method_Param_Reference) %>%
    summarise(RMSE = sqrt(mean(((cellProportion - truth) * 100)^2)), 
              .groups = "drop") 
  
  return(rmseDf)
}

rmseBySampleMedian <- function(df, pred_type, col) {
  
  col <- sym(col)
  
  dfPred <- df %>%
    filter(sampleType %in% c(pred_type)) %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion, 
           !!col)
  
  dfTruth <- df %>%
    filter(sampleType == "sn") %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion) %>%
    dplyr::rename(truth = cellProportion)
  
  # 7 single-cell samples, 10 bulk cell and 10 bulk nuc samples. 
  intersect <- intersect(dfPred$sampleNameEdit, dfTruth$sampleNameEdit)
  
  dfPredIntersect <- dfPred %>%
    filter(sampleNameEdit %in% c(intersect))
  
  dfTruthIntersect <- dfTruth %>%
    filter(sampleNameEdit %in% c(intersect))
  
  # Join and compute RMSE
  rmseDf <- dfPredIntersect %>%
    left_join(dfTruthIntersect, 
              by = c("sampleNameEdit", "manual.log.anot.fine")) %>%
    group_by(sampleNameEdit, 
             !!col) %>%
    summarise(RMSE = sqrt(median(((cellProportion - truth) * 100)^2)), 
              .groups = "drop") 
  
  return(rmseDf)
}

rmseCombinedCellNucSampleMean <- function(df, pred_type) {
  
  dfPred <- df %>%
    # filter(sampleType %in% c(pred_type)) %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion, 
           Method_Param_Reference) %>%
    mutate(sampleName)
  
  dfTruth <- df %>%
    # filter(sampleType == "sn") %>%
    select(sampleNameEdit, 
           sampleName, 
           manual.log.anot.fine, 
           cellProportion) %>%
    dplyr::rename(truth = cellProportion)
  
  # 7 single-cell samples, 10 bulk cell and 10 bulk nuc samples. 
  intersect <- intersect(dfPred$sampleNameEdit, dfTruth$sampleNameEdit)
  
  dfPredIntersect <- dfPred %>%
    filter(sampleNameEdit %in% c(intersect))
  
  dfTruthIntersect <- dfTruth %>%
    filter(sampleNameEdit %in% c(intersect))
  
  # Join and compute RMSE
  rmseDf <- dfPredIntersect %>%
    left_join(dfTruthIntersect, 
              by = c("sampleNameEdit", "manual.log.anot.fine")) %>%
    group_by(sampleNameEdit, 
             Method_Param_Reference) %>%
    summarise(RMSE = sqrt(mean(((cellProportion - truth) * 100)^2)), 
              .groups = "drop") 
  
  return(rmseDf)
}

plot_metric_tile <- function(data, 
                             estimate_type, 
                             metric_name, 
                             fill_low, 
                             fill_high, 
                             title,
                             reverse = FALSE) {
  
  df <- data %>%
    filter(EstimateType == estimate_type, Metric == metric_name)
  
  # Define the color scale limits based on the metric name
  color_limits <- if (metric_name %in% c("SCorr", "CCorr")) {
    c(-1, 1)  # for SCorr or CCorr
  } else {
    c(0, 1)   # for MAE or MAECorr
  }
  
  # Create the plot
  # p <- ggplot(df, aes(x = Metric, y = Method_Param_Reference, fill = Value)) +
  #   geom_tile(color = "white") +
  #   scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  #   scale_x_discrete(
  #   expand = expansion(mult = c(0, 0)),
  #   labels = function(x) str_wrap(x, width = 10)
  # ) +
  #   (if (reverse) 
  #       scale_fill_gradient(low = fill_high, high = fill_low, na.value = "white", limits = color_limits)
  #    else 
  #       scale_fill_gradient(low = fill_low, high = fill_high, na.value = "white", limits = color_limits)) +
  #   theme_minimal(base_size = 14) +
  #   labs(x = NULL, y = NULL, fill = metric_name) +
  #   ggtitle(title) +
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     strip.text = element_text(size = 14),
  #     plot.title = element_text(hjust = 0.5),
  #     legend.title = element_blank(),
  #     axis.text.y = element_text(angle = 45, size = 8, hjust = 1)  # Rotate y-axis labels
  #   )
  
  # With text 
  p <- ggplot(df, aes(x = Metric, y = Method_Param_Reference, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Value, 2)), size = 3) +  # <-- Add this line
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(
      expand = expansion(mult = c(0, 0)),
      labels = function(x) str_wrap(x, width = 10)
    ) +
    (if (reverse) 
      scale_fill_gradient(low = fill_high, high = fill_low, na.value = "white", limits = color_limits)
     else 
       scale_fill_gradient(low = fill_low, high = fill_high, na.value = "white", limits = color_limits)) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = NULL, fill = metric_name) +
    ggtitle(title) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      axis.text.y = element_text(angle = 45, size = 8, hjust = 1)
    )
  
  
  return(p)
}

plot_metric_dot <- function(data, 
                            estimate_type, 
                            metric_name, 
                            color_low, 
                            color_high, 
                            title,
                            reverse = FALSE) {
  
  # Filter just once for plotting
  df <- data %>%
    filter(EstimateType == estimate_type, Metric == metric_name)
  
  # Use full data for global Stdev range
  global_stdev_range <- range(data$Stdev, na.rm = TRUE)
  global_min_stdev <- global_stdev_range[1]
  global_max_stdev <- global_stdev_range[2]
  
  # Define color limits based on metric
  color_limits <- if (metric_name %in% c("SCorr", "CCorr")) {
    c(-1, 1)
  } else {
    c(0, 1)
  }
  
  # Create dot plot
  p <- ggplot(df, aes(x = Metric, y = Method_Param_Reference)) +
    geom_point(aes(color = Mean, size = Stdev), alpha = 0.8) +
    geom_text(aes(label = round(Mean, 2)), hjust = 0, size = 3) +
    scale_color_gradient(
      low = if (reverse) color_high else color_low,
      high = if (reverse) color_low else color_high,
      limits = color_limits,
      na.value = "white"
    ) +
    scale_size(
      range = c(4, 10),
      limits = c(global_min_stdev, global_max_stdev),
      name = "Stdev"
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = NULL, color = metric_name) +
    ggtitle(title) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )
  
  return(p)
}

