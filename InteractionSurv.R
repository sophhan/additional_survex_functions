#' Function to compute PDP & ICE Cuves for 1 numerical feature
#'
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model The explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param grid_length  One dimensional vector determining the number of quantile values for the numerical feature
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether partial dependence and or individual conditional expectation plots should be marginalized over time
#' @param plot_type A string to specify the type of plot: ’pdp’ for partial dependence plot, ’ice’ for individual conditional expectation curves, ’pdp + ice’ for partial dependence plot and ice curves within the same plot (for 2 features only pd plots can be generated)
#' @param center_at A numerical value at which the plot should be centered

calculate_pdp_ice_num <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  grid_length = NULL,
                                  times,
                                  marginalize_over_time,
                                  plot_type,
                                  center_at = NULL) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 1)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = TRUE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))
  assert_numeric(data_dt[,feature_name, with = FALSE][[1]])
  assert_choice(plot_type, c("pdp", "pdp+ice", "ice"))
  assert_numeric(center_at, null.ok = TRUE, max.len = 1)

  # Create grid values of feature of interest to calculate ice or pd curves
  if (!is.null(grid_length)){
    grid_dt <- seq(from = min(data_dt[,feature_name, with = FALSE][[1]]),
                  to = max(data_dt[,feature_name, with = FALSE][[1]]),
                  length.out = grid_length)
  }
  else if (is.null(grid_length)){
    grid_dt <- sort(unique(data_dt[,feature_name, with = FALSE][[1]]))
  }

  # Obtain 'cartesian product' of grid values and all other feature values not equal to feature of interest
  # New blown up dataset is obtained of size grid_length*number_of_instances
  expanded_dt <- as.data.table(merge(grid_dt, data_dt[,-feature_name, with = FALSE]))
  expanded_dt <- setnames(expanded_dt, "x", feature_name)

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times
  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt)
  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt, predictions)
  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Plot PDP/ICE for multiple separate time points
  if ((length(times) > 1) && (marginalize_over_time == FALSE)){

    # Center in case centering value is provided
    if (!is.null(center_at)){
      # Cartesian product of center_at value and original data.frame
      expanded_grid_center <- as.data.table(merge(center_at, data_dt[,-feature_name, with = FALSE]))
      expanded_grid_center <- setnames(expanded_grid_center, "x", feature_name)

      # Predict centering values
      predictions_center = predict_function(model = model,
                                            newdata = expanded_grid_center,
                                            times = times)

      # Rename columns to prediction times
      colnames(predictions_center) <- times
      # Save feature names from original/expanded dataset for aggregation
      original_col_names <- colnames(expanded_grid_center)
      # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
      center_dt <- cbind(expanded_grid_center, predictions_center)
      # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
      center_dt <- data.table::melt(center_dt,
                                    id.vars = original_col_names,
                                    variable.name = "time",
                                    value.name = "centering_val")

      # Add predictions to original data.table and keep centering value and id and time columns
      center_dt <- center_dt[, .(centering_val, time, id)]

      # Merge predictions data.table with uncentered ICE data.table on time and id value
      merged_dt <- ice_dt[center_dt, on = c("id","time")]

      # Center ICE values
      merged_dt <- merged_dt[, centered_predictions := predictions - centering_val]

      # Replace uncentered ICE data.table with centered ICE data.table
      ice_dt <- merged_dt[, predictions := NULL]
      ice_dt <- ice_dt[, predictions := centered_predictions]
    }

    # Average over time x feature values from expanded dataset (deletes unnecessary columns remaining from centering)
    ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                     by = c(original_col_names, "time")]

    # Average predictions over identical grid values to obtain the pd values
    feature_name_sym <- sym(feature_name)
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c(feature_name, "time")]
    #pdp_dt <- setnames(pdp_dt, "feature_name_sym", feature_name)

    # Obtain upper and lower limits for y axis
    y_floor_pd <- floor(min(pdp_dt[,pd])*10)/10
    y_ceiling_pd <- ceiling(max(pdp_dt[,pd])*10)/10
    y_floor_ice <- floor(min(ice_dt[,predictions])*10)/10
    y_ceiling_ice <- ceiling(max(ice_dt[,predictions])*10)/10

    # ICE
    if (plot_type == "ice") {
      plot <- ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                geom_line(alpha = 0.2, mapping = aes(group = interaction(id, time), color = time)) +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_ice, y_ceiling_ice)
      output <- list(result_ice = ice_dt, plot = plot)
    }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      plot <- ggplot() +
                geom_line(data = ice_dt, aes(x = !!feature_name_sym, y = predictions, group = interaction(id, time), color = time), alpha = 0.1) +
                geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, color = time), linewidth = 1.5, lineend = "round", linejoin = "round") +
                geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, group = time), color = "black", linewidth = 0.5, linetype = "dashed", lineend = "round", linejoin = "round") +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides="b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_ice, y_ceiling_ice)
      output <- list(result_pd = pdp_dt, result_ice = ice_dt, plot = plot)
    }

    # PDP
    else if (plot_type == "pdp") {
      plot <- ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                geom_line(aes(color = time)) +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_pd, y_ceiling_pd)
      output <- list(result_pd = pdp_dt, plot = plot)
    }
  }

  # Plot PDP/ICE marginalized/averaged over time if multiple time points are given for survival prediction or for one specific time point
  else if ((length(times) == 1) || (marginalize_over_time == TRUE)){
    # Average/marginalize over time if multiple time values are given
    ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                     by = c(original_col_names)]

    if (!is.null(center_at)){
      # Cartesian product of center_at value and original data.frame
      expanded_grid_center <- as.data.table(merge(center_at, data_dt[,-feature_name, with = FALSE]))
      expanded_grid_center <- setnames(expanded_grid_center, "x", feature_name)

      # Predict centering values
      predictions_center = predict_function(model = model,
                                            newdata = expanded_grid_center,
                                            times = times)

      # Rename columns to prediction times
      colnames(predictions_center) <- times
      # Save feature names from original/expanded dataset for aggregation
      original_col_names <- colnames(expanded_grid_center)
      # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
      center_dt <- cbind(expanded_grid_center, predictions_center)
      # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
      center_dt <- data.table::melt(center_dt,
                                    id.vars = original_col_names,
                                    variable.name = "time",
                                    value.name = "predictions")

      # Average/marginalize over time if multiple time values are given
      center_dt <- center_dt[, .(centering_val = mean(predictions)),
                             by = c(original_col_names)]

      # Keep only centering value and id column
      center_dt <- center_dt[, .(centering_val, id)]

      # Merge predictions data.table with uncentered ICE data.table
      merged_dt <- ice_dt[center_dt, on = "id"]

      # Center ICE values
      merged_dt <- merged_dt[, centered_predictions := predictions - centering_val]

      # Replace uncentered ICE data.table with centered ICE data.table
      ice_dt <- merged_dt[, predictions := NULL]
      ice_dt <- ice_dt[, predictions := centered_predictions]
    }

    # Average predictions over identical grid values to obtain the pd values
    feature_name_sym <- sym(feature_name)
    pdp_dt <- ice_dt[, .(pd = mean(predictions)), by = c(feature_name)]
    #pdp_dt <- setnames(pdp_dt, "feature_name_sym", feature_name)

    # Obtain upper and lower limits for y axis
    y_floor_pd <- floor(min(pdp_dt[,pd])*10)/10
    y_ceiling_pd <- ceiling(max(pdp_dt[,pd])*10)/10

    # ICE
    if (plot_type == "ice") {
      plot <- ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                 geom_line(alpha = 0.2, mapping = aes(group = id)) +
                 geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter")
      output <- list(result_ice = ice_dt, plot = plot)
    }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      plot <- ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
                geom_line(mapping = aes(group = id), alpha = 0.2) +
                geom_line(data = pdp_dt, aes(x = !!feature_name_sym, y = pd), linewidth = 2, color = "gold") +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter")
      output <- list(result_pd = pdp_dt, result_ice = ice_dt, plot = plot)
    }

    # PDP
    else if (plot_type == "pdp") {
      plot <- ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
                geom_line() +
                geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter") +
                ylim(y_floor_pd, y_ceiling_pd)
      output <- list(result_pd = pdp_dt, plot = plot)
    }

  }
}

#' Function to compute PDP & ICE Curves for 1 categorical feature
#'
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model The explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether partial dependence and or individual conditional expectation plots should be marginalized over time
#' @param plot_type A string to specify the type of plot: ’pdp’ for partial dependence plot, ’ice’ for individual conditional expectation curves, ’pdp + ice’ for partial dependence plot and ice curves within the same plot (for 2 features only pd plots can be generated)

calculate_pdp_ice_cat <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  times,
                                  marginalize_over_time,
                                  plot_type) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 1)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))
  assert_factor(data_dt[,feature_name, with = FALSE][[1]])
  assert_choice(plot_type, c("pdp", "pdp+ice", "ice"))

  # Grid values are equivalent to different factor levels
  level_names <- levels(data_dt[,feature_name, with = FALSE][[1]])
  count_dt <- data.table(count = as.numeric(table(data_dt[,feature_name, with = FALSE][[1]])),
                         levels = level_names)
  setnames(count_dt, "levels", feature_name)

  # Obtain 'cartesian product' of grid values and all other feature values not equal to feature of interest
  # New blown up dataset is obtained of size grid_length*number_of_instances
  expanded_dt <- as.data.table(merge(level_names, data_dt[,-feature_name, with = FALSE]))
  expanded_dt <- setnames(expanded_dt, "x", feature_name)

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times
  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt)
  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt, predictions)
  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Plot PDP/ICE for multiple separate time points
  if ((length(times) > 1) && (marginalize_over_time == FALSE)){

    # Add factor level counts to data.table for adding them as plot labels later
    ice_dt <- merge(ice_dt, count_dt)
    # Create new column of type: level name (level count)
    ice_dt <- ice_dt[, graph_label := paste(get(feature_name)," (", count, ")", sep = "")]
    # Create column name: feature_name (count) & convert to symbol & rename column
    feature_name_count <- paste(feature_name, "(count)", sep = " ")
    feature_name_count_sym = sym(feature_name_count)
    ice_dt <- setnames(ice_dt, "graph_label", feature_name_count)

    # Average predictions over identical grid values to obtain the pd values
    feature_name_sym <- sym(feature_name)
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c(feature_name, feature_name_count, "time")]
    # Rename columns (bug from symbol evaluation)
    #pdp_dt <- setnames(pdp_dt, c("feature_name_sym","feature_name_count_sym"),
    #                   c(feature_name,feature_name_count))

    # ICE
    if (plot_type == "ice") {
      plot <- ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
                geom_boxplot(alpha = 0.2, mapping = aes(color = time))
      output <- list(result_ice = ice_dt, plot = plot)
    }
    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      plot <- ggplot(mapping = aes(color = time))  +
                geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions), alpha = 0.2) +
                geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = time), linewidth = 0.6)
      output <- list(result_ice = ice_dt, result_pd = pdp_dt, plot = plot)
    }
    # PDP
    else if (plot_type == "pdp") {
      plot <- ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, fill = time)) +
                geom_bar(stat = "identity", width = 0.5, position = "dodge")
      output <- list(result_pd = pdp_dt, plot = plot)
    }
  }

  # Plot PDP/ICE for individual time point or marginalized over multiple time points
  else if ((length(times) == 1) || (marginalize_over_time == TRUE)){

    # Average/marginalize over time if multiple time values are given
    ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                     by = c(original_col_names)]

    # Add factor level counts to data.table for adding them as plot labels later
    ice_dt <- merge(ice_dt, count_dt)
    # Create new column of type: level name (level count)
    ice_dt <- ice_dt[, graph_label := paste(get(feature_name)," (", count, ")", sep = "")]
    # Create column name: feature_name (count) & convert to symbol & rename column
    feature_name_count <- paste(feature_name, "(count)", sep = " ")
    feature_name_count_sym = sym(feature_name_count)
    ice_dt <- setnames(ice_dt, "graph_label", feature_name_count)

    # Average predictions over identical grid values to obtain the pd values
    feature_name_sym <- sym(feature_name)
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c(feature_name, feature_name_count)]
    # Rename columns (bug from symbol evaluation)
    #pdp_dt <- setnames(pdp_dt, c("feature_name_sym","feature_name_count_sym"),
    #                   c(feature_name,feature_name_count))

    # ICE
    if (plot_type == "ice") {
      plot <- ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
                geom_boxplot(alpha = 0.2)
      output <- list(result_ice = ice_dt, plot = plot)
      }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      plot <- ggplot()  +
                geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions), alpha = 0.2) +
                geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = 1), linewidth = 2, color = "gold")
      output <- list(result_ice = ice_dt, result_pd = pdp_dt, plot = plot)
    }

    # PDP
    else if (plot_type == "pdp") {
      plot <- ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd),) +
                geom_bar(stat = "identity", width = 0.5)
      output <- list(result_pd = pdp_dt, plot = plot)
    }
  }

}

#' Function to compute PDP & ICE Cuves for 2 numerical features
#'
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model The explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param grid_length  One dimensional vector determining the number of quantile values for the numerical feature
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param plot_type A string to specify the type of plot pnly ’pdp’ for partial dependence plot possible

calculate_pdp_num_num <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  grid_length = NULL,
                                  times,
                                  marginalize_over_time,
                                  plot_type) {

  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 2)
  assert_numeric(grid_length, any.missing = FALSE, len = 2, null.ok = TRUE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))
  assert_numeric(data_dt[,feature_name[1], with = FALSE][[1]])
  assert_numeric(data_dt[,feature_name[2], with = FALSE][[1]])
  assert_choice(plot_type, c("pdp"))

  # Obtain grid values and prediction data.table for partial dependence plots
  if (!is.null(grid_length)){
    # Define grid values for both numerical features of interest using equally spaced grid
    grid_dt1 <- seq(from = min(data_dt[,feature_name[1], with = FALSE][[1]]),
                    to = max(data_dt[,feature_name[1], with = FALSE][[1]]),
                    length.out = grid_length[1])
    grid_dt2 <- seq(from = min(data_dt[,feature_name[2], with = FALSE][[1]]),
                    to = max(data_dt[,feature_name[2], with = FALSE][[1]]),
                    length.out = grid_length[2])

    # Obtain 'cartesian product'/cross join of both sets of grid values and all other feature values not equal to features of interest
    # New blown up dataset is obtained of size grid_length[1]*grid_length[2]*number_of_instances
    # First obtain cross join of feature 1's grid and original data
    expanded_dt0 <- cross_join(grid_dt1, data_dt[,-feature_name, with = FALSE])
    # Secondly obtain cross join of feature 2's grid and previous cross join
    expanded_dt <- cross_join(grid_dt2, expanded_dt0)
  }

  # Obtain grid values and prediction data.table for partial dependence values needed for H statistic calculation
  else if (is.null(grid_length)){
    # Define grid values for both numerical features of interest using observed values of feature combination (feature1_i, feature2_i) i=1,...,n
    grid_dt <- data_dt[,feature_name, with = FALSE]
    # Obtain cross join of feature1 x feature2 grid data.table and original data
    expanded_dt <- cross_join(grid_dt, data_dt[,-feature_name, with = FALSE])
  }

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times

  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt, predictions)

  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt)

  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Calculate PD for multiple separate time points
  if ((length(times) > 1) && (marginalize_over_time == FALSE)){
    # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
    # Convert feature name string to symbols
    feature_name_sym1 <- sym(feature_name[1])
    feature_name_sym2 <- sym(feature_name[2])
    # Calculate pd values
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c("time", feature_name)]
    # Rename symbol columns to true feature names
    #pdp_dt <- setnames(pdp_dt, c("feature_name_sym1","feature_name_sym2"),
    #                   c(feature_name[1], feature_name[2]))
    # Return pd output without plot
    output <- list(result_pd = pdp_dt)
  }
  # Plot & Calculate 2D PDP for one time point or marginalized over time
  else if ((length(times) == 1) || (marginalize_over_time == TRUE)){
    # Average/marginalize over time if multiple time values are given
    ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                     by = c(original_col_names)]

    # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
    # Convert feature name string to symbols
    feature_name_sym1 <- sym(feature_name[1])
    feature_name_sym2 <- sym(feature_name[2])
    # Calculate pd values
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c(feature_name)]
    # Rename symbol columns to true feature names
    #pdp_dt <- setnames(pdp_dt, c("feature_name_sym1","feature_name_sym2"),
    #                   c(feature_name[1], feature_name[2]))

    # 2D PDP
    plot <- ggplot(pdp_dt, aes(x = !!feature_name_sym1, y = !!feature_name_sym2)) +
      geom_tile(aes(fill = pd)) +
      scale_fill_gradient(low = "#ffff33",
                          high = "#990000",
                          guide = "colorbar") +
      stat_contour(aes(z = pd), bins = 15, colour = "#000000") +
      geom_rug(data = data_dt,
               aes(x = !!feature_name_sym1, y = !!feature_name_sym2),
               alpha = 0.8,
               position = "jitter")

    # Return output as list of plot and numeric results
    output <- list(result_pd = pdp_dt, plot = plot)
  }
}

#' Function to compute PDP & ICE Cuves for 1 numerical and 1 categorical feature
#'
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model The explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param grid_length  One dimensional vector determining the number of quantile values for the numerical feature
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param plot_type A string to specify the type of plot pnly ’pdp’ for partial dependence plot possible

calculate_pdp_num_cat <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  grid_length = NULL,
                                  times,
                                  marginalize_over_time,
                                  plot_type) {

  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 2)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = TRUE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))
  assert_choice(plot_type, c("pdp"))

  # Determine which feature is categorical and which numeric
  num_index <- ifelse(inherits(data_dt[, feature_name[1], with = FALSE][[1]],
                               "numeric"), 1, 2)
  cat_index <- ifelse(num_index == 1,2,1)

  # Check numeric and categorical argument
  assert_numeric(data_dt[,feature_name[num_index], with = FALSE][[1]])
  assert_factor(data_dt[,feature_name[cat_index], with = FALSE][[1]])

  # Define grid values for both numerical features of interest
  # Create grid values of feature of interest to calculate ice or pd curves
  if (!is.null(grid_length)){
    grid_dt <- as.data.table(seq(from = min(data_dt[,feature_name[num_index], with = FALSE][[1]]),
                                 to = max(data_dt[,feature_name[num_index], with = FALSE][[1]]),
                                 length.out = grid_length))
    setnames(grid_dt, "V1", feature_name[num_index])
  }
  else if (is.null(grid_length)){
    grid_dt <- data_dt[,feature_name[num_index], with = FALSE]
  }
  level_names <- unique(data_dt[,feature_name[cat_index], with = FALSE])
  # Create vector with elements "level name (level count)" for later plotting
  count <- table(data_dt[, feature_name[cat_index], with = FALSE][[1]])
  name_count <- paste(level_names[[1]], " (", count, ")", sep = "")

  # Obtain 'cartesian product' of both sets of grid values and all other feature values not equal to features of interest
  # New blown up dataset is obtained of size grid_length[1]*grid_length[2]*number_of_instances
  # First obtain Cartesian product of feature 1's grid and original data
  expanded_dt0 <- cross_join(grid_dt, data_dt[,-feature_name, with = FALSE])

  # Secondly obtain Cartesian product of expanded dataset 1 and feature 2's grid
  expanded_dt <- cross_join(level_names, expanded_dt0)


  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times
  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt)
  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt, predictions)
  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Calculate PD for multiple separate time points
  if ((length(times) > 1) && (marginalize_over_time == FALSE)){
    # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
    # Convert feature name string to symbols
    feature_name_num_sym <- sym(feature_name[num_index])
    feature_name_cat_sym <- sym(feature_name[cat_index])
    # Calculate pd values
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = c("time", feature_name[num_index], feature_name[cat_index])]
    # Convert factor column to numeric for plotting
    #pdp_dt <- pdp_dt[, feature_name_cat_sym := as.numeric(as.factor(feature_name_cat_sym))]
    # Rename symbol columns to true feature names
    #pdp_dt <- setnames(pdp_dt, c("feature_name_num_sym","feature_name_cat_sym"),
    #                   c(feature_name[num_index], feature_name[cat_index]))
    # Return pd output without plot
    output <- list(result_pd = pdp_dt)
  }
  # Plot & Calculate 2D PDP for one time point or marginalized over time
  else if ((length(times) == 1) || (marginalize_over_time == TRUE)){
    # Average/marginalize over time if multiple time values are given
    # Average over time x feature values from expanded dataset
    ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                     by = c(original_col_names)]

    # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
    # Convert feature name string to symbols
    feature_name_num_sym <- sym(feature_name[num_index])
    feature_name_cat_sym <- sym(feature_name[cat_index])
    # Calculate pd values
    pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                     by = .(feature_name[num_index], feature_name[cat_index])]
    # Convert factor column to numeric for plotting
    pdp_dt[, (feature_name[cat_index]) := as.numeric(as.factor(get(feature_name[cat_index])))]
    # Rename symbol columns to true feature names
    #pdp_dt <- setnames(pdp_dt, c("feature_name_num_sym","feature_name_cat_sym"),
    #                   c(feature_name[num_index], feature_name[cat_index]))

    # 2D PDP
    plot <- ggplot(pdp_dt, aes(x = !!feature_name_cat_sym, y = !!feature_name_num_sym, z = pd)) +
      geom_tile(aes(fill = pd)) +
      scale_fill_gradient(low = "#ffff33",
                          high = "#990000",
                          guide = "colorbar") +
      stat_contour(colour = "#000000") +
      scale_x_discrete(limits = name_count)

    # Return output as list of plot and numeric results
    output <- list(result_pd = pdp_dt, plot = plot)
  }
}

#' Function to compute PDP for all but one specific feature j
#'
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model The explained model
#' @param feature_name The column name of the feature for which to compute PDP_(-j)
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function

calculate_pdp_no_j <- function(data_dt,
                               predict_function,
                               model,
                               feature_name,
                               times){
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 1)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))

  # Obtain vector of column names
  col_names <- names(data_dt)
  # Obtain vector of column names without feature of interest
  col_names_no_j <- setdiff(col_names, feature_name)

  # Define grid values for both numerical features of interest
  # Create grid values of feature of interest to calculate ice or pd curves
  grid_list <- lapply(col_names, function(col) unique(data_dt[, col, with = FALSE][[1]]))
  names(grid_list) <- col_names

  # Obtain Cartesian product
  expanded_dt <- as.data.table(do.call(expand.grid, grid_list))

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times
  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt)
  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt, predictions)
  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Calculate PD for multiple separate time points
  # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
  # Calculate pd values
  pdp_dt <- ice_dt[, .(pd = mean(predictions)),
                   by = c("time", col_names_no_j)]

  # Output PD values as data.frame
  output <- pdp_dt
}


#' Function to compute H_j(k) statistics
#'
#' @param explainer An explainer object - model preprocessed by the `survex::explain()` function
#' @param feature_name One dimensional vector of strings containing the names of feature j for which H_j(k) statistic should be computed and plotted, can only be factor or numerical variable
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether H_j(k) plots should be marginalized over time

surv_interaction <- function(explainer,
                             feature_name = NULL,
                             times = NULL,
                             sample_percentage = NULL,
                             marginalize_over_time = FALSE){

  #### Extract and adjust necessary elements from explainer
  # Extract data as data.table from explainer
  data_dt <- as.data.table(explainer$data)
  # Sample from data
  if (!is.null(sample_percentage)) {
    # Calculate the number of rows to sample
    sample_size <- round(nrow(data_dt) * (sample_percentage / 100))
    # Generate random indices to sample rows
    sample_indices <- sample(seq_len(nrow(data_dt)), size = sample_size)
    # Sample the rows from data.table
    data_dt <- data_dt[sample_indices]
    }
  # Add id column to data_dt
  data_dt_id <- copy(data_dt)
  data_dt_id <- data_dt_id[, id := seq_len(.N)]
  # Extract predict_function from explainer
  predict_function <- explainer$predict_survival_function
  # Extract model from explainer
  model <- explainer$model
  # Set default of times vector if necessary
  if (missing(times) && marginalize_over_time == FALSE) {
    times = median(explainer$y[,1])
  }
  else if (missing(times) && marginalize_over_time == TRUE) {
    times = explainer$times
  }

  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, len = 1, null.ok = TRUE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))

  #### Obtain lists of names of numeric and factor columns
  col_names <- colnames(data_dt)
  # List of names of all categorical/factor columns
  factor_cols <- col_names[sapply(data_dt, is.factor)]
  # List of names of all numeric columns
  numeric_cols <- setdiff(col_names, factor_cols)
  # Sorted list of all column names (first all numeric columns then all categorical columns)
  col_names_sort <- c(numeric_cols, factor_cols)

  ### Logic to obtain H_j(k) if feature j is given
  if (!is.null(feature_name)){
    # Argument check
    assert_true(feature_name %in% colnames(data_dt))

    #### Obtain list of 1D PD values
    # Use calculate_pdp_ice_num function to obtain PD for all numeric features (list of data.tables as outcome)
    num_pd_list <- lapply(numeric_cols, function(feature) {
      res <- calculate_pdp_ice_num(data_dt = data_dt,
                                 predict_function = predict_function,
                                 model = model,
                                 feature_name = feature,
                                 grid_length = NULL,
                                 marginalize_over_time = FALSE,
                                 plot_type = "pdp",
                                 times = times)
      res2 <- res$result_pd
      res2$feature <- feature
      colnames(res2)[which(colnames(res2) == feature)] <- "vals"
      res2
    })

    #### Logic to obtain 1D & 2D PD values needed for H-statistic calculation if feature j is numeric
    if (is.numeric(data_dt[, feature_name, with = FALSE][[1]])){
      ### Obtain 1D PD values needed for H-statistic calculation
      # Use calculate_pdp_ice_cat function to obtain PD for all categorical features (list of data.tables as outcome)
      cat_pd_list <- lapply(factor_cols, function(feature) {
        res <- calculate_pdp_ice_cat(data_dt = data_dt,
                                     predict_function = predict_function,
                                     model = model,
                                     feature_name = feature,
                                     marginalize_over_time = FALSE,
                                     plot_type = "pdp",
                                     times = times)
        res2 <- res$result_pd
        res2$feature <- feature
        colnames(res2)[which(colnames(res2) == feature)] <- "vals"
        cols_to_keep <- c("vals", "time", "pd", "feature")  # Specify the columns to keep
        res2 <- res2[, ..cols_to_keep, with = FALSE]
        res2
      })
      # Append lists of numeric and categorical feature 1D PD values
      oned_pd_list <- c(num_pd_list, cat_pd_list)
      # Rbind list of numeric and categorical feature 1D PD value (features k) data.tables to one data.table excluding feature j
      oned_pd_dt <- do.call(rbind, oned_pd_list[-which(col_names_sort == feature_name)])

      ### Create separate data.table for feature j
      feature_j_dt <- oned_pd_list[which(col_names_sort == feature_name)][[1]]
      colnames(feature_j_dt)[which(colnames(feature_j_dt) == "vals")] <- feature_name
      feature_j_dt <- feature_j_dt[, -"feature"]

      ### Obtain 2D PD values needed for H-statistic calculation
      # Remove feature j from list of names of all numeric columns
      numeric_cols_noj <- setdiff(numeric_cols, feature_name)
      # Create list of tuples (feature_j, feature_k) for k = 1,...,K for k only numeric features
      numeric_cols_comb <- lapply(numeric_cols_noj, function(x) c(feature_name, x))
      # Use calculate_pdp_num_num function to obtain PD for tuples (feature_j, feature_k)
      num_num_pd_list <- lapply(seq(numeric_cols_comb), function(i) {
        res <- calculate_pdp_num_num(data_dt = data_dt,
                                     predict_function = predict_function,
                                     model = model,
                                     feature_name = numeric_cols_comb[[i]],
                                     grid_length = NULL,
                                     marginalize_over_time = FALSE,
                                     plot_type = "pdp",
                                     times = times)
        res2 <- res$result_pd
        res2$feature <- numeric_cols_comb[[i]][2]
        colnames(res2)[which(colnames(res2) == numeric_cols_comb[[i]][2])] <- "vals"
        res2
      })
      # Create list of tuples (feature_j,feature_k) for k = 1,...,K for k only categorical features
      factor_cols_comb <- lapply(factor_cols, function(x) c(feature_name, x))
      ### Calculate 2D PD for tuples (feature_j, feature_k)
      # Use calculate_pdp_num_cat function to obtain PD for tuples (feature_j, feature_k)
      num_cat_pd_list <- lapply(seq(factor_cols_comb), function(i) {
        res <- calculate_pdp_num_cat(data_dt = data_dt,
                                     predict_function = predict_function,
                                     model = model,
                                     feature_name = factor_cols_comb[[i]],
                                     grid_length = NULL,
                                     marginalize_over_time = FALSE,
                                     plot_type = "pdp",
                                     times = times)
        res2 <- res$result_pd
        res2$feature <- factor_cols_comb[[i]][2]
        colnames(res2)[which(colnames(res2) == factor_cols_comb[[i]][2])] <- "vals"
        res2
      })
      # Append lists of 2D PD values
      twod_pd_list <- c(num_num_pd_list, num_cat_pd_list)
      # Rbind list of numeric and categorical feature 2D PD value data.tables to one data.table
      twod_pd_dt <- do.call(rbind, twod_pd_list)
    }

    #### Logic to obtain 1D & 2D PD values needed for H-statistic calculation if feature j is categorical
    if (is.factor(data_dt[,feature_name, with = FALSE][[1]])){
      ### Obtain 1D PD values needed for H-statistic calculation
      # Use calculate_pdp_ice_cat function to obtain 1D PD for feature j
      res <- calculate_pdp_ice_cat(data_dt = data_dt,
                                   predict_function = predict_function,
                                   model = model,
                                   feature_name = feature_name,
                                   marginalize_over_time = FALSE,
                                   plot_type = "pdp",
                                   times = times)
      res2 <- res$result_pd
      cols_to_keep <- c("time", "pd", feature_name)
      feature_j_dt <- res2[, ..cols_to_keep, with = FALSE]
      # Rbind list of numeric feature 1D PD value (features k) data.tables to one data.table
      oned_pd_dt <- do.call(rbind, num_pd_list)

      ### Obtain 2D PD values needed for H-statistic calculation
      # Create list of tuples (feature_j,feature_k) for k = 1,...,K for k only categorical features
      numeric_cols_comb <- lapply(numeric_cols, function(x) c(feature_name, x))
      ### Calculate 2D PD for tuples (feature_j, feature_k)
      # Use calculate_pdp_num_ice function to obtain PD for tuples (feature_j, feature_k)
      cat_num_pd_list <- lapply(seq(numeric_cols_comb), function(i) {
        res <- calculate_pdp_num_cat(data_dt = data_dt,
                                     predict_function = predict_function,
                                     model = model,
                                     feature_name = numeric_cols_comb[[i]],
                                     grid_length = NULL,
                                     marginalize_over_time = FALSE,
                                     plot_type = "pdp",
                                     times = times)
        res2 <- res$result_pd
        res2$feature <- numeric_cols_comb[[i]][2]
        colnames(res2)[which(colnames(res2) == numeric_cols_comb[[i]][2])] <- "vals"
        res2
      })
      # Rbind list of numeric and categorical feature 2D PD value data.tables to one data.table
      twod_pd_dt <- do.call(rbind, cat_num_pd_list)
    }

    #### Merge all PD data.tables into one final table for H_j(k) statistic calculation
    # Merge pd values with data_dt by feature values to obtain data.table with
    # corresponding PD value for feature_j per id and time combination
    merge_feature_j_dt <- merge(feature_j_dt,
                                data_dt_id[, c(feature_name, "id"), with = FALSE],
                                by = c(feature_name),
                                allow.cartesian = TRUE)

    # Melt data_dt into long format equivalent to format of oned_pd_dt,
    # keeping only features for which PD values are recorded in oned_pd_dt
    # melt has to be separated by numeric vs factor columns to avoid coercion warning
    cols_to_keep <- c(unique(oned_pd_dt$feature), "id")
    oned_data_dt <- data_dt_id[, .SD, .SDcols = cols_to_keep]
    # Delete invalid elements for melting from numeric and factor column lists
    if (is.numeric(data_dt[, feature_name, with = FALSE][[1]])) {
      cols_to_remove <- c("id", feature_name)
      numeric_cols <- numeric_cols[!(numeric_cols %in% cols_to_remove)]
    }
    else if (is.factor(data_dt[, feature_name, with = FALSE][[1]])){
      cols_to_remove <- c("id", feature_name)
      factor_cols <- factor_cols[factor_cols != feature_name]
      numeric_cols <- numeric_cols[numeric_cols != "id"]
    }
    # Melt numeric variables
    if (length(intersect(numeric_cols, colnames(oned_data_dt))) != 0) {
      melt_oned_numeric <- melt.data.table(oned_data_dt,
                                           id.vars = c("id"),
                                           measure.vars = numeric_cols,
                                           variable.name = "feature",
                                           value.name = "vals")
    } else {melt_oned_numeric <- data.table(empty = 0)}
    # Melt factor variables
    if (length(intersect(factor_cols, colnames(oned_data_dt))) != 0) {
      melt_oned_factor <- melt.data.table(oned_data_dt,
                                          id.vars = c("id"),
                                          measure.vars = factor_cols,
                                          variable.name = "feature",
                                          value.name = "vals")
    } else {melt_oned_factor <- data.table(empty = 0)}

    # Merge melted results
    if ((ncol(melt_oned_factor) == 1) & (ncol(melt_oned_numeric) == 3)){
      melt_oned_data_dt <- melt_oned_numeric
    } else if ((ncol(melt_oned_factor) == 3) & (ncol(melt_oned_numeric) == 1)){
      melt_oned_data_dt <- melt_oned_factor
    } else if ((ncol(melt_oned_factor) == 3) & (ncol(melt_oned_numeric) == 3)){
      melt_oned_data_dt <- rbind(melt_oned_numeric, melt_oned_factor)
    }

    # Merge 1D PD values with melted data_dt by feature values to obtain data.table with
    # corresponding PD value for feature_k, k = 1,...,K per id and time combination
    merge_oned_pd_dt <- merge(oned_pd_dt,
                              melt_oned_data_dt,
                              by = c("feature","vals"),
                              allow.cartesian = TRUE)

    # Melt data_dt into long format equivalent to format of twod_pd_dt,
    # keeping only features for which pd values are recorded in twod_pd_dt
    # melt has to be separated by numeric vs factor columns to avoid coercion warning
    cols_to_keep <- c(unique(twod_pd_dt$feature), feature_name, "id")
    twod_data_dt <- data_dt_id[, .SD, .SDcols = cols_to_keep]
    # Melt numeric variables
    if (length(intersect(numeric_cols, colnames(twod_data_dt))) != 0) {
      melt_twod_numeric <- melt.data.table(twod_data_dt,
                                           id.vars = c("id", feature_name),
                                           measure.vars = numeric_cols,
                                           variable.name = "feature",
                                           value.name = "vals")
    } else {melt_twod_numeric <- data.table(empty = 0, empty2 = 0)}

    # Melt factor variables
    if (length(intersect(factor_cols, colnames(twod_data_dt))) != 0) {
      melt_twod_factor <- melt.data.table(twod_data_dt,
                                          id.vars = c("id", feature_name),
                                          measure.vars = factor_cols,
                                          variable.name = "feature",
                                          value.name = "vals")
    } else {melt_twod_factor <- data.table(empty = 0, empty2 = 0)}

    # Merge melted results
    if ((ncol(melt_twod_factor) == 2) & (ncol(melt_twod_numeric) == 4)){
      melt_twod_data_dt <- melt_twod_numeric
    } else if ((ncol(melt_twod_factor) == 4) & (ncol(melt_twod_numeric) == 2)){
      melt_twod_data_dt <- melt_twod_factor
    } else if ((ncol(melt_twod_factor) == 4) & (ncol(melt_twod_numeric) == 4)){
      melt_twod_data_dt <- rbind(melt_twod_numeric,melt_twod_factor)
    }

    # Merge 2D PD values with melted data_dt by feature values to obtain data.table with
    # corresponding PD value for feature_k, k = 1,...,K per id and time combination
    merge_twod_data_dt <- merge(twod_pd_dt,
                                melt_twod_data_dt,
                                by = c(feature_name, "feature", "vals"),
                                allow.cartesian = TRUE)

    # Final merges to have 3 PD values (PD_j, PD_k and PD_jk) in one data.table pd_dt
    pd_dt0 <- merge(merge_oned_pd_dt, merge_feature_j_dt,
                    by = c("id","time"),
                    suffixes = c("_k", "_j"))
    pd_dt <- merge(pd_dt0, merge_twod_data_dt,
                   by = c("id","time",feature_name,"feature","vals"))

    # Mean-center PD values for H_j(k) statistic calculation
    cols_to_scale <- c("pd_k", "pd_j", "pd")
    for (col in cols_to_scale) {
      pd_dt[, (col) := scale(get(col), center = TRUE, scale = FALSE),
            by = .(feature, time)]
    }
  }

  ### Logic to obtain H_j if no feature j is given
  else if (is.null(feature_name)){
    ### Obtain list of PD_(-j) values for every feature
    # Use calculate_pdp_no_j function to obtain PD_(-j) for all features (list of data.tables as outcome)
    pd_no_j_list <- lapply(seq(col_names), function(i) {
      res <- calculate_pdp_no_j(data_dt = data_dt,
                                predict_function = predict_function,
                                model = model,
                                feature_name = col_names[i],
                                times = times)
      res
    })

    # Merge data.tables in d_no_j_list list with original data.table to include id and feature j and only keep columns relevant for matching with other PD data.tables
    pd_no_j_list2 <- lapply(seq(col_names), function(i) {
      # Merge data.tables in pd_no_j_list with data_dt_id by all features except feature j
      res_merge <- merge.data.table(pd_no_j_list[[i]],
                                    data_dt_id,
                                    by = c(col_names[-i]),
                                    allow.cartesian = TRUE)
      # Keep only relevant columns for later matching
      cols_to_keep <- c("id", "time", col_names[i], "pd")
      res_subset <- res_merge[, .SD, .SDcols = cols_to_keep]
      # Melt data.table from wide to long format to match format for later merging of PD values
      res_melt <- data.table::melt(res_subset,
                                   id.vars = c("id","time","pd"),
                                   variable.name = "feature",
                                   value.name = "vals")
      res_melt
    })
    # Rbind list of 1D PD_(-j) value data.tables to one data.table
    pd_no_j_dt <- do.call(rbind, pd_no_j_list2)

    ### Obtain data.table of f(x) prediction values for every feature
    # Make predictions on original data
    pd = predict_function(model = model,
                          newdata = data_dt,
                          times = times)
    # Rename predictions according to different time points
    colnames(pd) <- times
    # Convert to data.table
    pred_data_dt <- as.data.table(cbind(data_dt_id, pd))
    # Convert data.table from wide to long format for later merging
    pred_long_dt <- melt(pred_data_dt,
                         id.vars = c(col_names, "id"),
                         variable.name = "time",
                         value.name = "pd")
    # Delete unnecessary columns from prediction data.table
    cols_to_keep <- c("id", "time", "pd")
    pred_long_dt <- pred_long_dt[, ..cols_to_keep, with = FALSE]

    #### Obtain list of 1D PD values
    # Use calculate_pdp_ice_num function to obtain PD for all numeric features (list of data.tables as outcome)
    num_pd_j_list <- lapply(numeric_cols, function(feature) {
      res <- calculate_pdp_ice_num(data_dt = data_dt,
                                   predict_function = predict_function,
                                   model = model,
                                   feature_name = feature,
                                   grid_length = NULL,
                                   marginalize_over_time = FALSE,
                                   plot_type = "pdp",
                                   times = times)
      res2 <- res$result_pd
      res2$feature <- feature
      colnames(res2)[which(colnames(res2) == feature)] <- "vals"
      res2
    })

    # Use calculate_pdp_ice_cat function to obtain PD for all numeric features (list of data.tables as outcome)
    cat_pd_j_list <- lapply(factor_cols, function(feature) {
      res <- calculate_pdp_ice_cat(data_dt = data_dt,
                                   predict_function = predict_function,
                                   model = model,
                                   feature_name = feature,
                                   marginalize_over_time = FALSE,
                                   plot_type = "pdp",
                                   times = times)
      res2 <- res$result_pd
      res2$feature <- feature
      colnames(res2)[which(colnames(res2) == feature)] <- "vals"
      cols_to_keep <- c("vals", "time", "pd", "feature")
      res2 <- res2[, ..cols_to_keep, with = FALSE]
      res2
    })

    # Append lists of 1D PD values
    pd_j_list <- c(num_pd_j_list, cat_pd_j_list)
    # Rbind list of numeric and categorical feature 1D PD value data.tables to one data.table
    pd_j_dt <- do.call(rbind, pd_j_list)

    #### Merge all PD data.tables into one final table for H_j statistic calculation
    # Convert original data into long format for merging
    # Melt numeric variables
    melt_numeric <- melt.data.table(data_dt_id,
                                    id.vars = c("id"),
                                    measure.vars = numeric_cols,
                                    variable.name = "feature",
                                    value.name = "vals")

    # Melt factor variables
    melt_factor <- melt.data.table(data_dt_id,
                                   id.vars = c("id"),
                                   measure.vars = factor_cols,
                                   variable.name = "feature",
                                   value.name = "vals")
    # Create 1D PD value data.table
    if ((ncol(melt_factor) == 1) & (ncol(melt_numeric) == 3)){
      melt_data_dt_id <- melt_numeric
    }
    else if ((ncol(melt_factor) == 3) & (ncol(melt_numeric) == 1)){
      melt_data_dt_id <- melt_factor
    }
    else if ((ncol(melt_factor) == 3) & (ncol(melt_numeric) == 3)){
      melt_data_dt_id <- rbind(melt_numeric, melt_factor)
    }
    # Merge original data in long format with with 1D PD values for each feature j
    merge_pd_j_dt <- merge(pd_j_dt,
                           melt_data_dt_id,
                           by = c("feature", "vals"),
                           allow.cartesian = TRUE)
    # Merge 1D PD values with prediction data.table
    pd_dt0 <- merge(merge_pd_j_dt,
                    pred_long_dt,
                    by = c("id", "time"),
                    suffixes = c("_j", ""),
                    allow.cartesian=TRUE)
    # Create final PD data.table by merging with 1D PD_(-j) values for
    pd_dt <- merge(pd_dt0,
                   pd_no_j_dt,
                   by = c("id", "time", "feature", "vals"),
                   suffixes = c("", "_k"),
                   allow.cartesian=TRUE)

    # Mean-center PD values for H_j statistic calculation
    cols_to_scale <- c("pd_j", "pd_k", "pd")
    for (col in cols_to_scale) {
      pd_dt[, (col) := scale(get(col), center = TRUE, scale = FALSE),
            by = .(feature, time)]
    }
  }

  #### Logic to obtain result_df and plot of H_j(k)/H_j statistic for separate time values (not marginalized over time)
  if (marginalize_over_time == FALSE){
    # Compute H_j(k)-statistic for every feature_k at each individual time point
    H_stat <- pd_dt[, .(H = sum((pd - pd_j - pd_k)^2)/sum(pd^2)),
                  by = .(feature, time)]
    # Convert data.table to data.frame for plotting
    H_stat_df <- as.data.frame(H_stat)
    # Convert time column to numeric ensuring that levels are time values
    H_stat_df$time <- as.numeric(levels(H_stat$time))[H_stat$time]
    # Create separate data.frame to create the times rug
    rug_times_df <- data.frame(time = explainer$y[,1])
    rug_times_df <- subset(rug_times_df,
                           time >= min(times) & time <= max(times))

    # Plot H_j(k)-statistic for each feature k over time
    plot <- ggplot(data = H_stat_df) +
      geom_rug(data = rug_times_df, aes(x = time), sides = "b") +
      geom_line(aes(x = time, y = H, color = feature, group = feature))

    # Output as list of result_df and plot
    output <- list(result_H_stat = H_stat_df, plot = plot)
  }

  #### Logic to obtain result_df and plot of H_j(k)/H_j statistic marginalized over time
  else if (marginalize_over_time == TRUE){
    # Compute H_j(k)-statistic for every feature k marginalized over time
    H_stat <- pd_dt[, .(H = sum((pd - pd_j - pd_k)^2)/sum(pd^2)),
                  by = .(feature)]
    # Sort the data.table by H_j(k)
    sorted_H_stat <- H_stat[order(H, decreasing = FALSE)]
    # Convert data.table to data.frame for plotting
    H_stat_df <- as.data.frame(sorted_H_stat)
    H_stat_df$feature <- factor(H_stat_df$feature, levels = unique(H_stat_df$feature))

    # Plot H_j(k)-statistic
    plot <- ggplot(H_stat_df, aes(y = feature, x = H)) +
      geom_point() +
      geom_segment(aes(yend = feature, x = 0, xend = H))

    # Output as list of result_df and plot
    output <- list(result_H_stat = H_stat_df, plot = plot)
  }
}


# Quick Test surv_interaction function

library(R6)
library(survival)
library(survex)
library(ranger)
library(data.table)
library(ggplot2)
library(checkmate)

set.seed(123)

# Create example ranger
rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran,
                             respect.unordered.factors = TRUE, num.trees = 100,
                             mtry = 3, max.depth = 5)
# Create explainer
rsf_ranger_exp <- survex::explain(rsf_ranger, data = veteran[, -c(3, 4)],
                                  y = Surv(veteran$time, veteran$status))


# H_jk statistic plot over time for "trt"
inter_trt <- surv_interaction(rsf_ranger_exp,
                                  feature_name = "trt",
                                  times = seq(1,400, by = 10),
                                  sample_percentage = NULL,
                                  marginalize_over_time = FALSE)

inter_trt$plot

# H_jk statistic plot over time for "diagtime"
inter_diagtime <- surv_interaction(rsf_ranger_exp,
                              feature_name = "diagtime",
                              times = 1:600,
                              sample_percentage = NULL,
                              marginalize_over_time = FALSE)

inter_diagtime$plot

# H_jk statistic plot over time for "celltype"
inter_celltype <- surv_interaction(rsf_ranger_exp,
                                   feature_name = "celltype",
                                   times = seq(1,400, by = 10),
                                   sample_percentage = NULL,
                                   marginalize_over_time = FALSE)

inter_celltype$plot

# H_jk statistic plot marginalized over time for "celltype"
inter_celltype <- surv_interaction(rsf_ranger_exp,
                                   feature_name = "celltype",
                                   sample_percentage = NULL,
                                   marginalize_over_time = TRUE)

inter_celltype$plot

# H_jk statistic plot marginalized over time for "trt"
inter_trt <- surv_interaction(rsf_ranger_exp,
                                   feature_name = "trt",
                                   sample_percentage = NULL,
                                   marginalize_over_time = TRUE)

inter_trt$plot

# H_j statistic plot over time
inter_total <- surv_interaction(rsf_ranger_exp,
                                times = 1:650,
                                sample_percentage = 50,
                                marginalize_over_time = FALSE)

inter_total$plot

# H_j statistic plot marginalized over time
inter_total <- surv_interaction(rsf_ranger_exp,
                                times = seq(1,400, by = 10),
                                sample_percentage = 50,
                                marginalize_over_time = TRUE)

inter_total$plot



