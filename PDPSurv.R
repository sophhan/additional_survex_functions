library(R6)
library(survival)
library(survex)
library(ranger)
library(data.table)
library(ggplot2)
library(checkmate)


#' Function to compute PD and ICE Plots
#'
#' @param explainer An explainer object - model preprocessed by the `survex::explain()` function
#' @param feature_name One or two dimensional vector of strings containing the names of the features for which partial dependence or accumulated local effects should be plotted, can only be factor or numerical variable, for two dimensional plot either two numerical or one factor and one numerical variable
#' @param grid_length  One or two dimensional vector determining the number of quantile values for numerical features
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether partial dependence and or individual conditional expectation plots should be marginalized over time
#' @param plot_type A string to specify the type of plot: ’pdp’ for partial dependence plot, ’ice’ for individual conditional expectation curves, ’pdp + ice’ for partial dependence plot and ice curves within the same plot (for 2 features only pd plots can be generated)
#' @param center_at A numerical value at which the plot should be centered, ignored in case of 2 features or one categorical feature

surv_pdp <- function(explainer,
                     feature_name,
                     grid_length = 100,
                     times = NULL,
                     marginalize_over_time = FALSE,
                     plot_type = "pdp+ice",
                     center_at = NULL){
  # Argument check
  assert_logical(marginalize_over_time)

  # Extract data as data.table from explainer
  data_dt <- as.data.table(explainer$data)
  # Add id column to data
  data_dt <- data_dt[, id := 1:nrow(data_dt)]
  # Extract predict_function from explainer
  predict_function <- explainer$predict_survival_function
  # Extract model from explainer
  model <- explainer$model

  # Set default of times vector if necessary
  # Check if marginalization setting matches length of times vector
  if (missing(times) && marginalize_over_time == FALSE) {
    times = median(explainer$y[,1])
  }
  else if (missing(times) && marginalize_over_time == TRUE) {
    times = explainer$times
  }
  else if (length(feature_name) > 1 && length(times) > 1 && marginalize_over_time == FALSE) {
    warning("Two dimensional PD Plots require marginalization over time for multiple provided time values. Marginalization over time will be performed.")
    marginalize_over_time == TRUE
  }

  # Determine type of PD/ICE Plot and corresponding function to use based on number of feature names provided and their data type
  # One-dimensional PD Plot for numerical feature
  if ((length(feature_name) == 1) && (inherits(data_dt[, feature_name, with = FALSE][[1]], "numeric"))) {
    calculate_pdp_ice_num(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          grid_length = grid_length,
                          times = times,
                          marginalize_over_time = marginalize_over_time,
                          plot_type = plot_type,
                          center_at = center_at)
  }

  # One-dimensional PD Plot for categorical feature
  else if ((length(feature_name) == 1) && (inherits(data_dt[, feature_name, with = FALSE][[1]], "factor"))) {
    calculate_pdp_ice_cat(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          times = times,
                          marginalize_over_time = marginalize_over_time,
                          plot_type = plot_type)
  }

  # Two-dimensional PD Plot for two numerical features
  else if ((length(feature_name) == 2) && (inherits(data_dt[, feature_name[1], with = FALSE][[1]], "numeric")) &&
           (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "numeric"))) {
    # Check adequacy of grid_length values provided
    if ((length(grid_length) == 1) && grid_length == 100) {
      # Set new two-dimensional default grid_length values
      grid_length = c(100,100)
    }
    else if ((length(grid_length) == 1) && (grid_length != 100)) {
      # Ensure that two grid_length values are set, if non-default grid_length value was chosen
      stop("Only one non-default grid_length value is chosen. Two are required for two-dimensional PD Plots with two numerical features.")
    }
    calculate_pdp_num_num(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          grid_length = grid_length,
                          plot_type = plot_type,
                          times = times,
                          marginalize_over_time == marginalize_over_time)
  }

  # Two-dimensional PD Plot for one numerical and one categorical feature
  else if ((length(feature_name) == 2) &&
           ((inherits(data_dt[, feature_name[1], with = FALSE][[1]], "factor")) & (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "numeric"))) |
           ((inherits(data_dt[, feature_name[1], with = FALSE][[1]], "numeric")) & (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "factor")))) {

    if (length(grid_length) >= 2) {
      stop("Multiple grid_length values are chosen. For a two-dimensional ALE Plot of one numerical and one categorical feature only one grid_length value is required for the numerical feature.")
    }

    calculate_pdp_num_cat(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          grid_length = grid_length,
                          plot_type = plot_type,
                          times = times,
                          marginalize_over_time == marginalize_over_time)
  }
  # Warning message if more than 2 features are provided
  else if (length(feature_name) > 2) {
    stop("More than two feature names are provided. Function can handle two-dimensional PD Plots at maximum.")
  }
  else {
    stop("An error occured choosing the correct PD plotting function. Check whether appropriate feature names and types provided.")
  }

}

#' Function to compute PDP & ICE Cuves for 1 categorical feature
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
                     by = .(eval(feature_name_sym), eval(feature_name_count_sym), time)]
    # Rename columns (bug from symbol evaluation)
    pdp_dt <- setnames(pdp_dt, c("feature_name_sym","feature_name_count_sym"),
                       c(feature_name,feature_name_count))

    # ICE
    if (plot_type == "ice") {
    ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
      geom_boxplot(alpha = 0.2, mapping = aes(color = time))
    }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
    ggplot(mapping = aes(color = time))  +
      geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions), alpha = 0.2) +
      geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = time), linewidth = 0.6)
    }
    # PDP
    else if (plot_type == "pdp") {
    ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, fill = time)) +
      geom_bar(stat = "identity", width = 0.5, position = "dodge")
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
                     by = .(eval(feature_name_sym), eval(feature_name_count_sym))]
    # Rename columns (bug from symbol evaluation)
    pdp_dt <- setnames(pdp_dt, c("feature_name_sym","feature_name_count_sym"),
                       c(feature_name,feature_name_count))

    # ICE
    if (plot_type == "ice") {
      ggplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions)) +
        geom_boxplot(alpha = 0.2)}

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      ggplot()  +
        geom_boxplot(data = ice_dt, aes(x = !!feature_name_count_sym, y = predictions), alpha = 0.2) +
        geom_line(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd, group = 1), linewidth = 2, color = "gold")
    }

    # PDP
    else if (plot_type == "pdp") {
      ggplot(data = pdp_dt, aes(x = !!feature_name_count_sym, y = pd),) +
        geom_bar(stat = "identity", width = 0.5)
    }
  }

}


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
                                  grid_length,
                                  times,
                                  marginalize_over_time,
                                  plot_type,
                                  center_at = NULL) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 1)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))
  assert_numeric(data_dt[,feature_name, with = FALSE][[1]])
  assert_choice(plot_type, c("pdp", "pdp+ice", "ice"))
  assert_numeric(center_at, null.ok = TRUE, max.len = 1)

  # Create grid values of feature of interest to calculate ice or pd curves
  grid_dt <- seq(from = min(data_dt[,feature_name, with = FALSE][[1]]),
                 to = max(data_dt[,feature_name, with = FALSE][[1]]),
                 length.out = grid_length)

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
                     by = .(eval(feature_name_sym), time)]
    pdp_dt <- setnames(pdp_dt, "feature_name_sym", feature_name)

    # Obtain upper and lower limits for y axis
    y_floor_pd <- floor(min(pdp_dt[,pd])*10)/10
    y_ceiling_pd <- ceiling(max(pdp_dt[,pd])*10)/10
    y_floor_ice <- floor(min(ice_dt[,predictions])*10)/10
    y_ceiling_ice <- ceiling(max(ice_dt[,predictions])*10)/10

    # ICE
    if (plot_type == "ice") {
      ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
        geom_line(alpha = 0.2, mapping = aes(group = interaction(id, time), color = time)) +
        geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides = "b", alpha = 0.8, position = "jitter") +
        ylim(y_floor_ice, y_ceiling_ice)
    }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
     ggplot() +
       geom_line(data = ice_dt, aes(x = !!feature_name_sym, y = predictions, group = interaction(id, time), color = time), alpha = 0.1) +
       geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, color = time), linewidth = 1.5, lineend = "round", linejoin = "round") +
       geom_path(data = pdp_dt, aes(x = !!feature_name_sym, y = pd, group = time), color = "black", linewidth = 0.5, linetype = "dashed", lineend = "round", linejoin = "round") +
       geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_ice), sides="b", alpha = 0.8, position = "jitter") +
       ylim(y_floor_ice, y_ceiling_ice)
    }

    # PDP
    else if (plot_type == "pdp") {
    ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
      geom_line(aes(color = time)) +
      geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter") +
      ylim(y_floor_pd, y_ceiling_pd)
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
    pdp_dt <- ice_dt[, .(pd = mean(predictions)), by = .(eval(feature_name_sym))]
    pdp_dt <- setnames(pdp_dt, "feature_name_sym", feature_name)

    # Obtain upper and lower limits for y axis
    y_floor_pd <- floor(min(pdp_dt[,pd])*10)/10
    y_ceiling_pd <- ceiling(max(pdp_dt[,pd])*10)/10

    # ICE
    if (plot_type == "ice") {
      ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
        geom_line(alpha = 0.2, mapping = aes(group = id)) +
        geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides = "b", alpha = 0.8, position = "jitter")
      }

    # PDP + ICE
    else if (plot_type == "pdp+ice") {
      ggplot(data = ice_dt, aes(x = !!feature_name_sym, y = predictions)) +
        geom_line(mapping = aes(group = id), alpha = 0.2) +
        geom_line(data = pdp_dt, aes(x = !!feature_name_sym, y = pd), linewidth = 2, color = "gold") +
        geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter")
      }

    # PDP
    else if (plot_type == "pdp") {
      ggplot(data = pdp_dt, aes(x = !!feature_name_sym, y = pd)) +
        geom_line() +
        geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_ceiling_pd), sides="b", alpha = 0.8, position = "jitter") +
        ylim(y_floor_pd, y_ceiling_pd)
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
                                  grid_length,
                                  times,
                                  marginalize_over_time,
                                  plot_type) {

  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 2)
  assert_numeric(grid_length, any.missing = FALSE, len = 2, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))
  assert_numeric(data_dt[,feature_name[1], with = FALSE][[1]])
  assert_numeric(data_dt[,feature_name[2], with = FALSE][[1]])
  assert_choice(plot_type, c("pdp"))
  assert_true(marginalize_over_time)

  # Define grid values for both numerical features of interest
  grid_dt1 <- seq(from = min(data_dt[,feature_name[1], with = FALSE][[1]]),
                  to = max(data_dt[,feature_name[1], with = FALSE][[1]]),
                  length.out = grid_length[1])
  grid_dt2 <- seq(from = min(data_dt[,feature_name[2], with = FALSE][[1]]),
                  to = max(data_dt[,feature_name[2], with = FALSE][[1]]),
                  length.out = grid_length[2])

  # Obtain 'cartesian product' of both sets of grid values and all other feature values not equal to features of interest
  # New blown up dataset is obtained of size grid_length[1]*grid_length[2]*number_of_instances
  # First obtain Cartesian product of feature 1's grid and original data
  expanded_dt <- merge(grid_dt1, data_dt[,-feature_name, with = FALSE])
  # Rename automatically named column to feature_name
  names(expanded_dt)[names(expanded_dt) == "x"] <- feature_name[1]
  # Secondly obtain Cartesian product of expanded dataset 1 and feature 2's grid
  expanded_dt2 <- as.data.table(merge(grid_dt2, expanded_dt))
  # Rename automatically named column to feature_name
  expanded_dt2 <- setnames(expanded_dt2, "x", feature_name[2])

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt2,
                                 times = times)

  # Marginalize/Average over time if multiple time points are given for survival prediction
  predictions = rowMeans(predictions)

  # Add predictions to expanded dataset
  ice_dt <- expanded_dt2[, predictions := predictions]

  # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
  # Convert feature name string to symbols
  feature_name_sym1 <- sym(feature_name[1])
  feature_name_sym2 <- sym(feature_name[2])
  # Calculate pd values
  pdp_dt <- ice_dt[, .(pd = mean(predictions)), by = .(eval(feature_name_sym1), eval(feature_name_sym2))]
  # Rename symbol columns to true feature names
  pdp_dt <- setnames(pdp_dt, c("feature_name_sym1","feature_name_sym2"), c(feature_name[1], feature_name[2]))

  # 2D PDP
  ggplot(pdp_dt, aes(x = !!feature_name_sym1, y = !!feature_name_sym2)) +
    geom_tile(aes(fill = pd)) +
    scale_fill_gradient(low = "#ffff33",
                        high = "#990000",
                        guide = "colorbar") +
    stat_contour(aes(z = pd), bins = 15, colour = "#000000") +
    geom_rug(data = data_dt,
             aes(x = !!feature_name_sym1, y = !!feature_name_sym2),
             alpha = 0.8,
             position = "jitter")
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
                                  grid_length,
                                  times,
                                  marginalize_over_time,
                                  plot_type) {

  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 2)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))
  assert_choice(plot_type, c("pdp"))
  assert_true(marginalize_over_time)

  # Determine which feature is categorical and which numeric
  num_index <- ifelse(inherits(data_dt[, feature_name[1], with = FALSE][[1]],
                               "numeric"), 1, 2)
  cat_index <- ifelse(num_index == 1,2,1)

  # Check numeric and categorical argument
  assert_numeric(data_dt[,feature_name[num_index], with = FALSE][[1]])
  assert_factor(data_dt[,feature_name[cat_index], with = FALSE][[1]])

  # Define grid values for both numerical features of interest
  grid_dt1 <- seq(from = min(data_dt[,feature_name[num_index], with = FALSE][[1]]),
                  to = max(data_dt[,feature_name[num_index], with = FALSE][[1]]),
                  length.out = grid_length)
  level_names <- levels(data_dt[,feature_name[cat_index], with = FALSE][[1]])
  # Create vector with elements "level name (level count)" for later plotting
  count <- table(data_dt[,feature_name[cat_index], with = FALSE][[1]])
  name_count <- paste(level_names, " (", count, ")", sep = "")

  # Obtain 'cartesian product' of both sets of grid values and all other feature values not equal to features of interest
  # New blown up dataset is obtained of size grid_length[1]*grid_length[2]*number_of_instances
  # First obtain Cartesian product of feature 1's grid and original data
  expanded_dt <- merge(grid_dt1, data_dt[,-feature_name, with = FALSE])
  # Rename automatically named column to feature_name
  names(expanded_dt)[names(expanded_dt) == "x"] <- feature_name[num_index]
  # Secondly obtain Cartesian product of expanded dataset 1 and feature 2's grid
  expanded_dt2 <- as.data.table(merge(level_names, expanded_dt))
  # Rename automatically named column to feature_name
  expanded_dt2 <- setnames(expanded_dt2, "x", feature_name[cat_index])

  # Obtain predictions for expanded dataset
  predictions = predict_function(model = model,
                                 newdata = expanded_dt2,
                                 times = times)

  # Rename columns to prediction times
  colnames(predictions) <- times
  # Save feature names from original/expanded dataset for aggregation
  original_col_names <- colnames(expanded_dt2)
  # Combine prediction dataset with expanded dataset to create dataset in wide format with features on the LHS and predictions at different time points on the RHS
  ice_dt <- cbind(expanded_dt2, predictions)
  # Transform wide dataset to long format (one column for predictions, one column to indicate prediction time point)
  ice_dt <- data.table::melt(ice_dt,
                             id.vars = original_col_names,
                             variable.name = "time",
                             value.name = "predictions")

  # Average over time x feature values from expanded dataset
  ice_dt <- ice_dt[, .(predictions = mean(predictions)),
                   by = c(original_col_names)]

  # Average predictions over identical grid value pairs (feature 1 x feature 2) to obtain the pd values
  # Convert feature name string to symbols
  feature_name_num_sym <- sym(feature_name[num_index])
  feature_name_cat_sym <- sym(feature_name[cat_index])
  # Calculate pd values
  pdp_dt <- ice_dt[, .(pd = mean(predictions)), by = .(eval(feature_name_num_sym),
                                                       eval(feature_name_cat_sym))]
  # Convert factor column to numeric for plotting
  pdp_dt <- pdp_dt[, feature_name_cat_sym := as.numeric(as.factor(feature_name_cat_sym))]
  # Rename symbol columns to true feature names
  pdp_dt <- setnames(pdp_dt, c("feature_name_num_sym","feature_name_cat_sym"),
                     c(feature_name[num_index], feature_name[cat_index]))

  # 2D PDP
  ggplot(pdp_dt, aes(x = !!feature_name_cat_sym, y = !!feature_name_num_sym, z = pd)) +
    geom_tile(aes(fill = pd)) +
    scale_fill_gradient(low = "#ffff33",
                        high = "#990000",
                        guide = "colorbar") +
    stat_contour(colour = "#000000") +
    scale_x_discrete(limits = name_count)
}


#' Demonstration
set.seed(123)

veteran <- survival::veteran
veteran$prior <- as.factor(veteran$prior)
levels(veteran$prior) <- c("no","yes")
veteran$trt <- as.factor(veteran$trt)
levels(veteran$trt) <- c("standard","test")
rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran,
                             respect.unordered.factors = TRUE, num.trees = 100,
                             mtry = 3, max.depth = 5)
rsf_ranger_exp <- survex::explain(rsf_ranger, data = veteran[, -c(3, 4)],
                                  y = Surv(veteran$time, veteran$status))

# PDP+ICE Plot numerical feature for median(time)
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp+ice")

# PDP+ICE Plots numerical feature for multiple time points without marginalization over time
theme_set(theme_minimal())
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "prior",
         grid_length = 100,
         times = c(1,10,100,200,300,400,500,600,700,800),
         marginalize_over_time = FALSE,
         plot_type = "pdp+ice")

surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "prior",
         grid_length = 100,
         times = 1:850,
         marginalize_over_time = FALSE,
         plot_type = "pdp") +
  theme(legend.position = "none")

surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "prior",
         grid_length = 100,
         marginalize_over_time = TRUE,
         plot_type = "pdp+ice")

# PDP+ICE Plots numerical feature for multiple time points with marginalization over time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp+ice",
         marginalize_over_time = TRUE)

########## PDP for numerical feature with marginalization over time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp",
         marginalize_over_time = TRUE)

# PDP for numerical feature without marginalization over time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp",
         times = c(1,10,100,200,300),
         marginalize_over_time = FALSE)

# ICE for numerical feature without marginalization over time for median(t)
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "ice")

# ICE for numerical feature without marginalization over time for multiple time points
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         times = c(1,10,100,200,300),
         grid_length = 100,
         plot_type = "ice")

# PDP+ICE Plot with centering for numerical feature for median(t)
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp+ice",
         center_at = 1)

# PDP+ICE Plot with centering for numerical feature for multiple time points without marginalization over time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         times = c(1,10,100,200,300),
         plot_type = "pdp+ice",
         center_at = 1)

# PDP+ICE Plot with centering for numerical feature for individual time point median(t)
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "pdp+ice",
         center_at = 1)

# PDP Plot with centering for numerical feature for multiple time points
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         times = c(1,10,100,200,300),
         plot_type = "pdp",
         center_at = 1)

# ICE Plot with centering for numerical feature for one time point with centering
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         plot_type = "ice",
         center_at = 1)

# ICE Plot with centering for numerical feature for multiple time points with centering
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "diagtime",
         grid_length = 100,
         times = c(1,10,100,200,300),
         plot_type = "ice",
         center_at = 1)

# PDP + ICE Plot for categorical feature without marginalization over time for median time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "celltype",
         plot_type = "pdp+ice")

# PDP + ICE Plot for categorical feature without marginalization over time for multiple time values
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "celltype",
         time = c(1,10,100,300,500),
         marginalize_over_time = FALSE,
         plot_type = "pdp+ice")


# PDP Plot for categorical feature without marginalization over time for median time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "celltype",
         time = c(1,10,100,300,500),
         marginalize_over_time = FALSE,
         plot_type = "pdp")

# ICE Plot for categorical feature without marginalization over time for median time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "celltype",
         plot_type = "ice")

# ICE Plot for categorical feature without marginalization over time for median time
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = "celltype",
         time = c(1,10,100,300,500),
         marginalize_over_time = FALSE,
         plot_type = "ice")

# 2D PDP+ICE plot
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = c("celltype","diagtime"),
         grid_length = 100,
         plot_type = "pdp")

# 2D PDP+ICE plot
surv_pdp(explainer = rsf_ranger_exp,
         feature_name = c("age","diagtime"),
         grid_length = c(100,100),
         plot_type = "pdp")

