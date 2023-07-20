#' Function to compute ALE Plots
#'
#' @param explainer An explainer object - model preprocessed by the `survex::explain()` function
#' @param feature_name One or two dimensional vector of strings containing the names of the features for which accumulated local effects should be plotted, can only be factor or numerical variable, for two dimensional plot either two numerical or one factor and one numerical variable
#' @param grid_length  One or two dimensional vector determining the number of quantile values for numerical features
#' @param times A vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether accumulated local effects plot should be marginalized over time


surv_ale <- function(explainer,
                     feature_name,
                     grid_length = 100,
                     times = NULL,
                     marginalize_over_time = FALSE)
{
  # Argument check
  assert_logical(marginalize_over_time)

  # Extract data as data.table from explainer
  data_dt <- as.data.table(explainer$data)
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
    warning("Two dimensional ALE Plots require marginalization over time for multiple provided time values. Marginalization over time will be performed.")
    marginalize_over_time == TRUE
  }

  # Determine type of ALE Plot and corresponding function to use based on number of feature names provided and their data type
  # One-dimensional ALE Plot for numerical feature
  if ((length(feature_name) == 1) && (inherits(data_dt[, feature_name, with = FALSE][[1]], "numeric"))) {
    calculate_ale_num(data_dt = data_dt,
                      predict_function = predict_function,
                      model = model,
                      feature_name = feature_name,
                      grid_length = grid_length,
                      times = times,
                      marginalize_over_time = marginalize_over_time)
  }
  # One-dimensional ALE Plot for categorical feature
  else if ((length(feature_name) == 1) && (inherits(data_dt[, feature_name, with = FALSE][[1]], "factor"))) {
    calculate_ale_cat(data_dt = data_dt,
                      predict_function = predict_function,
                      model = model,
                      feature_name = feature_name,
                      times = times,
                      marginalize_over_time = marginalize_over_time)
  }
  # Two-dimensional ALE Plot for two numerical features
  else if ((length(feature_name) == 2) && (inherits(data_dt[, feature_name[1], with = FALSE][[1]], "numeric")) &&
           (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "numeric"))) {
    # Check adequacy of grid_length values provided
    if (grid_length == 100) {
      # Set new two-dimensional default grid_length values
      grid_length = c(100,100)
    }
    else if ((length(grid_length) == 1) && (grid_length != 100)) {
      # Ensure that two grid_length values are set, if non-default grid_length value was chosen
      stop("Only one non-default grid_length value is chosen. Two are required for two-dimensional ALE Plots with two numerical features.")
    }
    calculate_ale_num_num(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          grid_length = grid_length,
                          times = times)
  }
  # Two-dimensional ALE Plot for one numerical and one categorical feature
  else if ((length(feature_name) == 2) &&
           ((inherits(data_dt[, feature_name[1], with = FALSE][[1]], "factor")) & (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "numeric"))) |
           ((inherits(data_dt[, feature_name[1], with = FALSE][[1]], "numeric")) & (inherits(data_dt[, feature_name[2], with = FALSE][[1]], "factor")))) {

    if (length(grid_length) >= 2) {
      stop("Multiple grid_length values are chosen. For a two-dimensional ALE Plot of one numerical and one categorical feature only one grid_length value is required for the numerical feature.")
    }

    calculate_ale_num_cat(data_dt = data_dt,
                          predict_function = predict_function,
                          model = model,
                          feature_name = feature_name,
                          grid_length = grid_length,
                          times = times)
  }
  # Warning message if more than 2 features are provided
  else if (length(feature_name) > 2) {
    stop("More than two feature names are provided. Function can handle two-dimensional ALE Plots at maximum.")
  }
  else {
    stop("An error occured choosing the correct ALE plotting function. Check whether appropriate feature names and types provided.")
  }

}


#' Function to compute ALE for 1 numerical feature
#'
#' @importFrom data.table melt as.data.table setnames
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model the explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param grid_length Numerical value determining the number of quantile values for numerical feature
#' @param times a vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether accumulated local effects plot should be marginalized over time
calculate_ale_num <- function(data_dt,
                              predict_function,
                              model,
                              feature_name,
                              grid_length,
                              times,
                              marginalize_over_time) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, unique = TRUE, len = 1)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))

  # Determine feature index/indices based on feature names provided
  feature_index <- match(feature_name, colnames(data_dt))

  # Number of quantile points for determined by grid length
  quantile_vals <- as.numeric(quantile(data_dt[, ..feature_index][[1]],
                                       seq(1/100, 1, length.out = grid_length),
                                       type = 1))
  # Quantile points vector of length grid_length + 1
  quantile_vec <- c(min(data_dt[, ..feature_index][[1]]), quantile_vals)
  # Convert quantile points to data.table, make sure quantile points are unique
  grid_dt <- data.table(quantile_vec = unique(quantile_vec))
  # Match feature instances to quantile intervals
  interval_index <- findInterval(data_dt[[feature_name]], grid_dt[[1]], left.open = TRUE)

  # Points in interval 0 should be in interval 1
  interval_index[interval_index == 0] <- 1
  # Prepare datasets with upper and lower interval limits replacing original feature values
  X_lower <- X_upper <- data_dt
  X_lower[, feature_name] <- grid_dt[interval_index,]
  X_upper[, feature_name] <- grid_dt[interval_index + 1,]
  # Get survival predictions for instances of upper and lower interval limits
  predictions_lower = predict_function(model = model,
                                       newdata = X_lower,
                                       times = times)
  predictions_upper = predict_function(model = model,
                                       newdata = X_upper,
                                       times = times)
  # First order finite differences
  prediction_deltas <- predictions_upper - predictions_lower
  # Rename columns to timepoints for which predictions were made
  colnames(prediction_deltas) <- times
  # Add predictions to data.table with corresponding feature values and interval values
  deltas <- cbind(X_lower[, feature_name, with = FALSE],
                  prediction_deltas,
                  data.frame("interval" = interval_index))
  # Reshape from long to wide format in case there predictions for multiple time points
  deltas <- data.table::melt(deltas,
                             id.vars = setdiff(colnames(deltas), times),
                             variable.name = "time",
                             value.name = "prediction_surv")

  if ((marginalize_over_time == TRUE) | length(times) == 1){
    # Average/marginalize over time if multiple time values are given
    deltas <- deltas[, .(prediction_surv = mean(prediction_surv)),
                     by = setdiff(colnames(deltas), c("time", "prediction_surv"))]
    # Reshape from long to wide format in case there are multiple columns of predictions types
    y_hat_names <- as.data.table(setdiff(colnames(deltas), c(colnames(data_dt), "interval")))
    y_hat_names = unlist(y_hat_names)
    deltas <- data.table::melt(deltas,
                               variable.name = "class",
                               value.name = "yhat_diff_surv",
                               measure.vars = y_hat_names)

    # Average over instances within each prediction type and interval
    setkeyv(deltas, c("class", "interval"))
    deltas <- deltas[,
                     list(yhat_diff_surv = mean(yhat_diff_surv)),
                     by = c("interval", "class")]
    # Accumulate local effects over intervals
    deltas_accumulated <- deltas[,
                                 list(y_hat_cumsum = cumsum_na(c(0, yhat_diff_surv))),
                                 by = "class"]

    # The mean effect used for centering is the weighted mean of the interval mid
    # point effects, weighted by the number of points in the interval
    interval_n <- as.numeric(table(interval_index))
    ale_means <- deltas_accumulated[,
                                    list(ale0 = sum(((y_hat_cumsum[1:(nrow(.SD) - 1)] +
                                                        y_hat_cumsum[2:nrow(.SD)]) / 2) *
                                                      interval_n) / sum(interval_n)),
                                    by = "class"
    ]
    # Centering the ALEs to obtain final ALE values
    deltas_accumulated <- merge(deltas_accumulated,
                                ale_means,
                                all.x = TRUE,
                                by = "class")
    ale_values <- deltas_accumulated[, list(ale = y_hat_cumsum - ale0,
                                            id = 1:nrow(.SD)),
                                     by = "class"
    ]
    # Add quantile points corresponding to ALE values to data.table
    grid_dt$id <- 1:nrow(grid_dt)
    ale_values <- merge(ale_values,
                        grid_dt,
                        by = "id")

    # Convert data.table to data.frame for plotting
    ale_values = data.frame(ale_values)
    colnames(ale_values)[which(colnames(ale_values) == "quantile_vec")] <- feature_name

    # ALE Plot
    # ylim values
    y_floor <- floor(min(ale_values[,"ale"])*100)/100
    y_ceiling <- ceiling(max(ale_values[,"ale"])*100)/100
    feature_name_sym <- sym(feature_name)

    ggplot() +
      geom_line(data = ale_values,
                aes(x = !!feature_name_sym, y = ale)) +
      geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_floor),
               sides = "b",
               alpha = 0.8,
               position = "jitter") +
      ylim(y_floor, y_ceiling)
  }

  else if ((marginalize_over_time == FALSE) & (length(times) > 1)){
    # Reshape from long to wide format in case there are multiple columns of predictions types
    y_hat_names <- as.data.table(setdiff(colnames(deltas), c(colnames(data_dt), "interval", "time")))
    y_hat_names = unlist(y_hat_names)
    deltas <- data.table::melt(deltas,
                               variable.name = "class",
                               value.name = "yhat_diff_surv",
                               measure.vars = y_hat_names)

    # Average over instances within each prediction type and interval
    setkeyv(deltas, c("class", "interval","time"))
    deltas <- deltas[,
                     list(yhat_diff_surv = mean(yhat_diff_surv)),
                     by = c("interval", "class", "time")]
    # Accumulate local effects over intervals
    deltas_accumulated <- deltas[,
                                 list(y_hat_cumsum = cumsum_na(c(0, yhat_diff_surv))),
                                 by = c("class","time")]

    # The mean effect used for centering is the weighted mean of the interval mid
    # point effects, weighted by the number of points in the interval
    interval_n <- as.numeric(table(interval_index))
    ale_means <- deltas_accumulated[,
                                    list(ale0 = sum(((y_hat_cumsum[1:(nrow(.SD) - 1)] +
                                                        y_hat_cumsum[2:nrow(.SD)]) / 2) *
                                                      interval_n) / sum(interval_n)),
                                    by = c("class","time")
    ]
    # Centering the ALEs to obtain final ALE values
    deltas_accumulated <- merge(deltas_accumulated,
                                ale_means,
                                all.x = TRUE,
                                by = c("class","time"))
    ale_values <- deltas_accumulated[, list(ale = y_hat_cumsum - ale0,
                                            id = 1:nrow(.SD)),
                                     by = c("class","time")
    ]
    # Add quantile points corresponding to ALE values to data.table
    grid_dt$id <- 1:nrow(grid_dt)
    ale_values <- merge(ale_values,
                        grid_dt,
                        by = "id")

    # Convert data.table to data.frame for plotting
    ale_values = data.frame(ale_values)
    colnames(ale_values)[which(colnames(ale_values) == "quantile_vec")] <- feature_name

    # ALE Plot
    # ylim values
    y_floor <- floor(min(ale_values[, "ale"])*100)/100
    y_ceiling <- ceiling(max(ale_values[, "ale"])*100)/100
    feature_name_sym <- sym(feature_name)

    ggplot() +
      geom_line(data = ale_values, alpha = 0.8, mapping = aes(x = !!feature_name_sym, y = ale, color = time)) +
      geom_rug(data = data_dt, aes(x = !!feature_name_sym, y = y_floor),
               sides = "b",
               alpha = 0.8,
               position = "jitter") +
      ylim(y_floor, y_ceiling)
  }
}

#' Compute ALE for 1 categorical feature
#' @importFrom data.table melt as.data.table setnames
#' @param data_dt data.table object with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model the explained model
#' @param feature_name The column name of the feature for which to compute ALE
#' @param times a vector of times, that are used for evaluation of survival function and cumulative hazard function
#' @param marginalize_over_time A logical vector that specifies whether accumulated local effects plot should be marginalized over time
calculate_ale_cat <- function(data_dt,
                              predict_function,
                              model,
                              feature_name,
                              times,
                              marginalize_over_time) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, len = 1, unique = TRUE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name %in% colnames(data_dt))

  # Get original level order and number of levels
  x_cat <- data_dt[, feature_name, with = FALSE][[1]]
  levels_original <- levels(droplevels(x_cat))
  levels_n <- nlevels(droplevels(x_cat))

  # Reorder levels according to distance to other features,
  # if levels are already ordered use preset order
  if (inherits(x_cat, "ordered")) {
    level_order <- 1:levels_n
  } else {
    level_order <- order_levels(data_dt, feature_name)
  }

  # The new order of the levels
  levels_ordered <- levels_original[level_order]
  # The feature with the levels in the new order
  x_ordered <- order(level_order)[as.numeric(droplevels(x_cat))]
  X_lower <- X_upper <- data_dt

  # Filter rows which are not already at maximum or minimum level values
  row_ind_increase <- (1:nrow(data_dt))[x_ordered < levels_n]
  row_ind_decrease <- (1:nrow(data_dt))[x_ordered > 1]

  # Increase / decrease the levels where not at minimum or maximum already
  X_lower[row_ind_decrease, feature_name] <-
    levels_ordered[x_ordered[row_ind_decrease] - 1]
  X_upper[row_ind_increase, feature_name] <-
    levels_ordered[x_ordered[row_ind_increase] + 1]

  # Make predictions for original levels
  y_hat <- predict_function(model = model,
                            newdata = data_dt,
                            times = times)
  # Make predictions for increased levels (excluding maximum levels)
  y_hat_increase <- predict_function(model = model,
                                     newdata = X_upper[row_ind_increase, ],
                                     times = times)
  # Make predictions for decreased levels (excluding minimum levels)
  y_hat_decrease <- predict_function(model = model,
                                     newdata = X_lower[row_ind_decrease, ],
                                     times = times)

  # To measure the difference or 'jump' in prediction between categories,
  # we come from both directions: We measure the difference in prediction
  # when we increase category k to k+1 for instances with category k
  # We also measure the difference in prediction when we decrease category
  # k+1 to k for instance with category k+1
  d_increase <- y_hat_increase - y_hat[row_ind_increase, ]
  d_decrease <- y_hat[row_ind_decrease, ] - y_hat_decrease
  # Combine increases and decreases rowwise
  deltas <- rbind(d_increase, d_decrease)
  # Rename columns to prediction time points
  colnames(deltas) <- times
  # Add level jump identifier variable
  deltas <- as.data.table(cbind(deltas,
                                "level_jump" = c(x_ordered[row_ind_increase],
                                                 x_ordered[row_ind_decrease] - 1)))
  # Reshape from long to wide format in case there are predictions for multiple time points
  deltas <- data.table::melt(deltas,
                             id.vars = setdiff(colnames(deltas), times),
                             variable.name = "time",
                             value.name = "prediction_surv")

  if ((marginalize_over_time == TRUE) | (length(times) == 1)){
    # Average/marginalize over time if multiple time values are given
    deltas <- deltas[, .(prediction_surv = mean(prediction_surv)),
                     by = setdiff(colnames(deltas),c("prediction_surv", "time"))]
    # Transform data.table from wide to long format, only relevant if there are
    # multiple types of predictions
    deltas <- data.table(data.table::melt(deltas,
                                          id.vars = c("level_jump"),
                                          value.name = "yhat_diff_surv",
                                          variable.name = "class"
    ))
    # All those difference are aggregated (averaged)
    # grouped by the jump between levels (and by class for multi output)
    setkeyv(deltas, c("class", "level_jump"))
    deltas <- deltas[, list(yhat_diff_surv = mean(yhat_diff_surv)),
                     by = c("class", "level_jump")]
    # We accumulate the jumps, starting at 0 for the first category
    deltas <- deltas[, list(ale = c(0, cumsum(yhat_diff_surv))), by = c("class")]
    # Extract frequency and proportion of instances in each factor level
    x_count <- as.numeric(table(x_cat))
    x_prob <- x_count / sum(x_count)
    # Center results by the mean interval effect,
    # weighted by number of instances in the interval
    deltas <- deltas[, list(ale = ale - sum(ale * x_prob[level_order]),
                            level = levels_ordered,
                            count = x_count[level_order]),
                     by = "class"]
    colnames(deltas) <- c("class", "ale", feature_name, "count")
    # Ensure that feature of interest is of type factor
    deltas[, feature_name] <- factor(deltas[, feature_name, with = FALSE][[1]],
                                     levels = levels_ordered)

    # Order rows by the new level order and convert to data.frame for plotting
    ale_values <- data.frame(deltas[order(deltas[[feature_name]]), ])
    # Add column for x axis labels of type: level name (count)
    vals_name_count <- paste(ale_values[,feature_name], " (", ale_values[,"count"], ")", sep = "")
    feature_name_count <- paste(feature_name, " (count)", sep = "")
    ale_values[,"col_name"] <- vals_name_count
    colnames(ale_values)[which(colnames(ale_values) == "col_name")] <- feature_name_count
    # Convert column name to symbol
    feature_name_count_sym <- sym(feature_name_count)

    # ALE Barplot
    ggplot(data = ale_values,
           aes(x = !!feature_name_count_sym, y = ale),) +
      geom_bar(stat = "identity", width = 0.5)
  }

  else if ((marginalize_over_time == FALSE) & (length(times) > 1)){
    # Transform data.table from wide to long format, only relevant if there are
    # multiple types of predictions
    deltas <- data.table(data.table::melt(deltas,
                                          id.vars = c("level_jump", "time"),
                                          value.name = "yhat_diff_surv",
                                          variable.name = c("class")
    ))
    # All those difference are aggregated (averaged)
    # grouped by the jump between levels (and by .class for multi output)
    setkeyv(deltas, c("class", "level_jump", "time"))
    deltas <- deltas[, list(yhat_diff_surv = mean(yhat_diff_surv)),
                     by = c("class", "level_jump","time")]
    # We accumulate the jumps, starting at 0 for the first category
    deltas <- deltas[, list(ale = c(0, cumsum(yhat_diff_surv))),
                     by = c("class","time")]
    # Extract frequency and proportion of instances in each factor level
    x_count <- as.numeric(table(x_cat))
    x_prob <- x_count / sum(x_count)
    # Center results by the mean interval effect,
    # weighted by number of instances in the interval
    deltas <- deltas[, list(ale = ale - sum(ale * x_prob[level_order]),
                            level = levels_ordered,
                            count = x_count[level_order]),
                     by = c("class", "time")]
    colnames(deltas) <- c("class", "time", "ale", feature_name, "count")
    # Convert feature to factor
    deltas[, feature_name] <- factor(deltas[, feature_name, with = FALSE][[1]],
                                     levels = levels_ordered)

    # Order rows by the new level order and convert to data.frame for plotting
    ale_values <- data.frame(deltas[order(deltas[[feature_name]]), ])
    # Add column for x axis labels of type: level name (count)
    vals_name_count <- paste(ale_values[,feature_name], " (", ale_values[,"count"], ")", sep = "")
    feature_name_count <- paste(feature_name, " (count)", sep = "")
    ale_values[,"col_name"] <- vals_name_count
    colnames(ale_values)[which(colnames(ale_values) == "col_name")] <- feature_name_count
    # Convert column name to symbol
    feature_name_count_sym <- sym(feature_name_count)

    # ALE Barplot
    ggplot() +
      geom_bar(data = ale_values, aes(x = !!feature_name_count_sym, y = ale, fill = time),
               stat = "identity", width = 0.5, position = "dodge")
  }
}


#' Compute ALE for 2 numerical features
#'
#' @param data_dt the data.frame with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model the explained model
#' @param feature_name The column names of the two numerical features for which to compute ALE as vector
#' @param grid_length Two numerical values determining the number of quantile values for each numerical feature
#' @param times a vector of times, that are used for evaluation of survival function and cumulative hazard function
calculate_ale_num_num <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  grid_length,
                                  times) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, len = 2, unique = TRUE)
  assert_numeric(grid_length, any.missing = FALSE, len = 2, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))

  # Determine feature index/indices based on feature names provided
  feature_index <- match(feature_name, colnames(data_dt))

  # Number of quantile points for determined by grid length
  quantile_vals1 <- as.numeric(quantile(data_dt[, feature_name[1], with = FALSE][[1]],
                                        seq(1/100, 1, length.out = grid_length[1]),
                                        type = 1))
  quantile_vals2 <- as.numeric(quantile(data_dt[, feature_name[2], with = FALSE][[1]],
                                        seq(1/100, 1, length.out = grid_length[2]),
                                        type = 1))
  # Quantile points vector of length grid_length + 1
  quantile_vec1 <- c(min(data_dt[, feature_name[1], with = FALSE][[1]]),
                     quantile_vals1)
  quantile_vec2 = c(min(data_dt[, feature_name[2], with = FALSE][[1]]),
                    quantile_vals2)
  # Convert to data.table, if variable is discrete quantile points may be non-unique
  grid_dt1 <- as.data.table(unique(quantile_vec1))
  grid_dt2 <- as.data.table(unique(quantile_vec2))

  # Remove data outside of boundaries
  data_dt <- data_dt[(data_dt[[feature_name[1]]] <= max(grid_dt1[[1]])) &
                       (data_dt[[feature_name[1]]] >= min(grid_dt1[[1]])) &
                       (data_dt[[feature_name[2]]] <= max(grid_dt2[[1]])) &
                       (data_dt[[feature_name[2]]] >= min(grid_dt2[[1]])),]
  # Matching instances to the grids of both features
  interval_index1 <- findInterval(data_dt[[feature_name[1]]],
                                  grid_dt1[[1]],
                                  left.open = TRUE)
  interval_index2 <- findInterval(data_dt[[feature_name[2]]],
                                  grid_dt2[[1]],
                                  left.open = TRUE)
  # Data point in the left most interval should be in interval 1, not zero
  interval_index1[interval_index1 == 0] <- 1
  interval_index2[interval_index2 == 0] <- 1

  # Preparation for predictions of cell corners
  X_low1_low2 <- X_up1_low2 <- X_low1_up2 <- X_up1_up2 <- copy(data_dt)
  X_low1_low2[, feature_name] <- data.table(grid_dt1[interval_index1, ],
                                            grid_dt2[interval_index2, ])
  X_up1_low2[, feature_name] <- data.table(grid_dt1[interval_index1 + 1, ],
                                           grid_dt2[interval_index2, ])
  X_low1_up2[, feature_name] <- data.table(grid_dt1[interval_index1, ],
                                           grid_dt2[interval_index2 + 1, ])
  X_up1_up2[, feature_name] <- data.table(grid_dt1[interval_index1 + 1, ],
                                          grid_dt2[interval_index2 + 1, ])

  # Getting all predictions from the model
  y_hat_11 <- predict_function(model = model,
                               newdata = X_low1_low2,
                               times = times)
  y_hat_21 <- predict_function(model = rsf_ranger_exp$model,
                               newdata = X_up1_low2,
                               times = times)
  y_hat_12 <- predict_function(model = model,
                               newdata = X_low1_up2,
                               times = times)
  y_hat_22 <- predict_function(model = rsf_ranger_exp$model,
                               newdata = X_up1_up2,
                               times = times)
  # Second order differences per cell and instance
  # (top right corner - top left corner) - (bottom right corner  - bottom left corner)
  y_hat <- (y_hat_22 - y_hat_21) - (y_hat_12 - y_hat_11)
  # Rename columns to corresponding prediction times
  colnames(y_hat) <- times
  # Add corresponding quantiles to each delta value
  deltas <- cbind(data_dt[, feature_name, with = FALSE],
                  y_hat,
                  data.frame("interval1" = interval_index1,
                             "interval2" = interval_index2))
  # Melt data.table from wide to long format, with time time column and
  # rows corresponding to y_hat values for different times
  deltas <- data.table::melt(deltas,
                             variable.name = "time",
                             value.name = "y_hat",
                             measure.vars = colnames(y_hat)
  )

  # Average/marginalize over time if multiple time values are given
  deltas <- deltas[, .(y_hat = mean(y_hat)),
                   by = setdiff(colnames(deltas), c("time", "y_hat"))]
  # Determine names of multi-dimensional predictions
  y_hat_names <- setdiff(colnames(deltas),
                         c(colnames(data_dt),
                           c("interval1", "interval2")))

  # instead of a matrix, we work with melted version of the data
  # This allows us to work with multi-dimensional predictions
  deltas <- data.table::melt(deltas,
                             variable.name = "class",
                             value.name = "y_hat",
                             measure.vars = y_hat_names
  )
  # Create data.frame of all possible quantile combinations,
  # to make sure all intervals are included
  interval_grid <- expand.grid(
    interval1 = unique(deltas$interval1),
    interval2 = unique(deltas$interval2),
    class = unique(deltas$class)
  )
  # Merge expanded interval grid with delta values computed (creates many NA values)
  deltas <- merge(deltas,
                  interval_grid,
                  on = c("interval1", "interval2", "class"),
                  all.y = TRUE)
  # Average over all instances within a cell (cell defined by interval & prediction class)
  setkeyv(deltas, c("class", "interval1", "interval2"))
  deltas <- deltas[, list(yhat_diff = mean(y_hat)),
                   by = c("interval1", "interval2", "class")]

  # Indicate NA status of deltas in "missing" column
  deltas_na_cell <- copy(deltas)
  deltas_na_cell$missing <- is.na(deltas_na_cell$yhat_diff)
  # Delete column of original deltas
  deltas_na_cell$yhat_diff <- NULL
  # replace the missing ones with the closest non-missing difference (measured in number of intervals)
  deltas <- rbindlist(future.apply::future_lapply(unique(deltas$class), function(cl) {
    impute_cells(deltas[class == cl, ],
                 grid1 = grid_dt1,
                 grid2 = grid_dt2,
                 x1.ind = "interval1",
                 x2.ind = "interval2"
    )
  }))

  # Accumulate the predictions from bottom left to top right
  # Accumulate over interval of feature 1
  ale <- deltas[, list(y_hat_cumsum = cumsum_na(c(0, yhat_diff)),
                       interval2 = c(0, interval2)),
                by = c("class", "interval1")]
  # Accumulate over interval of feature 2
  ale <- ale[, list(y_hat_cumsum = cumsum_na(c(0, y_hat_cumsum)),
                    interval1 = c(0, interval1)),
             by = c("class", "interval2")]

  # Number of cells (cell = interval combination)
  cell_counts <- as.matrix(table(interval_index1, interval_index2))
  # Data.table with cell counts for every interval combination
  cell_counts_m <- data.table::melt(as.data.table(cell_counts),
                                    measure.vars = "N")
  # Remove variable column
  cell_counts_m$variable <- NULL
  # Rename value column to .count
  colnames(cell_counts_m) <- c("interval1", "interval2", "count")
  # Convert character interval values to numeric
  cell_counts_m$interval1 <- as.numeric(as.character(cell_counts_m$interval1))
  cell_counts_m$interval2 <- as.numeric(as.character(cell_counts_m$interval2))
  # Add cell counts to ale data.table
  ale <- merge(ale, cell_counts_m, on = c("interval1", "interval2"), all.x = TRUE)

  # Computing the first-order effect of feature 1
  # First, take the differences across feature 1
  ale1 <- ale[, list(yhat_diff = y_hat_cumsum[2:nrow(.SD)] - y_hat_cumsum[1:(nrow(.SD) - 1)],
                     interval1 = interval1[2:(nrow(.SD))],
                     count = count[2:(nrow(.SD))]),
              by = c("class", "interval2")
  ]

  # Then take the prediction at the mid point of each interval,
  # which is the mean of the prediction at the end points and
  # use it calculate the mean, weighted by the number of data instances per cell
  ale1 <- ale1[, list(ale1 = sum(count[2:nrow(.SD)] *
                                   (yhat_diff[1:(nrow(.SD) - 1)] +
                                      yhat_diff[2:(nrow(.SD))]) / 2) /
                        sum(count[2:nrow(.SD)])),
               by = c("class", "interval1")]
  ale1 <- ale1[, list(ale1 = c(0, cumsum_na(ale1)),
                      interval1 = c(0, interval1)),
               by = c("class")]

  # Computing the first-order effect of feature 2
  # First, take the differences across feature 2
  ale2 <- ale[, list(yhat_diff = y_hat_cumsum[2:nrow(.SD)] - y_hat_cumsum[1:(nrow(.SD) - 1)],
                     interval2 = interval2[2:(nrow(.SD))],
                     count = count[2:(nrow(.SD))]),
              by = c("class", "interval1")
  ]
  # Then take the prediction at the mid point of each interval,
  # which is the mean of the prediction at the end points and
  # take calculate the mean, weighted by the number of data instances per cell
  ale2 <- ale2[, list(ale2 = sum(count[2:nrow(.SD)] *
                                   (yhat_diff[1:(nrow(.SD) - 1)] +
                                      yhat_diff[2:(nrow(.SD))]) / 2) /
                        sum(count[2:nrow(.SD)])),
               by = c("class", "interval2")]
  ale2 <- ale2[, list(ale2 = c(0, cumsum_na(ale2)),
                      interval2 = c(0, interval2)),
               by = "class"]

  # For each cell compute the average prediction through mean of the cell corners
  # then again a mean over all cells, weighted by the number of data points in the cell
  cls <- unique(ale$class)
  cls
  # Compute mean over all cells weighted by number of data points
  fJ0 <- unlist(lapply(cls, function(cl) {
    ale_cl <- ale[class == cl, ]
    ale1_cl <- ale1[class == cl, ]
    ale2_cl <- ale2[class == cl, ]
    dd <- data.table::dcast(ale_cl,
                            interval1 ~ interval2,
                            value.var = "y_hat_cumsum",
                            drop = FALSE)[, -1]
    dd <- dd - outer(ale1_cl$ale1, rep(1, nrow(ale2_cl))) -
      outer(rep(1, nrow(ale1_cl)), ale2_cl$ale2)
    sum(cell_counts * (dd[1:(nrow(dd) - 1), 1:(ncol(dd) - 1)] +
                         dd[1:(nrow(dd) - 1), 2:ncol(dd)] +
                         dd[2:nrow(dd), 1:(ncol(dd) - 1)] +
                         dd[2:nrow(dd), 2:ncol(dd)]) / 4, na.rm = TRUE) / sum(cell_counts)
  }))
  # Add computed weighted mean and main effects of two features to ale data.frame
  fJ0 <- data.frame("fJ0" = fJ0, class = cls)
  ale <- merge(ale, fJ0, by = c("class"))
  ale <- merge(ale, ale1, by = c("class", "interval1"))
  ale <- merge(ale, ale2, by = c("class", "interval2"))
  # Remove weighted mean prediction and feature main effects from ALE
  ale$ale <- ale$y_hat_cumsum - ale$ale1 - ale$ale2 - ale$fJ0

  # For later plotting, define the rectangles
  # These are not the same as the cells, which is a bit counterintuitive
  # each value in fJ is where the cells cross
  # instead of coloring the cells, we color a rectangle around each cell cross point
  # and for this we need to compute rectangles around these cross points
  # in image() this happens automatically (see ALEPlot::ALEPlot)
  # for the edges, simply use the grid value as the outer values
  # Compute interval distances based on feature 1
  interval_dists <- diff(grid_dt1[c(1, 1:nrow(grid_dt1), nrow(grid_dt1))][[1]])
  interval_dists <- 0.5 * interval_dists
  # Use interval distances to compute right and left rectangle points
  ale$right <- grid_dt1[ale$interval1 + 1, ] + interval_dists[ale$interval1 + 2]
  ale$left <- grid_dt1[ale$interval1 + 1, ] - interval_dists[ale$interval1 + 1]
  # Compute interval distances based on feature 2
  interval_dists2 <- diff(grid_dt2[c(1, 1:nrow(grid_dt2), nrow(grid_dt2))][[1]])
  interval_dists2 <- 0.5 * interval_dists2
  # Use interval distances to compute bottom and top rectangle points
  ale$bottom <- grid_dt2[ale$interval2 + 1, ] + interval_dists2[ale$interval2 + 2]
  ale$top <- grid_dt2[ale$interval2 + 1, ] - interval_dists2[ale$interval2 + 1]
  # Add feature values for x and y axis
  ale[, feature_name[1]] <- grid_dt1[ale$interval1 + 1, ]
  ale[, feature_name[2]] <- grid_dt2[ale$interval2 + 1, ]
  # Delete unnecessary columns
  ale <- ale[, setdiff(colnames(ale), c(
    "fJ0", "ale1", "ale2", "y_hat_cumsum", "count",
    "interval1", "interval2"
  )), with = FALSE]

  # Convert to data.frame for plotting
  ale_df <- data.frame(ale)
  # Convert feature names to symbols
  feature_name_sym1 <- sym(feature_name[1])
  feature_name_sym2 <- sym(feature_name[2])

  # ALE Plot
  ggplot(ale_df, aes(x = !!feature_name_sym1, y = !!feature_name_sym2)) +
    geom_rect(aes(ymin = bottom, ymax = top, fill = ale, xmin = left, xmax = right)) +
    scale_x_continuous(feature_name[1]) +
    scale_y_continuous(feature_name[2]) +
    stat_contour(aes(z = ale), bins = 15, colour = "#000000") +
    scale_fill_gradient(low = "#ffff33",
                        high = "#990000",
                        guide = "colorbar") +
    geom_rug(data = data_dt, aes(x = !!feature_name_sym1, y = !!feature_name_sym2),
             alpha = 0.8, position = "jitter")
}


#' Compute ALE for 2 features, one numerical, one categorical
#'
#' @param data_dt the data.frame with same columns as training data
#' @param predict_function Predict function of type: f(model, newdata, times)
#' @param model the explained model
#' @param feature_name The column names of the two features for which to compute ALE as vector
#' @param grid_length Numerical values determining the number of quantile values for numerical feature
#' @param times a vector of times, that are used for evaluation of survival function and cumulative hazard function
calculate_ale_num_cat <- function(data_dt,
                                  predict_function,
                                  model,
                                  feature_name,
                                  grid_length,
                                  times) {
  # Argument checks
  assert_data_table(data_dt)
  assert_character(feature_name, len = 2, unique = TRUE)
  assert_numeric(grid_length, any.missing = FALSE, len = 1, null.ok = FALSE)
  assert_numeric(times, any.missing = FALSE, min.len = 1, null.ok = FALSE)
  assert_function(predict_function, args = c("model","newdata","times"))
  assert_true(feature_name[1] %in% colnames(data_dt))
  assert_true(feature_name[2] %in% colnames(data_dt))

  # Figure out which feature is numeric and which categeorical
  feature_index <- match(feature_name, colnames(data_dt))
  x_num_index <- ifelse(inherits(data_dt[, feature_name[1], with = FALSE][[1]],
                                 "numeric"), 1, 2)
  x_num <- data_dt[, feature_name[x_num_index], with = FALSE][[1]]
  feature_index_num <- feature_index[x_num_index]

  # Number of quantile points for determined by grid length
  quantile_vals <- as.numeric(quantile(data_dt[, ..feature_index_num][[1]],
                                       seq(1/100, 1, length.out = grid_length),
                                       type = 1))
  # Quantile points vector of length grid_length + 1
  quantile_vec <- c(min(data_dt[, ..feature_index_num][[1]]), quantile_vals)
  # Convert quantile points to data.table, make sure quantile points are unique
  grid_dt <- data.table(unique(quantile_vec))

  # We can only compute ALE within min and max boundaries of given intervals
  # This part is only relevat for user-defined intervals
  #data <- data[which((x.num >= min(grid.dt[[1]])) &
  #                   (x.num <= max(grid.dt[[1]])))]

  # Figure out which feature is categorical
  x_cat_index <- setdiff(c(1, 2), x_num_index)
  x_cat <- data_dt[, feature_name[x_cat_index], with = FALSE][[1]]
  # Extract original level order
  levels_original <- levels(x_cat)

  # Determine level order according to distances to other features in the training
  # data, if already ordered, than use that
  if (inherits(x_cat, "ordered")) {
    level_order <- 1:nlevels(x_cat)
  } else {
    # reorders according to the distance of the levels in the other features
    level_order <- order_levels(data_dt, feature_name[x_cat_index])
  }

  # The new order of the levels
  levels_ordered <- levels_original[level_order]
  # The feature with the levels in the new order
  x_cat_ordered <- order(level_order)[as.numeric(x_cat)]
  # Filter rows for which the category can be increased/decreased (this
  # excludes minimum and maximum category rows)
  row_ind_increase <- (1:nrow(data_dt))[x_cat_ordered < nlevels(x_cat)]
  row_ind_decrease <- (1:nrow(data_dt))[x_cat_ordered > 1]

  # Obtain index values for which quantile interval each feature value is assigned to
  interval_index <- findInterval(data_dt[[feature_name[x_num_index]]],
                                 grid_dt[[1]], left.open = TRUE)
  # Data point in the left most interval should be in interval 1, not zero
  interval_index[interval_index == 0] <- 1


  # Since a categorical feature is involved, we do the second-order differences twice
  # Once for the instances that jump from k to k+1 and once for those that jump from category k+1 to k
  # Crete auxiliary data.tables for category increase predictions
  X_low1_low2 <- X_up1_low2 <- X_low1_up2 <- X_up1_up2 <- copy(data_dt)

  # Put category as first dimension without loss of generality
  X_low1_low2[row_ind_increase, feature_name[x_num_index]] <-
    grid_dt[interval_index, ][[1]][row_ind_increase]
  X_low1_up2[row_ind_increase, feature_name[x_num_index]] <-
    grid_dt[interval_index + 1, ][[1]][row_ind_increase]
  X_up1_low2[row_ind_increase, feature_name[x_cat_index]] <-
    levels_ordered[x_cat_ordered[row_ind_increase] + 1]
  X_up1_low2[row_ind_increase, feature_name[x_num_index]] <-
    grid_dt[interval_index, ][[1]][row_ind_increase]
  X_up1_up2[row_ind_increase, feature_name[x_cat_index]] <-
    levels_ordered[x_cat_ordered[row_ind_increase] + 1]
  X_up1_up2[row_ind_increase, feature_name[x_num_index]] <-
    grid_dt[interval_index + 1, 1][[1]][row_ind_increase]

  # Getting all predictions from the model for category increases
  y_hat_11 <- predict_function(model = model,
                               newdata = X_low1_low2[row_ind_increase, ],
                               times = times)
  y_hat_21 <- predict_function(model = model,
                               newdata = X_up1_low2[row_ind_increase, ],
                               times = times)
  y_hat_12 <- predict_function(model = model,
                               newdata = X_low1_up2[row_ind_increase, ],
                               times = times)
  y_hat_22 <- predict_function(model = model,
                               newdata = X_up1_up2[row_ind_increase, ],
                               times = times)

  # Second-order differences for instances jumping from category k to k+1
  d_increase <- as.data.table((y_hat_22 - y_hat_21) - (y_hat_12 - y_hat_11))
  # Average/Marginalize over time if survival predictions were made for multiple time points
  d_increase <- d_increase[, .(.mean = rowMeans(.SD))]

  # Crete auxiliary data.tables for category decrease prediction
  X_low1_low2 <- X_up1_low2 <- X_low1_up2 <- X_up1_up2 <- copy(data_dt)

  # Put category as first dimension without loss of generality
  X_low1_low2[row_ind_decrease, feature_name[x_num_index]] <-
    grid_dt[interval_index, ][[1]][row_ind_decrease]
  X_low1_low2[row_ind_decrease, feature_name[x_cat_index]] <-
    levels_ordered[x_cat_ordered[row_ind_decrease] - 1]
  X_low1_up2[row_ind_decrease, feature_name[x_cat_index]] <-
    levels_ordered[x_cat_ordered[row_ind_decrease] - 1]
  X_low1_up2[row_ind_decrease, feature_name[x_num_index]] <-
    grid_dt[interval_index + 1, ][[1]][row_ind_decrease]
  X_up1_low2[row_ind_decrease, feature_name[x_num_index]] <-
    grid_dt[interval_index, ][[1]][row_ind_decrease]
  X_up1_up2[row_ind_decrease, feature_name[x_num_index]] <-
    grid_dt[interval_index + 1, 1][[1]][row_ind_decrease]


  # Getting all predictions from the model for category decreases
  y_hat_11 <- predict_function(model = model,
                               newdata = X_low1_low2[row_ind_decrease, ],
                               times = times)
  y_hat_21 <- predict_function(model = model,
                               newdata = X_up1_low2[row_ind_decrease, ],
                               times = times)
  y_hat_12 <- predict_function(model = model,
                               newdata = X_low1_up2[row_ind_decrease, ],
                               times = times)
  y_hat_22 <- predict_function(model = model,
                               newdata = X_up1_up2[row_ind_decrease, ],
                               times = times)

  # Second-order differences for instances jumping from category k+1 to k
  d_decrease <- as.data.table((y_hat_22 - y_hat_21) - (y_hat_12 - y_hat_11))
  # Average/Marginalze over time in case survival predictions were made for multiple time points
  d_decrease <- d_decrease[, .(.mean = rowMeans(.SD))]

  # Combine the deltas over k to k+1 and k+1 to k to the same jump
  deltas <- rbind(d_increase, d_decrease)
  # Create delta data.frame with columns for interval values of the categorical and numerical features corresponding to specific delta value
  deltas <- cbind(.deltas = deltas[[1]],
                  level = c(x_cat_ordered[row_ind_increase], x_cat_ordered[row_ind_decrease] - 1),
                  num = interval_index[c(row_ind_increase, row_ind_decrease)]
  )
  deltas <- as.data.table(deltas)

  # Reshape long data.table to wide, only relevant if there are multidimensional/multi-type predictions
  deltas <- data.table(data.table::melt(deltas,
                                        measure.vars = ".deltas",
                                        variable.name = "class",
                                        value.name = ".yhat.diff"))

  # Ensure that all potential grid combinations of categorical and numerical feature are included
  interval_grid <- expand.grid(
    level = unique(deltas$level),
    num = unique(deltas$num),
    class = unique(deltas$class)
  )
  deltas <- merge(deltas, interval_grid,
                  on = c("level", "num", "class"),
                  all.y = TRUE)

  # Average results per cell (cell = interval combinations of 2 features)
  setkeyv(deltas, c("class", "level", "num"))
  deltas <- deltas[, list(.yhat.diff = mean(.yhat.diff)),
                   by = c("class", "level", "num")]

  # Fill empty cells with the closest neighbor cells value
  # Remember cell status (missing or non-missinf) for later
  deltas_na_cell <- copy(deltas)
  deltas_na_cell$missing <- is.na(deltas$.yhat.diff)
  deltas_na_cell$.yhat.diff <- NULL

  # Replace the missing ones with the closest non-missing difference (measured in number of intervals)
  deltas <- rbindlist(future.apply::future_lapply(unique(deltas$class), function(cl) {
    impute_cells(deltas[class == cl, ],
                 grid2 = grid_dt,
                 x1.ind = "level", x2.ind = "num"
    )
  }))

  # Accumulate the predictions from bottom left to top right for numerical and categorical feature
  deltas <- deltas[, list(ale = cumsum(c(0, .yhat.diff)), num = c(0, num)),
                   by = c("class", "level")]
  deltas <- deltas[, list(ale = cumsum(c(0, ale)), level = c(0, level)),
                   by = c("class", "num")]

  # Obtain number of cells needed for weighting
  # Frequency table for factor level (categorical feature) and quantile value (numerical feature) combination
  cell_counts <- table(x_cat_ordered, interval_index)
  cell_counts_m <- as.data.table(cell_counts)
  data.table::setnames(cell_counts_m, "N", "value")
  # Rename columns
  colnames(cell_counts_m) <- c("level", "num", "count")
  # Decrease category levels by 1
  cell_counts_m$level <- as.numeric(as.character(cell_counts_m$level)) - 1
  # Convert quantile values to numeric
  cell_counts_m$num <- as.numeric(as.character(cell_counts_m$num))
  # Merge delta data.table with cell counts
  deltas <- merge(deltas, cell_counts_m, on = c("level", "num"), all.x = TRUE)
  # Replace NA counts with 0
  deltas$count[is.na(deltas$count)] <- 0

  # Computing the first-order effect of feature 2, the numerical feature.
  # Take the differences across feature 2
  deltas2 <- deltas[, list(ale2 = ale[2:nrow(.SD)] - ale[1:(nrow(.SD) - 1)],
                           num = num[2:(nrow(.SD))],
                           count = count[2:(nrow(.SD))]),
                    by = c("class", "level")
  ]
  # Then take the weighted average over the categories to get the effect at each numerical grid point.
  deltas2 <- deltas2[, list(ale2 = sum(count * ale2) / sum(count)),
                     by = c("class", "num")]
  ale2 <- deltas2[, list(ale2 = c(0, cumsum(ale2)), num = c(0, num)),
                  by = c("class")]

  # Computing the first-order effect of feature 1, the categorical feature.
  # Take the differences across feature 1
  deltas1 <- deltas[, list(
    ale1 = ale[2:nrow(.SD)] - ale[1:(nrow(.SD) - 1)],
    level = level[2:(nrow(.SD))],
    count = (count[2:nrow(.SD)] + count[1:(nrow(.SD) - 1)]) / 2
  ), by = c("class", "num")]
  # Mid points between numerical intervals
  deltas1 <- deltas1[, list(
    ale1 = (ale1[2:nrow(.SD)] + ale1[1:(nrow(.SD) - 1)]) / 2,
    count = count[2:nrow(.SD)],
    num = num[2:nrow(.SD)]
  ), by = c("class", "level")]
  deltas1 <- deltas1[, list(ale1 = sum(ale1 * count) / sum(count)),
                     by = c("class", "level")]
  ale1 <- deltas1[, list(ale1 = c(0, cumsum_na(ale1)),
                         level = c(0, level)),
                  by = c("class")]

  # For each cell compute the average prediction through mean of the cell corners
  # then again a mean over all cells, weighted by the number of data points in the cell
  # Extract different prediction classes
  cls <- unique(deltas$class)
  # Compute mean over all cells separately for each prediction class weighted by number of data points
  fJ0 <- unlist(lapply(cls, function(cl) {
    deltas_cl <- deltas[class == cl, ]
    ale1_cl <- ale1[class == cl, ]
    ale2_cl <- ale2[class == cl, ]
    dd <- as.matrix(data.table::dcast(deltas_cl, level ~ num,
                                      value.var = "ale", drop = FALSE))[, -1]
    dd <- dd - outer(ale1_cl$ale1,
                     rep(1, nrow(ale2_cl))) - outer(rep(1, nrow(ale1_cl)),
                                                    ale2_cl$ale2)
    sum(cell_counts * (dd[, 1:(ncol(dd) - 1)] + dd[, 2:ncol(dd)]) / 2,
        na.rm = TRUE) / sum(cell_counts)
  }))

  # Add computed weighted mean and main effects of two features to ale dataframe
  fJ0 <- data.frame("fJ0" = fJ0, "class" = cls)
  deltas <- merge(deltas, fJ0, by = c("class"))
  deltas <- merge(deltas, ale1, by = c("class", "level"))
  deltas <- merge(deltas, ale2, by = c("class", "num"))
  # Remove weighted mean prediction to center globally and feature main effects (first order ALE effects) from ALE
  deltas$ale <- deltas$ale - deltas$ale1 - deltas$ale2 - deltas$fJ0

  # For later plotting, define the rectangles
  # These are not the same as the cells, which is a bit counterintuitive
  # each value in fJ is where the cells cross
  # instead of coloring the cells, we color a rectangle around each cell cross point
  # and for this we need to compute rectangles around these cross points
  # in image() this happens automatically (see ALEPlot::ALEPlot)
  # for the edges, simply use the grid value as the outer values
  # Right and left rectangle corners based on categorical variable
  deltas$right <- 1 + deltas$level + 0.5
  deltas$left <- 1 + deltas$level - 0.5
  # Bottom and top rectangle corners based on numerical variable
  interval_dists2 <- diff(grid_dt[c(1, 1:nrow(grid_dt), nrow(grid_dt))][[1]])
  interval_dists2 <- 0.5 * interval_dists2
  deltas$bottom <- grid_dt[deltas$num + 1, ] + interval_dists2[deltas$num + 2]
  deltas$top <- grid_dt[deltas$num + 1, ] - interval_dists2[deltas$num + 1]
  # Add corresponding numerical and categorical feature values from the grid for each ale value
  deltas[, feature_name[x_num_index]] <- grid_dt[deltas$num + 1, ]
  deltas[, feature_name[x_cat_index]] <- factor(levels_ordered[deltas$level + 1],
                                                levels = levels_ordered)
  deltas[, feature_name[x_cat_index]] <- as.numeric(deltas[, feature_name[x_cat_index],
                                                           with = FALSE][[1]])
  # Delete unnecessary columns
  deltas <- deltas[, setdiff(colnames(deltas), c(
    "fJ0", "ale1", "ale2", "count",
    "level", "num"
  )), with = FALSE]
  deltas$.type <- "ale"
  # Convert data.table to data.frame for plotting
  deltas_df <- data.frame(deltas)
  # Convert feature name string to symbol
  feature_name_sym1 <- sym(feature_name[x_cat_index])
  feature_name_sym2 <- sym(feature_name[x_num_index])

  # ALE Plot
  ggplot(deltas_df, aes(x = !!feature_name_sym1, y = !!feature_name_sym2)) +
    geom_rect(aes(ymin = bottom,
                  ymax = top,
                  fill = ale,
                  xmin = left,
                  xmax = right)) +
    scale_x_continuous(feature_name[x_cat_index]) +
    scale_y_continuous(feature_name[x_num_index]) +
    stat_contour(aes(z = ale), bins = 15, colour = "#000000") +
    scale_fill_gradient(low = "#ffff33",
                        high = "#990000",
                        guide = "colorbar") +
    scale_x_discrete(limits = levels_original) +
    geom_rug(data = data_dt, aes(x = min(as.numeric(x_cat)), y = !!feature_name_sym2),
             alpha = 0.8, position = "jitter", sides = 'l')
}


####' Quick test surv_ale function
library(survival)
library(survex)
library(ranger)
library(data.table)
library(ggplot2)
library(checkmate)

#' vet data + ranger test
set.seed(123)

vet <- survival::veteran

rsf_ranger <- ranger::ranger(survival::Surv(time, status) ~ ., data = veteran,
                             respect.unordered.factors = TRUE, num.trees = 100,
                             mtry = 3, max.depth = 5)
rsf_ranger_exp <- survex::explain(rsf_ranger, data = veteran[, -c(3, 4)],
                                  y = Surv(veteran$time, veteran$status))

surv_ale(explainer = rsf_ranger_exp, feature_name = "age")
surv_ale(explainer = rsf_ranger_exp, feature_name = "karno", marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "age", times = 87)
surv_ale(explainer = rsf_ranger_exp, feature_name = "age", times = 1:80, marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "age", times = c(1,10,50,100,200,300), marginalize_over_time = FALSE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "celltype")
surv_ale(explainer = rsf_ranger_exp, feature_name = "celltype", marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "celltype", times = 87)
surv_ale(explainer = rsf_ranger_exp, feature_name = "celltype", times = 1:80, marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "celltype", times = c(1,10,50,100,200,300), marginalize_over_time = FALSE)
surv_ale(explainer = rsf_ranger_exp, feature_name = c("diagtime", "age"), marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = c("diagtime", "age"))
surv_ale(explainer = rsf_ranger_exp, feature_name = c("celltype", "age"), marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = c("celltype", "age"))


#' DHS data + ranger test
rsf_ranger <- ranger::ranger(survival::Surv(survivaltime, censor) ~ ., data = analysis_df,
                             respect.unordered.factors = TRUE, num.trees = 100,
                             mtry = 3, max.depth = 5)
rsf_ranger_exp <- survex::explain(rsf_ranger, data = analysis_df[, -c(1, 2)],
                                  y = Surv(analysis_df$survivaltime, analysis_df$censor))

surv_ale(explainer = rsf_ranger_exp, feature_name = "mother_age")
surv_ale(explainer = rsf_ranger_exp, feature_name = "mother_age", marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "mother_age", times = 87)
surv_ale(explainer = rsf_ranger_exp, feature_name = "mother_age", times = 1:80, marginalize_over_time = TRUE)

surv_ale(explainer = rsf_ranger_exp, feature_name = "wealth_idx")
surv_ale(explainer = rsf_ranger_exp, feature_name = "wealth_idx", marginalize_over_time = TRUE)
surv_ale(explainer = rsf_ranger_exp, feature_name = "wealth_idx", times = 87)
surv_ale(explainer = rsf_ranger_exp, feature_name = "wealth_idx", times = 1:80, marginalize_over_time = TRUE)

surv_ale(explainer = rsf_ranger_exp, feature_name = c("wealth_idx", "mother_age"))
surv_ale(explainer = rsf_ranger_exp, feature_name = c("wealth_idx", "mother_age"), marginalize_over_time = TRUE)







