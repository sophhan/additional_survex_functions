#' @title Order levels of a categorical features
#'
#' @description
#' Orders the levels by their similarity in other features. Computes per feature
#' the distance, sums up all distances and does multi-dimensional scaling
#'
#' @details
#' Goal: Compute the distances between two categories.
#' Input: Instances from category 1 and 2
#'
#' 1. For all features, do (excluding the categorical feature for which we are computing the order):
#'  - If the feature is numerical: Take instances from category 1, calculate the
#'  empirical cumulative probability distribution function (ecdf) of the
#'  feature. The ecdf is a function that tells us for a given feature value, how
#'  many values are smaller. Do the same for category 2. The distance is the
#'  absolute maximum point-wise distance of the two ecdf. Practically, this
#'  value is high when the distribution from one category is strongly shifted
#'  far away from the other. This measure is also known as the
#'  Kolmogorov-Smirnov distance
#'  (\url{https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test}).
#'  - If the feature is categorical: Take instances from category 1 and
#'  calculate a table with the relative frequency of each category of the other
#'  feature. Do the same for instances from category 2. The distance is the sum
#'  of the absolute difference of both relative frequency tables.
#' 2. Sum up the distances over all features
#'
#' This algorithm we run for all pairs of categories.
#' Then we have a k times k matrix, when k is the number of categories, where
#' each entry is the distance between two categories. Still not enough to have a
#' single order, because, a (dis)similarity tells you the pair-wise distances,
#' but does not give you a one-dimensional ordering of the classes. To kind of
#' force this thing into a single dimension, we have to use a dimension
#' reduction trick called multi-dimensional scaling. This can be solved using
#' multi-dimensional scaling, which takes in a distance matrix and returns a
#' distance matrix with reduced dimension. In our case, we only want 1 dimension
#' left, so that we have a single ordering of the categories and can compute the
#' accumulated local effects. After reducing it to a single ordering, we are
#' done and can use this ordering to compute ALE. This is not the Holy Grail how
#' to order the factors, but one possibility.
#'
#' @param data_dt data.table with the training data
#' @param feature_name the name of the categorical feature
#' @return the order of the levels (not levels itself)

order_levels <- function(data_dt, feature_name) {
  # Perform argument tests
  assert_data_frame(data_dt)
  assert_character(feature_name)
  assert_true(feature_name %in% names(data_dt))
  assert_factor(data_dt[, feature_name, with = FALSE][[1]])

  data_dt[, feature_name] <- droplevels(data_dt[, feature_name, with = FALSE])
  feature <- data_dt[, feature_name, with = FALSE][[1]]
  x.count <- as.numeric(table(data_dt[, feature_name, with = FALSE]))
  x.prob <- x.count / sum(x.count)
  K <- nlevels(data_dt[, feature_name, with = FALSE])

  dists <- lapply(setdiff(colnames(data_dt), feature_name), function(x) {
    feature.x <- data_dt[, x, with = FALSE][[1]]
    dists <- expand.grid(levels(feature), levels(feature))
    colnames(dists) <- c("from.level", "to.level")
    if (inherits(feature.x, "factor")) {
      A <- table(feature, feature.x) / x.count
      dists$dist <- rowSums(abs(A[dists[, "from.level"], ] - A[dists[, "to.level"], ])) / 2
    } else {
      quants <- quantile(feature.x, probs = seq(0, 1, length.out = 100), na.rm = TRUE, names = FALSE)
      ecdfs <- data.frame(lapply(levels(feature), function(lev) {
        x.ecdf <- ecdf(feature.x[feature == lev])(quants)
      }))
      colnames(ecdfs) <- levels(feature)
      ecdf.dists.all <- abs(ecdfs[, dists$from.level] - ecdfs[, dists$to.level])
      dists$dist <- apply(ecdf.dists.all, 2, max)
    }
    dists
  })

  dists.cumulated.long <- as.data.table(Reduce(function(d1, d2) {
    d1$dist <- d1$dist + d2$dist
    d1
  }, dists))
  dists.cumulated <- data.table::dcast(dists.cumulated.long, from.level ~ to.level, value.var = "dist")[, -1]
  diag(dists.cumulated) <- 0
  scaled <- cmdscale(dists.cumulated, k = 1)
  order(scaled)
}

#' Impute missing cells of grid
#'
#' by default assumes first column of cell.dat is x1 and second is x2
#' leave grid1 NULL if feature x1 is a factor
#' the difference variable has to be named .yhat.diff
#'
#' @param cell.dat data.table with at least 4 columns: .yhat.diff and the two interval indices.
#' Make sure that empty cells are also included and cell.dat is not the sparse representation.
#' @param grid1 data.frame where each row is the actual value for a given interval index for feature 1.
#' If empty impute_cells  assumes that the feature is categorical (factor).
#' @param grid2 data.frame where each row is the actual value for a given interval index for feature 2
#' @param x1.ind column number or name of cell.dat for feature 1. If one feature is categorical, has to be x1
#' @param x2.ind column number or name of cell.dat for feature 2
#' @keywords internal
impute_cells <- function(cell.dat, grid1 = NULL, grid2, x1.ind = 1, x2.ind = 2) {
  # Perform argument tests
  assert_data_table(cell.dat)
  assert_data_table(grid1, null.ok = TRUE)
  assert_data_table(grid2)

  if (!requireNamespace("yaImpute")) {
    stop("Please install package yaImpute")
  }

  # Making sure cell.dat contains all possible cells
  stopifnot(nrow(cell.dat) == length(unique(cell.dat[, x1.ind, with = FALSE][[1]])) * length(unique(cell.dat[, x2.ind, with = FALSE][[1]])))
  d.miss.ind <- is.na(cell.dat$.yhat.diff)

  # We don't have to impute anything when all cells are missing.
  if (!any(d.miss.ind)) {
    return(cell.dat)
  }

  # Distinguishes for normalization between categorical and numerical feature x1
  if (is.null(grid1)) {
    # For categorical feature, the range is the number of levels, and
    # scaled to range 1/n.categories to 1
    range.x1 <- length(unique(cell.dat[[x1.ind]]))
    x1.normalized <- cell.dat[[x1.ind]] / range.x1
  } else {
    # For numerical feature range is from min to max
    # normalizing means taking the mean of neighbours and dividing by range
    range.x1 <- max(grid1[[1]]) - min(grid1[[1]])
    x1.normalized <- (grid1[cell.dat[[x1.ind]], ][[1]] + grid1[cell.dat[[x1.ind]] + 1, ][[1]]) / (2 * range.x1)
  }
  # same for the second numerical feature
  range.x2 <- max(grid2[[1]]) - min(grid2[[1]])
  x2.normalized <- (grid2[cell.dat[[x2.ind]], ][[1]] + grid2[cell.dat[[x2.ind]] + 1, ][[1]]) / (2 * range.x2)
  # preparation for yaImpute::ann function
  z.na <- cbind(x1.normalized[d.miss.ind], x2.normalized[d.miss.ind])
  z.non.na <- cbind(x1.normalized[!d.miss.ind], x2.normalized[!d.miss.ind])
  # matches closest non-missing neighbour
  nbrs <- yaImpute::ann(as.matrix(z.non.na), as.matrix(z.na), k = 1, verbose = FALSE)$knnIndexDist[, 1]
  cell.dat[d.miss.ind, ".yhat.diff"] <- cell.dat[which(!d.miss.ind)[nbrs], .yhat.diff]
  cell.dat
}

# accumulate over the intervals
cumsum_na <- function(values) {
  values[is.na(values)] <- 0
  cumsum(values)
}

# function to cross join two data.tables with multiple rows
cross_join <- function(dt1,dt2){
  cj = CJ(1:nrow(dt1),1:nrow(dt2))
  cbind(dt1[cj[[1]],],dt2[cj[[2]],])
}
