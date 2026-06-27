#' @importFrom dplyr bind_rows
#' @importFrom randomForest randomForest
#' @importFrom e1071 svm
#' @importFrom nnet multinom nnet class.ind
#' @importFrom MASS mvrnorm
#' @importFrom stats predict quantile pbinom sd setNames rnorm runif
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics abline axis boxplot grid layout legend mtext par plot.new plot.window points segments strwidth text
#'
#' @title Classes Mapping function for HNP Algorithm
#' Map ordered class labels to internal levels "1", ..., "T"
#'
#' @description Validate the class column and re-label the provided class names
#'   to internal factor levels "1", "2", ..., "T". The input order of the
#'   labels in `...` is treated as the importance order: the first label maps to
#'   level "1" (highest priority), and the last label maps to the
#'   level "T" (lowest priority). Supports any number of classes `T >= 2`.
#' @param data A data.frame or data.table containing the dataset.
#' @param class_col Character scalar. Name of the class/label column in `data`.
#' @param ... Two or more original class labels given in importance order. The
#'   first label maps to "1", the second to "2", and so on. Labels must be
#'   non-empty, non-NA, and unique.
#' @return A data.frame with the input `data` added with `class_col` converted to a factor whose levels
#'   are the internal labels c("1","2",...,"T").
#' @examples
#' df <- data.frame(y = c("low", "mid", "high", "mid"), x1 = rnorm(4))
#' df2 <- hnp_map_classes(df, class_col = "y", "high", "mid", "low")
#' @export
hnp_map_classes <- function(data, class_col, ...) {
  # ---------------------------
  # Supports ordered class mapping for any T >= 2.
  # The input order is treated as importance_order:
  # the 1st label -> "1" (highest priority), and the T-th label -> "T" (lowest priority).
  # ---------------------------
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  if (!class_col %in% colnames(data)) {
    stop("error: the specified class column '", class_col, "' does not exist in the data")
  }
  
  classes <- unlist(list(...), use.names = FALSE)
  classes <- as.character(classes)
  
  if (length(classes) < 2) {
    stop("error: please provide at least 2 class labels to map, e.g., hnp_map_classes(data, 'y', 'A', 'B').")
  }
  
  if (any(is.na(classes)) || any(classes == "")) {
    stop("error: class labels in '...' must be non-empty and non-NA.")
  }
  
  if (any(duplicated(classes))) {
    stop("error: duplicated class labels in mapping: ",
         paste(unique(classes[duplicated(classes)]), collapse = ", "))
  }
  
  unique_classes <- unique(as.character(data[[class_col]]))
  missing_in_data <- setdiff(classes, unique_classes)
  
  if (length(missing_in_data) > 0) {
    warning("warning: some specified classes do not exist in the data: ",
            paste(missing_in_data, collapse = ", "))
    print("classes in the data:")
    print(unique_classes)
    print("specified classes in importance_order:")
    print(classes)
  }
  
  Tn <- length(classes)
  target_levels <- as.character(seq_len(Tn))
  
  data[[class_col]] <- factor(
    as.character(data[[class_col]]),
    levels = classes,
    labels = target_levels
  )
  
  cat("finished mapping:\n")
  for (i in seq_len(Tn)) {
    cat("  ", classes[i], " -> ", i, "\n", sep = "")
  }
  
  data
}

hnp_delta_search <- function(n, level, delta) {
  # ---------------------------
  # Find the largest k such that
  #   P{Binomial(n, level) <= k - 1} <= delta.
  # If the sample size is insufficient, return 0.
  # Use pbinom() to avoid underflow or overflow from choose(n, k) in large samples.
  # ---------------------------
  n <- as.integer(n)
  level <- as.numeric(level)
  delta <- as.numeric(delta)
  
  if (length(n) != 1 || is.na(n) || n <= 0) {
    stop("n must be a positive integer.")
  }
  
  if (length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie in (0, 1).")
  }
  
  if (length(delta) != 1 || !is.finite(delta) || delta <= 0 || delta >= 1) {
    stop("delta must lie in (0, 1).")
  }
  
  if (n < log(delta) / log(1 - level)) {
    return(0L)
  }
  
  best_k <- 0L
  
  for (k in seq_len(n)) {
    v_k <- stats::pbinom(k - 1, size = n, prob = level)
    
    if (is.finite(v_k) && v_k <= delta) {
      best_k <- k
    } else {
      break
    }
  }
  
  best_k + 1
}

hnp_upper_bound <- function(S_it, level, delta_i, score_functions, thresholds, i) {
  S_it <- as.data.frame(S_it)
  n_i <- nrow(S_it)
  
  if (n_i == 0) {
    stop("S_it is empty.")
  }
  
  T_i <- score_functions[[i]]
  Tau <- as.numeric(T_i(S_it))
  
  if (length(Tau) != n_i) {
    stop("score_functions[[", i, "]] must return a vector of length nrow(S_it).")
  }
  
  if (any(!is.finite(Tau)) || any(is.na(Tau))) {
    stop("score_functions[[", i, "]] returned NA, NaN, or non-finite values.")
  }
  
  Tau <- sort(Tau)
  k_i <- hnp_delta_search(n_i, level, delta_i)
  
  if (k_i == 0) {
    stop("Exceed minimum sample size required.")
  }
  
  if (k_i > length(Tau)) {
    stop("Invalid order statistic index in hnp_upper_bound().")
  }
  
  t_i_bar <- Tau[k_i]
  
  if (i <= 1) {
    return(t_i_bar)
  }
  
  if (is.null(thresholds) || length(thresholds) < (i - 1)) {
    stop("thresholds must contain previous thresholds for classes 1,...,i-1.")
  }
  
  filtered_indices <- rep(TRUE, n_i)
  
  for (j in seq_len(i - 1)) {
    score_j <- as.numeric(score_functions[[j]](S_it))
    
    if (length(score_j) != n_i) {
      stop("score_functions[[", j, "]] must return a vector of length nrow(S_it).")
    }
    
    if (any(!is.finite(score_j)) || any(is.na(score_j))) {
      stop("score_functions[[", j, "]] returned NA, NaN, or non-finite values.")
    }
    
    filtered_indices <- filtered_indices & (score_j < thresholds[j])
  }
  
  S_it_prime <- S_it[filtered_indices, , drop = FALSE]
  n_i_prime <- nrow(S_it_prime)
  
  if (n_i_prime == 0) {
    return(t_i_bar)
  }
  
  Tau_i_prime <- as.numeric(T_i(S_it_prime))
  
  if (length(Tau_i_prime) != n_i_prime) {
    stop("score_functions[[", i, "]] must return a vector of length nrow(S_it_prime).")
  }
  
  if (any(!is.finite(Tau_i_prime)) || any(is.na(Tau_i_prime))) {
    stop("score_functions[[", i, "]] returned NA, NaN, or non-finite values on the filtered set.")
  }
  
  Tau_i_prime <- sort(Tau_i_prime)
  
  p_i_hat <- n_i_prime / n_i
  cn_i <- 2 / sqrt(n_i)
  p_i <- p_i_hat + cn_i
  
  if (!is.finite(p_i) || p_i <= 0) {
    return(t_i_bar)
  }
  
  a_i_prime <- level / p_i
  delta_i_prime <- delta_i - exp(-2 * n_i * cn_i^2)
  
  if (!is.finite(delta_i_prime) || delta_i_prime <= 0) {
    return(t_i_bar)
  }
  
  if (!is.finite(a_i_prime) || a_i_prime <= 0 || a_i_prime >= 1) {
    return(t_i_bar)
  }
  
  min_required <- log(delta_i_prime) / log(1 - a_i_prime)
  
  if (is.finite(min_required) && n_i_prime >= min_required) {
    k_i_prime <- hnp_delta_search(n_i_prime, a_i_prime, delta_i_prime)
    
    if ((k_i_prime != 0) && k_i_prime <= length(Tau_i_prime)) {
      t_i_bar <- Tau_i_prime[k_i_prime]
    }
  }
  
  t_i_bar
}

base_function <- function(x, y, method = "randomforest") {
  # ---------------------------
  # Built-in base learners: randomforest / svm / logistic.
  # Training labels should already be internal labels "1", ..., "T".
  # ---------------------------
  method_choices <- c("randomforest", "svm", "logistic")
  method <- match.arg(method, method_choices)
  
  x <- as.data.frame(x)
  y <- as.factor(y)
  
  if (method == "randomforest") {
    model <- randomForest::randomForest(
      x = x,
      y = y,
      ntree = 100,
      mtry = min(2, ncol(x))
    )
  } else if (method == "svm") {
    dataset <- data.frame(x, y = y)
    colnames(dataset)[ncol(dataset)] <- "y"
    model <- e1071::svm(y ~ ., data = dataset, probability = TRUE, scale = TRUE)
  } else if (method == "logistic") {
    dataset <- data.frame(x, y = y)
    colnames(dataset)[ncol(dataset)] <- "y"
    model <- nnet::multinom(y ~ ., data = dataset, trace = FALSE)
  }
  
  model
}

hnp_summary_get_prediction <- function(classifier,
                                       X,
                                       required_levels,
                                       original_levels) {
  X <- as.data.frame(X)
  n <- nrow(X)
  Tn <- length(required_levels)
  
  # ------------------------------------------------------------
  # Convert classifier output to a unified predicted-label vector
  # Supported outputs:
  #   1. A vector of class labels
  #   2. n x T probability matrix
  #   3. n x T score matrix
  # ------------------------------------------------------------
  output_to_label <- function(out) {
    
    if (is.null(out)) {
      return(NULL)
    }
    
    # For a single observation, some predict() methods return a probability vector of length T.
    if (!is.matrix(out) && !is.data.frame(out) && is.numeric(out)) {
      if (n == 1 && length(out) == Tn) {
        out_mat <- matrix(out, nrow = 1)
        
        if (!is.null(names(out))) {
          colnames(out_mat) <- names(out)
        } else {
          colnames(out_mat) <- original_levels
        }
        
        out <- out_mat
      }
    }
    
    # Case A: the output is a matrix/data.frame and is interpreted as class scores or class probabilities
    if (is.matrix(out) || is.data.frame(out)) {
      out <- as.data.frame(out)
      
      if (nrow(out) != n) {
        return(NULL)
      }
      
      if (ncol(out) != Tn) {
        return(NULL)
      }
      
      if (is.null(colnames(out))) {
        colnames(out) <- original_levels
      }
      
      cols <- colnames(out)
      
      # A1: column names are original class labels, e.g., C, A, B
      if (all(original_levels %in% cols)) {
        out <- out[, original_levels, drop = FALSE]
        score_mat <- as.matrix(out)
        storage.mode(score_mat) <- "numeric"
        
        pred <- original_levels[
          max.col(score_mat, ties.method = "first")
        ]
        
        return(as.character(pred))
      }
      
      # A2: column names are internal class labels, e.g., 1, 2, 3
      if (all(required_levels %in% cols)) {
        out <- out[, required_levels, drop = FALSE]
        score_mat <- as.matrix(out)
        storage.mode(score_mat) <- "numeric"
        
        pred <- required_levels[
          max.col(score_mat, ties.method = "first")
        ]
        
        return(as.character(pred))
      }
      
      # A3: column names do not match, but the column count is correct; interpret columns in original_levels order
      colnames(out) <- original_levels
      score_mat <- as.matrix(out)
      storage.mode(score_mat) <- "numeric"
      
      pred <- original_levels[
        max.col(score_mat, ties.method = "first")
      ]
      
      return(as.character(pred))
    }
    
    # Case B: the output is a vector and is interpreted as class labels
    out <- as.character(out)
    
    if (length(out) != n) {
      return(NULL)
    }
    
    out
  }
  
  # ------------------------------------------------------------
  # 1) classifier is a function
  #    Supported outputs:
  #      function(X) returns class labels
  #      function(X) returns an n x T score/probability matrix
  # ------------------------------------------------------------
  if (is.function(classifier)) {
    out <- tryCatch(
      classifier(X),
      error = function(e) NULL
    )
    
    pred <- output_to_label(out)
    
    if (!is.null(pred)) {
      return(pred)
    }
  }
  
  # ------------------------------------------------------------
  # 2) classifier is a fitted model object
  #    Keep the logic consistent with pretrained_model in hnp_umbrella
  # ------------------------------------------------------------
  
  # randomForest and similar models
  out <- tryCatch(
    predict(classifier, newdata = X, type = "prob"),
    error = function(e) NULL
  )
  
  pred <- output_to_label(out)
  if (!is.null(pred)) {
    return(pred)
  }
  
  # nnet::multinom
  out <- tryCatch(
    predict(classifier, newdata = X, type = "probs"),
    error = function(e) NULL
  )
  
  pred <- output_to_label(out)
  if (!is.null(pred)) {
    return(pred)
  }
  
  # nnet::nnet softmax
  out <- tryCatch(
    predict(classifier, newdata = X, type = "raw"),
    error = function(e) NULL
  )
  
  pred <- output_to_label(out)
  if (!is.null(pred)) {
    return(pred)
  }
  
  # e1071::svm(probability = TRUE)
  out <- tryCatch({
    svm_pred <- predict(classifier, newdata = X, probability = TRUE)
    attr(svm_pred, "probabilities")
  }, error = function(e) NULL)
  
  pred <- output_to_label(out)
  if (!is.null(pred)) {
    return(pred)
  }
  
  # Plain predict(), usually returning class labels directly
  out <- tryCatch(
    predict(classifier, newdata = X),
    error = function(e) NULL
  )
  
  pred <- output_to_label(out)
  if (!is.null(pred)) {
    return(pred)
  }
  
  stop(
    "classifier must be a callable function or a fitted model object ",
    "that can produce class labels, class probabilities, or class scores."
  )
}


#' @title HNP Classifier's Performance Summary
#'
#' @description Evaluate a `T`-class classifier produced by
#'   the HNP pipeline (or any compatible model/function). Given features and
#'   true labels, it computes the confusion matrix, class-wise false
#'   positive/negative rates, under- and over-classification errors, overall
#'   accuracy, and a row-normalized error table. The class priority order is
#'   resolved from `importance_order`, from attributes attached to the
#'   `classifier`, or inferred from the labels.
#' @param classifier Either a function `function(X) ...` returning class labels
#'   or an `n x T` score/probability matrix, or a fitted model object (e.g.
#'   randomForest, support vector machine, multinomial logistic regression) from which predictions can be derived.
#' @param X A data.frame of features (or a score matrix) to predict on.
#' @param Y A vector of true class labels of length `nrow(X)`.
#' @param importance_order Optional character vector giving the class priority
#'   order (most severe first). If `NULL`, the order is taken from the
#'   classifier's attributes or inferred from `Y`.
#' @return A list with components of the HNP experiment results including `confusion_matrix`,
#'   `false_positive_rate`, `false_negative_rate`, `overall_accuracy`,
#'   `predictions`, `predictions_sample`, `under_classification_error`,
#'   `over_classification_error`, `total_under_classification_error`,
#'   `total_over_classification_error`, `error_table`, `class_levels`,
#'   `internal_levels`, `label_mapping`, `importance_order`, and `order_source`.
#' @examples
#' set.seed(123)
#' n <- 300
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' Y <- sample(c("A", "B", "C"), n, replace = TRUE)
#' clf <- hnp_umbrella(
#'   X = X, Y = Y,
#'   levels = c(0.1, 0.1),
#'   tolerances = c(0.1, 0.1),
#'   importance_order = c("A", "B", "C"),
#'   method = "logistic"
#' )
#' result <- hnp_summary(clf, X = X, Y = Y, importance_order = c("A", "B", "C"))
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_summary <- function(classifier,
                        X,
                        Y,
                        importance_order = NULL) {
  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  X <- as.data.frame(X)
  
  if (is.null(Y)) {
    stop("Y must be provided.")
  }
  
  if (nrow(X) != length(Y)) {
    stop("X and Y have inconsistent numbers of rows.")
  }
  
  true_raw <- as.character(Y)
  
  # ---------------------------
  # 1) Read existing attributes from classifier
  # ---------------------------
  clf_class_levels     <- attr(classifier, "class_levels")
  clf_label_mapping    <- attr(classifier, "label_mapping")
  clf_importance_order <- attr(classifier, "importance_order")
  clf_severity_order   <- attr(classifier, "severity_order")
  
  # ---------------------------
  # 2) If importance_order is not explicitly provided, infer the default order from labels
  # ---------------------------
  infer_order_from_labels <- function(labels) {
    labels <- unique(as.character(labels))
    labels <- labels[!is.na(labels) & labels != ""]
    
    numeric_labels <- suppressWarnings(as.numeric(labels))
    
    warning_message <- paste0(
      "importance_order is not specified; the default class order will be used. ",
      "Please provide importance_order if a custom priority order is desired."
    )
    
    if (length(labels) > 0 && all(!is.na(numeric_labels))) {
      ordered_labels <- labels[order(numeric_labels)]
      
      return(list(
        ordered_labels = ordered_labels,
        order_source = "default numerical class order",
        warning_message = warning_message
      ))
    } else {
      ordered_labels <- sort(labels)
      
      return(list(
        ordered_labels = ordered_labels,
        order_source = "default character class order",
        warning_message = warning_message
      ))
    }
  }
  
  # ---------------------------
  # 3) Determine the original class order
  # Priority:
  # user-provided importance_order >
  # classifier importance_order attribute >
  # classifier legacy severity_order attribute >
  # names of label_mapping >
  # automatic ordering from Y
  # ---------------------------
  issue_order_warning <- FALSE
  order_warning_message <- NULL
  
  if (!is.null(importance_order)) {
    ordered_original <- as.character(importance_order)
    order_source <- "user-provided importance_order"
  } else if (!is.null(clf_importance_order)) {
    ordered_original <- as.character(clf_importance_order)
    order_source <- "classifier importance_order attribute"
  } else if (!is.null(clf_severity_order)) {
    ordered_original <- as.character(clf_severity_order)
    order_source <- "classifier severity_order attribute"
  } else if (!is.null(clf_label_mapping) && !is.null(names(clf_label_mapping))) {
    ordered_original <- as.character(names(clf_label_mapping))
    order_source <- "names of classifier label_mapping"
  } else {
    inferred <- infer_order_from_labels(true_raw)
    ordered_original <- inferred$ordered_labels
    order_source <- inferred$order_source
    issue_order_warning <- TRUE
    order_warning_message <- inferred$warning_message
  }
  
  if (issue_order_warning) {
    warning(order_warning_message, call. = FALSE)
  }
  
  if (length(ordered_original) < 2) {
    stop("At least two classes are required.")
  }
  
  if (any(is.na(ordered_original)) || any(ordered_original == "")) {
    stop("importance_order contains empty or NA class labels.")
  }
  
  if (any(duplicated(ordered_original))) {
    stop("importance_order contains duplicated class labels.")
  }
  
  # ---------------------------
  # 4) Construct the original -> internal mapping
  # ---------------------------
  if (!is.null(clf_label_mapping)) {
    label_mapping <- as.character(clf_label_mapping)
    names(label_mapping) <- names(clf_label_mapping)
  } else {
    label_mapping <- stats::setNames(
      as.character(seq_along(ordered_original)),
      ordered_original
    )
  }
  
  # ---------------------------
  # 5) Determine internal_levels
  # ---------------------------
  if (!is.null(label_mapping)) {
    internal_levels <- as.character(unname(label_mapping[ordered_original]))
  } else if (!is.null(clf_class_levels)) {
    internal_levels <- as.character(clf_class_levels)
  } else {
    internal_levels <- as.character(seq_along(ordered_original))
  }
  
  if (any(is.na(internal_levels))) {
    stop("importance_order is inconsistent with classifier label_mapping.")
  }
  
  internal_levels <- unique(as.character(internal_levels))
  class_number <- length(internal_levels)
  
  if (class_number != length(ordered_original)) {
    stop("The number of internal class levels is inconsistent with importance_order.")
  }
  
  display_levels <- ordered_original
  
  # ---------------------------
  # 6) Obtain predictions in a unified way
  # ---------------------------
  pred_all <- hnp_summary_get_prediction(
    classifier = classifier,
    X = X,
    required_levels = internal_levels,
    original_levels = display_levels
  )
  
  pred_raw <- as.character(pred_all)
  
  # ---------------------------
  # 7) Map true labels and predicted labels to internal labels
  # ---------------------------
  map_to_internal <- function(v, mapping = NULL) {
    v <- as.character(v)
    out <- v
    
    if (!is.null(mapping) && !is.null(names(mapping))) {
      idx <- match(v, names(mapping))
      ok <- !is.na(idx)
      out[ok] <- unname(mapping[idx[ok]])
    }
    
    as.character(out)
  }
  
  true_internal <- map_to_internal(true_raw, label_mapping)
  pred_internal <- map_to_internal(pred_raw, label_mapping)
  
  # ---------------------------
  # 8) Check unknown classes
  # ---------------------------
  unknown_true <- setdiff(unique(true_internal), internal_levels)
  unknown_pred <- setdiff(unique(pred_internal), internal_levels)
  
  if (length(unknown_true) > 0) {
    stop(
      "Y contains classes not covered by importance_order or classifier label mapping: ",
      paste(unknown_true, collapse = ", ")
    )
  }
  
  if (length(unknown_pred) > 0) {
    stop(
      "classifier produced classes not covered by importance_order or classifier label mapping: ",
      paste(unknown_pred, collapse = ", ")
    )
  }
  
  # ---------------------------
  # 9) Construct predictions and the confusion matrix
  # ---------------------------
  true_factor <- factor(true_internal, levels = internal_levels, labels = display_levels)
  pred_factor <- factor(pred_internal, levels = internal_levels, labels = display_levels)
  
  predictions <- data.frame(
    true_class = true_factor,
    predicted_class = pred_factor,
    stringsAsFactors = FALSE
  )
  
  sample_size <- min(5, nrow(predictions))
  
  if (sample_size > 0) {
    sample_index <- sample(seq_len(nrow(predictions)), size = sample_size)
    predictions_sample <- predictions[sample_index, , drop = FALSE]
  } else {
    predictions_sample <- predictions
  }
  
  conf_matrix <- table(True = true_factor, Predicted = pred_factor)
  
  # ---------------------------
  # 10) Compute class-wise error metrics
  # ---------------------------
  false_positive_rate <- rep(NA_real_, class_number)
  false_negative_rate <- rep(NA_real_, class_number)
  under_classification_error <- rep(NA_real_, class_number)
  over_classification_error  <- rep(NA_real_, class_number)
  
  for (class in seq_len(class_number)) {
    TP <- conf_matrix[class, class]
    FN <- sum(conf_matrix[class, ]) - TP
    FP <- sum(conf_matrix[, class]) - TP
    TN <- sum(conf_matrix) - TP - FN - FP
    
    denom_fpr <- FP + TN
    denom_fnr <- TP + FN
    denom_cls <- sum(conf_matrix[class, ])
    
    false_positive_rate[class] <- if (denom_fpr > 0) FP / denom_fpr else NA_real_
    false_negative_rate[class] <- if (denom_fnr > 0) FN / denom_fnr else NA_real_
    
    if (denom_cls > 0) {
      if (class < class_number) {
        under_classification_error[class] <-
          sum(conf_matrix[class, (class + 1):class_number, drop = FALSE]) / denom_cls
      } else {
        under_classification_error[class] <- 0
      }
      
      if (class > 1) {
        over_classification_error[class] <-
          sum(conf_matrix[class, 1:(class - 1), drop = FALSE]) / denom_cls
      } else {
        over_classification_error[class] <- 0
      }
    }
  }
  
  weight_proportion <- rowSums(conf_matrix) / sum(conf_matrix)
  
  total_under_classification_error <- sum(
    under_classification_error * weight_proportion,
    na.rm = TRUE
  )
  
  total_over_classification_error <- sum(
    over_classification_error * weight_proportion,
    na.rm = TRUE
  )
  
  correct_predictions <- sum(diag(conf_matrix))
  total_samples <- sum(conf_matrix)
  overall_accuracy <- if (total_samples > 0) correct_predictions / total_samples else NA_real_
  
  # ---------------------------
  # 11) Construct the row-normalized error_table
  # ---------------------------
  row_totals <- rowSums(conf_matrix)
  error_table <- conf_matrix
  
  for (i in seq_len(nrow(error_table))) {
    if (row_totals[i] > 0) {
      error_table[i, ] <- error_table[i, ] / row_totals[i]
    } else {
      error_table[i, ] <- NA_real_
    }
  }
  
  error_table <- cbind(error_table, false_negative_rate)
  colnames(error_table) <- c(
    paste0("Pred_", display_levels),
    "total_classification_error"
  )
  rownames(error_table) <- paste0("True_", display_levels)
  
  # ---------------------------
  # 12) Return results
  # ---------------------------
  list(
    confusion_matrix                 = conf_matrix,
    false_positive_rate              = false_positive_rate,
    false_negative_rate              = false_negative_rate,
    overall_accuracy                 = overall_accuracy,
    predictions                      = predictions,
    predictions_sample               = predictions_sample,
    under_classification_error       = under_classification_error,
    over_classification_error        = over_classification_error,
    remaining_error                  = total_over_classification_error,
    total_under_classification_error = total_under_classification_error,
    total_over_classification_error  = total_over_classification_error,
    error_table                      = error_table,
    class_levels                     = display_levels,
    internal_levels                  = internal_levels,
    label_mapping                    = label_mapping,
    importance_order                 = ordered_original,
    order_source                     = order_source
  )
}

#' @title HNP Under-classification Error Box Plot
#'
#' @description Summarize and visualize the under-classification errors and the
#'   overall misclassification error across repeated runs for a `T`-class
#'   problem. Given a list of confusion matrices (one per run), it draws grouped
#'   boxplots of the class-wise under-classification errors and the overall
#'   error. When a second list of confusion matrices is supplied, the two
#'   methods (e.g. a classical classifier vs. H-NP) are compared side by side.
#'   Optional control levels and tolerances add reference lines and
#'   `(1 - delta)` quantile markers.
#' @param conf_1 A list of `T x T` confusion matrices (one per run) for the
#'   first method.
#' @param conf_2 Optional list of `T x T` confusion matrices (one per run) for
#'   the second method, with the same number of runs as `conf_1`. If `NULL`,
#'   only `conf_1` is plotted (single mode).
#' @param levels Optional numeric vector of length `T - 1`. Control levels
#'   (alpha) per class; drawn as dashed reference lines. Each must lie in (0, 1).
#' @param tolerances Optional numeric vector of length `T - 1`. Tolerances
#'   (delta) per class; used to mark the `(1 - delta)` quantile as red points.
#'   Each must lie in (0, 1).
#' @param name_1 Character. Legend label for the first method.
#' @param name_2 Character. Legend label for the second method (used only when
#'   `conf_2` is provided).
#' @return Invisibly draws the boxplot and returns a list with two components:
#'   `classwise` (a data.frame of per-class control level, tolerance,
#'   under-classification error mean/sd, and violation rate) and `overall`
#'   (a data.frame of the overall misclassification error mean and sd).
#' @examples
#' set.seed(123)
#' make_cm <- function() {
#'   m <- matrix(sample(5:20, 9, replace = TRUE), 3, 3)
#'   dimnames(m) <- list(c("1", "2", "3"), c("1", "2", "3"))
#'   m
#' }
#' conf_classical <- replicate(20, make_cm(), simplify = FALSE)
#' conf_hnp <- replicate(20, make_cm(), simplify = FALSE)
#' tab <- hnp_boxplot(
#'   conf_1 = conf_classical,
#'   conf_2 = conf_hnp,
#'   levels = c(0.1, 0.1),
#'   tolerances = c(0.1, 0.1)
#' )
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_boxplot <- function(conf_1, 
                        conf_2 = NULL,
                        levels = NULL,
                        tolerances = NULL,   
                        name_1 = "Classical",
                        name_2 = "H-NP") {
  
  conf_before <- conf_1
  conf_after  <- conf_2
  single_mode <- is.null(conf_2)
  
  if (!single_mode && length(conf_before) != length(conf_after)) {
    stop("conf_1 and conf_2 must have the same number of runs.")
  }
  
  if (length(conf_before) == 0) {
    stop("conf_1 cannot be empty.")
  }
  
  T <- nrow(as.matrix(conf_before[[1]]))
  
  if (is.null(T) || T < 2) {
    stop("Cannot infer the number of classes from conf_1[[1]].")
  }
  
  if (ncol(as.matrix(conf_before[[1]])) != T) {
    stop("Each confusion matrix must be square.")
  }
  
  if (!is.null(levels)) {
    if (length(levels) != T - 1) {
      stop("levels must have length T - 1.")
    }
    if (any(!is.finite(levels)) || any(levels <= 0) || any(levels >= 1)) {
      stop("Each element of levels must lie in (0, 1).")
    }
  }
  
  if (!is.null(tolerances)) {
    if (length(tolerances) != T - 1) {
      stop("tolerances must have length T - 1.")
    }
    if (any(!is.finite(tolerances)) || any(tolerances <= 0) || any(tolerances >= 1)) {
      stop("Each element of tolerances must lie in (0, 1).")
    }
  }
  
  one_run_metrics <- function(CM) {
    CM <- as.matrix(CM)
    
    if (nrow(CM) != T || ncol(CM) != T) {
      stop("All confusion matrices must have the same T by T dimension.")
    }
    
    N <- sum(CM)
    
    if (N <= 0) {
      stop("Each confusion matrix must contain at least one observation.")
    }
    
    err <- 1 - sum(diag(CM)) / N
    under <- numeric(T)
    
    for (k in seq_len(T)) {
      nk <- sum(CM[k, ])
      
      under[k] <- if (nk > 0 && k < T) {
        sum(CM[k, (k + 1):T]) / nk
      } else {
        0
      }
    }
    
    list(
      err = err,
      under = under
    )
  }
  
  n_runs <- length(conf_before)
  
  err_before <- numeric(n_runs)
  under_before <- matrix(NA_real_, n_runs, T)
  
  if (!single_mode) {
    err_after <- numeric(n_runs)
    under_after <- matrix(NA_real_, n_runs, T)
  }
  
  for (r in seq_len(n_runs)) {
    mb <- one_run_metrics(conf_before[[r]])
    
    err_before[r] <- mb$err
    under_before[r, ] <- mb$under
    
    if (!single_mode) {
      ma <- one_run_metrics(conf_after[[r]])
      
      err_after[r] <- ma$err
      under_after[r, ] <- ma$under
    }
  }
  
  colnames(under_before) <- paste0("Class_", seq_len(T))
  
  if (!single_mode) {
    colnames(under_after) <- paste0("Class_", seq_len(T))
  }
  
  col_before <- "#edf8fb"
  col_after <- "#41ae76"
  
  get_single_col <- function() {
    if (grepl("hnp|h-np", name_1, ignore.case = TRUE)) {
      col_after
    } else {
      col_before
    }
  }
  
  draw_err_under_combo <- function() {
    
    n_under <- T - 1L
    m <- n_under + 1L
    
    metric_names <- c(as.character(seq_len(n_under)), "Overall Error")
    
    if (single_mode) {
      vals <- vector("list", length = m)
      names(vals) <- rep("", m)
      
      for (k in seq_len(n_under)) {
        vals[[k]] <- under_before[, k]
      }
      vals[[m]] <- err_before
      
      at_pos <- seq_len(m)
      group_centers <- at_pos
      box_cols <- rep(get_single_col(), m)
      xlim_combo <- c(0.5, m + 0.35)
      
    } else {
      vals <- vector("list", length = 2L * m)
      names(vals) <- rep("", 2L * m)
      
      idx <- 0L
      
      for (k in seq_len(n_under)) {
        idx <- idx + 1L
        vals[[idx]] <- under_before[, k]
        
        idx <- idx + 1L
        vals[[idx]] <- under_after[, k]
      }
      
      idx <- idx + 1L
      vals[[idx]] <- err_before
      
      idx <- idx + 1L
      vals[[idx]] <- err_after
      
      at_pos <- seq_len(2L * m)
      group_centers <- seq(1.5, 2 * m - 0.5, by = 2)
      box_cols <- rep(c(col_before, col_after), m)
      xlim_combo <- c(0.5, 2 * m + 0.35)
    }
    
    extra_y <- numeric(0)
    
    if (!is.null(levels)) {
      extra_y <- c(extra_y, levels)
    }
    
    if (!is.null(tolerances)) {
      qvals <- vapply(
        seq_len(n_under),
        function(k) {
          if (single_mode) {
            quantile(under_before[, k], probs = 1 - tolerances[k], na.rm = TRUE)
          } else {
            quantile(under_after[, k], probs = 1 - tolerances[k], na.rm = TRUE)
          }
        },
        numeric(1)
      )
      extra_y <- c(extra_y, qvals)
    }
    
    ylim_all <- range(c(unlist(vals, use.names = FALSE), extra_y), na.rm = TRUE)
    
    if (!is.finite(ylim_all[1]) || !is.finite(ylim_all[2])) {
      ylim_all <- c(0, 1)
    }
    
    boxplot(
      vals,
      at = at_pos,
      names = rep("", length(at_pos)),
      xlab = "",
      ylab = "Error",
      col = box_cols,
      border = "black",
      outline = FALSE,
      xaxt = "n",
      ylim = ylim_all,
      xlim = xlim_combo,
      cex.lab = 2.01,
      cex.axis = 2.01
    )
    
    grid(nx = NA, ny = NULL, lty = 3, col = "gray85")
    
    if (!single_mode && m >= 2) {
      abline(
        v = seq(2.5, 2 * m - 0.5, by = 2),
        lty = 3,
        col = "gray85",
        lwd = 1
      )
    }
    
    if (single_mode && m >= 2) {
      abline(
        v = seq(1.5, m - 0.5, by = 1),
        lty = 3,
        col = "gray85",
        lwd = 1
      )
    }
    
    abline(
      v = if (single_mode) n_under + 0.5 else 2 * n_under + 0.5,
      lty = 1,
      col = "black",
      lwd = 2.5
    )
    
    axis(
      1,
      at = group_centers,
      labels = metric_names,
      tick = FALSE,
      cex.axis = 2.01,
      line = 0
    )
    
    usr <- par("usr")
    under_centers <- if (single_mode) {
      seq_len(n_under)
    } else {
      seq(1.5, 2 * n_under - 0.5, by = 2)
    }
    
    text(
      x = mean(under_centers),
      y = usr[3] - 0.11 * (usr[4] - usr[3]),
      labels = "Under-classification Error",
      cex = 2.01,
      xpd = NA
    )
    
    for (k in seq_len(n_under)) {
      if (single_mode) {
        xC <- k
        
        if (!is.null(levels)) {
          segments(
            x0 = xC - 0.35,
            y0 = levels[k],
            x1 = xC + 0.35,
            y1 = levels[k],
            lty = 2,
            col = "gray55",
            lwd = 3
          )
        }
        
        if (!is.null(tolerances)) {
          q_val <- quantile(
            under_before[, k],
            probs = 1 - tolerances[k],
            na.rm = TRUE
          )
          
          points(xC, q_val, pch = 19, col = "red", cex = 1.4)
        }
        
      } else {
        xC <- 2 * k - 1
        xH <- 2 * k
        
        if (!is.null(levels)) {
          segments(
            x0 = xC - 0.45,
            y0 = levels[k],
            x1 = xH + 0.45,
            y1 = levels[k],
            lty = 2,
            col = "gray55",
            lwd = 3
          )
        }
        
        if (!is.null(tolerances)) {
          q_val_before <- quantile(
            under_before[, k],
            probs = 1 - tolerances[k],
            na.rm = TRUE
          )
          
          q_val_after <- quantile(
            under_after[, k],
            probs = 1 - tolerances[k],
            na.rm = TRUE
          )
          
          points(xC, q_val_before, pch = 19, col = "red", cex = 1.4)
          points(xH, q_val_after,  pch = 19, col = "red", cex = 1.4)
        }
      }
    }
    
    invisible(NULL)
  }
  
  draw_right_legend <- function() {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    par(mar = c(0, 0, 0, 0), xpd = NA)
    
    if (single_mode) {
      legend_names <- name_1
      legend_cols <- get_single_col()
    } else {
      legend_names <- c(name_1, name_2)
      legend_cols <- c(col_before, col_after)
    }
    
    legend(
      x = -0.1,
      y = 0.55,
      legend = legend_names,
      fill = legend_cols,
      border = "black",
      bty = "n",
      cex = 2.01,
      pt.cex = 3.2,
      x.intersp = 0.8,
      y.intersp = 1.25,
      text.width = max(strwidth(legend_names, cex = 2.01)),
      xjust = 0,
      yjust = 0.5
    )
    
    invisible(NULL)
  }
  
  compute_violation_rate <- function(mat) {
    if (is.null(levels)) {
      return(rep(NA_real_, T - 1))
    }
    
    out <- vapply(
      seq_len(T - 1),
      function(k) {
        mean(mat[, k] > levels[k], na.rm = TRUE)
      },
      numeric(1)
    )
    
    names(out) <- paste0("Class_", seq_len(T - 1))
    out
  }
  
  level_row <- if (!is.null(levels)) levels else rep(NA_real_, T - 1)
  tolerance_row <- if (!is.null(tolerances)) tolerances else rep(NA_real_, T - 1)
  
  if (single_mode) {
    classwise_table <- rbind(
      `control level` = level_row,
      tolerance = tolerance_row,
      `under-classification error mean` = colMeans(
        under_before[, seq_len(T - 1), drop = FALSE],
        na.rm = TRUE
      ),
      `under-classification error sd` = apply(
        under_before[, seq_len(T - 1), drop = FALSE],
        2,
        stats::sd,
        na.rm = TRUE
      ),
      `violation rate` = compute_violation_rate(under_before)
    )
    
    colnames(classwise_table) <- paste0(name_1, "_Class_", seq_len(T - 1))
    
    overall_table <- rbind(
      `overall misclassification error mean` = mean(err_before, na.rm = TRUE),
      `overall misclassification error sd` = stats::sd(err_before, na.rm = TRUE)
    )
    
    colnames(overall_table) <- name_1
    
  } else {
    classwise_table <- rbind(
      `control level` = c(level_row, level_row),
      tolerance = c(tolerance_row, tolerance_row),
      `under-classification error mean` = c(
        colMeans(
          under_before[, seq_len(T - 1), drop = FALSE],
          na.rm = TRUE
        ),
        colMeans(
          under_after[, seq_len(T - 1), drop = FALSE],
          na.rm = TRUE
        )
      ),
      `under-classification error sd` = c(
        apply(
          under_before[, seq_len(T - 1), drop = FALSE],
          2,
          stats::sd,
          na.rm = TRUE
        ),
        apply(
          under_after[, seq_len(T - 1), drop = FALSE],
          2,
          stats::sd,
          na.rm = TRUE
        )
      ),
      `violation rate` = c(
        compute_violation_rate(under_before),
        compute_violation_rate(under_after)
      )
    )
    
    colnames(classwise_table) <- c(
      paste0(name_1, "_Class_", seq_len(T - 1)),
      paste0(name_2, "_Class_", seq_len(T - 1))
    )
    
    overall_table <- rbind(
      `overall misclassification error mean` = c(
        mean(err_before, na.rm = TRUE),
        mean(err_after, na.rm = TRUE)
      ),
      `overall misclassification error sd` = c(
        stats::sd(err_before, na.rm = TRUE),
        stats::sd(err_after, na.rm = TRUE)
      )
    )
    
    colnames(overall_table) <- c(name_1, name_2)
  }
  
  classwise_table <- as.data.frame(classwise_table, check.names = FALSE)
  overall_table <- as.data.frame(overall_table, check.names = FALSE)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    layout(1)
    par(oldpar)
  }, add = TRUE)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.0, 1.1))
  
  par(mar = c(6.2, 4.4, 4.8, 0.6), xpd = FALSE)
  draw_err_under_combo()
  
  par(mar = c(6.2, 0.2, 4.8, 1.8), xpd = NA)
  draw_right_legend()
  
  out <- list(
    classwise = classwise_table,
    overall = overall_table
  )
  
  return(out)
}




gen_class_disk <- function(n, k, cfg) {
  cx <- cfg$centers[k, 1]
  cy <- cfg$centers[k, 2]
  R  <- cfg$radii[k]
  
  u <- runif(n)
  v <- runif(n)
  r <- R * sqrt(u)
  th <- 2 * pi * v
  
  x1 <- cx + r * cos(th)
  x2 <- cy + r * sin(th)
  
  if (cfg$noise_sd > 0) {
    x1 <- x1 + rnorm(n, 0, cfg$noise_sd)
    x2 <- x2 + rnorm(n, 0, cfg$noise_sd)
  }
  
  data.frame(x1 = x1, x2 = x2, y = as.character(k), stringsAsFactors = FALSE)
}


gen_data <- function(n_per_class, cfg) {
  S <- do.call(
    rbind,
    lapply(seq_len(cfg$T), function(k) gen_class_disk(n_per_class[k], k, cfg))
  )
  S$y <- factor(S$y, levels = as.character(seq_len(cfg$T)))
  S
}

#' @title Neural Network Classifier Helper Function
#'
#' @description Fit a single-hidden-layer softmax neural network
#'   (`nnet::nnet`) for multi-class classification and return the trained model
#'   together with the metadata needed to produce class scores. This can be used 
#'   to supply a pretrained scoring model to the HNP Umbrella pipeline.
#' @param X A data.frame or matrix of predictors/features.
#' @param Y A vector or factor of class labels of length `nrow(X)`, with at
#'   least two classes.
#' @return A list with components: `model` (the fitted `nnet` object),
#'   `class_levels` (the ordered class labels), and `feature_names` (the column
#'   names of `X`).
#' @examples
#' set.seed(123)
#' X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
#' Y <- factor(sample(c("1", "2", "3"), 60, replace = TRUE))
#' fit <- train_nn_and_get_scores(X, Y)
#' @export
train_nn_and_get_scores <- function(X, Y) {
  X <- as.data.frame(X)
  
  if (missing(Y) || is.null(Y)) {
    stop("Y must be provided.")
  }
  
  if (nrow(X) != length(Y)) {
    stop("X and Y have inconsistent numbers of rows.")
  }
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  
  Y <- as.factor(Y)
  y_levels <- levels(Y)
  Tn <- length(y_levels)
  
  if (Tn < 2) {
    stop("Need at least 2 classes to train the neural network.")
  }
  
  Y_mat <- nnet::class.ind(Y)
  colnames(Y_mat) <- y_levels
  
  nn <- nnet::nnet(
    x = as.matrix(X),
    y = Y_mat,
    size = 8,
    softmax = TRUE,
    maxit = 200,
    decay = 5e-4,
    trace = FALSE
  )
  
  list(
    model = nn,
    class_levels = y_levels,
    feature_names = colnames(X)
  )
}


gen_ball <- function(n, center, R) {
  Z <- matrix(rnorm(n * 3), n, 3)
  U <- Z / sqrt(rowSums(Z^2))
  r <- R * runif(n)^(1/3)
  sweep(sweep(U, 1, r, `*`), 2, center, `+`)
}


#' @title Generate Synthetic Three-class Ball Data
#'
#' @description Generate a synthetic three-class dataset in which each class is
#'   sampled uniformly from a 3-dimensional ball with a given center and radius.
#'   Useful for demonstrating and testing the HNP Umbrella pipeline.
#' @param n Integer. Number of samples to draw for each class.
#' @param centers A list of three numeric vectors of length 3, giving the
#'   center coordinates of the balls for classes "A", "B", and "C".
#' @param radii Numeric vector of length 3 giving the radius of each class ball.
#' @return A data.frame with feature columns `x1`, `x2`, `x3` and a factor
#'   column `y` with levels c("A","B","C"), containing `3 * n` rows.
#' @examples
#' set.seed(123)
#' data <- generate_ball_data(
#'   n = 100,
#'   centers = list(c(0, 0, 0), c(2, 2, 2), c(4, 0, 0)),
#'   radii = c(1, 1, 1)
#' )
#' @export
generate_ball_data <- function(n, centers, radii) {
  X1 <- gen_ball(n, centers[[1]], radii[1])
  X2 <- gen_ball(n, centers[[2]], radii[2])
  X3 <- gen_ball(n, centers[[3]], radii[3])
  data <- as.data.frame(rbind(X1, X2, X3))
  colnames(data) <- paste0("x", 1:3)
  data$y <- factor(c(rep("A", n), rep("B", n), rep("C", n)),
                   levels = c("A", "B", "C"))
  return(data)
}


gen_class_normal <- function(n, k, cfg) {
  X <- MASS::mvrnorm(n = n, mu = cfg$means[k, ], Sigma = cfg$Sigma_list[[k]])
  X <- as.data.frame(X)
  colnames(X) <- c("x1", "x2")
  X$y <- as.character(k)
  X
}


gen_normal_data <- function(n_per_class, cfg) {
  S <- do.call(
    rbind,
    lapply(seq_len(cfg$T), function(k) gen_class_normal(n_per_class[k], k, cfg))
  )
  S$y <- factor(S$y, levels = as.character(seq_len(cfg$T)))
  S
}

hnp_read_split_value <- function(cfg, name, default = 0) {
  if (!is.null(cfg) && !is.null(names(cfg)) && (name %in% names(cfg))) {
    v <- as.numeric(cfg[[name]])
    if (length(v) == 1 && !is.na(v)) {
      return(v)
    }
  }
  default
}


hnp_split_one_class <- function(Si, cfg, class_i) {
  ni <- nrow(Si)
  if (ni <= 1) {
    stop("Class ", class_i, " has too few samples: ", ni)
  }
  
  idx <- sample.int(ni, ni)
  
  r_train <- hnp_read_split_value(cfg, "train", 0)
  r_thre  <- hnp_read_split_value(cfg, "threshold", 0)
  r_err   <- hnp_read_split_value(cfg, "error", NA_real_)
  
  if (r_train + r_thre > 1) {
    stop("For each class, train + threshold must be <= 1.")
  }
  
  if (!is.na(r_err) && r_train + r_thre + r_err > 1) {
    stop("For each class, train + threshold + error must be <= 1.")
  }
  
  r_train <- max(0, min(1, r_train))
  r_thre  <- max(0, min(1, r_thre))
  
  if (is.na(r_err)) {
    r_err <- max(0, 1 - r_train - r_thre)
  } else {
    r_err <- max(0, min(1, r_err))
  }
  
  if (r_train + r_thre > 1) {
    r_thre <- max(0, 1 - r_train)
  }
  
  n_train <- floor(r_train * ni)
  n_thre_end <- floor((r_train + r_thre) * ni)
  
  n_train <- max(0, min(ni, n_train))
  n_thre_end <- max(n_train, min(ni, n_thre_end))
  
  train_idx <- if (n_train >= 1) idx[1:n_train] else integer(0)
  thre_idx  <- if (n_thre_end > n_train) idx[(n_train + 1):n_thre_end] else integer(0)
  err_idx   <- if (ni > n_thre_end) idx[(n_thre_end + 1):ni] else integer(0)
  
  list(
    Ss = Si[train_idx, , drop = FALSE],
    St = Si[thre_idx,  , drop = FALSE],
    Se = Si[err_idx,   , drop = FALSE]
  )
}


hnp_predict_with_thresholds <- function(X, t_vec, score_functions, feature_names, Tn) {
  X <- as.data.frame(X)
  X <- X[, feature_names, drop = FALSE]
  
  if (nrow(X) == 0) {
    return(integer(0))
  }
  
  scores <- vector("list", Tn - 1)
  
  for (i in 1:(Tn - 1)) {
    s <- as.numeric(score_functions[[i]](X))
    s[is.na(s)] <- -Inf
    scores[[i]] <- s
  }
  
  decision <- rep.int(Tn, nrow(X))
  
  for (i in seq.int(Tn - 1, 1)) {
    idx <- which(scores[[i]] >= t_vec[i])
    if (length(idx) > 0) {
      decision[idx] <- i
    }
  }
  
  as.integer(decision)
}


hnp_eval_objective_from_thresholds <- function(t_vec, by_class, feature_names, Tn,
                                               pi_hat, score_functions) {
  obj <- 0
  
  for (i in 2:Tn) {
    Se_i <- by_class[[i]]$Se
    if (nrow(Se_i) == 0) {
      next
    }
    
    preds <- hnp_predict_with_thresholds(
      X = Se_i[, feature_names, drop = FALSE],
      t_vec = t_vec,
      score_functions = score_functions,
      feature_names = feature_names,
      Tn = Tn
    )
    
    e_i <- mean(preds < i)
    obj <- obj + pi_hat[i] * e_i
  }
  
  obj
}


hnp_search_thresholds_rec <- function(i, prev_t, Tn, by_class, feature_names,
                                      levels, tolerances, score_functions,
                                      grid_set, max_grid, max_combinations,
                                      pi_hat, verbose, state_env) {
  if (state_env$comb_count >= max_combinations) {
    return(invisible(NULL))
  }
  
  # Last layer: directly use hnp_upper_bound to compute t_{T-1}
  if (i == (Tn - 1)) {
    St_last <- by_class[[Tn - 1]]$St
    if (nrow(St_last) == 0) {
      stop("Threshold set for class ", Tn - 1, " is empty.")
    }
    
    t_last <- hnp_upper_bound(
      S_it = St_last[, feature_names, drop = FALSE],
      level = levels[Tn - 1],
      delta_i = tolerances[Tn - 1],
      score_functions = score_functions,
      thresholds = prev_t,
      i = (Tn - 1)
    )
    
    t_vec <- c(prev_t, t_last)
    
    obj <- hnp_eval_objective_from_thresholds(
      t_vec = t_vec,
      by_class = by_class,
      feature_names = feature_names,
      Tn = Tn,
      pi_hat = pi_hat,
      score_functions = score_functions
    )
    
    if (!is.finite(obj)) {
      return(invisible(NULL))
    }
    
    if (obj < state_env$best_obj) {
      state_env$best_obj <- obj
      state_env$best_thresholds <- t_vec
      
      if (isTRUE(verbose)) {
        cat(sprintf("[update best] obj=%.6f, thresholds=%s\n",
                    obj, paste(round(t_vec, 6), collapse = ",")))
      }
    }
    
    return(invisible(NULL))
  }
  
  # Non-last layers: first compute the upper bound, then enumerate candidates
  St_i <- by_class[[i]]$St
  if (nrow(St_i) == 0) {
    stop("Threshold set for class ", i, " is empty.")
  }
  
  tbar_i <- hnp_upper_bound(
    S_it = St_i[, feature_names, drop = FALSE],
    level = levels[i],
    delta_i = tolerances[i],
    score_functions = score_functions,
    thresholds = if (i == 1) NULL else prev_t,
    i = i
  )
  
  cand_i <- hnp_get_candidates(
    i = i,
    tbar = tbar_i,
    grid_set = grid_set,
    by_class = by_class,
    score_functions = score_functions,
    feature_names = feature_names,
    max_grid = max_grid
  )
  
  if (length(cand_i) == 0) {
    return(invisible(NULL))
  }
  
  for (k in seq_along(cand_i)) {
    if (state_env$comb_count >= max_combinations) {
      break
    }
    
    ti <- cand_i[k]
    state_env$comb_count <- state_env$comb_count + 1L
    
    hnp_search_thresholds_rec(
      i = i + 1,
      prev_t = c(prev_t, ti),
      Tn = Tn,
      by_class = by_class,
      feature_names = feature_names,
      levels = levels,
      tolerances = tolerances,
      score_functions = score_functions,
      grid_set = grid_set,
      max_grid = max_grid,
      max_combinations = max_combinations,
      pi_hat = pi_hat,
      verbose = verbose,
      state_env = state_env
    )
  }
  
  invisible(NULL)
}


hnp_map_labels_by_importance <- function(Y, importance_order) {
  if (is.null(Y)) {
    stop("Y must be provided in the current implementation.")
  }
  
  y_char <- as.character(Y)
  importance_order <- as.character(importance_order)
  
  if (any(is.na(y_char))) {
    stop("Y contains NA. Please handle it first.")
  }
  
  if (length(importance_order) < 2) {
    stop("importance_order must contain at least two classes.")
  }
  
  if (any(is.na(importance_order)) || any(importance_order == "")) {
    stop("importance_order must contain non-empty class labels.")
  }
  
  if (anyDuplicated(importance_order)) {
    stop("importance_order contains duplicated class labels: ",
         paste(unique(importance_order[duplicated(importance_order)]), collapse = ", "))
  }
  
  missing_in_order <- setdiff(unique(y_char), importance_order)
  if (length(missing_in_order) > 0) {
    stop("The following labels in Y are not covered by importance_order: ",
         paste(missing_in_order, collapse = ", "))
  }
  
  missing_in_y <- setdiff(importance_order, unique(y_char))
  if (length(missing_in_y) > 0) {
    stop("The following classes in importance_order do not appear in Y: ",
         paste(missing_in_y, collapse = ", "))
  }
  
  mapped <- match(y_char, importance_order)
  
  factor(
    as.character(mapped),
    levels = as.character(seq_along(importance_order))
  )
}

hnp_align_score_input <- function(X, importance_order) {
  X <- as.data.frame(X)
  importance_order <- as.character(importance_order)
  Tn <- length(importance_order)
  
  if (ncol(X) != Tn) {
    stop("When input_is_score = TRUE, X must have exactly T columns. ",
         "Need ", Tn, " columns, got ", ncol(X), ".")
  }
  
  # Case 1: the score matrix has original class-name columns, e.g., A, B, C.
  # Reorder columns automatically according to importance_order.
  if (!is.null(colnames(X)) && all(importance_order %in% colnames(X))) {
    X <- X[, importance_order, drop = FALSE]
  } else {
    # Case 2: the score matrix has no usable class-name columns.
    # In this case, importance_order may be integer column indices, e.g., c(3, 1, 2).
    idx <- suppressWarnings(as.integer(importance_order))
    
    if (any(is.na(idx))) {
      stop(
        "When input_is_score = TRUE and X does not contain class-name columns, ",
        "importance_order must be coercible to integer column indices."
      )
    }
    
    if (!all(idx %in% seq_len(ncol(X)))) {
      stop("Some entries of importance_order are invalid column indices for X.")
    }
    
    X <- X[, idx, drop = FALSE]
  }
  
  X <- as.data.frame(lapply(X, as.numeric))
  
  if (any(vapply(X, function(z) any(is.na(z)), logical(1)))) {
    stop("Score matrix X contains NA after alignment.")
  }
  
  colnames(X) <- as.character(seq_len(Tn))
  X
}

hnp_align_probability_matrix <- function(prob,
                                         required_levels,
                                         original_levels,
                                         use_pretrained,
                                         context = "pretrained_model") {
  ## Convert model output to an internal probability matrix:
  ## column order is required_levels = "1", ..., "T";
  ## for user-supplied pretrained models, prefer matching columns by original_levels = importance_order.
  
  if (is.null(dim(prob))) {
    prob <- t(as.matrix(prob))
  } else {
    prob <- as.matrix(prob)
  }
  
  required_levels <- as.character(required_levels)
  original_levels <- as.character(original_levels)
  Tn <- length(required_levels)
  
  if (ncol(prob) < Tn) {
    stop(context, " output has too few columns: need ", Tn,
         ", got ", ncol(prob), ".")
  }
  
  cn <- colnames(prob)
  
  if (isTRUE(use_pretrained)) {
    ## User-defined model: first expect columns named by original labels, e.g., A, B, C.
    if (!is.null(cn) && all(original_levels %in% cn)) {
      prob <- prob[, original_levels, drop = FALSE]
      colnames(prob) <- required_levels
      
    } else if (!is.null(cn) && all(required_levels %in% cn)) {
      ## Compatibility for advanced users: internal column names 1, 2, 3 are also allowed.
      prob <- prob[, required_levels, drop = FALSE]
      
    } else if (is.null(cn) && ncol(prob) == Tn) {
      ## If there are no column names, assume columns are already ordered by importance_order.
      prob <- prob[, seq_len(Tn), drop = FALSE]
      colnames(prob) <- required_levels
      
    } else {
      stop(
        context, " output columns cannot be matched to class labels.\n",
        "Expected columns named by original labels in importance_order: ",
        paste(original_levels, collapse = ", "), "\n",
        "or internal labels: ", paste(required_levels, collapse = ", "), ".\n",
        "Observed columns: ",
        if (is.null(cn)) "NULL" else paste(cn, collapse = ", ")
      )
    }
    
  } else {
    ## Internal base learner: training labels are already 1, ..., T.
    if (is.null(cn)) {
      if (ncol(prob) == Tn) {
        colnames(prob) <- required_levels
      } else {
        stop(context, " output has no column names and has wrong number of columns.")
      }
    }
    
    cn <- colnames(prob)
    
    if (all(required_levels %in% cn)) {
      prob <- prob[, required_levels, drop = FALSE]
    } else {
      stop(
        context, " output columns cannot be matched to internal class labels.\n",
        "Expected columns: ", paste(required_levels, collapse = ", "), "\n",
        "Observed columns: ", paste(cn, collapse = ", ")
      )
    }
  }
  
  prob <- as.matrix(prob)
  storage.mode(prob) <- "numeric"
  prob[is.na(prob)] <- 0
  
  if (any(!is.finite(prob))) {
    prob[!is.finite(prob)] <- 0
  }
  
  prob
}


hnp_predict_proba <- function(model, X, feature_names, required_levels,
                              use_pretrained, method,
                              original_levels = NULL) {
  X <- as.data.frame(X)
  X <- X[, feature_names, drop = FALSE]
  
  required_levels <- as.character(required_levels)
  
  if (is.null(original_levels)) {
    original_levels <- required_levels
  } else {
    original_levels <- as.character(original_levels)
  }
  
  prob <- NULL
  
  if (isTRUE(use_pretrained)) {
    ## Case 1: pretrained_model is a user-defined function.
    ## It should take raw feature data as input and return an n x T score/probability matrix.
    if (is.function(model)) {
      prob <- model(X)
      
    } else {
      ## Case 2: pretrained_model is a fitted model object.
      ## Try common predict interfaces in sequence and strictly check whether the output is an n x T numeric matrix.
      
      check_probability_matrix <- function(prob, source_name) {
        
        if (is.null(prob)) {
          return(NULL)
        }
        
        ## Some predict() methods return a list; here we only accept common probability fields.
        if (is.list(prob) && !is.data.frame(prob)) {
          if (!is.null(prob$prob)) {
            prob <- prob$prob
          } else if (!is.null(prob$probs)) {
            prob <- prob$probs
          } else if (!is.null(prob$probabilities)) {
            prob <- prob$probabilities
          } else {
            return(NULL)
          }
        }
        
        ## Only matrix or data.frame outputs are accepted.
        if (!(is.matrix(prob) || is.data.frame(prob))) {
          return(NULL)
        }
        
        prob <- as.data.frame(prob, check.names = FALSE)
        
        n_required <- nrow(X)
        T_required <- length(required_levels)
        
        ## The number of rows must equal the input sample size n.
        if (nrow(prob) != n_required) {
          return(NULL)
        }
        
        prob_names <- colnames(prob)
        
        ## Case 1: column names exist and contain all required_levels.
        if (!is.null(prob_names) && all(original_levels %in% prob_names)) {
          prob <- prob[, original_levels, drop = FALSE]
          colnames(prob) <- required_levels
          
        } else if (!is.null(prob_names) && all(required_levels %in% prob_names)) {
          prob <- prob[, required_levels, drop = FALSE]
        } else if ((is.null(prob_names) || all(is.na(prob_names)) || all(prob_names == "")) &&
                   ncol(prob) == T_required) {
          colnames(prob) <- required_levels
          
          ## All other cases are treated as invalid outputs.
        } else {
          return(NULL)
        }
        
        ## After alignment, the output must be exactly n x T.
        if (nrow(prob) != n_required || ncol(prob) != T_required) {
          return(NULL)
        }
        
        ## The output must be convertible to a numeric matrix.
        prob_mat <- suppressWarnings(as.matrix(prob))
        
        if (!is.numeric(prob_mat)) {
          return(NULL)
        }
        
        storage.mode(prob_mat) <- "double"
        
        ## NA, NaN, and Inf are not allowed.
        if (anyNA(prob_mat) || any(!is.finite(prob_mat))) {
          return(NULL)
        }
        
        ## Negative values are not allowed.
        if (any(prob_mat < 0)) {
          return(NULL)
        }
        
        ## Do not check row sums or normalize rows.
        prob <- as.data.frame(prob_mat, check.names = FALSE)
        colnames(prob) <- required_levels
        
        prob
      }
      
      
      try_probability_output <- function(expr, source_name) {
        prob_try <- tryCatch(
          expr,
          error = function(e) NULL
        )
        
        check_probability_matrix(prob_try, source_name)
      }
      
      
      prob <- try_probability_output(
        predict(model, newdata = X, type = "prob"),
        "predict(type = 'prob')"
      )
      
      if (is.null(prob)) {
        prob <- try_probability_output(
          predict(model, newdata = X, type = "probs"),
          "predict(type = 'probs')"
        )
      }
      
      if (is.null(prob)) {
        prob <- try_probability_output(
          predict(model, newdata = X, type = "raw"),
          "predict(type = 'raw')"
        )
      }
      
      if (is.null(prob)) {
        prob <- try_probability_output({
          pred <- predict(model, newdata = X, probability = TRUE)
          attr(pred, "probabilities")
        }, "predict(probability = TRUE)")
      }
      
      if (is.null(prob)) {
        stop(
          "pretrained_model cannot produce a valid n x T numeric output matrix. ",
          "Expected n = ", nrow(X), ", T = ", length(required_levels), ". ",
          "Required columns are: ",
          paste(required_levels, collapse = ", "), ". ",
          "Tried predict(type = 'prob'), predict(type = 'probs'), ",
          "predict(type = 'raw'), and predict(probability = TRUE).",
          call. = FALSE
        )
      }
    }
    
    if (is.null(prob)) {
      stop("pretrained_model cannot produce class probability scores.")
    }
    
    prob <- hnp_align_probability_matrix(
      prob = prob,
      required_levels = required_levels,
      original_levels = original_levels,
      use_pretrained = TRUE,
      context = "pretrained_model"
    )
    
  } else {
    ## Internally trained base learner.
    if (method == "svm") {
      pred <- predict(model, newdata = X, probability = TRUE)
      prob <- attr(pred, "probabilities")
    } else if (method == "randomforest") {
      prob <- predict(model, newdata = X, type = "prob")
    } else if (method == "logistic") {
      prob <- predict(model, newdata = X, type = "probs")
    } else {
      stop("Unknown method: ", method)
    }
    
    prob <- hnp_align_probability_matrix(
      prob = prob,
      required_levels = required_levels,
      original_levels = original_levels,
      use_pretrained = FALSE,
      context = paste0("base learner: ", method)
    )
  }
  
  if (nrow(prob) != nrow(X)) {
    stop("Probability matrix row count mismatch.")
  }
  
  prob
}


hnp_make_score_function <- function(score_index,
                                    trained_model,
                                    feature_names,
                                    required_levels,
                                    original_levels,
                                    use_pretrained,
                                    method,
                                    Tn,
                                    eps = .Machine$double.eps) {
  force(score_index)
  force(trained_model)
  force(feature_names)
  force(required_levels)
  force(original_levels)
  force(use_pretrained)
  force(method)
  force(Tn)
  force(eps)
  
  ## First layer: T_1(x) = score_{class 1}(x)
  if (score_index == 1) {
    return(function(X) {
      prob <- hnp_predict_proba(
        model = trained_model,
        X = X,
        feature_names = feature_names,
        required_levels = required_levels,
        original_levels = original_levels,
        use_pretrained = use_pretrained,
        method = method
      )
      
      as.numeric(prob[, "1"])
    })
  }
  
  ## Later layers:
  ## T_i(x) = score_{class i}(x) / sum_{j=i+1}^T score_{class j}(x)
  ii <- score_index
  
  function(X) {
    prob <- hnp_predict_proba(
      model = trained_model,
      X = X,
      feature_names = feature_names,
      required_levels = required_levels,
      original_levels = original_levels,
      use_pretrained = use_pretrained,
      method = method
    )
    
    p_num <- as.numeric(prob[, as.character(ii)])
    
    denom_cols <- as.character((ii + 1):Tn)
    p_den <- rowSums(prob[, denom_cols, drop = FALSE])
    
    p_den[!is.finite(p_den) | p_den <= 0 | is.na(p_den)] <- eps
    
    out <- p_num / p_den
    out[!is.finite(out)] <- 0
    
    as.numeric(out)
  }
}


hnp_build_score_functions <- function(Tn,
                                      trained_model,
                                      feature_names,
                                      required_levels,
                                      original_levels,
                                      use_pretrained,
                                      method,
                                      eps = .Machine$double.eps) {
  score_functions <- vector("list", Tn - 1)
  
  for (i in 1:(Tn - 1)) {
    score_functions[[i]] <- hnp_make_score_function(
      score_index = i,
      trained_model = trained_model,
      feature_names = feature_names,
      required_levels = required_levels,
      original_levels = original_levels,
      use_pretrained = use_pretrained,
      method = method,
      Tn = Tn,
      eps = eps
    )
  }
  
  score_functions
}


hnp_build_score_functions_from_scores <- function(Tn,
                                                  feature_names = as.character(seq_len(Tn)),
                                                  eps = .Machine$double.eps) {
  ## When input_is_score = TRUE, X has already been processed by hnp_align_score_input().
  ## It has been reordered by importance_order and renamed to internal column names 1, ..., T.
  
  score_functions <- vector("list", Tn - 1)
  
  score_functions[[1]] <- function(X) {
    X <- as.data.frame(X)
    X <- X[, feature_names, drop = FALSE]
    as.numeric(X[, "1", drop = TRUE])
  }
  
  if (Tn >= 3) {
    for (i in 2:(Tn - 1)) {
      score_functions[[i]] <- local({
        ii <- i
        function(X) {
          X <- as.data.frame(X)
          X <- X[, feature_names, drop = FALSE]
          
          p_i <- as.numeric(X[, as.character(ii), drop = TRUE])
          denom_cols <- as.character((ii + 1):Tn)
          denom <- rowSums(as.matrix(X[, denom_cols, drop = FALSE]))
          
          denom[!is.finite(denom) | denom <= 0 | is.na(denom)] <- eps
          
          out <- p_i / denom
          out[!is.finite(out)] <- 0
          
          as.numeric(out)
        }
      })
    }
  }
  
  score_functions
}


hnp_build_output_classifier <- function(best_thresholds,
                                        score_functions,
                                        feature_names,
                                        Tn,
                                        input_is_score = FALSE,
                                        importance_order = NULL) {
  force(best_thresholds)
  force(score_functions)
  force(feature_names)
  force(Tn)
  force(input_is_score)
  force(importance_order)
  
  function(new_data) {
    new_data <- as.data.frame(new_data)
    
    if (isTRUE(input_is_score)) {
      new_data <- hnp_align_score_input(new_data, importance_order)
    } else {
      new_data <- new_data[, feature_names, drop = FALSE]
    }
    
    hnp_predict_with_thresholds(
      X = new_data,
      t_vec = best_thresholds,
      score_functions = score_functions,
      feature_names = feature_names,
      Tn = Tn
    )
  }
}


hnp_wrap_classifier_labels <- function(internal_clf, label_mapping) {
  ## label_mapping has the form:
  ##   C   A   B
  ##  "1" "2" "3"
  inverse_mapping <- stats::setNames(
    names(label_mapping),
    as.character(label_mapping)
  )
  
  out_clf <- function(newX, output_internal = FALSE) {
    pred_internal <- internal_clf(newX)
    pred_internal_chr <- as.character(pred_internal)
    
    if (isTRUE(output_internal)) {
      return(pred_internal_chr)
    }
    
    pred_original <- unname(inverse_mapping[pred_internal_chr])
    
    if (any(is.na(pred_original))) {
      warning("Some internal predicted labels cannot be mapped back to original labels.")
    }
    
    pred_original
  }
  
  attr(out_clf, "internal_classifier") <- internal_clf
  attr(out_clf, "label_mapping") <- label_mapping
  attr(out_clf, "inverse_label_mapping") <- inverse_mapping
  
  out_clf
}


hnp_get_candidates <- function(i, tbar, grid_set, by_class, score_functions,
                               feature_names, max_grid) {
  Ai <- NULL
  
  if (!is.null(grid_set) && length(grid_set) >= i) {
    Ai <- grid_set[[i]]
  }
  
  St_i <- by_class[[i]]$St
  
  if (nrow(St_i) == 0) {
    stop("Threshold set for class ", i, " is empty; cannot compute thresholds.")
  }
  
  ## If the user does not provide A_i, use the score values on the threshold set as the default candidate grid.
  if (is.null(Ai)) {
    sc <- as.numeric(score_functions[[i]](St_i[, feature_names, drop = FALSE]))
    cand <- sc
  } else {
    cand <- as.numeric(Ai)
  }
  
  ## feasible set: F_i = A_i cap (-Inf, tbar]
  cand <- cand[is.finite(cand) & !is.na(cand) & cand <= tbar]
  cand <- sort(unique(cand))
  
  ## If the feasible set is too large, randomly sample max_grid candidate points.
  if (length(cand) > max_grid) {
    cand <- sort(sample(cand, size = max_grid, replace = FALSE))
  }
  
  cand
}


#' @title HNP Umbrella Algorithm
#'
#' @description Implementation of the Hierarchical Neyman-Pearson (HNP) Umbrella
#'   algorithm for general `T`-class classification (`T >= 2`). Classes are
#'   ordered by `importance_order` and re-labelled internally to the classes "1", ..., "T".
#'   The algorithm trains (or reuses) a scoring-type classifier (e.g. random forest, support 
#'   vector machine, or logistic regression) and selects thresholds that control the
#'   under-classification error of the more severe classes at the requested
#'   levels with the given tolerances, either by treating upper bounds as thresholds or by
#'   a grid search that minimizes the weighted remaining error.
#' @param X A data.frame of features, or an `n x T` score matrix when
#'   `input_is_score = TRUE`.
#' @param Y A vector of class labels of length `nrow(X)`. All labels must be
#'   covered by `importance_order`.
#' @param levels Numeric vector of length `T - 1`. Control levels (alpha) for
#'   the under-classification error of classes 1, ..., T-1. Each must lie in (0, 1).
#' @param tolerances Numeric vector of length `T - 1`. Tolerance (delta) values
#'   for the confidence bounds. Each must lie in (0, 1).
#' @param importance_order Character vector giving the class priority order,
#'   most severe first. Defines the mapping to internal levels "1", ..., "T".
#' @param method Character string: one of 'randomforest', 'svm', or 'logistic'.
#'   The base method used when training a scoring-type classifier internally.
#' @param pretrained_model Optional. A user-supplied scoring model: a fitted
#'   model object, a function returning an `n x T` score/probability matrix, or
#'   a list containing `score_functions`. When provided, no internal model is
#'   trained.
#' @param input_is_score Logical. If `TRUE`, `X` is treated as a precomputed
#'   `n x T` score matrix and `pretrained_model` is ignored.
#' @param grid_search Logical. If `TRUE`, search candidate thresholds to
#'   minimize the weighted remaining error; if `FALSE`, use the upper-bound as
#'   thresholds directly.
#' @param grid_set Optional list of length `T - 1` of candidate threshold
#'   vectors per class. If `NULL`, candidates are derived from the threshold set.
#' @param max_grid Integer. Maximum number of candidate thresholds considered
#'   per class in the grid search.
#' @param max_combinations Integer. Maximum number of threshold combinations
#'   explored in the grid search.
#' @param hnp_split Optional list of length `T` specifying the class-wise data
#'   split ratios (train / threshold / error). If `NULL`, sensible defaults are
#'   used based on `grid_search` and whether scores are already available.
#' @param verbose Logical. If `TRUE`, print progress information.
#' @return A HNP Umbrella classifier function `function(new_data, output_internal = FALSE)`
#'   that returns predicted class labels (original labels by default, internal
#'   labels if `output_internal = TRUE`), with attributes recording the selected
#'   thresholds, objective value, `importance_order`, label mapping, and other
#'   metadata. Returns `NULL` if no feasible classifier is found.
#' @examples
#' set.seed(123)
#' n <- 500
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' Y <- sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.2, 0.3, 0.5))
#' clf <- hnp_umbrella(
#'   X = X, Y = Y,
#'   levels = c(0.1, 0.1),
#'   tolerances = c(0.1, 0.1),
#'   importance_order = c("A", "B", "C"),
#'   method = "logistic"
#' )
#' preds <- clf(X)
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_umbrella <- function(X, Y, levels, tolerances,
                         importance_order,
                         method = "logistic",
                         pretrained_model = NULL,
                         input_is_score = FALSE,
                         grid_search = TRUE,
                         grid_set = NULL,
                         max_grid = 15,
                         max_combinations = 2000,
                         hnp_split = NULL,
                         verbose = FALSE) {
  
  # ---------------------------
  # 0) Basic checks
  # ---------------------------
  method_choices <- c("randomforest", "svm", "logistic")
  method <- match.arg(method, method_choices)
  
  if (missing(importance_order) || is.null(importance_order)) {
    stop("importance_order must be provided.")
  }
  
  if (isTRUE(input_is_score) && !is.null(pretrained_model)) {
    warning("scores in X are used directly and pretrained_model is ignored")
    pretrained_model <- NULL
  }
  
  importance_order <- as.character(importance_order)
  Tn <- length(importance_order)
  required_levels <- as.character(seq_len(Tn))
  original_levels <- importance_order
  
  if (Tn < 2) {
    stop("importance_order must contain at least two classes.")
  }
  
  if (length(levels) != (Tn - 1)) {
    stop("levels must have length T-1.")
  }
  
  if (length(tolerances) != (Tn - 1)) {
    stop("tolerances must have length T-1.")
  }
  
  if (any(!is.finite(levels)) || any(levels <= 0) || any(levels >= 1)) {
    stop("Each element of levels must lie in (0, 1).")
  }
  
  if (any(!is.finite(tolerances)) || any(tolerances <= 0) || any(tolerances >= 1)) {
    stop("Each element of tolerances must lie in (0, 1).")
  }
  
  if (!is.finite(max_grid) || max_grid <= 0) {
    stop("max_grid must be a positive integer.")
  }
  
  if (!is.finite(max_combinations) || max_combinations <= 0) {
    stop("max_combinations must be a positive integer.")
  }
  
  max_grid <- as.integer(max_grid)
  max_combinations <- as.integer(max_combinations)
  
  
  if (!is.numeric(max_grid) || length(max_grid) != 1 ||
      !is.finite(max_grid) || max_grid < 1 || max_grid != floor(max_grid)) {
    stop("max_grid must be a positive integer.")
  }
  
  if (!is.numeric(max_combinations) || length(max_combinations) != 1 ||
      !is.finite(max_combinations) || max_combinations < 1 ||
      max_combinations != floor(max_combinations)) {
    stop("max_combinations must be a positive integer.")
  }
  
  # ---------------------------
  # 0.1) Grid size control
  #      Used only when grid_search = TRUE
  # ---------------------------
  if (isTRUE(grid_search) && Tn > 2) {
    effective_grid_cap <- floor(max_combinations^(1 / (Tn - 2)))
    effective_grid_cap <- max(1L, as.integer(effective_grid_cap))
    max_grid_effective <- min(max_grid, effective_grid_cap)
  } else {
    max_grid_effective <- max_grid
  }
  
  if (isTRUE(verbose)) {
    cat("[grid control] max_grid =", max_grid,
        ", max_combinations =", max_combinations,
        ", effective max_grid =", max_grid_effective, "\n")
  }
  
  # ---------------------------
  # 1) Construct the internal data frame S and internal label .y_internal
  # ---------------------------
  y_internal <- hnp_map_labels_by_importance(Y, importance_order)
  
  label_mapping <- stats::setNames(
    required_levels,
    importance_order
  )
  
  class_levels <- importance_order
  
  if (!isTRUE(input_is_score)) {
    X <- as.data.frame(X)
    
    if (nrow(X) != length(y_internal)) {
      stop("X and Y have inconsistent numbers of rows.")
    }
    
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("x", seq_len(ncol(X)))
    }
    
    feature_names <- colnames(X)
    X_internal <- X
    
  } else {
    X_internal <- hnp_align_score_input(X, importance_order)
    
    if (nrow(X_internal) != length(y_internal)) {
      stop("X and Y have inconsistent numbers of rows.")
    }
    
    feature_names <- colnames(X_internal)
  }
  
  S <- X_internal
  S$.y_internal <- y_internal
  class_col <- ".y_internal"
  y_raw <- as.character(S[[class_col]])
  
  # ---------------------------
  # 2) Split configuration
  # ---------------------------
  use_pretrained <- !is.null(pretrained_model)
  score_already_available <- isTRUE(input_is_score) || use_pretrained
  
  if (is.null(hnp_split)) {
    hnp_split <- vector("list", Tn)
    
    if (!score_already_available) {
      
      if (isTRUE(grid_search)) {
        hnp_split[[1]] <- c(train = 0.5, threshold = 0.5, error = 0)
        
        if (Tn > 2) {
          for (i in 2:(Tn - 1)) {
            hnp_split[[i]] <- c(train = 0.45, threshold = 0.5, error = 0.05)
          }
        }
        
        hnp_split[[Tn]] <- c(train = 0.95, threshold = 0, error = 0.05)
        
      } else {
        hnp_split[[1]] <- c(train = 0.5, threshold = 0.5, error = 0)
        
        if (Tn > 2) {
          for (i in 2:(Tn - 1)) {
            hnp_split[[i]] <- c(train = 0.5, threshold = 0.5, error = 0)
          }
        }
        
        hnp_split[[Tn]] <- c(train = 1, threshold = 0, error = 0)
      }
      
    } else {
      
      if (isTRUE(grid_search)) {
        hnp_split[[1]] <- c(train = 0, threshold = 1, error = 0)
        
        if (Tn > 2) {
          for (i in 2:(Tn - 1)) {
            hnp_split[[i]] <- c(train = 0, threshold = 0.95, error = 0.05)
          }
        }
        
        hnp_split[[Tn]] <- c(train = 0, threshold = 0, error = 1)
        
      } else {
        hnp_split[[1]] <- c(train = 0, threshold = 1, error = 0)
        
        if (Tn > 2) {
          for (i in 2:(Tn - 1)) {
            hnp_split[[i]] <- c(train = 0, threshold = 1, error = 0)
          }
        }
        
        hnp_split[[Tn]] <- c(train = 0, threshold = 0, error = 1)
      }
    }
    
  } else {
    if (length(hnp_split) != Tn) {
      stop("hnp_split must be a list of length T.")
    }
    
    # When grid_search = FALSE, the error subset is not needed
    if (!isTRUE(grid_search)) {
      hnp_split2 <- vector("list", Tn)
      
      for (i in 1:Tn) {
        cfg_i <- hnp_split[[i]]
        
        r_train <- hnp_read_split_value(cfg_i, "train", if (i == Tn) 1 else 0.5)
        r_thre  <- hnp_read_split_value(cfg_i, "threshold", if (i == Tn) 0 else 0.5)
        
        r_train <- max(0, min(1, r_train))
        r_thre  <- max(0, min(1, r_thre))
        
        if (i < Tn) {
          s <- r_train + r_thre
          
          if (s <= 0) {
            r_train <- 0.5
            r_thre <- 0.5
          } else {
            r_train <- r_train / s
            r_thre  <- r_thre / s
          }
          
          hnp_split2[[i]] <- c(train = r_train, threshold = r_thre, error = 0)
        } else {
          hnp_split2[[i]] <- c(train = r_train, threshold = 0, error = 0)
        }
      }
      
      hnp_split <- hnp_split2
    }
  }
  
  # ---------------------------
  # 3) Split by class: Ss / St / Se
  # ---------------------------
  by_class <- vector("list", Tn)
  
  for (i in 1:Tn) {
    Si <- S[y_raw == as.character(i), , drop = FALSE]
    
    if (nrow(Si) == 0) {
      stop("Class ", i, " has no samples.")
    }
    
    by_class[[i]] <- hnp_split_one_class(Si, hnp_split[[i]], i)
    
    if (!isTRUE(grid_search)) {
      by_class[[i]]$Se <- NULL
    }
    
    if (isTRUE(verbose)) {
      cat(sprintf(
        "[split] class %d: score=%d, threshold=%d, error=%s\n",
        i,
        nrow(by_class[[i]]$Ss),
        nrow(by_class[[i]]$St),
        if (is.null(by_class[[i]]$Se)) "NULL" else as.character(nrow(by_class[[i]]$Se))
      ))
    }
  }
  
  # ---------------------------
  # 4) Construct score_functions
  # ---------------------------
  if (isTRUE(input_is_score)) {
    
    trained_model <- NULL
    
    score_functions <- hnp_build_score_functions_from_scores(
      Tn = Tn,
      feature_names = feature_names,
      eps = .Machine$double.eps
    )
    
  } else {
    
    direct_score_functions <- FALSE
    
    if (use_pretrained) {
      
      # Case A: the user directly passes list(score_functions = ...)
      if (is.list(pretrained_model) && !is.null(pretrained_model$score_functions)) {
        score_functions <- pretrained_model$score_functions
        
        trained_model <- if (!is.null(pretrained_model$model)) {
          pretrained_model$model
        } else {
          pretrained_model
        }
        
        direct_score_functions <- TRUE
        
      } else {
        # ------------------------------------------------------------
        # Key modification:
        # Do not forcibly wrap the fitted object as predict(..., type = "raw").
        # Keep the original pretrained_model object directly.
        # Later, hnp_predict_proba() will try type = "prob",
        # type = "probs", type = "raw", and the SVM probability attribute in sequence.
        # ------------------------------------------------------------
        trained_model <- pretrained_model
      }
      
    } else {
      
      ss_parts <- vector("list", Tn)
      
      for (i in 1:Tn) {
        ss_parts[[i]] <- by_class[[i]]$Ss
      }
      
      S_s <- do.call(rbind, ss_parts)
      S_s <- as.data.frame(S_s)
      
      if (nrow(S_s) == 0) {
        stop("The internal score-training set is empty.")
      }
      
      x_train <- S_s[, feature_names, drop = FALSE]
      y_train <- factor(S_s[[class_col]], levels = required_levels)
      
      trained_model <- base_function(x_train, y_train, method = method)
      
      if (is.null(trained_model)) {
        stop("base_function training failed, method = ", method)
      }
    }
    
    if (!direct_score_functions) {
      score_functions <- hnp_build_score_functions(
        Tn = Tn,
        trained_model = trained_model,
        feature_names = feature_names,
        required_levels = required_levels,
        original_levels = original_levels,
        use_pretrained = use_pretrained,
        method = method,
        eps = .Machine$double.eps
      )
    }
  }
  
  if (length(score_functions) != (Tn - 1)) {
    stop("score_functions must have length T-1.")
  }
  
  # ---------------------------
  # 5) Estimate class prior proportions pi_hat
  #    Used for the objective when grid_search = TRUE;
  #    retained for attribute consistency when grid_search = FALSE, but not used for selection.
  # ---------------------------
  n_total <- nrow(S)
  pi_hat <- numeric(Tn)
  
  for (i in 1:Tn) {
    pi_hat[i] <- sum(y_raw == as.character(i)) / n_total
  }
  
  # ---------------------------
  # 6) grid_search = FALSE
  #
  # According to the LaTeX description:
  # If the user chooses to skip the grid search,
  # the algorithm sets A_i = {bar_t_i}.
  #
  # Therefore, hnp_get_candidates() is no longer called here,
  # and sample(cand_i, 1) is no longer used.
  # Each layer directly takes t_i = bar_t_i.
  # ---------------------------
  if (!isTRUE(grid_search)) {
    
    t_vec <- numeric(0)
    
    for (i in 1:(Tn - 1)) {
      St_i <- by_class[[i]]$St
      
      if (nrow(St_i) == 0) {
        stop("Threshold set for class ", i, " is empty; cannot compute thresholds.")
      }
      
      tbar_i <- hnp_upper_bound(
        S_it = St_i[, feature_names, drop = FALSE],
        level = levels[i],
        delta_i = tolerances[i],
        score_functions = score_functions,
        thresholds = if (i == 1) NULL else t_vec,
        i = i
      )
      
      if (!is.finite(tbar_i) || is.na(tbar_i)) {
        stop("Invalid upper-bound threshold for class ", i, ".")
      }
      
      t_vec <- c(t_vec, tbar_i)
      
      if (isTRUE(verbose)) {
        cat(sprintf("[upper bound only] i=%d, tbar=%.6f\n", i, tbar_i))
      }
    }
    
    out_clf_internal <- hnp_build_output_classifier(
      best_thresholds = t_vec,
      score_functions = score_functions,
      feature_names = feature_names,
      Tn = Tn,
      input_is_score = input_is_score,
      importance_order = importance_order
    )
    
    out_clf <- hnp_wrap_classifier_labels(
      internal_clf = out_clf_internal,
      label_mapping = label_mapping
    )
    
    attr(out_clf, "thresholds") <- t_vec
    attr(out_clf, "objective") <- NA_real_
    attr(out_clf, "T") <- Tn
    attr(out_clf, "method") <- if (isTRUE(input_is_score)) {
      "score_input"
    } else if (use_pretrained) {
      "pretrained"
    } else {
      method
    }
    attr(out_clf, "grid_search") <- FALSE
    attr(out_clf, "selection") <- "upperbound_only"
    attr(out_clf, "importance_order") <- importance_order
    attr(out_clf, "class_levels") <- class_levels
    attr(out_clf, "internal_levels") <- required_levels
    attr(out_clf, "label_mapping") <- label_mapping
    attr(out_clf, "inverse_label_mapping") <- stats::setNames(
      names(label_mapping),
      as.character(label_mapping)
    )
    attr(out_clf, "max_grid") <- max_grid
    attr(out_clf, "max_grid_effective") <- max_grid_effective
    attr(out_clf, "max_combinations") <- max_combinations
    
    return(out_clf)
  }
  
  # ---------------------------
  # 7) grid_search = TRUE
  # ---------------------------
  state_env <- new.env(parent = emptyenv())
  state_env$best_obj <- Inf
  state_env$best_thresholds <- NULL
  state_env$comb_count <- 0L
  
  hnp_search_thresholds_rec(
    i = 1,
    prev_t = numeric(0),
    Tn = Tn,
    by_class = by_class,
    feature_names = feature_names,
    levels = levels,
    tolerances = tolerances,
    score_functions = score_functions,
    grid_set = grid_set,
    max_grid = max_grid_effective,
    max_combinations = max_combinations,
    pi_hat = pi_hat,
    verbose = verbose,
    state_env = state_env
  )
  
  if (is.null(state_env$best_thresholds)) {
    warning("No feasible multi-class HNP Umbrella classifier found ",
            "(empty candidate grid or insufficient samples).")
    return(NULL)
  }
  
  out_clf_internal <- hnp_build_output_classifier(
    best_thresholds = state_env$best_thresholds,
    score_functions = score_functions,
    feature_names = feature_names,
    Tn = Tn,
    input_is_score = input_is_score,
    importance_order = importance_order
  )
  
  out_clf <- hnp_wrap_classifier_labels(
    internal_clf = out_clf_internal,
    label_mapping = label_mapping
  )
  
  attr(out_clf, "thresholds") <- state_env$best_thresholds
  attr(out_clf, "objective") <- state_env$best_obj
  attr(out_clf, "T") <- Tn
  attr(out_clf, "method") <- if (isTRUE(input_is_score)) {
    "score_input"
  } else if (use_pretrained) {
    "pretrained"
  } else {
    method
  }
  attr(out_clf, "grid_search") <- TRUE
  attr(out_clf, "selection") <- "grid_objective_min"
  attr(out_clf, "importance_order") <- importance_order
  attr(out_clf, "class_levels") <- class_levels
  attr(out_clf, "internal_levels") <- required_levels
  attr(out_clf, "label_mapping") <- label_mapping
  attr(out_clf, "inverse_label_mapping") <- stats::setNames(
    names(label_mapping),
    as.character(label_mapping)
  )
  attr(out_clf, "max_grid") <- max_grid
  attr(out_clf, "max_grid_effective") <- max_grid_effective
  attr(out_clf, "max_combinations") <- max_combinations
  attr(out_clf, "comb_count") <- state_env$comb_count
  
  out_clf
}

