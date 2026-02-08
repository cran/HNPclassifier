#' @importFrom dplyr bind_rows
#' @importFrom randomForest randomForest
#' @importFrom e1071 svm
#' @importFrom nnet multinom
#' @importFrom stats predict quantile
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics abline boxplot grid mtext par points
#'
#' @title Classes Mapping function for HNP Algorithm
#' Map class labels to canonical levels "1", "2", "3"
#'
#' @description Validate the class column and re-label provided class names to
#'   canonical factor levels "1", "2", and "3". Useful for preparing datasets
#'   before training and evaluation in the HNP Umbrella pipeline.
#' @param data A data.frame or data.table containing the dataset.
#' @param class_col Character scalar. Name of the class/label column in `data`.
#' @param class_1 Character. Original label that should map to level "1" (most severe with most attentions).
#' @param class_2 Character. Original label that should map to level "2" (median severe).
#' @param class_3 Character. Original label that should map to level "3" (normal or less important).
#' @return The input `data` with `class_col` converted to a factor with
#'   levels c("1","2","3").
#' @examples
#' df <- data.frame(y = c("low","mid","high","mid"), x1 = rnorm(4))
#' df2 <- hnp_map_classes(df, class_col = "y", class_1 = "low", class_2 = "mid", class_3 = "high")
#' table(df2$y)
#' @export
hnp_map_classes <- function(data, class_col, class_1, class_2, class_3) {
  # Testing input parameters
  if (!class_col %in% colnames(data)) {
    stop("error: the specified class column '", class_col, "' does not exist in the data")
  }

  # Checking whether the classes exist
  unique_classes <- unique(data[[class_col]])
  if (!all(c(class_1, class_2, class_3) %in% unique_classes)) {
    warning("warning: some specified classes do not exist in the data")
    message("classes in the data:")
    message(paste(unique_classes, collapse = ", "))
    message("specified classes:")
    message(paste(c(class_1, class_2, class_3), collapse = ", "))
  }
  
  # Create class mapping
  class_mapping <- c(class_1, class_2, class_3)
  names(class_mapping) <- c("1", "2", "3")
  
  # re-label
  data[[class_col]] <- factor(data[[class_col]], 
                              levels = class_mapping,
                              labels = c("1", "2", "3"))
  
  return(data)
}


############################
# Algorithm 1
############################
#' @title Delta search 
#' 
#' @description Calculate the order k of the statistic that 
#'   satisfies the given confidence requirements for 
#'   determining classification thresholds.
#' 
#' @param n Integer specifying the cardinality of the 
#'   grid set Tau (size of `S_it`).
#' @param level Numeric between 0 and 1 representing the 
#' desired control level (alpha) for the ith under-classification error.
#' @param delta Numeric tolerance parameter for the confidence bound.
#'
#' @return An integer `k` representing the order of the 
#'   statistic that meets the 
#'   confidence requirements. Returns `NA` if no valid solution exists.
#'
#' @examples
#' k <- hnp_delta_search(n = 100, level = 0.05, delta = 0.01)
#'
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_delta_search <- function(n, level, delta) {
  k <- 0
  v_k <- 0
  
  # Warning of minimum sample size
  if (n < (log(delta) / log(1 - level))) {
    return(k)
  }
  else{
        while (v_k <= delta) {
        v_k <- v_k + choose(n, k) * level^k * (1 - level)^(n - k)
        k <- k + 1
    }
  }
  return(k)
}


############################
# Algorithm 2
############################
#' @title Upper Bound of the ith Threshold (Optimal ith Threshold)
#' 
#' @description Compute the optimal threshold for class i using score 
#'   functions and confidence bounds, given tolerance and under classification 
#'   error level.
#'
#' @param S_it The left-out class-i samples.
#' @param level (alpha) desired control level for the ith under 
#'   classification error.
#' @param delta_i ith tolerance parameter.
#' @param score_functions A list of score functions (T_1, ..., T_i).
#' @param thresholds Numeric vector of length `i - 1` with thresholds for
#'   previously evaluated classes; ignored when `i == 1`.
#' @param i Class-i.
#' 
#' @examples
#' set.seed(123)
#' n <- 200
#' S_it <- data.frame(
#'   feature1 = rnorm(n, mean = 2, sd = 1),
#'   feature2 = runif(n, min = 0, max = 5)
#' )
#' level <- 0.05
#' delta_i <- 0.01
#' score_functions <- list(
#'   function(data) runif(nrow(data)),
#'   function(data) runif(nrow(data))
#' )
#' thresholds <- c(2.5, NA)
#' i <- 1
#' t_i_bar <- hnp_upper_bound(S_it, level, delta_i, score_functions, thresholds, i)
#'
#' @return t_i_bar Optimal ith threshold.
#' 
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_upper_bound <- function(S_it, level, delta_i, score_functions, thresholds, i){
  n_i <- nrow(S_it)
  T_i <- score_functions[[i]]

  #Tau
  Tau <- vector("list", n_i)
  for (ii in 1:n_i) {
    Tau[[ii]] <- T_i(S_it[ii, , drop = FALSE])
  }
  if (!is.atomic(Tau)) Tau <- unlist(Tau) 
  if (any(is.na(Tau))) warning("NA values detected in Tau")
  Tau <- sort(Tau, na.last = TRUE)

  k_i = hnp_delta_search(n_i, level, delta_i)

  if(k_i == 0){
    stop("Exceed minimum sample size required")
  }
  t_i_bar <- Tau[k_i]

  if(i>1){
    filter_condition <- function(x_row_df) {
      x_row_df <- as.data.frame(x_row_df)
      all(sapply(1:(i-1), function(j) {
        score_functions[[j]](x_row_df) < thresholds[j]
      }))
    }
    
    filtered_indices <- vapply(seq_len(nrow(S_it)), function(r) {
      filter_condition(S_it[r, , drop = FALSE])
    }, logical(1))
    
    S_it_prime <- S_it[filtered_indices, , drop = FALSE]
    n_i_prime <- nrow(S_it_prime)
    if(n_i_prime == 0){
      return(t_i_bar)
    }
    
    Tau_i_prime <- sort(T_i(S_it_prime), na.last = TRUE)
    
    if(n_i_prime>0){
      #correct the calculations of cn_i and delta_i_prime
      p_i_hat <- n_i_prime/n_i
      cn_i <- 2 / sqrt(n_i)
      p_i <- p_i_hat + cn_i
      a_i_prime <- level/p_i
      delta_i_prime <- delta_i - exp(-2 * n_i * cn_i^2)
      
      if( a_i_prime<1  && n_i_prime>=(log(delta_i_prime)/log(1-a_i_prime))){
        k_i_prime <- hnp_delta_search(n_i_prime, a_i_prime, delta_i_prime)
        if ((k_i_prime != 0) && k_i_prime <= length(Tau_i_prime)) {
          #Using adjusted threshold based on filtered subset
          t_i_bar <- Tau_i_prime[k_i_prime]
        }
      }
    }
  }
 
  return(t_i_bar)
}


############################
# Algorithm 3
############################

#' @title Base Classifier Training function
#' Train a base multi-class model (RF / SVM / Multinomial Logistic)
#'
#' @description Fit one of the supported classifiers for ternary classification:
#'   Random Forest, SVM (with probabilities), or multinomial logistic regression
#'   via `nnet::multinom`.
#' @param x A data.frame of predictors/features.
#' @param y A factor response with levels "1","2","3".
#' @param method Character string: one of 'randomforest', 'svm', or 'logistic'.
#' @return A trained model object compatible with the downstream scoring
#'   functions.
#' @examples
#' set.seed(123)
#' x <- data.frame(a = rnorm(20), b = rnorm(20))
#' y <- factor(sample(c("1","2","3"), 20, TRUE))
#' model <- base_function(x, y, method = 'randomforest')
#' @export
base_function <- function(x, y, method = 'randomforest') {
  # Define available methods for ternary classification
  method_choices = c("randomforest", "svm", "logistic")
  method = match.arg(method, method_choices)
  
  x = as.data.frame(x)
  y = as.factor(y)
  
  # Train base algorithm based on method
  if(method == 'randomforest') {
    # Random Forest - good for ternary classification
    model = randomForest::randomForest(x = x, y = y, ntree = 100, mtry = 2)
  } 
  else if(method == 'svm') {
    # Support Vector Machine - supports multiple classes
    dataset = data.frame(x, y = y)
    colnames(dataset)[ncol(dataset)] = "y"  
    model = e1071::svm(y ~ ., data = dataset, probability = TRUE, scale = TRUE)
  }
  else if(method == 'logistic') {
    # Multinomial logistic regression via nnet::multinom
    dataset = data.frame(x, y = y)
    colnames(dataset)[ncol(dataset)] = "y"
    model = nnet::multinom(y ~ ., data = dataset, trace = FALSE)
  }
  else{
    return(NULL)
  }
  return(model)
}


#' @title T1 Calculation
#' Create T1 scoring function from a fitted model
#'
#' @description Return a function that takes new data and outputs the score
#'   for class 1, typically the predicted probability P(Y=1|X). Works with the
#'   supported `method`s used by `base_function`.
#' @param model A fitted model returned by `base_function` or equivalent.
#' @param method Character string specifying the model family used: one of
#'   'svm', 'randomforest', or 'logistic'.
#' @return A function of the form `function(X) numeric`, where `X` is a
#'   data.frame of features and the returned numeric vector are scores for class 1.
#' @examples
#' set.seed(123)
#' x <- data.frame(a = rnorm(20), b = rnorm(20))
#' y <- factor(sample(c("1","2","3"), 20, TRUE))
#' model <- base_function(x, y, method = 'randomforest')
#' T1 <- probability_to_score_1(model, method = 'randomforest')
#' newx <- data.frame(a = rnorm(5), b = rnorm(5))
#' scores <- T1(newx)
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
probability_to_score_1 <- function(model, method) {
  function(X) {
    X <- as.data.frame(X)
    if (method == "svm") {
      # SVM
      svm_predictions <- predict(model, newdata = X, probability = TRUE, decision.values = TRUE)
      prob <- attr(svm_predictions, "probabilities")
      if (is.matrix(prob) || is.data.frame(prob)) {
        if ("1" %in% colnames(prob)) {
          prob[, "1"]
        } else if (ncol(prob) >= 1) {
          prob[, 1]  # first column
        } else {
          stop("Cannot extract probability for class 1 from SVM")
        }
      } else {
        stop("SVM probability prediction failed")
      }
    }

    else if (method == "randomforest") {
      # Random Forest
      prob <- predict(model, newdata = X, type = "prob")
      if (is.matrix(prob) || is.data.frame(prob)) {
        if ("1" %in% colnames(prob)) {
          prob[, "1"]
        } else if (ncol(prob) >= 1) {
          prob[, 1]
        } else {
          stop("Cannot extract probability for class 1 from Random Forest")
        }
      } else {
        stop("Random Forest probability prediction failed")
      }
    }

    else if(method == 'logistic') {
      # Logistic: use multinomial probabilities
      probs <- predict(model, newdata = X, type = "probs")
      if (is.null(dim(probs))) {
        prob <- as.data.frame(t(as.matrix(probs)))
      } else {
        prob <- as.data.frame(probs)
      }


      wanted <- c("1","2","3")
      if (is.null(colnames(prob))) {
        #If there are no column names, assign names directly 
        #(assuming the labels have been mapped to factor levels "1", "2", "3" during training)
        colnames(prob) <- wanted[seq_len(ncol(prob))]
      }
      missing <- setdiff(wanted, colnames(prob))
      if (length(missing) > 0) {
        # zeros padding
        for (m in missing) prob[[m]] <- 0
      }
      prob <- prob[, wanted, drop = FALSE]

      if (is.matrix(prob) || is.data.frame(prob)) {
        if (!is.null(colnames(prob)) && ("1" %in% colnames(prob))) {
          prob[, "1"]
        } else if (ncol(prob) >= 1) {
          prob[, 1]
        } else {
          stop("Cannot extract probability for class 1 from logistic")
        }
      } else {
        stop("Logistic probability prediction failed")
      }
    }

    else {
      stop("Unknown method: ", method)
    }
  }
}

#' @title T2 Calculation
#' Create T2 scoring function as ratio P(class 2)/P(class 3)
#'
#' @description Return a function that produces the ratio of predicted
#'   probabilities P(Y=2|X) / P(Y=3|X), with safeguards for zeros/NA and
#'   infinite values. Works with the supported `method`s used by `base_function`.
#' @param model A fitted model returned by `base_function` or equivalent.
#' @param method Character string specifying the model family used: one of
#'   'svm', 'randomforest', or 'logistic'.
#' @return A function of the form `function(X) numeric`, where `X` is a
#'   data.frame of features and the returned numeric vector are T2 scores.
#' @examples
#' set.seed(123)
#' x <- data.frame(a = rnorm(20), b = rnorm(20))
#' y <- factor(sample(c("1","2","3"), 20, TRUE))
#' model <- base_function(x, y, method = 'randomforest')
#' T2 <- probability_to_score_2(model, method = 'randomforest')
#' newx <- data.frame(a = rnorm(5), b = rnorm(5))
#' scores <- T2(newx)
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
probability_to_score_2 <- function(model, method) {
  function(X) {
    X <- as.data.frame(X)
    if (method == "svm") {
      # SVM
      pred <- predict(model, newdata = X, probability = TRUE)
      prob <- attr(pred, "probabilities")
      if (is.matrix(prob) || is.data.frame(prob)) {
        if ("2" %in% colnames(prob) && "3" %in% colnames(prob)) {
          # protect from na
          p2 <- prob[, "2"]
          p3 <- prob[, "3"]
        } else if (ncol(prob) >= 3) {
          p2 <- prob[, 2]
          p3 <- prob[, 3]
        } else {
          stop("Cannot extract probability for classes 2 and 3 from SVM")
        }
        # Safely handle denominator: replace 0 or very small values with epsilon
        p3[p3 < .Machine$double.eps] <- .Machine$double.eps
        ratio <- p2 / p3
        ratio[!is.finite(ratio)] <- 0
        ratio
      } else {
        stop("SVM probability prediction failed")
      }
    } 
    else if (method == "randomforest") {
      # Random Forest
      # protect from na
      prob <- predict(model, newdata = X, type = "prob")
      if (is.matrix(prob) || is.data.frame(prob)) {
        if ("2" %in% colnames(prob) && "3" %in% colnames(prob)) {
          p2 <- prob[, "2"]
          p3 <- prob[, "3"]
        } else if (ncol(prob) >= 3) {
          p2 <- prob[, 2]
          p3 <- prob[, 3]
        } else {
          stop("Cannot extract probability for classes 2 and 3 from Random Forest")
        }
        # Safely handle denominator: replace 0 or very small values with epsilon
        p3[p3 < .Machine$double.eps] <- .Machine$double.eps
        ratio <- p2 / p3
        ratio[!is.finite(ratio)] <- 0
        ratio
      } else {
        stop("Random Forest probability prediction failed")
      }
    }
    else if(method == 'logistic') {
      probs <- predict(model, newdata = X, type = "probs")
      if (is.null(dim(probs))) {
        prob <- as.data.frame(t(as.matrix(probs)))
      } else {
        prob <- as.data.frame(probs)
      }

      # rename the factor to "1","2","3"
      wanted <- c("1","2","3")
      if (is.null(colnames(prob))) {
        #If there are no column names, assign names directly 
        #(assuming the labels have been mapped to factor levels "1", "2", "3" during training)
        colnames(prob) <- wanted[seq_len(ncol(prob))]
      }
      missing <- setdiff(wanted, colnames(prob))
      if (length(missing) > 0) {
        # set the missing columns as 0s
        for (m in missing) prob[[m]] <- 0
      }
      prob <- prob[, wanted, drop = FALSE]

      if (is.matrix(prob) || is.data.frame(prob)) {
        cols <- colnames(prob)
        if (all(c("2", "3") %in% cols)) {
          p2 <- prob[, "2"]
          p3 <- prob[, "3"]
        } else if (ncol(prob) >= 3) {
          p2 <- prob[, 2]
          p3 <- prob[, 3]
        } else {
          stop("Cannot extract probability for classes 2 and 3 from logistic")
        }
        # Safely handle denominator: replace 0 or very small values with epsilon
        p3[p3 < .Machine$double.eps] <- .Machine$double.eps
        ratio <- p2 / p3
        ratio[!is.finite(ratio)] <- 0
        ratio
      } else {
        stop("Logistic probability prediction failed")
      }
    }
    else {
      stop("Unknown method: ", method)
    }
  }
}

#' @title HNP Umbrella Algorithm
#' @description Implementation of the HNP Umbrella algorithm for ternary classification
#' @param S Training dataset
#' @param levels Confidence levels (alpha) for each class
#' @param tolerances Tolerance parameters (delta) for each class
#' @param A1 Candidate thresholds for class 1
#' @param method Classification method to use ('randomforest', 'svm', 'logistic')
#' @param hnp_split Data splitting ratios for each class
#' @param class_col Character scalar. Name of the class column in the dataset
#'   (must be mapped to levels "1","2","3").
#' @return A classifier function that takes new data and classifies it into a class with controlled 
#' type-one error rate
#' @examples
#' set.seed(123)
#' n <- 500
#' features <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' y <- factor(sample(c("1", "2", "3"), n, replace = TRUE, prob = c(0.2, 0.3, 0.5)))
#' data <- cbind(features, y)
#' clf <- hnp_umbrella(
#'   S = data,
#'   levels = c(0.1, 0.1),
#'   tolerances = c(0.1, 0.1),
#'   class_col = "y",
#'   method = "randomforest"
#' )
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_umbrella <- function(S, levels, tolerances, A1 = NULL,
                         method = 'randomforest',
                         hnp_split = NULL,
                         class_col) {
  # Convert data.table to data.frame to avoid ..feature_names issues
  S <- as.data.frame(S)

  # Split ratio by default
  if (is.null(hnp_split)) {
    hnp_split <- list(
      c(train = 0.5, threshold = 0.5),                # Class 1
      c(train = 0.45, threshold = 0.5, error = 0.05),  # Class 2
      c(train = 0.95, error = 0.05)               # Class 3
    )
  }

  ########## if (error ratio==0): use t1_bar for threshold, skip search
  
  # Class 1
  S_1 <- S[S[[class_col]] == "1", ]
  n_S_1 <- nrow(S_1)
  index_1 <- sample(1:n_S_1, n_S_1)
  train_end_1 <- floor(hnp_split[[1]]["train"] * n_S_1)
  S_1s <- S_1[index_1[1:train_end_1], ]
  S_1t <- S_1[index_1[(train_end_1+1):n_S_1], ]

  # Class 2
  S_2 <- S[S[[class_col]] == "2", ]
  n_S_2 <- nrow(S_2)
  index_2 <- sample(1:n_S_2, n_S_2)
  train_end_2 <- floor(hnp_split[[2]]["train"] * n_S_2)
  threshold_end_2 <- floor((hnp_split[[2]]["train"] + hnp_split[[2]]["threshold"]) * n_S_2)
  S_2s <- S_2[index_2[1:train_end_2], ]
  S_2t <- S_2[index_2[(train_end_2+1):threshold_end_2], ]
  S_2e <- S_2[index_2[(threshold_end_2+1):n_S_2], ]

  # Class 3
  S_3 <- S[S[[class_col]] == "3", ]
  n_S_3 <- nrow(S_3)
  index_3 <- sample(1:n_S_3, n_S_3)
  train_end_3 <- floor(hnp_split[[3]]["train"] * n_S_3)
  S_3s <- S_3[index_3[1:train_end_3], ]
  S_3e <- S_3[index_3[(train_end_3+1):n_S_3], ]


  # Class probabilities
  pie_2_hat <- nrow(S_2) / nrow(S)
  pie_3_hat <- nrow(S_3) / nrow(S)
  
  S_s <- dplyr::bind_rows(S_1s, S_2s, S_3s)

  # Get scoring functions using base_function
  feature_names <- setdiff(names(S_s), class_col)
  x_train <- S_s[, feature_names, drop = FALSE]
  y_train <- S_s[[class_col]]

  # Train model using base_function
  trained_model <- base_function(x_train, y_train, method = method)
  if(is.null(trained_model)){
    stop("Unknown method: ", method)
    return(NULL)
  }

  # Create scoring functions
  T1 <- probability_to_score_1(trained_model, method = method)
  T2 <- probability_to_score_2(trained_model, method = method)

  score_functions <- list(T1, T2)

  t1_bar <- hnp_upper_bound(S_1t[, feature_names, drop = FALSE], levels[1], tolerances[1], score_functions, NULL, 1)
  Total_error <- 1
  best_classifier <- NULL


  # If error validation sets are empty (error ratio is 0), skip search and use t1_bar directly (for small dataset)
  if (nrow(S_2e) == 0 && nrow(S_3e) == 0) {
     t_sample <- c(t1_bar)
  } else {
    if (is.null(A1)) {
      valid_A1 <- T1(S_1t[, feature_names, drop = FALSE])
    }
    else{
      valid_A1 <- A1
    }
  
    valid_A1 <- valid_A1[!is.na(valid_A1) & valid_A1 <= t1_bar]
  
    if (length(valid_A1) > 0 && all(valid_A1 <= t1_bar, na.rm = TRUE)) {
      t_sample <- sample(valid_A1, min(15, length(valid_A1)), replace = FALSE)
    } else {
      warning("Not valid candidate A1 values, use S1t instead")
      valid_A1 <- T1(S_1t[, feature_names, drop = FALSE])
      valid_A1 <- valid_A1[!is.na(valid_A1) & valid_A1 <= t1_bar]
      t_sample <- sample(valid_A1, min(15, length(valid_A1)), replace = FALSE)
  
    }
  }

    for(t1 in t_sample) {
    
    t2 <- hnp_upper_bound(S_2t[, feature_names, drop = FALSE], levels[2], tolerances[2], score_functions, c(t1), 2)

    classifier <- function(X) {
      X <- as.data.frame(X)
      X <- X[, feature_names, drop = FALSE]
      if (nrow(X) == 0) {
        return(integer(0))
      }
      s1 <- as.numeric(T1(X))
      s2 <- as.numeric(T2(X))
      s1[is.na(s1)] <- -Inf
      s2[is.na(s2)] <- -Inf
      set1 <- which(s1 >= t1)
      set2 <- which(s2 >= t2)
      decision <- rep.int(3L, nrow(X))
      decision[set2] <- 2L
      decision[set1] <- 1L
      decision
    }

    # Empirical errors
    e21 <- if (nrow(S_2e) > 0) {
      preds_2e <- classifier(S_2e[, feature_names, drop = FALSE])
      mean(preds_2e == 1L)
    } 
    else {
      0
    }
    e3  <- if (nrow(S_3e) > 0) {
      preds_3e <- classifier(S_3e[, feature_names, drop = FALSE])
      mean(preds_3e %in% c(1L, 2L))
    } 
    else {
      0
    }
      
    # Weighted error
    current_error <- pie_2_hat * e21 + pie_3_hat * e3
    
    # Update best classifier
    if(current_error < Total_error) {
      Total_error <- current_error
      best_classifier <- classifier
    }
  }

  best_classifier

  if (!is.null(best_classifier)) {
    hnp_umbrella_classifier <- function(new_data) {
      new_data <- as.data.frame(new_data) 
      new_data <- new_data[, feature_names, drop = FALSE]

      if (nrow(new_data) == 0) {
        return(integer(0))
      }

      return(as.integer(best_classifier(new_data)))
    }

    return(hnp_umbrella_classifier)
  } else {
    warning("No valid classifier found")
    return(NULL)
  }
}


#' HNP Umbrella (flex): use custom score functions and pre-split data
#'
#' @description Flexible variant of the HNP Umbrella algorithm that accepts
#'   user-provided scoring functions and explicit data splits for thresholding
#'   and error estimation. This bypasses model training inside and focuses on
#'   threshold selection with confidence controls.
#' @param score_data A data.frame for fitting/deriving scoring behavior.
#' @param threshold_data A data.frame used to compute thresholds.
#' @param error_data A data.frame used to estimate empirical errors.
#' @param levels Numeric vector of length 2. Confidence levels (alpha) for
#'   class 1 and class 2 under-classification controls.
#' @param tolerances Numeric vector of length 2. Tolerance (delta) values for
#'   the corresponding classes.
#' @param A1 Optional numeric vector of candidate thresholds for class 1.
#' @param score_functions A list with at least two functions: `T1`, `T2`.
#'   Each must accept a data.frame and return numeric scores.
#' @param class_col Character scalar. Name of the class column in the data.
#' @return A classifier function `function(new_data) data.frame(result=...)`,
#'   or `NULL` if no valid classifier is found.
#' @examples
#' set.seed(123)
#' n <- 500
#' score_data <- data.frame(x=rnorm(n), y=factor(sample(1:3, n, replace=TRUE)))
#' threshold_data <- data.frame(x=rnorm(n), y=factor(sample(1:3, n, replace=TRUE)))
#' error_data <- data.frame(x=rnorm(n), y=factor(sample(1:3, n, replace=TRUE)))
#' T1 <- function(d) as.numeric(d$x > 0)
#' T2 <- function(d) as.numeric(d$x > 0.5)
#' clf <- hnp_umbrella_flex(score_data, threshold_data, error_data,
#'                          levels = c(0.05, 0.05), tolerances = c(0.01, 0.01),
#'                          score_functions = list(T1, T2), class_col = 'y')
#' preds <- clf(score_data)
#' @references Lijia Wang, Y. X. Rachel Wang, Jingyi Jessica Li, and Xin Tong (2024).
#' "Hierarchical Neyman-Pearson Classification for Prioritizing Severe Disease
#' Categories in COVID-19 Patient Data."
#' \emph{Journal of the American Statistical Association}, 119(545), 39-51.
#' \doi{10.1080/01621459.2023.2270657}
#' @export
hnp_umbrella_flex <- function(score_data, threshold_data, error_data,
                              levels, tolerances, A1 = NULL,
                              score_functions = NULL,
                              class_col) {
  if (is.null(score_functions) || length(score_functions) < 2) {
    stop("score_functions must be a list of at least two functions")
  }


  score_data <- as.data.frame(score_data)
  threshold_data <- as.data.frame(threshold_data)
  error_data <- as.data.frame(error_data)

  datasets <- list(score = score_data, threshold = threshold_data, error = error_data)
  for (name in names(datasets)) {
    data_i <- datasets[[name]]
    if (nrow(data_i) == 0) next
    if (!(class_col %in% colnames(data_i))) {
      stop("class column '", class_col, "' not found in ", name, "_data")
    }
  }

  feature_names <- setdiff(colnames(score_data), class_col)
  if (length(feature_names) == 0) {
    stop("score_data must contain at least one feature column apart from the class column")
  }

  T1 <- score_functions[[1]]
  T2 <- score_functions[[2]]
  if (!is.function(T1) || !is.function(T2)) {
    stop("score_functions[[1]] and score_functions[[2]] must be functions")
  }

  S <- dplyr::bind_rows(score_data, threshold_data, error_data)
  if (nrow(S) == 0) {
    stop("Combined dataset is empty; check the provided data splits")
  }

  S_1t <- threshold_data[threshold_data[[class_col]] == "1", , drop = FALSE]
  S_2t <- threshold_data[threshold_data[[class_col]] == "2", , drop = FALSE]
  S_3t <- threshold_data[threshold_data[[class_col]] == "3", , drop = FALSE]
  S_1e <- error_data[error_data[[class_col]] == "1", , drop = FALSE]
  S_2e <- error_data[error_data[[class_col]] == "2", , drop = FALSE]
  S_3e <- error_data[error_data[[class_col]] == "3", , drop = FALSE]

  S_1 <- S[S[[class_col]] == "1", , drop = FALSE]
  S_2 <- S[S[[class_col]] == "2", , drop = FALSE]
  S_3 <- S[S[[class_col]] == "3", , drop = FALSE]

  pie_2_hat <- if (nrow(S) > 0) nrow(S_2) / nrow(S) else 0
  pie_3_hat <- if (nrow(S) > 0) nrow(S_3) / nrow(S) else 0

  if (nrow(S_1t) == 0) {
    stop("threshold_data must contain class 1 samples")
  }
  if (nrow(S_2t) == 0) {
    stop("threshold_data must contain class 2 samples")
  }

  t1_bar <- hnp_upper_bound(S_1t[, feature_names, drop = FALSE], levels[1], tolerances[1], score_functions, NULL, 1)
  Total_error <- Inf
  best_classifier <- NULL

  if (is.null(A1)) {
    valid_A1 <- as.numeric(T1(S_1t[, feature_names, drop = FALSE]))
  } else {
    valid_A1 <- A1
  }

  valid_A1 <- valid_A1[!is.na(valid_A1) & valid_A1 <= t1_bar]

  if (length(valid_A1) > 0 && all(valid_A1 < t1_bar, na.rm = TRUE)) {
    t_sample <- c(t1_bar)
  } else if (length(valid_A1) > 0) {
    t_sample <- sample(valid_A1, min(5, length(valid_A1)), replace = FALSE)
    t_sample <- c(t1_bar, t_sample)
  } else {
    warning("No valid candidate A1 values; use t1_bar instead")
    t_sample <- c(t1_bar)
  }

  for (t1 in t_sample) {
    t2 <- hnp_upper_bound(S_2t[, feature_names, drop = FALSE], levels[2], tolerances[2], score_functions, c(t1), 2)

    classifier <- function(X) {
      X <- as.data.frame(X)
      X <- X[, feature_names, drop = FALSE]
      if (nrow(X) == 0) {
        return(integer(0))
      }
      s1 <- as.numeric(T1(X))
      s2 <- as.numeric(T2(X))
      s1[is.na(s1)] <- -Inf
      s2[is.na(s2)] <- -Inf
      set1 <- which(s1 >= t1)
      set2 <- which(s2 >= t2)
      decision <- rep.int(3L, nrow(X))
      decision[set2] <- 2L
      decision[set1] <- 1L
      decision
    }

    e21 <- if (nrow(S_2e) > 0) {
      preds_2e <- classifier(S_2e[, feature_names, drop = FALSE])
      mean(preds_2e == 1L)
    } else {
      0
    }

    e3 <- if (nrow(S_3e) > 0) {
      preds_3e <- classifier(S_3e[, feature_names, drop = FALSE])
      mean(preds_3e %in% c(1L, 2L))
    } else {
      0
    }

    current_error <- pie_2_hat * e21 + pie_3_hat * e3

    if (current_error < Total_error) {
      Total_error <- current_error
      best_classifier <- classifier
    }
  }

  if (!is.null(best_classifier)) {
    hnp_umbrella_classifier <- function(new_data) {
      new_data <- as.data.frame(new_data)
      new_data <- new_data[, feature_names, drop = FALSE]

      if (nrow(new_data) == 0) {
        return(integer(0))
      }

      return(as.integer(best_classifier(new_data)))
    }

    return(hnp_umbrella_classifier)
  } else {
    warning("No valid classifier found")
    return(NULL)
  }
}



############################
# hnp_summary
############################

#' @title hnp_summary
#' Summarize a ternary classifier's performance
#'
#' @description Compute confusion matrix, class-wise false positive/negative
#'   rates, over- and under-classification errors, overall accuracy, and a
#'   normalized error table for a ternary classifier produced by the HNP
#'   pipeline.
#' @param classifier A function `function(X) { ... }` that returns class labels
#'   1/2/3 for a single-row data.frame or vectorized over rows.
#' @param data A data.frame containing features and the true class column.
#' @param class_col Character scalar. Name of the true class/label column.
#' @param class_number Optional integer. Number of classes; if `NULL`, inferred
#'   from the data.
#' @return A list with components: `confusion_matrix`, `false_positive_rate`,
#'   `false_negative_rate`, `overall_accuracy`, `predictions`,
#'   `under_classification_error`, `over_classification_error`,
#'   `total_over_classification_error`, `total_under_classification_error`, and
#'   `error_table`.
#' @examples
#' set.seed(123)
#' n <- 50
#' x <- data.frame(a = rnorm(n), b = rnorm(n))
#' y <- factor(sample(c("1","2","3"), n, TRUE))
#' df <- cbind(x, y)
#' clf <- function(X) sample(c(1,2,3), nrow(X), replace=TRUE)
#' res <- hnp_summary(clf, data = df, class_col = "y")
#' @export
hnp_summary <- function(classifier, data, class_col, class_number = NULL) {
  data <- as.data.frame(data)
  feature_names <- setdiff(names(data), class_col)
  
 
  X_all <- data[, feature_names, drop = FALSE]
  pred_all <- classifier(X_all)
  
  if (length(pred_all) != nrow(data)) {
    stop("The classifier must return a vector of length nrow(data).")
  }
  
  if (is.null(class_number)) {
 
    class_number <- max(as.numeric(as.character(data[[class_col]])), na.rm = TRUE)
  }
  
  predictions <- data.frame(
    true_class      = data[[class_col]],
    predicted_class = pred_all
  )
  
  predictions$true_class <- factor(
    predictions$true_class,
    levels = as.character(1:class_number)
  )
  predictions$predicted_class <- factor(
    as.character(predictions$predicted_class),
    levels = as.character(1:class_number)
  )
  
  # confusion matrix
  conf_matrix <- table(
    True      = predictions$true_class,
    Predicted = predictions$predicted_class
  )
  
  # False Positive Rate, False Negative Rate
  false_positive_rate <- numeric(class_number)
  false_negative_rate <- numeric(class_number)
  under_classification_error <- numeric(class_number)
  under_classification_error[class_number] <- 0
  
  over_classification_error  <- numeric(class_number)
  over_classification_error[1] <- 0
  
  for (class in 1:class_number) {
    true_negative  <- sum(conf_matrix[-class, -class])
    false_positive <- sum(conf_matrix[-class, class])
    false_positive_rate[class] <- false_positive / (false_positive + true_negative)
    
    true_positive  <- conf_matrix[class, class]
    false_negative <- sum(conf_matrix[class, -class])
    false_negative_rate[class] <- false_negative / (true_positive + false_negative)
    
    if (class != class_number) {
      classify_to_low <- conf_matrix[class, which(1:class_number > class)]
      under_classification_error[class] <-
        sum(classify_to_low) / sum(conf_matrix[class, ])
    }
    
    if (class != 1) {
      classify_to_high <- conf_matrix[class, which(1:class_number < class)]
      over_classification_error[class] <-
        sum(classify_to_high) / sum(conf_matrix[class, ])
    }
  }
  
  weight_proportion  <- rowSums(conf_matrix) / sum(conf_matrix)
  total_over_classification_error  <- sum(over_classification_error  * weight_proportion)
  total_under_classification_error <- sum(under_classification_error * weight_proportion)
  
  correct_predictions <- sum(diag(conf_matrix))
  total_samples       <- sum(conf_matrix)
  overall_accuracy    <- correct_predictions / total_samples
  
  error_table <- conf_matrix / rowSums(conf_matrix)
  error_table <- cbind(error_table, false_negative_rate)
  colnames(error_table) <- c(
    paste0("Pred_", 1:class_number),
    "total_classification_error"
  )
  rownames(error_table) <- paste0("True_", 1:class_number)
  
  list(
    confusion_matrix                 = conf_matrix,
    false_positive_rate              = false_positive_rate,
    false_negative_rate              = false_negative_rate,
    overall_accuracy                 = overall_accuracy,
    predictions                      = predictions,
    under_classification_error       = under_classification_error,
    over_classification_error        = over_classification_error,
    total_over_classification_error  = total_over_classification_error,
    total_under_classification_error = total_under_classification_error,
    error_table                      = error_table
  )
}

#' @title HNP Box Plot Experiment
#'
#' @description Runs multiple iterations of HNP experiment on a dataset (with random 7:3 splits)
#' and generates a PDF with 15 boxplots comparing Before vs After NP performance.
#'
#' @param data A data.frame containing features and class label.
#' @param class_col Character. Name of the class column (must be mapped to "1","2","3").
#' @param method Character. Base classifier method ('randomforest', 'svm', 'logistic').
#' @param n_runs Integer. Number of iterations to run.
#' @param levels Numeric vector. Alpha levels (constraints) for classes (e.g., c(0.05, 0.1)).
#' @param tolerances Numeric vector. Delta tolerances for classes (e.g., c(0.01, 0.02)).
#' @param output_file Character. Path to save the PDF output.
#' @param hnp_split List. Split configuration for HNP internal validation.
#' @param split_ratio Numeric vector. Ratio of data used for training and testing (e.g., c(0.7, 0.3)).
#' @return No return value, called for side effects.
#' @examples
#' set.seed(123)
#' n <- 2000
#' features <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' y <- factor(sample(c("1", "2", "3"), n, replace = TRUE, prob = c(0.2, 0.3, 0.5)))
#' data <- cbind(features, y)
#' hnp_box_plot(
#'   data = data,
#'   class_col = "y",
#'   method = "logistic", 
#'   n_runs = 2,
#'   levels = c(0.05, 0.05),
#'   tolerances = c(0.05, 0.05),
#'   output_file = tempfile(fileext = ".pdf")
#' )
#' @export
hnp_box_plot <- function(data, class_col, method = 'logistic', 
                         n_runs = 100,
                         levels = c(0.05, 0.05),
                         tolerances = c(0.05, 0.05),
                         output_file = NULL,
                         hnp_split = NULL,
                         split_ratio = c(0.7, 0.3)) {
  
  # Use the first element as train ratio
  train_ratio <- split_ratio[1]

  # If no output file path is provided, put it into tempfile()
  if (is.null(output_file)) {
    output_file <- tempfile(fileext = ".pdf")
  }

  # Ensure data is data.frame (fix for data.table inputs)
  data <- as.data.frame(data)
  
  if (!class_col %in% names(data)) stop(sprintf("Label column '%s' not found in data", class_col))
  
  feature_names <- setdiff(names(data), class_col)
  
  # Initialize results storage
  results_list <- vector("list", n_runs)
  
  # Experiment Loop
  for (r in 1:n_runs) {
    
    # Random Split using split_ratio
    n_total <- nrow(data)
    n_train <- floor(train_ratio * n_total)
    train_idx <- sample(seq_len(n_total), n_train)
    
    Train <- data[train_idx, ]
    Test  <- data[-train_idx, ]
    
    # Train Base Model
    fit_before <- base_function(
      x = Train[, feature_names, drop = FALSE],
      y = Train[[class_col]],
      method = method
    )
    
    # baseline prediction (batch-wise)
    baseline_clf <- function(Xnew) {
      Xnew <- as.data.frame(Xnew)
      # Ensure correct columns
      missing_cols <- setdiff(feature_names, names(Xnew))
      if(length(missing_cols) > 0) stop("Missing columns in prediction data")
      Xnew <- Xnew[, feature_names, drop = FALSE]
      
      prob <- NULL
      
      if (method == "svm") {
        pred <- predict(fit_before, newdata = Xnew, probability = TRUE)
        prob <- attr(pred, "probabilities")
      } else if (method == "randomforest") {
        prob <- predict(fit_before, newdata = Xnew, type = "prob")
      } else if (method == "logistic") {
        prob <- predict(fit_before, newdata = Xnew, type = "probs")
        if (!is.matrix(prob)) {
           prob <- t(as.matrix(prob))
        }
      }
      
      # Ensure prob is a matrix
      if (!is.matrix(prob)) prob <- as.matrix(prob)
      
      if (is.null(colnames(prob))) colnames(prob) <- c("1","2","3")[seq_len(ncol(prob))]
      
      # Find max probability class for each row
      max_indices <- max.col(prob, ties.method = "first")
      cls <- colnames(prob)[max_indices]
      
      as.integer(cls)
    }
    
    # Evaluate Before NP
    out_before <- hnp_summary(baseline_clf, Test, class_col = class_col)
    
    # Train HNP Umbrella
    clf_hnp <- hnp_umbrella(
      S           = Train,
      levels      = levels,
      tolerances  = tolerances,
      method      = method,
      hnp_split   = hnp_split,
      class_col   = class_col
    )
    
    if (is.null(clf_hnp)) {
      warning(sprintf("Run %d failed to train HNP model", r))
      next
    }
    
    # Evaluate After NP
    out_after <- hnp_summary(clf_hnp, Test, class_col = class_col)
    
    # Collect Results
    results_list[[r]] <- data.frame(
      run = r,
      # Row 1
      acc_before = out_before$overall_accuracy,
      acc_after  = out_after$overall_accuracy,
      tot_under_before = out_before$total_under_classification_error,
      tot_under_after  = out_after$total_under_classification_error,
      tot_over_before  = out_before$total_over_classification_error,
      tot_over_after   = out_after$total_over_classification_error,
      
      # Row 2 (FPR)
      fpr1_before = out_before$false_positive_rate[1],
      fpr2_before = out_before$false_positive_rate[2],
      fpr3_before = out_before$false_positive_rate[3],
      fpr1_after  = out_after$false_positive_rate[1],
      fpr2_after  = out_after$false_positive_rate[2],
      fpr3_after  = out_after$false_positive_rate[3],
      
      # Row 3 (FNR)
      fnr1_before = out_before$false_negative_rate[1],
      fnr2_before = out_before$false_negative_rate[2],
      fnr3_before = out_before$false_negative_rate[3],
      fnr1_after  = out_after$false_negative_rate[1],
      fnr2_after  = out_after$false_negative_rate[2],
      fnr3_after  = out_after$false_negative_rate[3],
      
      # Row 4 (Under Error)
      under1_before = out_before$under_classification_error[1],
      under2_before = out_before$under_classification_error[2],
      under3_before = out_before$under_classification_error[3],
      under1_after  = out_after$under_classification_error[1],
      under2_after  = out_after$under_classification_error[2],
      under3_after  = out_after$under_classification_error[3],
      
      # Row 5 (Over Error)
      over1_before  = out_before$over_classification_error[1],
      over2_before  = out_before$over_classification_error[2],
      over3_before  = out_before$over_classification_error[3],
      over1_after   = out_after$over_classification_error[1],
      over2_after   = out_after$over_classification_error[2],
      over3_after   = out_after$over_classification_error[3]
    )
  }
  
  results <- dplyr::bind_rows(results_list)
  if (nrow(results) == 0) stop("No results collected.")
  
  # Plotting
  
  # Helper to draw one boxplot panel
  draw_panel <- function(before_vec, after_vec, title, ylab, 
                         h_alpha = NULL, q_delta = NULL) {
    # Prepare data list
    plot_data <- list(`Before NP` = before_vec, `After NP` = after_vec)
    
    boxplot(plot_data,
            main = title, ylab = ylab,
            col = "white", border = "black",
            outline = FALSE) # Remove outliers as per common preference for clean plots
    
    grid(nx = NA, ny = NULL, lty = 3, col = "gray85")
    
    # Add alpha line (dashed)
    if (!is.null(h_alpha)) {
      abline(h = h_alpha, lty = 2, col = "gray40", lwd = 1.5)
    }
    
    # Add 1-delta quantile point (red dot)
    if (!is.null(q_delta)) {
      # Calculate the (1 - delta) quantile of the 'After' distribution
      q_val <- quantile(after_vec, probs = 1 - q_delta, na.rm = TRUE)
      points(2, q_val, pch = 19, col = "red", cex = 1.4)
    }
  }
  
  # Setup PDF
  pdf(output_file, width = 10, height = 16)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(5, 3), mar = c(3, 4, 3, 1))
  
  # Row 1: Overall Metrics
  draw_panel(results$acc_before, results$acc_after, 
             "Overall Accuracy", "Accuracy")
  draw_panel(results$tot_under_before, results$tot_under_after, 
             "Total Under-classification", "Error")
  draw_panel(results$tot_over_before, results$tot_over_after, 
             "Total Over-classification", "Error")
  
  # Row 2: FPR (Class 1, 2, 3)
  for (k in 1:3) {
    draw_panel(results[[paste0("fpr", k, "_before")]], 
               results[[paste0("fpr", k, "_after")]],
               sprintf("FPR (Class %d)", k), "FPR")
  }
  
  # Row 3: FNR (Class 1, 2, 3)
  for (k in 1:3) {
    draw_panel(results[[paste0("fnr", k, "_before")]], 
               results[[paste0("fnr", k, "_after")]],
               sprintf("FNR (Class %d)", k), "FNR")
  }
  
  # Row 4: Under Error (Class 1, 2, 3)
  # With Alpha Lines and Red Dots (1-delta quantile)
  for (k in 1:3) {
    # Determine alpha and delta for this class
    a_k <- if (k <= length(levels)) levels[k] else NULL
    d_k <- if (k <= length(tolerances)) tolerances[k] else NULL
    
    draw_panel(results[[paste0("under", k, "_before")]],
               results[[paste0("under", k, "_after")]],
               sprintf("Under Error (Class %d)", k), "Under Error",
               h_alpha = a_k, q_delta = d_k)
  }
  
  # Row 5: Over Error (Class 1, 2, 3)
  for (k in 1:3) {
    draw_panel(results[[paste0("over", k, "_before")]], 
               results[[paste0("over", k, "_after")]],
               sprintf("Over Error (Class %d)", k), "Over Error")
  }
  
  mtext("HNP Before vs After Performance (15 Metrics)", outer = TRUE, line = -1.5, cex = 1.1)
  
  message(sprintf("Experiment finished. Boxplots saved to %s", output_file))
}