# HNPclassifier 0.2.0

## Breaking changes

* `hnp_umbrella()` has a new interface. It now takes separate `X` and `Y`
  arguments together with `importance_order` instead of a single data frame
  `S` and a `class_col`. It supports an arbitrary number of ordered classes
  (`T >= 2`) rather than only ternary classification.
* `hnp_summary()` now takes `classifier`, `X`, `Y` and `importance_order`
  instead of `data` and `class_col`, and accepts classifiers that return class
  labels, probability matrices or score matrices, as well as fitted model
  objects.
* `hnp_map_classes()` now accepts a variable number of class labels via `...`
  (in decreasing priority order) instead of the fixed `class_1`, `class_2`,
  `class_3` arguments.
* Removed `probability_to_score_1()`, `probability_to_score_2()`,
  `hnp_umbrella_flex()` and `hnp_box_plot()`.

## New features

* Generalized the H-NP umbrella algorithm to any number of ordered classes
  (`T >= 2`).
* Added support for user-supplied pretrained models and pre-computed score
  matrices in `hnp_umbrella()` via the `pretrained_model` and `input_is_score`
  arguments.
* Added grid search over candidate thresholds (`grid_search`, `grid_set`,
  `max_grid`, `max_combinations`) to minimize the weighted misclassification
  objective, with a recursive multi-class threshold search.
* Added `hnp_boxplot()` to visualize and summarize under-classification and
  overall error from confusion matrices, supporting single- and two-method
  comparisons.
* Added data-generation helpers for examples and simulation
  (`gen_data()`, `gen_normal_data()`, `generate_ball_data()`, etc.) and a
  neural-network scoring helper (`train_nn_and_get_scores()`).

## Improvements and bug fixes

* `hnp_delta_search()` now uses `stats::pbinom()` instead of an explicit
  `choose()`-based summation, avoiding numerical overflow and underflow for
  large samples.
* `hnp_upper_bound()` is now vectorized over the score functions and adds
  validation for non-finite, `NA` and length-mismatched score outputs.
* Added stricter input validation across the public functions (class labels,
  control levels, tolerances and split ratios).
* Added internal label mapping so that original class labels are preserved in
  predictions and summaries.
