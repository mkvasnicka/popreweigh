#' Creates the solver model for reweighing
#'
#' @param F numeric feature matrix (n x k) where n is the number of observations
#'   and k is the number of features where n >> p.
#' @param w numeric vector of original weights (length n).
#' @param p numeric vector of desired prevalence (length k); if the features are
#'   binary, p must be in [0, 1]; if the features are continuous, p may be any
#'   real number; in this case, it is not prevalence but weighted mean.
#' @param lambda penalty on errors (scalar or vector of length k).
#' @param lb_weights lower bound for the weights (scalar of vector of length n).
#' @param ub_weights upper bound for the weights (scalar of vector of length n).
#' @param lb_errors lower bound for the errors (scalar or vector of length k).
#' @param ub_errors upper bound for the errors (scalar or vector of length k).
#' @param scaling scaling method for weights and features; one of "none",
#'   "weights", or "auto"; see details.
#' @param weight_scale scaling factor for the weights (scalar).
#'
#' @return model specification for cohort reweighing.
#'
#' @details
#' - scaling
#'   - equal to "none" means no scaling;
#'   - equal to "weights" means that weights are scaled so that their
#'     average is equal to weight_scale; features are not scaled;
#'   - equal to "auto" means that weights are scaled so that the average
#'     weight is equal to weight_scale if given or to 1 if n <= 1000 or to
#'     sqrt(n / 1000) if n > 1000; features are standardized to have mean 0 and
#'     standard deviation 1; p and bounds are scaled accordingly
#' - if a feature is continuous instead of binary, its corresponding prevalence
#'   is not a prevalence but its weighted mean; increased bounds may be necessary
#' - BEWARE: many continuous features may require big changes in weights;
#'   some weights may become huge, i.e., a great part of the new cohort may be
#'   represented by few individuals and the error in prevalence may still be huge
#'
#' @section Typical use:
#' - for binary features:
#'   `solver_model(F, W, p)` or
#'   `solver_model(F, W, p, scaling = "weights")` when n is large
#' - for mix of binary and continuous features or continuous features with very
#'   different means/ranges:
#'   `solver_model(F, W, p, scaling = "auto")`
#' - with many continuous features, lower lb_errors and increase ub_errors but
#'   be careful
#' - if several individuals have too high weights, lower ub_weights
#' - if no solution is found, lower lb_errors and increase ub_errors
#' @keywords internal
solver_model <- function(
  F, w, p,
  lambda = 1e3L,
  lb_weights = 0,
  ub_weights = 1,
  lb_errors = -1,
  ub_errors = 1,
  scaling = c("none", "weights", "auto"),
  weight_scale = NULL
) {
  # problem size
  n <- length(w)
  k <- ncol(F)
  stopifnot(nrow(F) == n)
  stopifnot(length(p) == k)
  # if weight scale is not given, it is set so that mean weight is 1
  # if the number of individuals is is max 1000; then it's multiplied
  # by increased slower
  stopifnot(
    is.null(weight_scale) || (
      (is.numeric(weight_scale) && length(weight_scale) == 1L &&
        weight_scale > 0)
    )
  )
  # stopifnot(!(scaling == "none" && is.null(weight_scale)))
  weight_scale <- if (is.null(weight_scale)) {
    min(c(n, sqrt(n * 1000)))
  } else {
    weight_scale
  }
  # feature scaling
  fidentity <- function(x, f) x
  scaling <- match.arg(scaling)
  scaling <- switch(
    scaling,
    none = list(a = 1, b = identity, c = fidentity),
    # weights ... weights are multiplied by weight_scale, F_js not scaled
    weights = list(
      a = weight_scale,
      b = identity,
      c = fidentity
    ),
    # auto ... weights are multiplied by weight_scale, F_js are standardized
    auto = list(
      a = weight_scale,
      b = function(x) (x - mean(x)) / sd(x),
      c = function(x, f) (x - mean(f)) / sd(f)
    )
  )
  # lambdas
  stopifnot(
    is.numeric(lambda) &&
      ((length(lambda) == 1L || length(lambda) == k)) &&
      all(lambda > 0)
  )
  lambda <- rep(lambda, length.out = k)
  # bounds
  stopifnot(length(lb_weights) == 1L || length(lb_weights) == n)
  stopifnot(is.numeric(lb_weights) && all(lb_weights >= 0))
  stopifnot(length(ub_weights) == 1L || length(ub_weights) == n)
  stopifnot(
    is.numeric(ub_weights) && all(ub_weights >= lb_weights) &&
      all(ub_weights <= 1)
  )
  lb_weights <- rep(lb_weights, length.out = n)
  ub_weights <- rep(ub_weights, length.out = n)
  stopifnot(length(lb_errors) == 1L || length(lb_errors) == k)
  stopifnot(is.numeric(lb_errors))
  stopifnot(length(ub_errors) == 1L || length(ub_errors) == k)
  stopifnot(is.numeric(ub_errors) && all(ub_errors >= lb_errors))
  lb_errors <- rep(lb_errors, length.out = k)
  ub_errors <- rep(ub_errors, length.out = k)

  # scaling
  p <- purrr::imap_dbl(p, ~scaling$c(.x, F[, .y])) * scaling$a
  F <- apply(F, 2, scaling$b)
  w <- w * scaling$a
  lb_weights <- lb_weights * scaling$a
  ub_weights <- ub_weights * scaling$a
  lb_errors <- lb_errors * scaling$a
  ub_errors <- ub_errors * scaling$a
  # objective function
  # - quadratic part (sparse matrix)
  Q <- as(Matrix::Diagonal(n + k), "generalMatrix")
  diag(Q) <- c(rep(1L, n), lambda)
  # - linear part
  obj <- c(-2 * w, rep(0, k))
  # - offset
  alpha <- sum(w^2)
  # constraints
  # - linear equality
  A <- t(rbind(cbind(F, 1), cbind(diag(k), 0)))
  rhs <- c(p, scaling$a)
  # model
  list(
    Q = Q,
    obj = obj,
    alpha = alpha,
    A = A,
    rhs = rhs,
    sense = rep("=", length(rhs)),
    lb = c(lb_weights, lb_errors),
    ub = c(ub_weights, ub_errors),
    .setting = list(n = n, k = k, weight_scale = scaling$a)
  )
}

#' Reweighs with Gurobi
#'
#' @param model The model created by `solver_model`.
#' @param verbose If `TRUE`, solver output is printed.
#' @param ... Additional parameters for the Gurobi solver.
#'
#' @return A list with the optimization results.
#' @keywords internal
reweigh_gurobi <- function(model, verbose, ...) {
  params <- list(...)
  if (!verbose) params$OutputFlag <- 0
  result <- gurobi::gurobi(model, params = params)
  list(
    status = tolower(result$status),
    weights = result$x[1:model$.setting$n],
    prevalence_error = result$x[- (1:model$.setting$n)],
    objective_value = result$objval,
    solver = result
  )
}

#' Reweighs with HiGHS
#'
#' @param model The model created by `solver_model`.
#' @param verbose If `TRUE`, solver output is printed.
#' @param ... Additional parameters for the HiGHS solver.
#'
#' @return A list with the optimization results.
#' @keywords internal
reweigh_highs <- function(model, verbose, ...) {
  cntrl <- list(...)
  if (verbose) cntrl$log_to_console <- TRUE
  cntrl$solver <- "choose"  # necessary to use QP solver
  result <- highs::highs_solve(
    Q = model$Q,
    L = model$obj,
    A = model$A,
    lhs = model$rhs,
    rhs = model$rhs,
    lower = model$lb,
    upper = model$ub,
    offset = model$alpha,
    control = cntrl
  )
  list(
    status = tolower(result$status_message),
    weights = result$primal_solution[1:model$.setting$n],
    prevalence_error = result$primal_solution[- (1:model$.setting$n)],
    objective_value = result$objective_value,
    solver = result
  )
}

#' Finds new weights that are as close as possible to the original weights to
#' achieve the desired prevalence
#'
#' @param F A numeric feature matrix (n x k) or a data.frame. If a data.frame,
#'   it may contain the weight column specified by `w`.
#' @param w A numeric vector of original weights (length n). If `F` is a
#'   data.frame, `w` can be the name of the weight column in `F`.
#' @param p A numeric vector of desired prevalences (length k).
#' @param lambda The penalty on errors (a scalar or a vector of length k).
#' @param lb_weights The lower bound for the weights (a scalar or a vector of
#'   length n).
#' @param ub_weights The upper bound for the weights (a scalar or a vector of
#'   length n).
#' @param lb_errors The lower bound for the errors (a scalar or a vector of
#'   length k).
#' @param ub_errors The upper bound for the errors (a scalar or a vector of
#'   length k).
#' @param scaling The scaling method for weights and features. One of "none",
#'   "weights", or "auto". See [solver_model()] for details.
#' @param weight_scale The scaling factor for the weights.
#' @param solver The solver to use, either "gurobi" or "highs". The solver must
#'   be installed.
#' @param verbose If `TRUE`, solver output is printed. Default is `FALSE`.
#' @param ... Additional parameters passed to the solver.
#'
#' @return A list with the following elements:
#'   - `status`: The optimization status.
#'   - `weights`: A numeric vector of new weights (length n).
#'   - `prevalence_error`: A numeric vector of prevalence errors (length k).
#'   - `objective_value`: The numeric value of the objective function.
#'   - `solver`: The original solver result.
#' @export
reweigh <- function(
  F, w, p,
  lambda = 1e3L,
  lb_weights = 0,
  ub_weights = 1,
  lb_errors = -1,
  ub_errors = 1,
  scaling = NULL,
  weight_scale = NULL,
  solver = c("gurobi", "highs"),
  verbose = FALSE,
  ...
) {
  solver <- match.arg(solver)
  if (inherits(F, "data.frame")) {
    if (is.character(w)) {
      wn <- w
      w <- F[[w]]
      F <- F |> select(-all_of(wn))
    }
    F <- F |>
      mutate(across(everything(), as.numeric)) |>
      as.matrix()
    stopifnot(!any(is.na(F)))
  }
  model <- solver_model(
    F, w, p,
    lambda = lambda,
    lb_weights = lb_weights,
    ub_weights = ub_weights,
    lb_errors = lb_errors,
    ub_errors = ub_errors,
    scaling = scaling,
    weight_scale = weight_scale
  )
  result <- switch(
    solver,
    gurobi = reweigh_gurobi(model, verbose, ...),
    highs = reweigh_highs(model, verbose, ...)
  )
  if (result$status != "optimal") {
    warning("Optimization failed")
  }
  # de-scale the weights
  result$weights <- result$weights / model$.setting$weight_scale
  result$prevalence_error <-
    result$prevalence_error / model$.setting$weight_scale
  # # de-scale the objective value
  result$objective_value <-
    result$objective_value / model$.setting$weight_scale^2
  # return
  result
}


#' Reweigh a donor pool to match target prevalences
#'
#' This function reweighs a donor pool to match target prevalences for a set of
#' variables. The target prevalences can be specified directly or as trends
#' applied to original prevalences.
#'
#' @param donors A data.frame representing the donor pool. It must contain the
#'   `weight_var` and the variables for which prevalences are to be matched.
#' @param weight_var A character string specifying the name of the column in
#'   `donors` that contains the weights. Defaults to "weight".
#' @param original_prevalences A named numeric vector of original prevalences.
#'   The names must correspond to columns in `donors`. If `NULL` (the default),
#'   the original prevalences are calculated from the `donors` data.frame.
#' @param target_prevalences A named numeric vector of target prevalences. The
#'   names must correspond to columns in `donors`. Either `target_prevalences`
#'   or `trends` must be provided.
#' @param trends A named numeric vector of trends. Trends are multipliers
#'   applied to the `original_prevalences` to obtain the `target_prevalences`.
#'   Either `target_prevalences` or `trends` must be provided.
#' @param population The target population size. Default is 1.
#' @param solver A character string specifying the solver to be used for the
#'   reweighing optimization. Defaults to "gurobi".
#' @param verbose A logical value indicating whether to print the prevalence
#'   errors after reweighing. Defaults to `TRUE`.
#' @param ... Additional arguments passed to the `reweigh()` function.
#'
#' @return A data.frame identical to the input `donors` but with the
#'   `weight_var` column updated with the new weights.
#' @export
# TODO: not working with HiGHS solver currently
# TODO: improve printed diagnostics
reweigh_fem <- function(
  donors,
  population = 1,
  weight_var = "weight",
  original_prevalences = NULL,
  target_prevalences = NULL,
  trends = NULL,
  solver = c("gurobi", "highs"),
  verbose = TRUE,
  ...
) {
  # get solver
  solver <- match.arg(solver)
  # check we have donors
  stopifnot(
    inherits(donors, "data.frame"),
    nrow(donors) > 0
  )
  # check trends original and target prevalences
  is_null_or_named_numeric_in_donors <- function(x, donors) {
    is.null(x) || (
      is.numeric(x) &&
        !is.null(names(x)) &&
        all(names(x) %in% names(donors))
    )
  }
  stopifnot(
    is_null_or_named_numeric_in_donors(original_prevalences, donors),
    is_null_or_named_numeric_in_donors(target_prevalences, donors),
    is_null_or_named_numeric_in_donors(trends, donors)
  )
  # check that trends and original and target prevalences have the same
  # variables
  same_variables_in <- function(a, b) {
    if (is.null(a) || is.null(b)) {
      return(TRUE)
    }
    all(sort(names(a)) == sort(names(b)))
  }
  stopifnot(
    same_variables_in(original_prevalences, trends),
    same_variables_in(target_prevalences, trends),
    same_variables_in(original_prevalences, target_prevalences)
  )
  # check population
  stopifnot(is.numeric(population), length(population) == 1, population > 0)
  # check weight variable
  stopifnot(
    is.character(weight_var),
    length(weight_var) == 1,
    weight_var %in% names(donors)
  )
  # get trended variables
  trend_vars <- if (!is.null(trends)) {
    names(trends)
  } else if (!is.null(target_prevalences)) {
    names(target_prevalences)
  } else {
    names(original_prevalences)
  }
  # get F matrix
  F <- as.matrix(donors[ , trend_vars])
  # get weights
  w <- donors[[weight_var]] / sum(donors[[weight_var]])
  # get original prevalences
  p_ori <- if (is.null(original_prevalences)) {
    setNames(as.double(w %*% F), trend_vars)
  }  else {
    original_prevalences[trend_vars]
  }
  # check that original prevalences are valid
  # get target prevalences
  if (!is.null(target_prevalences)) {
    p <- target_prevalences[trend_vars]
  } else {
    p <- setNames(p_ori * trends, trend_vars)
  }
  # reweigh
  res <- reweigh(F = F, w = w, p = p, solver = solver, verbose = verbose, ...)
  #
  if (res$status != "optimal") {
    stop("Reweighing did not succeed: ", res$status)
  }
  if (verbose) {
    cat("Prevalence errors:\n")
    cat(paste0(trend_vars, ": ", res$prevalence_error, collapse = "\n"))
  }
  donors[[weight_var]] <- res$weights * population
  donors
}
