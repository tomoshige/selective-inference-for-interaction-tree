# Fixed-penalty GLS interaction pruning.
# S=q'Y has variance sigma^2; its line direction is v/d, not q.
source("gls_tree.R")

.pr_children <- function(name) {
  if (name == "root") c("L", "R") else paste0(name, c("L", "R"))
}

.pr_values <- function(values, keys) {
  answer <- numeric(length(keys))
  positions <- match(keys, names(values))
  present <- !is.na(positions)
  answer[present] <- unname(values[positions[present]])
  answer
}

.pr_check_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0) stop("lambda must be one finite nonnegative number")
  invisible(TRUE)
}

prune_tree <- function(tree, Y, lambda) {
  .pr_check_lambda(lambda)
  if (any(!is.finite(Y))) stop("Nonfinite outcomes in pruning")
  if (!length(tree)) {
    empty <- setNames(numeric(), character())
    return(list(kept = character(), F = empty, V = empty, gains = empty))
  }
  labels <- names(tree)
  if (is.null(labels) || anyDuplicated(labels) || !("root" %in% labels)) {
    stop("The grown tree must have unique node addresses and a root")
  }
  gains <- vapply(tree, function(node) {
    if (length(node$u) != length(Y)) stop("Contrast/outcome dimensions differ")
    sum(node$u * Y)^2
  }, numeric(1))
  if (any(!is.finite(gains))) stop("Nonfinite selected gains")
  depths <- vapply(tree, function(node) as.integer(node$depth), integer(1))
  F <- V <- setNames(numeric(length(labels)), labels)
  for (i in order(depths, decreasing = TRUE)) {
    name <- labels[i]
    V[name] <- gains[name] - lambda + sum(.pr_values(F, .pr_children(name)))
    F[name] <- max(0, V[name])
  }
  kept <- character()
  visit <- function(name) {
    if (!(name %in% labels) || V[name] <= 0) return(invisible(NULL))
    kept <<- c(kept, name)
    for (child in .pr_children(name)) visit(child)
    invisible(NULL)
  }
  visit("root")
  list(kept = kept, F = F, V = V, gains = gains)
}

# A supplied pruning object must use this same grown tree, Y, and lambda.
pruning_delta <- function(tree, name, Y, lambda, prune = NULL) {
  .pr_check_lambda(lambda)
  target <- tree[[name]]
  if (is.null(target)) stop("Target is not an internal grown node")
  if (is.null(prune)) prune <- prune_tree(tree, Y, lambda)
  if (!all(names(tree) %in% names(prune$F)) ||
      !all(names(tree) %in% names(prune$gains))) {
    stop("Pruning object does not contain every grown node")
  }
  a <- -lambda + sum(.pr_values(prune$F, .pr_children(name)))
  delta <- -a
  ancestors <- rev(head(target$path, -1L))
  for (ancestor in ancestors) {
    parent <- tree[[ancestor]]
    if (is.null(parent)) stop("Missing ancestor in target path")
    direction <- substr(name, parent$depth + 1L, parent$depth + 1L)
    if (!(direction %in% c("L", "R"))) stop("Invalid target address")
    sibling <- .pr_children(ancestor)[if (direction == "L") 2L else 1L]
    a <- a + unname(prune$gains[ancestor]) - lambda +
         .pr_values(prune$F, sibling)
    delta <- max(delta, -a)
  }
  if (!is.finite(delta)) stop("Nonfinite pruning threshold")
  unname(delta)
}

# Cached grow_support may come from a deeper grown tree whose target path is
# unchanged. It must correspond to this same target contrast and outcome data.
selection_support_pruned <- function(tree, name, Y, lambda, mode = "path",
                                     prune = NULL, grow_support = NULL) {
  mode <- match.arg(mode, c("path", "full_tree"))
  if (is.null(tree[[name]])) stop("Target is not an internal grown node")
  if (is.null(prune)) prune <- prune_tree(tree, Y, lambda)
  if (!(name %in% prune$kept)) stop("Target split was collapsed by pruning")
  delta <- pruning_delta(tree, name, Y, lambda, prune)
  observed <- sum(tree[[name]]$u * Y)
  if (is.null(grow_support)) {
    grown <- selection_support(tree, name, Y, mode)
  } else {
    grown <- grow_support
    if (!is.list(grown) || length(grown$t) != 1L || !is.finite(grown$t) ||
        abs(grown$t - observed) > 1e-10 * max(1, abs(observed)) ||
        !is.matrix(grown$support) || ncol(grown$support) != 2L ||
        !nrow(grown$support) || anyNA(grown$support) ||
        any(grown$support[, 1L] >= grown$support[, 2L])) {
      stop("Invalid or inconsistent cached growth support")
    }
  }
  # Strict endpoints are retained as interval boundaries; Gaussian masses do
  # not depend on endpoint conventions. The observed split must strictly pass.
  if (!(observed^2 > delta)) {
    stop("Retained split is inconsistent with its strict pruning threshold")
  }
  gate <- if (delta < 0) {
    whole_line()
  } else {
    threshold <- sqrt(delta)
    rbind(c(-Inf, -threshold), c(threshold, Inf))
  }
  support <- intersect_intervals(grown$support, gate)
  if (!nrow(support)) stop("Empty retained-path selection support")
  if (!any(observed >= support[, 1L] - 1e-9 &
           observed <= support[, 2L] + 1e-9)) {
    stop("Observed statistic is outside retained-path support")
  }
  endpoints <- support[is.finite(support)]
  list(t = observed, support = support,
       boundary_distance = if (length(endpoints)) min(abs(observed - endpoints)) else Inf,
       delta = delta)
}

infer_pruned <- function(tree, name, Y, mu, lambda, sigma = 1, invert = FALSE,
                         prune = NULL, grow_support = NULL) {
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be one finite positive number")
  }
  if (length(mu) != length(Y) || any(!is.finite(mu))) {
    stop("Invalid conditional mean vector")
  }
  node <- tree[[name]]
  ss <- selection_support_pruned(tree, name, Y, lambda, "path", prune, grow_support)
  truth_unit <- sum(node$u * mu)
  truth <- truth_unit / node$norm
  null_odds <- tn_logodds(ss$t, 0, sigma, ss$support)
  true_odds <- tn_logodds(ss$t, truth_unit, sigma, ss$support)
  ci <- c(NA_real_, NA_real_)
  ci_status <- "not_requested"
  ci_error <- ""
  if (invert) {
    ci <- tryCatch(tn_interval(ss$t, sigma, ss$support), error = function(e) {
      ci_error <<- conditionMessage(e)
      c(NA_real_, NA_real_)
    })
    ci_status <- if (all(is.finite(ci))) "ok" else "failed"
    ci <- ci / node$norm
  }
  data.frame(method = "proposed", node = name, available = TRUE,
    split = node$spec$id, variable = node$spec$j, threshold = node$spec$cut,
    estimate = ss$t / node$norm, truth = truth,
    p_value = 2 * plogis(-abs(null_odds)),
    pivot_truth = plogis(true_odds), covered = abs(true_odds) <= qlogis(.975),
    lower = ci[1L], upper = ci[2L], width = ci[2L] - ci[1L],
    endpoint_covered = if (all(is.finite(ci))) ci[1L] <= truth && truth <= ci[2L] else NA,
    ci_status = ci_status, ci_error = ci_error,
    components = nrow(ss$support), boundary_distance = ss$boundary_distance,
    inference_status = "ok", inference_error = "", delta = ss$delta,
    stringsAsFactors = FALSE)
}
