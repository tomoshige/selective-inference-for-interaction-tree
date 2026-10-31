# GLS interaction trees with exact path selection geometry.
# Y|D ~ N(mu, sigma^2 diag(variance_shape)); shape is fixed, known, and positive.
# U/q are score loadings Wv/d, while direction is v/d. These must not be confused.

empty_intervals <- function() matrix(numeric(), ncol = 2L)
whole_line <- function() matrix(c(-Inf, Inf), ncol = 2L)

intersect_intervals <- function(x, y) {
  if (!nrow(x) || !nrow(y)) return(empty_intervals())
  out <- vector("list", nrow(x) * nrow(y)); k <- 0L
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(y))) {
    lo <- max(x[i, 1], y[j, 1]); hi <- min(x[i, 2], y[j, 2])
    if (lo < hi) { k <- k + 1L; out[[k]] <- c(lo, hi) }
  }
  if (!k) return(empty_intervals())
  ans <- do.call(rbind, out[seq_len(k)])
  ans[order(ans[, 1]), , drop = FALSE]
}

# Solve (a*t+b)*(c*t+d)>=0, retaining all positive-length components.
factor_intervals <- function(a, b, c, d) {
  tol <- 64 * .Machine$double.eps
  # Slopes below rounding error in inner products are treated as zero.
  if (abs(a) < tol) a <- 0
  if (abs(c) < tol) c <- 0
  if (a == 0 && c == 0) {
    return(if (b * d >= 0) whole_line() else empty_intervals())
  }
  if (a == 0 || c == 0) {
    slope <- if (a == 0) b * c else a * d
    offset <- b * d
    if (slope == 0) return(if (offset >= 0) whole_line() else empty_intervals())
    r <- -offset / slope
    return(if (slope > 0) matrix(c(r, Inf), ncol = 2L) else
      matrix(c(-Inf, r), ncol = 2L))
  }
  rr <- sort(c(-b/a, -d/c))
  if (rr[1] == rr[2]) return(if (a * c > 0) whole_line() else empty_intervals())
  if (a * c > 0) rbind(c(-Inf, rr[1]), c(rr[2], Inf)) else
    matrix(rr, ncol = 2L)
}

.gls_shape <- function(variance_shape, n) {
  if (!is.numeric(variance_shape) || length(variance_shape) != n ||
      any(!is.finite(variance_shape)) || any(variance_shape <= 0)) {
    stop("variance_shape must contain one known positive value per participant")
  }
  as.numeric(variance_shape)
}

make_candidates <- function(X, A, B, idx, cuts = c(-0.5, 0, 0.5), min_cell = 8L,
                            variance_shape = rep(1, nrow(X))) {
  n <- nrow(X)
  variance_shape <- .gls_shape(variance_shape, n)
  if (length(idx) < 4L * min_cell) return(NULL)
  unit_shape <- all(variance_shape[idx] == 1)
  sw <- 1 / sqrt(variance_shape[idx])
  D <- cbind(B[idx, , drop = FALSE], A[idx])
  qrD <- qr(if (unit_shape) D else D * sw)
  specs <- list(); vectors <- list(); directions <- list(); norms <- numeric()
  for (j in seq_len(ncol(X))) for (cut in cuts) {
    h <- as.numeric(X[idx, j] <= cut)
    n11 <- sum(A[idx] * h); n10 <- sum(A[idx] * (1 - h))
    n01 <- sum((1 - A[idx]) * h); n00 <- sum((1 - A[idx]) * (1 - h))
    if (min(n11, n10, n01, n00) < min_cell) next
    ht <- qr.resid(qrD, if (unit_shape) h else h * sw)
    gt <- qr.resid(qrD, if (unit_shape) A[idx] * h else A[idx] * h * sw)
    hh <- sum(ht^2)
    v <- if (hh > 1e-12) gt - ht * sum(ht * gt) / hh else gt
    nv <- sqrt(sum(v^2))
    if (nv < 1e-8) next
    u <- direction <- numeric(n)
    u[idx] <- if (unit_shape) v / nv else sw * v / nv
    direction[idx] <- if (unit_shape) v / nv else v / (sw * nv)
    if (length(vectors) && any(vapply(vectors, function(vv) {
      abs(sum(vv * variance_shape * u)) > 1 - 1e-12
    }, logical(1)))) next
    k <- length(vectors) + 1L
    vectors[[k]] <- u; directions[[k]] <- direction; norms[k] <- nv
    specs[[k]] <- list(j = j, cut = cut, id = sprintf("x%d_le_%g", j, cut),
                       left = idx[h == 1], right = idx[h == 0])
  }
  if (!length(vectors)) return(NULL)
  list(U = do.call(cbind, vectors), directions = do.call(cbind, directions),
       norms = norms, specs = specs, idx = idx)
}

fit_tree <- function(X, A, B, Y, max_depth = 2L, cuts = c(-0.5, 0, 0.5),
                     min_cell = 8L, cache = NULL,
                     variance_shape = rep(1, nrow(X))) {
  variance_shape <- .gls_shape(variance_shape, nrow(X))
  # A supplied cache must be restricted to this same fixed design and shape.
  nodes <- list()
  visit <- function(idx, name, depth, ancestors) {
    if (depth >= max_depth) return(invisible(NULL))
    key <- paste(idx, collapse = ",")
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
      cc <- get(key, envir = cache, inherits = FALSE)
    } else {
      cc <- make_candidates(X, A, B, idx, cuts, min_cell, variance_shape)
      if (!is.null(cache)) assign(key, cc, envir = cache)
    }
    if (is.null(cc)) return(invisible(NULL))
    w <- which.max(drop(crossprod(cc$U, Y))^2)
    nm <- if (name == "") "root" else name
    nodes[[nm]] <<- list(name = nm, depth = depth, candidates = cc, winner = w,
                         spec = cc$specs[[w]], u = cc$U[, w],
                         direction = cc$directions[, w], norm = cc$norms[w],
                         path = c(ancestors, nm))
    next_name <- if (name == "") "" else name
    visit(cc$specs[[w]]$left, paste0(next_name, "L"), depth + 1L, c(ancestors, nm))
    visit(cc$specs[[w]]$right, paste0(next_name, "R"), depth + 1L, c(ancestors, nm))
    invisible(NULL)
  }
  visit(seq_len(nrow(X)), "", 0L, character())
  nodes
}

selection_support <- function(tree, target_name, Y, mode = "path") {
  target <- tree[[target_name]]
  u <- target$u; direction <- target$direction
  if (is.null(direction)) stop("GLS node has no perturbation direction")
  t <- sum(u * Y); z <- Y - direction * t
  keys <- switch(mode, path = target$path, node_only = target_name,
                 full_tree = names(tree), stop("Unknown conditioning mode"))
  support <- whole_line()
  for (key in keys) {
    node <- tree[[key]]; U <- node$candidates$U; w <- node$winner
    aa <- drop(crossprod(U, direction)); bb <- drop(crossprod(U, z))
    # Enforce proven algebraic zero slopes rather than treating roundoff-sized
    # nonzero slopes as genuine boundaries at distances around 1e14.
    if (!(key %in% target$path)) {
      aa[] <- 0
      bb <- drop(crossprod(U, Y))
    } else if (key != target_name) {
      target_idx <- target$candidates$idx
      zero <- vapply(node$candidates$specs, function(ss) {
        all(target_idx %in% ss$left) || all(target_idx %in% ss$right)
      }, logical(1))
      # The selected ancestor always contains the target in one child. A
      # competing ancestor split is zeroed only when the same containment holds.
      if (!zero[w]) stop("Ancestor/target regions are inconsistent")
      aa[zero] <- 0
      bb[zero] <- drop(crossprod(U[, zero, drop = FALSE], Y))
    } else {
      aa[w] <- 1
      bb[w] <- 0
    }
    for (j in setdiff(seq_along(aa), w)) {
      ss <- factor_intervals(aa[w] - aa[j], bb[w] - bb[j],
                             aa[w] + aa[j], bb[w] + bb[j])
      support <- intersect_intervals(support, ss)
      if (!nrow(support)) stop("Empty selection support")
    }
  }
  if (!any(t >= support[, 1] - 1e-9 & t <= support[, 2] + 1e-9)) {
    stop("Observed statistic outside selection support")
  }
  list(t = t, support = support,
       boundary_distance = min(abs(t - support[is.finite(support)])))
}

log_sub <- function(a, b) {
  # log(exp(a)-exp(b)), a>=b. -expm1 is accurate near equality.
  out <- a + log(-expm1(pmin(0, b - a)))
  out[is.infinite(b) & b < 0] <- a[is.infinite(b) & b < 0]
  out
}
log_sum <- function(x) {
  if (!length(x) || all(x == -Inf)) return(-Inf)
  m <- max(x); m + log(sum(exp(x - m)))
}
normal_log_mass <- function(lo, hi, mean = 0, sd = 1) {
  if (sd <= 0) stop("sd must be positive")
  a <- (lo - mean)/sd; b <- (hi - mean)/sd
  out <- rep(-Inf, length(a)); ok <- a < b
  right <- ok & a >= 0
  left <- ok & !right
  if (any(right)) out[right] <- log_sub(
    pnorm(a[right], lower.tail = FALSE, log.p = TRUE),
    pnorm(b[right], lower.tail = FALSE, log.p = TRUE))
  if (any(left)) out[left] <- log_sub(pnorm(b[left], log.p = TRUE),
                                    pnorm(a[left], log.p = TRUE))
  out
}

# Log odds of the truncated CDF; avoids subtracting a CDF from one.
tn_logodds <- function(t, mean, sd, support) {
  lower <- normal_log_mass(support[, 1], pmin(t, support[, 2]), mean, sd)
  upper <- normal_log_mass(pmax(t, support[, 1]), support[, 2], mean, sd)
  ans <- log_sum(lower) - log_sum(upper)
  if (is.nan(ans)) stop("Undefined truncated-normal probability")
  ans
}

tn_interval <- function(t, sd, support, alpha = 0.05) {
  solve_mean <- function(prob) {
    target <- qlogis(prob)
    f <- function(m) tn_logodds(t, m, sd, support) - target
    lo <- t - sd; hi <- t + sd
    for (k in 0:35) {
      flo <- f(lo); fhi <- f(hi)
      if (is.finite(flo) && is.finite(fhi) && flo >= 0 && fhi <= 0) break
      radius <- sd * 2^(k + 1)
      lo <- t - radius; hi <- t + radius
      if (k == 35L) stop("Confidence endpoint could not be bracketed")
    }
    uniroot(f, c(lo, hi), tol = 1e-7 * sd, maxiter = 150L)$root
  }
  c(lower = solve_mean(1 - alpha/2), upper = solve_mean(alpha/2))
}

infer_node <- function(tree, name, Y, mu, mode, sigma = 1, invert = TRUE) {
  node <- tree[[name]]
  t <- sum(node$u * Y); theta_u <- sum(node$u * mu)
  theta <- theta_u/node$norm
  if (mode == "naive") {
    support <- whole_line(); boundary <- Inf
  } else {
    ss <- selection_support(tree, name, Y, mode)
    support <- ss$support; boundary <- ss$boundary_distance
  }
  l0 <- tn_logodds(t, 0, sigma, support)
  lt <- tn_logodds(t, theta_u, sigma, support)
  p <- 2 * plogis(-abs(l0))
  cover <- abs(lt) <= qlogis(0.975)
  ci <- c(NA_real_, NA_real_); ci_status <- "not_requested"; ci_error <- ""
  if (invert) {
    ci <- tryCatch(tn_interval(t, sigma, support), error = function(e) {
      ci_error <<- conditionMessage(e); c(NA_real_, NA_real_)
    })
    ci_status <- if (all(is.finite(ci))) "ok" else "failed"
    ci <- ci/node$norm
  }
  data.frame(method = mode, node = name, available = TRUE,
    split = node$spec$id, variable = node$spec$j, threshold = node$spec$cut,
    estimate = t/node$norm, truth = theta, p_value = p,
    pivot_truth = plogis(lt), covered = cover,
    lower = ci[1], upper = ci[2], width = ci[2] - ci[1],
    endpoint_covered = if (all(is.finite(ci))) ci[1] <= theta && theta <= ci[2] else NA,
    ci_status = ci_status, ci_error = ci_error,
    components = nrow(support), boundary_distance = boundary,
    inference_status = "ok", inference_error = "", stringsAsFactors = FALSE)
}

unavailable_row <- function(method, name, reason, status = "unavailable") {
  data.frame(method = method, node = name, available = FALSE,
    split = NA_character_, variable = NA_integer_, threshold = NA_real_,
    estimate = NA_real_, truth = NA_real_, p_value = NA_real_,
    pivot_truth = NA_real_, covered = NA,
    lower = NA_real_, upper = NA_real_, width = NA_real_, endpoint_covered = NA,
    ci_status = "not_requested", ci_error = "", components = NA_integer_,
    boundary_distance = NA_real_, inference_status = status,
    inference_error = reason, stringsAsFactors = FALSE)
}

holdout_node <- function(training_tree, name, X, A, B, min_cell = 1L,
                          variance_shape = rep(1, nrow(X))) {
  variance_shape <- .gls_shape(variance_shape, nrow(X))
  if (is.null(training_tree[[name]])) return(NULL)
  idx <- seq_len(nrow(X))
  for (ancestor in head(training_tree[[name]]$path, -1)) {
    ss <- training_tree[[ancestor]]$spec
    # Child direction is the character immediately after the ancestor's address.
    depth <- training_tree[[ancestor]]$depth
    direction <- substr(name, depth + 1L, depth + 1L)
    inside_left <- X[idx, ss$j] <= ss$cut
    idx <- idx[if (direction == "L") inside_left else !inside_left]
  }
  ss <- training_tree[[name]]$spec
  h <- as.numeric(X[idx, ss$j] <= ss$cut)
  if (length(idx) < ncol(B) + 4L ||
      min(table(factor(A[idx], 0:1), factor(h, 0:1))) < min_cell) return(NULL)
  C <- cbind(B[idx, , drop = FALSE], A[idx], h)
  unit_shape <- all(variance_shape[idx] == 1)
  sw <- 1 / sqrt(variance_shape[idx])
  v <- if (unit_shape) qr.resid(qr(C), A[idx] * h) else
    qr.resid(qr(C * sw), A[idx] * h * sw)
  nv <- sqrt(sum(v^2)); if (nv < 1e-8) return(NULL)
  u <- direction <- numeric(nrow(X))
  u[idx] <- if (unit_shape) v/nv else sw*v/nv
  direction[idx] <- if (unit_shape) v/nv else v/(sw*nv)
  list(u = u, direction = direction, norm = nv, spec = ss)
}
