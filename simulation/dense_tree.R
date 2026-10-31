# Dense, design-only candidate search for the 2026-10-31 study.
# Source gls_pruning.R first, then this file. The original routines are kept for
# independent comparisons. Both score and perturbation matrices retain the API.
if (!exists(".dense_original_make_candidates", inherits = FALSE)) {
  .dense_original_make_candidates <- make_candidates
  .dense_original_fit_tree <- fit_tree
  .dense_original_selection_support <- selection_support
}

.dense_cutpoints <- function(x) {
  z <- sort(unique(x))
  if (length(z) < 2L) return(numeric())
  lo <- head(z, -1L); hi <- tail(z, -1L)
  ans <- lo / 2 + hi / 2
  # A representable midpoint can round to the upper value. Retain the same
  # order-statistic partition by using its lower value in that case.
  ans[ans >= hi] <- lo[ans >= hi]
  ans[ans < lo] <- lo[ans < lo]
  ans
}

# Screen near-parallel residual directions without forming their K x K Gram
# matrix. If |q_i'q_j| > 1-1e-12, unit q_i and +/-q_j are within ~1.42e-6.
# A unit-vector projection therefore provides a necessary 2e-6 window. Actual
# correlations within that window use the original duplicate criterion, in the
# original feature/cut priority order. The wider window accommodates roundoff.
.dense_unique_directions <- function(Q) {
  k <- ncol(Q)
  if (k < 2L) return(seq_len(k))
  rr <- sin(seq_len(nrow(Q)) * sqrt(2)) + cos(seq_len(nrow(Q)) * sqrt(3))
  rr <- rr / sqrt(sum(rr^2))
  pp <- abs(drop(crossprod(rr, Q)))
  oo <- order(pp, seq_len(k)); sorted <- pp[oo]
  lower <- findInterval(pp - 2e-6, sorted) + 1L
  upper <- findInterval(pp + 2e-6, sorted)
  keep <- rep(TRUE, k)
  for (j in seq_len(k)) {
    if (lower[j] > upper[j]) next
    previous <- oo[seq.int(lower[j], upper[j])]
    previous <- previous[previous < j & keep[previous]]
    if (!length(previous)) next
    cc <- drop(crossprod(Q[, j], Q[, previous, drop = FALSE]))
    if (any(abs(cc) > 1 - 1e-12)) keep[j] <- FALSE
  }
  which(keep)
}

make_candidates <- function(X, A, B, idx, cuts = NULL, min_cell = 8L,
                            variance_shape = rep(1, nrow(X))) {
  n <- nrow(X); m <- length(idx)
  variance_shape <- .gls_shape(variance_shape, n)
  if (m < 4L * min_cell) return(NULL)
  if (any(!is.finite(X[idx, , drop = FALSE])) || anyNA(A[idx]) ||
      any(!(A[idx] %in% c(0, 1))) || any(!is.finite(B[idx, , drop = FALSE]))) {
    stop("Nonfinite design or nonbinary treatment")
  }
  aa <- A[idx]; sw <- 1 / sqrt(variance_shape[idx])
  hh <- list(); jj <- list(); cut_values <- list(); ranks <- list()
  for (j in seq_len(ncol(X))) {
    x <- X[idx, j]
    cc <- if (is.null(cuts)) .dense_cutpoints(x) else sort(unique(cuts))
    if (!length(cc)) next
    H <- outer(x, cc, "<=") * 1
    n1 <- drop(crossprod(aa, H)); n0 <- colSums(H) - n1
    valid <- pmin(n1, n0, sum(aa) - n1, m - sum(aa) - n0) >= min_cell
    if (!any(valid)) next
    hh[[j]] <- H[, valid, drop = FALSE]
    jj[[j]] <- rep.int(j, sum(valid)); cut_values[[j]] <- cc[valid]
    ranks[[j]] <- as.integer(colSums(H[, valid, drop = FALSE]))
  }
  if (!length(hh)) return(NULL)
  H <- do.call(cbind, hh)
  if (is.null(H) || !ncol(H)) return(NULL)
  feature <- unlist(jj, use.names = FALSE)
  cut <- unlist(cut_values, use.names = FALSE)
  rank <- unlist(ranks, use.names = FALSE)
  D <- cbind(B[idx, , drop = FALSE], aa) * sw
  qrD <- qr(D)
  Q <- qr.Q(qrD)[, seq_len(qrD$rank), drop = FALSE]
  ht <- H * sw; gt <- ht * aa
  ht <- ht - Q %*% crossprod(Q, ht)
  gt <- gt - Q %*% crossprod(Q, gt)
  hnorm2 <- colSums(ht^2)
  coefficient <- numeric(ncol(H))
  nonzero <- hnorm2 > 1e-12
  coefficient[nonzero] <- colSums(ht[, nonzero, drop = FALSE] *
    gt[, nonzero, drop = FALSE]) / hnorm2[nonzero]
  v <- gt - sweep(ht, 2L, coefficient, "*")
  norms <- sqrt(colSums(v^2))
  identifiable <- norms >= 1e-8
  if (!any(identifiable)) return(NULL)
  H <- H[, identifiable, drop = FALSE]
  v <- v[, identifiable, drop = FALSE]
  norms <- norms[identifiable]
  feature <- feature[identifiable]; cut <- cut[identifiable]; rank <- rank[identifiable]
  whitened <- sweep(v, 2L, norms, "/")
  keep <- .dense_unique_directions(whitened)
  H <- H[, keep, drop = FALSE]; whitened <- whitened[, keep, drop = FALSE]
  norms <- norms[keep]; feature <- feature[keep]; cut <- cut[keep]; rank <- rank[keep]
  k <- length(keep)
  U <- directions <- matrix(0, n, k)
  U[idx, ] <- whitened * sw
  directions[idx, ] <- whitened / sw
  specs <- lapply(seq_len(k), function(i) {
    list(j = feature[i], cut = cut[i],
         id = sprintf("x%d_r%d_le_%.17g", feature[i], rank[i], cut[i]),
         left = idx[H[, i] == 1], right = idx[H[, i] == 0])
  })
  list(U = U, directions = directions, norms = norms, specs = specs, idx = idx,
       raw_candidates = sum(vapply(hh, function(z) if (is.null(z)) 0L else ncol(z), integer(1))),
       identifiable_candidates = ncol(v), candidate_policy = "all_midpoints",
       feature = feature, cuts = cut, design_X = X)
}

.dense_materialized_make_candidates <- make_candidates

# Prefix representation of all candidate residual directions. It avoids the
# n x K matrices in production; materialization below is for validation only.
make_candidates <- function(X, A, B, idx, cuts = NULL, min_cell = 8L,
                            variance_shape = rep(1, nrow(X))) {
  n <- nrow(X); m <- length(idx)
  variance_shape <- .gls_shape(variance_shape, n)
  if (m < 4L * min_cell) return(NULL)
  if (any(!is.finite(X[idx, , drop = FALSE])) || anyNA(A[idx]) ||
      any(!(A[idx] %in% c(0, 1))) || any(!is.finite(B[idx, , drop = FALSE]))) {
    stop("Nonfinite design or nonbinary treatment")
  }
  aa <- A[idx]; sw <- 1 / sqrt(variance_shape[idx]); weights <- sw^2
  qrD <- qr(cbind(B[idx, , drop = FALSE], aa) * sw)
  Q <- qr.Q(qrD)[, seq_len(qrD$rank), drop = FALSE]
  orders <- vector("list", ncol(X)); feature <- rank <- numeric()
  cut <- sumw <- sumwa <- numeric(); qh <- qz <- list(); part <- 0L
  for (j in seq_len(ncol(X))) {
    oo <- order(X[idx, j], seq_len(m)); orders[[j]] <- oo
    xx <- X[idx[oo], j]
    cc <- if (is.null(cuts)) .dense_cutpoints(xx) else sort(unique(cuts))
    if (!length(cc)) next
    rr <- findInterval(cc, xx)
    take <- rr > 0L & rr < m
    cc <- cc[take]; rr <- rr[take]
    if (!length(cc)) next
    n1 <- cumsum(aa[oo])[rr]; n0 <- rr - n1
    valid <- pmin(n1, n0, sum(aa) - n1, m - sum(aa) - n0) >= min_cell
    rr <- rr[valid]; cc <- cc[valid]
    if (!length(cc)) next
    cw <- cumsum(weights[oo]); cwa <- cumsum(weights[oo] * aa[oo])
    cqh <- apply(Q[oo, , drop = FALSE] * sw[oo], 2L, cumsum)
    cqz <- apply(Q[oo, , drop = FALSE] * (sw[oo] * aa[oo]), 2L, cumsum)
    if (is.null(dim(cqh))) cqh <- matrix(cqh, ncol = 1L)
    if (is.null(dim(cqz))) cqz <- matrix(cqz, ncol = 1L)
    part <- part + 1L
    qh[[part]] <- t(cqh[rr, , drop = FALSE]); qz[[part]] <- t(cqz[rr, , drop = FALSE])
    feature <- c(feature, rep.int(j, length(rr))); rank <- c(rank, rr)
    cut <- c(cut, cc); sumw <- c(sumw, cw[rr]); sumwa <- c(sumwa, cwa[rr])
  }
  if (!length(feature)) return(NULL)
  qh <- do.call(cbind, qh); qz <- do.call(cbind, qz)
  hh <- sumw - colSums(qh^2)
  hz <- sumwa - colSums(qh * qz)
  gg <- sumwa - colSums(qz^2)
  coefficient <- numeric(length(feature)); nonzero <- hh > 1e-12
  coefficient[nonzero] <- hz[nonzero] / hh[nonzero]
  norm2 <- gg - coefficient * hz
  qv <- qz - sweep(qh, 2L, coefficient, "*")
  # Near zero, squared-norm subtraction can lose relative accuracy. Recompute
  # such candidates directly, including cases where h lies in the base span.
  uncertain <- which(norm2 < 1e-10 * pmax(1, sumwa))
  for (j in uncertain) {
    h <- numeric(m); h[orders[[feature[j]]][seq_len(rank[j])]] <- 1
    ht <- qr.resid(qrD, h * sw); gt <- qr.resid(qrD, aa * h * sw)
    h2 <- sum(ht^2)
    coefficient[j] <- if (h2 > 1e-12) sum(ht * gt) / h2 else 0
    vv <- gt - coefficient[j] * ht
    norm2[j] <- sum(vv^2)
    qv[, j] <- qz[, j] - coefficient[j] * qh[, j]
  }
  identifiable <- which(norm2 >= 1e-16)
  if (!length(identifiable)) return(NULL)
  cc <- list(idx = idx, n = n, sw = sw, weights = weights, A_local = aa, Q = Q,
    feature = as.integer(feature[identifiable]), ranks = as.integer(rank[identifiable]),
    cuts = cut[identifiable], coefficient = coefficient[identifiable],
    qv = qv[, identifiable, drop = FALSE], norms = sqrt(norm2[identifiable]),
    orders = orders, design_X = X, candidate_policy = "all_midpoints",
    raw_candidates = length(feature), identifiable_candidates = length(identifiable))
  # A projection window is a necessary condition for the original correlation
  # duplicate rule; correlations are evaluated from weighted intersections.
  rr <- sin(seq_len(m) * sqrt(2)) + cos(seq_len(m) * sqrt(3))
  rr <- rr / sqrt(sum(rr^2))
  projection_y <- numeric(n); projection_y[idx] <- rr / sw
  pp <- abs(.dense_candidate_scores(cc, projection_y))
  oo <- order(pp, seq_along(pp)); sorted <- pp[oo]
  lower <- findInterval(pp - 2e-6, sorted) + 1L
  upper <- findInterval(pp + 2e-6, sorted)
  keep <- rep(TRUE, length(pp))
  for (j in seq_along(pp)) {
    if (lower[j] > upper[j]) next
    previous <- oo[seq.int(lower[j], upper[j])]
    previous <- previous[previous < j & keep[previous]]
    if (!length(previous)) next
    hj <- orders[[cc$feature[j]]][seq_len(cc$ranks[j])]
    for (i in previous) {
      hi <- orders[[cc$feature[i]]][seq_len(cc$ranks[i])]
      both <- intersect(hi, hj)
      covariance <- (1 - cc$coefficient[i] - cc$coefficient[j]) * sum(weights[both] * aa[both]) +
        cc$coefficient[i] * cc$coefficient[j] * sum(weights[both]) -
        sum(cc$qv[, i] * cc$qv[, j])
      correlation <- covariance / (cc$norms[i] * cc$norms[j])
      if (abs(correlation) > 1 - 1e-12) { keep[j] <- FALSE; break }
    }
  }
  keep <- which(keep)
  for (field in c("feature", "ranks", "cuts", "coefficient", "norms")) cc[[field]] <- cc[[field]][keep]
  cc$qv <- cc$qv[, keep, drop = FALSE]
  ids <- sprintf("x%d_r%d_le_%.17g", cc$feature, cc$ranks, cc$cuts)
  cc$specs <- lapply(seq_along(keep), function(i) {
    list(j = cc$feature[i], cut = cc$cuts[i], id = ids[i])
  })
  cc
}

.dense_candidate_scores <- function(cc, Y) {
  if (!is.null(cc$U)) return(drop(crossprod(cc$U, Y)))
  yr <- cc$sw * Y[cc$idx]
  yr <- yr - drop(cc$Q %*% crossprod(cc$Q, yr))
  weighted <- cc$sw * yr
  answer <- numeric(length(cc$norms))
  for (j in unique(cc$feature)) {
    positions <- which(cc$feature == j); oo <- cc$orders[[j]]
    sy <- cumsum(weighted[oo]); say <- cumsum(weighted[oo] * cc$A_local[oo])
    rr <- cc$ranks[positions]
    answer[positions] <- (say[rr] - cc$coefficient[positions] * sy[rr]) / cc$norms[positions]
  }
  answer
}

.dense_selected_vector <- function(cc, w) {
  if (!is.null(cc$U)) return(list(u = cc$U[, w], direction = cc$directions[, w]))
  h <- numeric(length(cc$idx)); h[cc$orders[[cc$feature[w]]][seq_len(cc$ranks[w])]] <- 1
  vv <- cc$sw * (cc$A_local - cc$coefficient[w]) * h - drop(cc$Q %*% cc$qv[, w])
  u <- direction <- numeric(cc$n)
  u[cc$idx] <- cc$sw * vv / cc$norms[w]
  direction[cc$idx] <- vv / (cc$sw * cc$norms[w])
  list(u = u, direction = direction)
}

.dense_candidate_spec <- function(cc, w) {
  ss <- cc$specs[[w]]
  if (is.null(ss$left)) {
    left <- cc$design_X[cc$idx, ss$j] <= ss$cut
    ss$left <- cc$idx[left]; ss$right <- cc$idx[!left]
  }
  ss
}

dense_materialize_candidates <- function(cc) {
  if (is.null(cc) || !is.null(cc$U)) return(cc)
  vectors <- lapply(seq_along(cc$norms), function(w) .dense_selected_vector(cc, w))
  cc$U <- do.call(cbind, lapply(vectors, `[[`, "u"))
  cc$directions <- do.call(cbind, lapply(vectors, `[[`, "direction"))
  cc$specs <- lapply(seq_along(cc$specs), function(w) .dense_candidate_spec(cc, w))
  cc
}

dense_materialize_tree <- function(tree) {
  lapply(tree, function(node) { node$candidates <- dense_materialize_candidates(node$candidates); node })
}

fit_tree <- function(X, A, B, Y, max_depth = 2L, cuts = NULL,
                     min_cell = 8L, cache = NULL,
                     variance_shape = rep(1, nrow(X))) {
  variance_shape <- .gls_shape(variance_shape, nrow(X))
  if (length(Y) != nrow(X) || any(!is.finite(Y))) stop("Invalid outcome vector")
  if (!is.null(cache)) {
    signature <- list(X = X, A = A, B = B, shape = variance_shape,
                      cuts = cuts, min_cell = min_cell)
    signature_key <- ".__dense_design_signature__"
    if (exists(signature_key, envir = cache, inherits = FALSE)) {
      if (!identical(signature, get(signature_key, envir = cache, inherits = FALSE))) {
        stop("Candidate cache reused with a different design or variance shape")
      }
    } else assign(signature_key, signature, envir = cache)
  }
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
    score <- .dense_candidate_scores(cc, Y)
    w <- which.max(score^2)
    selected <- .dense_selected_vector(cc, w)
    selected_spec <- .dense_candidate_spec(cc, w)
    cc$specs[[w]] <- selected_spec
    nm <- if (name == "") "root" else name
    nodes[[nm]] <<- list(name = nm, depth = depth, candidates = cc, winner = w,
      spec = selected_spec, u = selected$u, direction = selected$direction,
      norm = cc$norms[w], path = c(ancestors, nm))
    visit(selected_spec$left, paste0(name, "L"), depth + 1L, c(ancestors, nm))
    visit(selected_spec$right, paste0(name, "R"), depth + 1L, c(ancestors, nm))
    invisible(NULL)
  }
  visit(seq_len(nrow(X)), "", 0L, character())
  nodes
}

# The path geometry is one outer interval intersected with the complement of
# one central interval. This vectorization is algebraically equivalent to the
# original pairwise factor_intervals/intersect_intervals loops.
selection_support <- function(tree, target_name, Y, mode = "path") {
  target <- tree[[target_name]]
  if (is.null(target)) stop("Unknown target")
  direction <- target$direction
  if (is.null(direction)) stop("GLS node has no perturbation direction")
  t <- sum(target$u * Y); z <- Y - direction * t
  keys <- switch(mode, path = target$path, node_only = target_name,
    full_tree = names(tree), stop("Unknown conditioning mode"))
  outer_lo <- -Inf; outer_hi <- Inf
  central_lo <- Inf; central_hi <- -Inf
  slope_tolerance <- 64 * .Machine$double.eps
  for (key in keys) {
    node <- tree[[key]]; cc <- node$candidates; w <- node$winner
    aa <- .dense_candidate_scores(cc, direction); bb <- .dense_candidate_scores(cc, z)
    if (!(key %in% target$path)) {
      # Orthogonality makes every comparison constant on the target line.
      aa[] <- 0; bb <- .dense_candidate_scores(cc, Y)
      if (any(abs(bb) > abs(bb[w]) + 1e-10 * max(1, abs(bb[w])))) {
        stop("Nonancestor winner inconsistent with the observed outcome")
      }
      next
    }
    if (key != target_name) {
      target_idx <- target$candidates$idx
      if (!is.null(node$candidates$design_X)) {
        target_x <- node$candidates$design_X[target_idx, , drop = FALSE]
        minimum <- apply(target_x, 2L, min); maximum <- apply(target_x, 2L, max)
        feature <- node$candidates$feature; cuts <- node$candidates$cuts
        zero <- maximum[feature] <= cuts | minimum[feature] > cuts
      } else zero <- vapply(node$candidates$specs, function(ss) {
        all(target_idx %in% ss$left) || all(target_idx %in% ss$right)
      }, logical(1))
      if (!zero[w]) stop("Ancestor/target regions are inconsistent")
      aa[zero] <- 0
      bb[zero] <- .dense_candidate_scores(cc, Y)[zero]
      aa[abs(aa) < slope_tolerance] <- 0
      competitors <- setdiff(seq_along(aa), w)
      constant <- competitors[aa[competitors] == 0]
      if (length(constant) && any((bb[w] - bb[constant]) * (bb[w] + bb[constant]) < 0)) {
        stop("Constant ancestor comparison failed")
      }
      varying <- competitors[aa[competitors] != 0]
      if (length(varying)) {
        low <- (-abs(bb[w]) - bb[varying]) / aa[varying]
        high <- (abs(bb[w]) - bb[varying]) / aa[varying]
        outer_lo <- max(outer_lo, pmin(low, high))
        outer_hi <- min(outer_hi, pmax(low, high))
      }
    } else {
      aa[w] <- 1; bb[w] <- 0
      competitors <- setdiff(seq_along(aa), w)
      # Extremely close residual directions use the general legacy solver.
      if (any(abs(aa[competitors]) >= 1 - 128 * .Machine$double.eps)) {
        return(.dense_original_selection_support(dense_materialize_tree(tree), target_name, Y, mode))
      }
      if (length(competitors)) {
        r1 <- -bb[competitors] / (1 + aa[competitors])
        r2 <- bb[competitors] / (1 - aa[competitors])
        central_lo <- min(r1, r2); central_hi <- max(r1, r2)
      }
    }
  }
  if (!(outer_lo < outer_hi)) stop("Empty outer selection interval")
  outer <- matrix(c(outer_lo, outer_hi), ncol = 2L)
  if (central_lo < central_hi) {
    rays <- rbind(c(-Inf, central_lo), c(central_hi, Inf))
    support <- intersect_intervals(outer, rays)
  } else support <- outer
  if (!nrow(support)) stop("Empty selection support")
  if (!any(t >= support[, 1] - 1e-9 & t <= support[, 2] + 1e-9)) {
    stop("Observed statistic outside selection support")
  }
  finite <- support[is.finite(support)]
  list(t = t, support = support,
       boundary_distance = if (length(finite)) min(abs(t - finite)) else Inf)
}
