# ACTG 175 all-midpoint analysis. Run from the repository root.
# Numerical construction and inference are unchanged from the manuscript analysis.
project_root <- normalizePath(getwd())
local({ old <- getwd(); setwd(file.path(project_root, "actg"))
  on.exit(setwd(old)); source("gls_pruning.R", local=.GlobalEnv)
  source("stable_tn.R", local=.GlobalEnv) })
out_dir <- file.path(project_root, "results", "actg")
.analysis_cache <- new.env(parent=emptyenv())

# Given the fixed basis and working weights, candidates use X and A only.
# Rank decisions match the prior
# application (tol=1e-9), including nuisance rank deficiency. No merely
# near-parallel residual contrasts are merged.
.midpoints <- function(x) {
  z <- sort(unique(x)); if(length(z)<2L)return(numeric())
  lo<-head(z,-1L);hi<-tail(z,-1L);v<-lo/2+hi/2
  v[v>=hi]<-lo[v>=hi];v[v<lo]<-lo[v<lo];v
}
make_candidates <- function(X,A,B,idx,cuts=NULL,min_cell=20L,
                            variance_shape=rep(1,nrow(X))) {
  n<-nrow(X);variance_shape<-.gls_shape(variance_shape,n)
  if(length(idx)<4L*min_cell)return(NULL)
  sw<-1/sqrt(variance_shape[idx]); aa<-A[idx]
  specs<-vectors<-directions<-aliases<-partitions<-list();norms<-numeric()
  keys<-new.env(hash=TRUE,parent=emptyenv());priority<-0L
  for(vj in seq_len(ncol(X))) {
    cc<-.midpoints(X[idx,vj])
    for(cut in cc) {
      priority<-priority+1L;h<-as.numeric(X[idx,vj]<=cut)
      n1<-sum(aa*h);n0<-sum((1-aa)*h)
      if(min(n1,n0,sum(aa)-n1,length(idx)-sum(aa)-n0)<min_cell)next
      id<-sprintf("%s_le_%.17g_r%d",colnames(X)[vj],cut,sum(h))
      canonical<-if(h[1L]==0)h else 1-h
      key<-paste0(canonical,collapse="")
      if(exists(key,keys,inherits=FALSE)) {
        first<-get(key,keys,inherits=FALSE)
        aliases[[length(aliases)+1L]]<-data.frame(alias=id,
          representative=specs[[first]]$id,complementary=identical(h,1-partitions[[first]]),
          stringsAsFactors=FALSE)
        next
      }
      C<-cbind(B[idx,,drop=FALSE],A=aa,h=h);g<-aa*h
      CW<-C*sw;gw<-g*sw;qc<-qr(CW,tol=1e-9)
      if(qr(cbind(CW,gw),tol=1e-9)$rank!=qc$rank+1L)next
      v<-qr.resid(qc,gw);nv<-sqrt(sum(v^2))
      if(!is.finite(nv)||nv<1e-8)next
      u<-direction<-numeric(n);u[idx]<-sw*v/nv;direction[idx]<-v/(sw*nv)
      k<-length(vectors)+1L;assign(key,k,keys)
      partitions[[k]]<-h;vectors[[k]]<-u;directions[[k]]<-direction;norms[k]<-nv
      specs[[k]]<-list(j=vj,variable=colnames(X)[vj],cut=cut,id=id,priority=priority,
        left=idx[h==1],right=idx[h==0],nuisance_rank=qc$rank)
    }
  }
  if(!length(vectors))return(NULL)
  list(U=do.call(cbind,vectors),directions=do.call(cbind,directions),norms=norms,
       specs=specs,idx=idx,aliases=if(length(aliases))do.call(rbind,aliases)else data.frame())
}

make_settings <- function() {
  data.frame(sd0=1,sd1=1,setting=c("main","lambda1","lambda2","lambda3","lambda4"),
             depth=3L,lambda=0:4,drop_zero=FALSE)
}

read_application_data <- function(drop_zero = FALSE) {
  d <- read.table(file.path(project_root, "data", "ACTG175.txt"), header=TRUE)
  p <- d[d$arms %in% c(1,2), ]; p <- p[order(p$pidnum), ]
  stopifnot(nrow(p)==1046L, sum(p$arms==1)==522L, sum(p$arms==2)==524L)
  # This seed and sorted-ID sampling reproduce the manuscript's fixed523/523 split.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261028)
  training_ids <- p$pidnum[sort(sample.int(nrow(p), nrow(p)%/%2L))]
  holdout_ids <- p$pidnum[!p$pidnum %in% training_ids]
  if (drop_zero) p <- p[p$cd40 != 0, ]
  B <- cbind(intercept=1, as.matrix(p[c("cd40", "cd80", "age", "wtkg",
    "karnof", "gender", "symptom", "preanti")]),
    strat2=as.integer(p$strat==2), strat3=as.integer(p$strat==3))
  for (j in 2:ncol(B)) B[,j] <- as.numeric(scale(B[,j]))
  list(ids=p$pidnum, A=as.integer(p$arms==1), B=B,
    X=as.matrix(p[c("cd40", "cd80", "age", "preanti", "str2", "symptom")]),
    Y=p$cd420, training=which(p$pidnum %in% training_ids),
    holdout=which(p$pidnum %in% holdout_ids))
}

pilot_variances <- function(dat, index = seq_along(dat$Y)) {
  do.call(rbind, lapply(0:1, function(a) {
    ix <- index[dat$A[index] == a]; q <- qr(dat$B[ix, , drop = FALSE], tol = 1e-9)
    rss <- sum(qr.resid(q, dat$Y[ix])^2); df <- length(ix) - q$rank
    stopifnot(df > 0, rss > 0)
    data.frame(arm = a, n = length(ix), rank_B = q$rank, df = df,
               rss = rss, variance = rss/df, sd = sqrt(rss/df))
  }))
}

subset_data <- function(dat, ix) {
  list(ids = dat$ids[ix], X = dat$X[ix, , drop = FALSE], A = dat$A[ix],
       B = dat$B[ix, , drop = FALSE], Y = dat$Y[ix])
}

holdout_contrast <- function(tree, name, hold, variance_shape) {
  ix <- seq_along(hold$Y); target <- tree[[name]]
  for (anc in head(target$path, -1L)) {
    ss <- tree[[anc]]$spec
    go_left <- substr(name, tree[[anc]]$depth + 1L, tree[[anc]]$depth + 1L) == "L"
    inside <- hold$X[ix, ss$j] <= ss$cut
    ix <- ix[if (go_left) inside else !inside]
  }
  h <- as.numeric(hold$X[ix, target$spec$j] <= target$spec$cut)
  if (length(ix) < 15L) return(list(node = NULL, reason = "parent_too_small", idx = ix))
  if (min(table(factor(hold$A[ix], 0:1), factor(h, 0:1))) < 1L)
    return(list(node = NULL, reason = "empty_treatment_child", idx = ix))
  C <- cbind(hold$B[ix, , drop = FALSE], A = hold$A[ix], h = h)
  sw <- 1/sqrt(variance_shape[ix]); gw <- hold$A[ix] * h * sw
  qc <- qr(C * sw, tol = 1e-9)
  if (qr(cbind(C * sw, gw), tol = 1e-9)$rank != qc$rank + 1L)
    return(list(node = NULL, reason = "unidentifiable", idx = ix))
  v <- qr.resid(qc, gw); nv <- sqrt(sum(v^2)); u <- direction <- numeric(length(hold$Y))
  u[ix] <- sw*v/nv; direction[ix] <- v/(sw*nv)
  nd <- list(u = u, direction = direction, norm = nv, spec = target$spec,
    candidates = list(idx = ix), path = target$path, depth = target$depth)
  nd$spec$left <- ix[h == 1]; nd$spec$right <- ix[h == 0]; nd$spec$nuisance_rank <- qc$rank
  list(node = nd, reason = "ok", idx = ix)
}

contrast_metadata <- function(node, dat, variance_shape) {
  ix <- node$candidates$idx; ss <- node$spec
  eta <- node$u/node$norm; w <- eta*dat$A
  left <- seq_along(dat$Y) %in% ss$left; right <- seq_along(dat$Y) %in% ss$right
  data.frame(parent_n = length(ix),
    n0left = sum(dat$A[ss$left] == 0), n1left = sum(dat$A[ss$left] == 1),
    n0right = sum(dat$A[ss$right] == 0), n1right = sum(dat$A[ss$right] == 1),
    candidate_count = if (is.null(node$candidates$U)) NA_integer_ else ncol(node$candidates$U),
    nuisance_rank = ss$nuisance_rank, contrast_norm = node$norm,
    standard_error = 1/node$norm, eta_l1 = sum(abs(eta)), eta_max_abs = max(abs(eta)),
    causal_weight_l1 = sum(abs(w)), causal_left_negative = sum(w[left] < -1e-10),
    causal_right_positive = sum(w[right] > 1e-10),
    normalization_error = abs(sum(node$u^2 * variance_shape) - 1),
    stringsAsFactors = FALSE)
}

evaluate_contrast <- function(node, name, dat, support, method, setting,
                              variance_shape, delta = NA_real_, boundary = Inf) {
  t <- sum(node$u * dat$Y); est <- t/node$norm
  p <- 2 * plogis(-abs(tn_logodds(t, 0, 1, support)))
  ci <- if (method %in% c("naive", "honest")) t + c(-1, 1)*qnorm(.975) else
    tn_interval(t, 1, support)
  endpoint_error <- max(abs(c(tn_logodds(t, ci[1], 1, support) - qlogis(.975),
                              tn_logodds(t, ci[2], 1, support) - qlogis(.025))))
  ci <- ci/node$norm
  cbind(data.frame(setting = setting, method = method, node = name,
    depth = node$depth, variable = node$spec$variable, threshold = node$spec$cut,
    split = node$spec$id, available = TRUE, estimate = est,
    p_value = p, lower = ci[1], upper = ci[2], width = diff(ci),
    components = nrow(support), delta = delta, score = t,
    boundary_distance = boundary, endpoint_logodds_error = endpoint_error,
    inference_status = "ok", inference_error = "", stringsAsFactors = FALSE),
    contrast_metadata(node, dat, variance_shape))
}

empty_results <- function() {
  data.frame(setting=character(),method=character(),node=character(),depth=integer(),
    variable=character(),threshold=numeric(),split=character(),available=logical(),
    estimate=numeric(),p_value=numeric(),lower=numeric(),upper=numeric(),width=numeric(),
    components=integer(),delta=numeric(),score=numeric(),boundary_distance=numeric(),
    endpoint_logodds_error=numeric(),inference_status=character(),inference_error=character(),
    parent_n=integer(),n0left=integer(),n1left=integer(),n0right=integer(),n1right=integer(),
    candidate_count=integer(),nuisance_rank=integer(),contrast_norm=numeric(),
    standard_error=numeric(),eta_l1=numeric(),eta_max_abs=numeric(),causal_weight_l1=numeric(),
    causal_left_negative=integer(),causal_right_positive=integer(),normalization_error=numeric(),
    stringsAsFactors=FALSE)
}

unavailable_result <- function(node,name,setting,reason,status="unavailable") {
  z <- empty_results()[1,,drop=FALSE]
  z$setting<-setting;z$method<-"honest";z$node<-name;z$depth<-node$depth
  z$variable<-node$spec$variable;z$threshold<-node$spec$cut;z$split<-node$spec$id
  z$available<-FALSE;z$inference_status<-status;z$inference_error<-reason
  z
}

fit_setting <- function(sc) {
  dat <- read_application_data(sc$drop_zero); full_pilot <- pilot_variances(dat)
  train_pilot <- pilot_variances(dat, dat$training)
  multipliers <- c(sc$sd0, sc$sd1)
  full_variances <- full_pilot$variance * multipliers^2
  honest_variances <- train_pilot$variance * multipliers^2
  vfull <- full_variances[dat$A + 1L]
  train <- subset_data(dat, dat$training); hold <- subset_data(dat, dat$holdout)
  vtrain <- honest_variances[train$A + 1L]; vhold <- honest_variances[hold$A + 1L]
  tree_key<-paste(sc$depth,sc$sd0,sc$sd1,sc$drop_zero,sep="|")
  if(exists(tree_key,.analysis_cache,inherits=FALSE)) {
    cached<-get(tree_key,.analysis_cache,inherits=FALSE)
    full_tree<-cached$full;training_tree<-cached$training
  } else {
    full_tree <- fit_tree(dat$X, dat$A, dat$B, dat$Y, max_depth=sc$depth,
                          min_cell=20L,variance_shape=vfull)
    training_tree <- fit_tree(train$X,train$A,train$B,train$Y,max_depth=sc$depth,
                             min_cell=20L,variance_shape=vtrain)
    assign(tree_key,list(full=full_tree,training=training_tree),.analysis_cache)
  }
  full_prune <- prune_tree(full_tree, dat$Y, sc$lambda)
  training_prune <- prune_tree(training_tree, train$Y, sc$lambda)
  supports <- holdout <- list(); rows <- list(); k <- 0L
  add <- function(z) { k <<- k + 1L; rows[[k]] <<- z }
  for (nm in full_prune$kept) {
    node <- full_tree[[nm]]
    growth <- selection_support(full_tree, nm, dat$Y, "path")
    selected <- selection_support_pruned(full_tree, nm, dat$Y, sc$lambda,
      "path", full_prune, growth)
    all_growth <- selection_support(full_tree, nm, dat$Y, "full_tree")
    all_selected <- selection_support_pruned(full_tree, nm, dat$Y, sc$lambda,
      "full_tree", full_prune)
    supports[[nm]] <- list(growth = growth, proposed = selected,
                           full_growth = all_growth, full_proposed = all_selected)
    add(evaluate_contrast(node, nm, dat, whole_line(), "naive", sc$setting, vfull))
    add(evaluate_contrast(node, nm, dat, growth$support, "growth_only", sc$setting,
      vfull, boundary = growth$boundary_distance))
    add(evaluate_contrast(node, nm, dat, selected$support, "proposed", sc$setting,
      vfull, selected$delta, selected$boundary_distance))
  }
  for (nm in training_prune$kept) {
    hc <- holdout_contrast(training_tree, nm, hold, vhold); holdout[[nm]] <- hc
    if (!is.null(hc$node)) add(evaluate_contrast(hc$node, nm, hold, whole_line(),
      "honest", sc$setting, vhold))
    else add(unavailable_result(training_tree[[nm]],nm,sc$setting,hc$reason))
  }
  results <- if (length(rows)) do.call(rbind, rows) else empty_results()
  results$lambda <- rep(sc$lambda,nrow(results))
  pilots <- rbind(transform(full_pilot, pilot = "full"),
                  transform(train_pilot, pilot = "Honest_training"))
  pilots$setting <- sc$setting
  pilots$sd_multiplier <- multipliers[pilots$arm + 1L]
  pilots$working_variance <- pilots$variance * pilots$sd_multiplier^2
  pilots$working_sd <- sqrt(pilots$working_variance)
  snapshot <- list(setting = sc, data = dat, train = train, hold = hold,
    full_pilot = full_pilot, training_pilot = train_pilot,
    full_variance = vfull, training_variance = vtrain, holdout_variance = vhold,
    full_tree = full_tree, training_tree = training_tree,
    full_prune = full_prune, training_prune = training_prune,
    supports = supports, holdout = holdout, results = results, pilots = pilots)
  snapshot
}
