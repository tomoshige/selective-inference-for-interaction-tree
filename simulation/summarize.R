# Summarize generated full-midpoint simulations. Independent trial bundles are the
# Monte Carlo units; nodes and reused rho/depth/penalty versions are not.
suppressPackageStartupMessages(library(data.table))
script <- grep("^--file=", commandArgs(), value=TRUE)
if (length(script)) setwd(dirname(normalizePath(sub("^--file=", "", script[1L]))))
root <- normalizePath("..")
args <- commandArgs(trailingOnly=TRUE)
opt <- list(stage="main", input=NULL, output=NULL, allow_partial=FALSE)
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--allow-partial") opt$allow_partial <- TRUE else {
    if (!a %in% c("--stage", "--input", "--output") || i == length(args))
      stop("Usage: summarize.R [--stage main|quick] [--input DIR] [--output DIR] [--allow-partial]")
    i <- i+1L; opt[[sub("^--", "", a)]] <- args[i]
  }
  i <- i+1L
}
stopifnot(opt$stage %in% c("main","quick"))
if (is.null(opt$input)) opt$input <- file.path(root,"results","simulation",opt$stage,"chunks")
if (is.null(opt$output)) opt$output <- file.path(root,"results","simulation",opt$stage,"summaries")
dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
files <- sort(list.files(opt$input, pattern="^scenario_[0-9]+_chunk_[0-9]+\\.rds$", full.names=TRUE))
if (!length(files)) stop("No result RDS files: ", opt$input)

wilson <- function(k, n) {
  if (!n) return(c(lo=NA_real_, hi=NA_real_))
  p <- k/n; z <- qnorm(.975); den <- 1+z*z/n
  mid <- (p+z*z/(2*n))/den
  half <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
  c(lo=max(0, mid-half), hi=min(1, mid+half))
}
binary <- function(x) {
  if (anyNA(x)) stop("Missing binary outcome in declared denominator")
  if (!length(x)) return(c(rate=NA_real_,mcse=NA_real_,lo=NA_real_,hi=NA_real_))
  p <- mean(x); c(rate=p,mcse=sqrt(p*(1-p)/length(x)),wilson(sum(x),length(x)))
}
cluster_ratio <- function(k, den) {
  N <- length(k); total <- sum(den)
  if (!total) return(c(rate=NA_real_,mcse=NA_real_,lo=NA_real_,hi=NA_real_))
  p <- sum(k)/total
  se <- if (N>1L) sqrt(N/(N-1)*sum((k-p*den)^2))/total else NA_real_
  ci <- if (p>0 && p<1 && sum(den>0)>=30L && is.finite(se)) p+c(-1,1)*qnorm(.975)*se else c(NA_real_,NA_real_)
  c(rate=p,mcse=se,lo=max(0,ci[1]),hi=min(1,ci[2]))
}
fields <- function(prefix, x) setNames(as.list(unname(x)),
  ifelse(names(x)=="rate",prefix,paste0(prefix,"_",names(x))))
q7 <- function(x, p) if (length(x)) unname(quantile(x,p,type=7L,names=FALSE)) else NA_real_

width_summary <- function(z) {
  ciok <- z$ci_status=="ok"
  w <- z$width[ciok]
  if (anyNA(w) || any(w<0)) stop("Invalid interval width marked successful")
  ep <- if ("endpoint_covered" %in% names(z)) z$endpoint_covered[ciok] else rep(NA,length(w))
  out <- list(ci_requested_total=sum(z$ci_planned),
    ci_attempted=sum(z$ci_status %in% c("ok","failed")),ci_n=sum(ciok),
    ci_failures=sum(z$ci_status=="failed"),width_small_cell=sum(ciok)<30L,
    width_median=q7(w,.5),width_p90=q7(w,.9),width_p99=q7(w,.99),
    width_max=if(length(w))max(w) else NA_real_,width_infinite=sum(is.infinite(w)),
    endpoint_coverage=if(length(ep)&&all(!is.na(ep)))mean(ep) else NA_real_,
    endpoint_pivot_disagreements=if(length(ep)&&all(!is.na(ep)))sum(ep!=z$covered[ciok]) else NA_integer_)
  for (cut in c(10,100,1000)) {
    out[[paste0("width_gt_",cut,"_n")]] <- sum(w>cut)
    out[[paste0("width_gt_",cut,"_rate")]] <- if(length(w))mean(w>cut) else NA_real_
  }
  out
}

count_summary <- function(z) {
  ok <- z$inference_status=="ok"
  list(attempted=nrow(z),grown=sum(z$grown),retained=sum(z$retained),
    eligible=sum(z$eligible),available=sum(ok),unavailable=sum(!z$eligible),
    not_grown=sum(!z$grown),pruned=sum(z$grown&!z$retained),
    holdout_ineligible=sum(z$retained&!z$eligible),
    failures=sum(z$inference_status=="failed"),
    reached=if("reached"%in%names(z))sum(z$reached)else NA_integer_,
    reached_without_candidate=if(all(c("reached","candidate_count")%in%names(z)))
      sum(z$reached&!is.na(z$candidate_count)&z$candidate_count==0)else NA_integer_)
}
fixed_node <- function(z) {
  if (anyDuplicated(z$replication)) stop("Duplicate fixed-address trial")
  ok <- z$inference_status=="ok"
  c(list(trials=nrow(z)),count_summary(z),
    fields("grown_rate",binary(z$grown)),
    fields("retained_per_possible",binary(z$retained)),
    fields("retained_given_grown",binary(z$retained[z$grown])),
    fields("availability_rate",binary(z$eligible)),
    fields("eligible_given_retained",binary(z$eligible[z$retained])),
    fields("reject_rate",binary(z$p_value[ok]<.05)),
    fields("coverage",binary(z$covered[ok])),
    fields("discovery_per_trial",binary(ifelse(ok,z$p_value<.05,FALSE))),
    list(small_cell=sum(ok)<100L,mc_interval_type="Wilson"),width_summary(z))
}
cluster_group <- function(z, tree=FALSE) {
  ok <- z$inference_status=="ok"
  M <- cbind(possible=1,grown=as.integer(z$grown),retained=as.integer(z$retained),
    eligible=as.integer(z$eligible),available=as.integer(ok),
    reject=as.integer(ifelse(ok,z$p_value<.05,FALSE)),
    covered=as.integer(ifelse(ok,z$covered,FALSE)),
    failures=as.integer(z$inference_status=="failed"))
  ct <- rowsum(M,z$replication,reorder=TRUE); N <- nrow(ct)
  out <- c(list(trials=N,possible_per_trial=nrow(z)/N),count_summary(z))
  ratios <- list(grown_rate=c("grown","possible"),
    retained_per_possible=c("retained","possible"),
    retained_given_grown=c("retained","grown"),
    availability_rate=c("eligible","possible"),
    eligible_given_retained=c("eligible","retained"),
    reject_rate=c("reject","available"),coverage=c("covered","available"),
    discovery_per_possible=c("reject","possible"))
  for (nm in names(ratios)) {
    key <- ratios[[nm]]; out <- c(out,fields(nm,cluster_ratio(ct[,key[1]],ct[,key[2]])))
  }
  out <- c(out,list(positive_available_trials=sum(ct[,"available"]>0),
    small_cell=sum(ct[,"available"])<100L||sum(ct[,"available"]>0)<30L,
    discovery_per_trial=mean(ct[,"reject"]),
    discovery_per_trial_mcse=if(N>1L)sd(ct[,"reject"])/sqrt(N) else NA_real_,
    mc_interval_type="trial-cluster normal; omitted at boundaries or <30 positive clusters"))
  if (tree) {
    for (nm in c("grown","retained","eligible","available")) {
      out[[paste0("mean_",nm,"_nodes")]] <- mean(ct[,nm])
      out[[paste0("mean_",nm,"_nodes_mcse")]] <- if(N>1L)sd(ct[,nm])/sqrt(N) else NA_real_
    }
    out$mean_retained_leaves <- mean(ct[,"retained"])+1
    for (nm in c("retained","available","rejection")) {
      col <- if(nm=="rejection")"reject" else nm
      out <- c(out,fields(paste0("any_",nm,"_rate"),binary(ct[,col]>0)))
    }
    K <- ct[as.character(z$replication),"retained"]
    bonf <- ifelse(ok,z$p_value<.05/pmax(K,1),FALSE)
    zero <- ifelse(ok,abs(z$truth)<1e-10,FALSE)
    false_reject <- ifelse(ok,z$p_value<.05,FALSE)&zero
    for (nm in c("any_bonferroni_rejection_rate","any_false_rejection_rate","any_bonferroni_false_rejection_rate")) {
      b <- switch(nm,any_bonferroni_rejection_rate=bonf,
        any_false_rejection_rate=false_reject,any_bonferroni_false_rejection_rate=bonf&zero)
      out <- c(out,fields(nm,binary(rowsum(as.integer(b),z$replication,reorder=TRUE)[,1]>0)))
    }
    out$trials_with_numerical_failure <- sum(ct[,"failures"]>0)
    out$zero_target_n <- sum(zero); out$nonzero_target_n <- sum(ok&!zero)
    out$zero_target_rejections <- sum(false_reject)
    out$nonzero_target_rejections <- sum(ifelse(ok,z$p_value<.05,FALSE)&!zero)
    out$zero_target_reject_rate <- if(sum(zero))sum(false_reject)/sum(zero) else NA_real_
    out$nonzero_target_reject_rate <- if(sum(ok&!zero))out$nonzero_target_rejections/sum(ok&!zero) else NA_real_
    out$max_retained_node_depth_median <- median(vapply(split(seq_len(nrow(z)),z$replication),function(ii) {
      if(any(z$retained[ii]))max(z$node_depth[ii][z$retained[ii]]) else -1
    },numeric(1)))
  }
  c(out,width_summary(z))
}

normalize_records <- function(z, file) {
  if (!is.data.frame(z)) stop("Expected records data.frame: ",file)
  z <- as.data.table(z)
  aliases <- c(rho="variance_ratio",effect="shape",D="max_depth",rep="replication",status="inference_status")
  for (old in names(aliases)) if(old%in%names(z)&&!aliases[[old]]%in%names(z))setnames(z,old,aliases[[old]])
  required <- c("n","shape","variance_ratio","replication","max_depth","lambda","method","node",
    "grown","retained","eligible","inference_status","p_value","truth","covered","ci_status","width")
  if(!all(required%in%names(z)))stop("Missing input fields in ",file,": ",paste(setdiff(required,names(z)),collapse=","))
  defaults <- list(variance_mode="known",node_depth=ifelse(z$node=="root",0L,nchar(z$node)),
    ci_planned=z$ci_status!="not_requested",endpoint_covered=NA,seed=NA_integer_)
  for(nm in names(defaults))if(!nm%in%names(z))z[,(nm):=defaults[[nm]]]
  if(!"base_group"%in%names(z))z[,base_group:=paste(n,shape,sep="_")]
  if(!"scenario"%in%names(z))z[,scenario:=paste(n,shape,variance_ratio,sep="_")]
  if(!"gamma"%in%names(z))z[,gamma:=unname(c(null=0,g05=.5,g10=1,g20=2,weak_a=1)[shape])]
  if(!"a"%in%names(z))z[,a:=ifelse(shape=="weak_a",.5,1)]
  if(anyNA(z[,.(scenario,n,shape,variance_ratio,replication,max_depth,lambda,method,node,
      grown,retained,eligible,inference_status,ci_status,variance_mode,node_depth)]))stop("Missing identifiers/status: ",file)
  if(any(z$retained&!z$grown)||any(z$eligible&!z$retained))stop("Invalid eligibility hierarchy: ",file)
  if(any(!z$inference_status%in%c("ok","unavailable","failed")))stop("Unknown inference status: ",file)
  if(any(!z$ci_status%in%c("ok","not_requested","failed")))stop("Unknown CI status: ",file)
  ok <- z$inference_status=="ok"
  if(any((ok|z$inference_status=="failed")!=z$eligible))stop("Eligibility/status mismatch: ",file)
  if(anyNA(z$p_value[ok])||any(!is.finite(z$p_value[ok]))||any(z$p_value[ok]<0|z$p_value[ok]>1)||
     anyNA(z$covered[ok])||any(!is.finite(z$truth[ok])))stop("Invalid successful inference record: ",file)
  z
}
id_cols <- c("scenario","base_group","n","shape","gamma","a","variance_ratio","max_depth","lambda","variance_mode","method")
raw_key <- c("scenario","replication","max_depth","lambda","variance_mode","method","node")
extract_object <- function(obj,file) {
  if(is.data.frame(obj))return(list(records=normalize_records(obj,file)))
  if(is.list(obj)) {
    key <- intersect(c("records","rows","data","results"),names(obj))
    if(length(key)&&is.data.frame(obj[[key[1]]])) {
      obj$records <- normalize_records(obj[[key[1]]],file);return(obj)
    }
  }
  stop("Unsupported RDS object in ",file)
}
paired_summary <- function(z, left="proposed",right="growth_only",by_mode=FALSE) {
  keep <- c("replication","node","inference_status","truth","covered","p_value","ci_status","width")
  x <- z[z$method==left,..keep];y <- z[z$method==right,..keep]
  if(by_mode) {
    x <- z[z$variance_mode==left,..keep]; y <- z[z$variance_mode==right,..keep]
  }
  if(!nrow(x)||!nrow(y))return(NULL)
  m <- merge(x,y,by=c("replication","node"),suffixes=c("_left","_right"),all=TRUE)
  if(anyNA(m$inference_status_left)||anyNA(m$inference_status_right))stop("Incomplete same-target pairs")
  ok <- m$inference_status_left=="ok"&m$inference_status_right=="ok"
  if(any(abs(m$truth_left[ok]-m$truth_right[ok])>1e-8*(1+abs(m$truth_left[ok]))))stop("Paired targets differ")
  M <- cbind(available=as.integer(ok),covered=ifelse(ok,as.integer(m$covered_left)-as.integer(m$covered_right),0L),
    reject=ifelse(ok,as.integer(m$p_value_left<.05)-as.integer(m$p_value_right<.05),0L))
  ct <- rowsum(M,m$replication,reorder=TRUE); N <- nrow(ct);den<-sum(ct[,"available"])
  out <- list(left=left,right=right,trials=N,n_pairs=den,positive_available_trials=sum(ct[,"available"]>0))
  for(key in c("covered","reject")) {
    p <- if(den)sum(ct[,key])/den else NA_real_
    se <- if(den&&N>1L)sqrt(N/(N-1)*sum((ct[,key]-p*ct[,"available"])^2))/den else NA_real_
    ci <- if(sum(ct[,"available"]>0)>=30L&&is.finite(se))p+c(-1,1)*qnorm(.975)*se else c(NA,NA)
    nm <- if(key=="covered")"coverage_difference" else "reject_difference"
    out[[nm]]<-p;out[[paste0(nm,"_mcse")]]<-se
    out[[paste0(nm,"_lo")]]<-max(-1,ci[1]);out[[paste0(nm,"_hi")]]<-min(1,ci[2])
  }
  iw <- ok&m$ci_status_left=="ok"&m$ci_status_right=="ok"
  out$ci_pairs<-sum(iw)
  wr<-m$width_left[iw]/m$width_right[iw]
  out$paired_width_ratio_undefined<-sum(is.na(wr))
  out$median_paired_width_ratio<-q7(wr[!is.na(wr)],.5)
  out
}

# Filenames should keep scenario chunks together. Read one scenario into memory,
# never the full production result grid. A single file may hold many scenarios.
file_groups <- split(files,sub("([_-](chunk|part|batch)[_-]?[0-9]+)?\\.rds$","",basename(files),ignore.case=TRUE))
tables <- list(node=list(),depth=list(),tree=list(),paired=list(),sigma_paired=list(),qq=list(),recovery=list(),timing=list(),sigma=list())
inventory <- list(); violations <- character(); seen_keys <- new.env(hash=TRUE,parent=emptyenv())
for(group in file_groups) {
  objs <- lapply(group,function(f)extract_object(readRDS(f),f))
  d <- rbindlist(lapply(objs,`[[`,"records"),fill=TRUE)
  if(anyDuplicated(d[,..raw_key]))stop("Duplicate result key in input group: ",paste(basename(group),collapse=","))
  # Configurations from different files must not silently overlap.
  ck <- unique(d[,.(scenario,replication)])
  signatures <- do.call(paste,c(as.list(ck),sep="|"))
  if(any(vapply(signatures,exists,logical(1),envir=seen_keys,inherits=FALSE)))
    stop("Scenario/trial/configuration repeated across separate file groups; use *_chunk_*.rds naming")
  for(k in signatures)assign(k,TRUE,envir=seen_keys)
  counts <- d[,.N,by=.(scenario,replication,max_depth,lambda,variance_mode,method)]
  if(any(counts$N!=2^counts$max_depth-1L))stop("Incomplete potential-node rows in ",paste(basename(group),collapse=","))
  if(any(d$node_depth>=d$max_depth))stop("Node deeper than growth limit")
  node <- d[,fixed_node(.SD),by=c(id_cols,"node","node_depth")]
  depth <- d[,cluster_group(.SD),by=c(id_cols,"node_depth")]
  tree <- d[,cluster_group(.SD,TRUE),by=id_cols]
  tables$node[[length(tables$node)+1L]]<-node
  tables$depth[[length(tables$depth)+1L]]<-depth
  tables$tree[[length(tables$tree)+1L]]<-tree
  pcols <- setdiff(id_cols,"method")
  for(comp in c("naive","growth_only")) {
    pp <- d[,paired_summary(.SD,"proposed",comp),by=pcols]
    if(nrow(pp))tables$paired[[length(tables$paired)+1L]]<-pp
  }
  scols <- setdiff(id_cols,"variance_mode")
  if(all(c("known","plugin")%in%d$variance_mode)) {
    pp <- d[,paired_summary(.SD,"plugin","known",TRUE),by=scols]
    if(nrow(pp))tables$sigma_paired[[length(tables$sigma_paired)+1L]]<-pp
  }
  qq <- d[inference_status=="ok", {
    pp <- c(.001,.005,.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99,.995,.999)
    list(probability=pp,p_value_quantile=as.numeric(quantile(p_value,pp,type=7L)),p_values=.N)
  },by=c(id_cols,"node_depth")]
  tables$qq[[length(tables$qq)+1L]]<-qq
  if("sigma_hat"%in%names(d)) {
    sig <- unique(d[,intersect(c("scenario","base_group","n","shape","gamma","a","variance_ratio","variance_mode","method","replication","sigma_hat","sigma_df"),names(d)),with=FALSE])
    sg <- setdiff(names(sig),c("replication","sigma_hat","sigma_df"))
    sig <- sig[is.finite(sigma_hat),.(trials=.N,sigma_mean=mean(sigma_hat),sigma_median=median(sigma_hat),
      sigma_p10=q7(sigma_hat,.1),sigma_p90=q7(sigma_hat,.9),sigma2_ratio_mean=mean(sigma_hat^2)),by=sg]
    tables$sigma[[length(tables$sigma)+1L]]<-sig
  }
  for(obj in objs) {
    if(is.data.frame(obj$recovery)&&nrow(obj$recovery))tables$recovery[[length(tables$recovery)+1L]]<-as.data.table(obj$recovery)
    if(is.data.frame(obj$trials)&&nrow(obj$trials))tables$timing[[length(tables$timing)+1L]]<-as.data.table(obj$trials)
  }
  if(!any(vapply(objs,function(x)is.data.frame(x$trials),logical(1)))&&"trial_seconds"%in%names(d)) {
    tm <- unique(d[,.(scenario,base_group,n,shape,gamma,a,variance_ratio,replication,trial_seconds)])
    if(anyDuplicated(tm[,.(scenario,replication)]))violations<-c(violations,"Multiple trial_seconds values per outcome version; timing CSV retained without aggregation.")
    tables$timing[[length(tables$timing)+1L]]<-tm
  }
  inventory[[length(inventory)+1L]]<-data.table(files=paste(basename(group),collapse=";"),records=nrow(d),
    outcome_versions=nrow(unique(d[,.(scenario,replication)])),
    independent_bundles=nrow(unique(d[,.(base_group,replication)])))
  cat(sprintf("Summarized %s: %s rows\n",paste(basename(group),collapse=","),format(nrow(d),big.mark=",")))
  rm(d,objs,node,depth,tree);gc(FALSE)
}
combined <- lapply(tables,function(x)if(length(x))rbindlist(x,fill=TRUE)else data.table())
for(nm in c("node","depth","tree","paired","sigma_paired","qq","sigma")) {
  z <- combined[[nm]]
  if(nrow(z)) {
    keys<-intersect(c(id_cols,"node_depth","node","left","right","probability"),names(z))
    setorderv(z,keys); fwrite(z,file.path(opt$output,paste0(nm,"_summary.csv")))
  }
}
if(nrow(combined$timing)) {
  tm <- unique(combined$timing)
  fwrite(tm,file.path(opt$output,"trial_timing.csv"))
  if("trial_seconds"%in%names(tm)) {
    cols<-intersect(c("n","shape","gamma","a","variance_ratio"),names(tm))
    ts<-tm[,.(trials=.N,seconds_mean=mean(trial_seconds),seconds_median=median(trial_seconds),
      seconds_p90=q7(trial_seconds,.9),seconds_max=max(trial_seconds)),by=cols]
    fwrite(ts,file.path(opt$output,"timing_summary.csv"))
    measurement_cols<-intersect(c("trial_seconds","grow_seconds","infer_seconds","full_root_candidates",
      "honest_root_candidates","full_tree_bytes","training_tree_bytes","full_sigma_hat","honest_sigma_hat"),names(tm))
    for(mc in measurement_cols)set(tm,j=mc,value=as.numeric(tm[[mc]]))
    meas<-melt(tm,id.vars=c(cols,"replication"),measure.vars=measurement_cols,
      variable.name="measurement",value.name="value")
    ms<-meas[!is.na(value),.(trials=.N,mean=mean(value),median=median(value),p90=q7(value,.9),max=max(value)),
      by=c(cols,"measurement")]
    fwrite(ms,file.path(opt$output,"computation_summary.csv"))
    sigcols<-intersect(c("full_sigma_hat","honest_sigma_hat"),names(tm))
    if(length(sigcols)) {
      sig<-melt(tm,id.vars=c(cols,"replication"),measure.vars=sigcols,
        variable.name="estimate_group",value.name="sigma_hat")
      ss<-sig[is.finite(sigma_hat),.(trials=.N,sigma_mean=mean(sigma_hat),sigma_median=median(sigma_hat),
        sigma_p10=q7(sigma_hat,.1),sigma_p90=q7(sigma_hat,.9),sigma2_ratio_mean=mean(sigma_hat^2)),
        by=c(cols,"estimate_group")]
      fwrite(ss,file.path(opt$output,"sigma_estimation_summary.csv"))
    }
  }
}
if(nrow(combined$recovery)) {
  rr <- combined$recovery
  fwrite(rr,file.path(opt$output,"recovery_trials.csv"))
  if(all(c("recovered","rejected")%in%names(rr))) {
    rc <- intersect(c("scenario","n","shape","gamma","a","variance_ratio","max_depth","lambda","variance_mode","method","true_split","true_node","true_depth"),names(rr))
    rs <- rr[,c(list(trials=.N),fields("recovered_rate",binary(recovered)),
      fields("recovered_and_rejected_rate",binary(recovered&rejected))),by=rc]
    fwrite(rs,file.path(opt$output,"recovery_summary.csv"))
  }
}
fwrite(rbindlist(inventory),file.path(opt$output,"input_inventory.csv"))
tree <- combined$tree
configs <- unique(tree[,intersect(c("n","shape","gamma","a","variance_ratio","max_depth","lambda","variance_mode","method"),names(tree)),with=FALSE])
fwrite(configs,file.path(opt$output,"observed_configurations.csv"))
manifest <- c(file.path(opt$input,"scenarios.csv"),"results/scenarios.csv","scenarios.csv")
manifest <- manifest[file.exists(manifest)]
complete <- FALSE
if(length(manifest)) {
  sp <- fread(manifest[1])
  if(all(c("id","repetitions")%in%names(sp))&&"scenario"%in%names(tree)) {
    expected <- sp[,.(scenario=as.character(id),expected_trials=repetitions)]
    observed <- unique(tree[,.(scenario=as.character(scenario),trials)])
    check <- merge(expected,observed,by="scenario",all=TRUE)
    complete <- !anyNA(check)&&all(check$expected_trials==check$trials)
    fwrite(check,file.path(opt$output,"completeness_check.csv"))
  }
}
if(!complete)violations<-c(violations,"Production completeness was not established from scenarios.csv. These summaries describe observed input only.")
if(!complete&&!opt$allow_partial)warning(tail(violations,1))
null <- tree[shape=="null"&variance_mode=="known"]
fmt_range <- function(x) {x<-x[is.finite(x)];if(length(x))sprintf("%.4f to %.4f",min(x),max(x))else"unavailable"}
report <- c("FULL-MIDPOINT INTERACTION-TREE SIMULATION ANALYSIS",
  paste("Stage:",opt$stage),paste("Production completeness established:",complete),
  paste("Input files:",length(files)),paste("Observed configuration rows:",nrow(configs)),"",
  "Independent trial bundles are the Monte Carlo sampling units. Reused variance ratios, depths, penalties, and nodes are correlated.",
  "Fixed-address intervals are Wilson; depth/tree summaries use trial-cluster ratio MCSEs and normal intervals (omitted at boundaries or <30 positive clusters).",
  "Rates condition on successfully computed eligible targets unless explicitly per-trial/per-possible. Structural unavailability and numerical failures remain separate.",
  "All potential-node rows and empty-tree trials are retained. Any-rejection rates use all trials. Failed computations do not count as observed rejections; inspect failure counts.",
  "Bonferroni denominators count every retained split, including unavailable honest targets. Plugin estimates have no finite-sample exactness claim.",
  "Widths are descriptive, use every successfully computed prespecified endpoint result, include infinite widths, and are never capped or trimmed.",
  "Naive, growth-only, and selective use shared full-data targets. Honest generally has a different target and selected population; generic width differences are not common-target efficiency comparisons.",
  "Alternative-model rejection/discovery rates are not automatically common-parameter power. Null calibration uses gamma=0; incidental zero targets use abs(truth)<1e-10 for diagnostics.","",
  paste("Total inference failures:",sum(tree$failures)),paste("Total CI failures:",sum(tree$ci_failures)),
  paste("Known-scale null selective pooled reject-rate range:",fmt_range(null[method=="proposed"]$reject_rate)),
  paste("Known-scale null naive pooled reject-rate range:",fmt_range(null[method=="naive"]$reject_rate)),
  paste("Known-scale null growth-only pooled reject-rate range:",fmt_range(null[method=="growth_only"]$reject_rate)),
  paste("Known-scale null honest pooled reject-rate range:",fmt_range(null[method=="honest"]$reject_rate)),"",
  "Cross-configuration ranges include sparse cells and are descriptive, not simultaneous calibration evidence.",
  if(length(violations))c("NOTES:",unique(violations))else"No additional analysis notes.")
writeLines(report,file.path(opt$output,"report.txt"))
writeLines(c(paste("stage",opt$stage),paste("files",length(files)),paste("complete",complete)),file.path(opt$output,"ANALYSIS_STATUS.txt"))
cat(paste(report,collapse="\n"),"\n")
