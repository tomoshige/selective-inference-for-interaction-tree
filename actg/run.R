# Recompute all five ACTG penalty settings. Invoke: Rscript actg/run.R
source("actg/analysis.R")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
settings <- make_settings()
results <- pilots <- nodes <- counts <- list()
for (i in seq_len(nrow(settings))) {
  sc <- settings[i,,drop=FALSE]
  cat("ACTG ",sc$setting," (lambda=",sc$lambda,")\n",sep="");flush.console()
  s <- fit_setting(sc)
  results[[i]] <- s$results; pilots[[i]] <- s$pilots
  rows <- list()
  for (kind in c("full","Honest_training")) {
    tr <- if(kind=="full")s$full_tree else s$training_tree
    pr <- if(kind=="full")s$full_prune else s$training_prune
    dat <- if(kind=="full")s$data else s$train
    vv <- if(kind=="full")s$full_variance else s$training_variance
    for (nm in names(tr)) {
      nd <- tr[[nm]]
      rows[[length(rows)+1L]] <- cbind(data.frame(setting=sc$setting,tree=kind,
        node=nm,depth=nd$depth,split=nd$spec$id,variable=nd$spec$variable,
        threshold=nd$spec$cut,retained=nm %in% pr$kept,score=sum(nd$u*dat$Y),
        gain=pr$gains[nm],F=pr$F[nm],V=pr$V[nm],path=paste(nd$path,collapse=";")),
        contrast_metadata(nd,dat,vv))
    }
  }
  nodes[[i]] <- do.call(rbind,rows)
  counts[[i]] <- cbind(sc,data.frame(n=length(s$data$Y),train_n=length(s$train$Y),
    holdout_n=length(s$hold$Y),full_grown=length(s$full_tree),
    full_retained=length(s$full_prune$kept),training_grown=length(s$training_tree),
    training_retained=length(s$training_prune$kept),
    honest_available=length(Filter(function(z)!is.null(z$node),s$holdout))))
}
contrasts <- do.call(rbind,results)
stopifnot(nrow(contrasts)==100L,all(contrasts$inference_status=="ok"),
          all(contrasts$available),all(contrasts$width>0))
write.csv(contrasts,file.path(out_dir,"contrasts.csv"),row.names=FALSE)
write.csv(do.call(rbind,pilots),file.path(out_dir,"pilots.csv"),row.names=FALSE)
write.csv(do.call(rbind,nodes),file.path(out_dir,"nodes.csv"),row.names=FALSE)
write.csv(do.call(rbind,counts),file.path(out_dir,"settings.csv"),row.names=FALSE)
cat("ACTG complete:100 records in results/actg/.\n")
