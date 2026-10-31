args<-commandArgs(trailingOnly=TRUE)
filearg<-grep("^--file=",commandArgs(),value=TRUE)
if(length(filearg))setwd(dirname(normalizePath(sub("^--file=","",filearg[1]))))
value<-function(name,default) {
  z<-grep(paste0("^--",name,"="),args,value=TRUE)
  if(length(z))sub(paste0("^--",name,"="),"",tail(z,1L)) else default
}
stage<-value("stage","main");stopifnot(stage%in%c("quick","main"))
workers<-as.integer(value("workers","1"));chunk_size<-as.integer(value("chunk","50"))
stopifnot(workers>=1L,chunk_size>=1L)
if(.Platform$OS.type=="windows"&&workers>1L)stop("Use --workers=1 on Windows")
root<-normalizePath("..")
source("study_dense.R")
if(stage=="quick") {
  nrep<-as.integer(value("reps","2"));width<-as.integer(value("width","2"))
  stopifnot(nrep>=1L,width>=0L)
  scenarios<-make_scenarios(nrep,nrep,width)
} else {
  scenarios<-make_scenarios(null_reps=5000L,alternative_reps=1000L,width_reps=500L)
}
lambdas<-c(0,4,8)
simulation_run_dir<-file.path(root,"results","simulation",stage)
folder<-file.path(simulation_run_dir,"chunks");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
source_files<-c("gls_tree.R","gls_pruning.R","dense_tree.R","stable_tn.R","study_dense.R","run.R")
stopifnot(all(file.exists(source_files)))
config<-list(stage=stage,scenarios=scenarios,lambdas=lambdas,depths=dense_depths,
 methods=dense_methods,plugin=list(shapes=c("null","g10"),lambda=4),
 width_rule="fixed first width_reps trial indices, every method; failures retained",
 source_md5=as.list(tools::md5sum(source_files)),chunk_size=chunk_size,
 seed_main=900000000L,seed_pilot=700000000L,evaluation_seed=6091031L,RNGkind=RNGkind())
cfgfile<-file.path(folder,"frozen_config.rds")
if(file.exists(cfgfile))stopifnot(identical(readRDS(cfgfile),config)) else saveRDS(config,cfgfile)
write.csv(scenarios,file.path(folder,"scenarios.csv"),row.names=FALSE)
save_atomic<-function(x,path) {
  tmp<-paste0(path,".tmp_",Sys.getpid());saveRDS(x,tmp,compress="gzip")
  if(!file.rename(tmp,path))stop("Could not commit ",path)
}
tasks<-list()
for(bg in unique(scenarios$base_group)) {
  sc<-scenarios[scenarios$base_group==bg,,drop=FALSE]
  groups<-split(seq_len(sc$repetitions[1]),ceiling(seq_len(sc$repetitions[1])/chunk_size))
  for(k in seq_along(groups))tasks[[length(tasks)+1L]]<-list(bg=bg,chunk=k,reps=as.integer(groups[[k]]))
}
check_chunk<-function(obj,sc,reps) {
  z<-obj$records
  stopifnot(is.data.frame(z),setequal(z$replication,reps),all(z$scenario==sc$id),
    all(z$n==sc$n),all(z$shape==sc$shape),all(z$variance_ratio==sc$variance_ratio),
    all(z$available==z$eligible),all(!z$eligible|z$retained),all(!z$retained|z$grown),
    !anyDuplicated(z[c("replication","max_depth","lambda","method","variance_mode","node")]))
  expected<-length(reps)*sum(2^dense_depths-1L)*length(lambdas)*4L
  if(sc$shape%in%c("null","g10")&&4%in%lambdas)
    expected<-expected+length(reps)*sum(2^dense_depths-1L)*4L
  stopifnot(nrow(z)==expected,nrow(obj$trials)==length(reps))
  invisible(TRUE)
}
run_task<-function(task) {
  scs<-scenarios[scenarios$base_group==task$bg,,drop=FALSE]
  files<-file.path(folder,sprintf("scenario_%02d_chunk_%03d.rds",scs$id,task$chunk))
  if(all(file.exists(files))) {
    for(j in 1:3)check_chunk(readRDS(files[j]),scs[j,,drop=FALSE],task$reps)
    return(files)
  }
  result<-vector("list",3L)
  for(j in 1:3)result[[j]]<-vector("list",length(task$reps))
  for(k in seq_along(task$reps)) {
    r<-task$reps[k];dat<-generate_base(scs[1,,drop=FALSE],r,stage)
    for(j in 1:3)result[[j]][[k]]<-one_version(dat,scs[j,,drop=FALSE],r,stage,
                           lambdas,r<=scs$inversion_repetitions[j])
  }
  fail<-ci_fail<-0L
  for(j in 1:3) {
    obj<-lapply(c("records","trials","recovery"),function(nm)
      do.call(rbind,lapply(result[[j]],`[[`,nm)))
    names(obj)<-c("records","trials","recovery");check_chunk(obj,scs[j,,drop=FALSE],task$reps)
    fail<-fail+sum(obj$records$inference_status=="failed")
    ci_fail<-ci_fail+sum(obj$records$ci_status=="failed")
    if(file.exists(files[j]))check_chunk(readRDS(files[j]),scs[j,,drop=FALSE],task$reps)
    else save_atomic(obj,files[j])
  }
  cat(sprintf("%s bg=%02d n=%d %s reps=%d-%d inference_fail=%d ci_fail=%d\n",
              stage,task$bg,scs$n[1],scs$shape[1],min(task$reps),max(task$reps),fail,ci_fail))
  flush.console();files
}
t0<-proc.time()[3]
writeLines(sprintf("%s START %s workers=%d tasks=%d",Sys.time(),stage,workers,length(tasks)),
           file.path(folder,"status.txt"))
files<-parallel::mclapply(tasks,run_task,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
if(any(vapply(files,inherits,logical(1),"try-error")))stop("Worker failed; checkpoints retained")
files<-unlist(files,use.names=FALSE)
stats<-vapply(files,function(path) {
  x<-readRDS(path)
  c(rows=nrow(x$records),versions=nrow(x$trials),
    inference_fail=sum(x$records$inference_status=="failed"),ci_fail=sum(x$records$ci_status=="failed"))
},numeric(4))
status<-c(sprintf("%s COMPLETE %s",Sys.time(),stage),
 sprintf("Independent bundles: %d",sum(stats["versions",])/3L),
 sprintf("Outcome versions: %d",sum(stats["versions",])),sprintf("Rows: %d",sum(stats["rows",])),
 sprintf("Elapsed wall seconds: %.3f",unname(proc.time()[3]-t0)),
 sprintf("Inference failures: %d",sum(stats["inference_fail",])),sprintf("CI failures: %d",sum(stats["ci_fail",])))
writeLines(status,file.path(folder,"status.txt"));cat(paste(status,collapse="\n"),"\n")
capture.output(sessionInfo(),file=file.path(folder,"session.txt"))
