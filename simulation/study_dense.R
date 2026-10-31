# Data generation and interaction-tree inference for the manuscript simulations.
source("gls_pruning.R")
source("dense_tree.R")
if (file.exists("stable_tn.R")) source("stable_tn.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
dense_depths <- 1:3
dense_methods <- c("naive", "growth_only", "proposed", "honest")
dense_effects <- data.frame(shape=c("null","g05","g10","g20","weak_a"),
                            gamma=c(0,.5,1,2,1),a=c(1,1,1,1,.5))
make_scenarios <- function(null_reps=5000L, alternative_reps=1000L,
                           width_reps=500L) {
  base <- merge(data.frame(n=c(100L,200L,400L,800L)),dense_effects,sort=FALSE)
  base <- base[order(base$n,match(base$shape,dense_effects$shape)),]
  base$base_group <- seq_len(nrow(base))
  out <- do.call(rbind,lapply(seq_len(nrow(base)),function(i) {
    x<-base[rep(i,3L),];x$variance_ratio<-c(.25,1,4);x
  }))
  out$id<-seq_len(nrow(out));out$p<-.5
  out$repetitions<-ifelse(out$gamma==0,null_reps,alternative_reps)
  out$inversion_repetitions<-pmin(out$repetitions,width_reps)
  rownames(out)<-NULL;out
}
study_seed <- function(base_group,r,stage) {
  as.integer((if(stage=="main")900000000L else 700000000L)+
               100000L*base_group+r)
}
generate_base <- function(sc,r,stage="main") {
  set.seed(study_seed(sc$base_group,r,stage))
  X<-matrix(runif(sc$n*5L,-1,1),sc$n,5L);A<-rbinom(sc$n,1,.5)
  B<-cbind(1,X);eps<-rnorm(sc$n)
  tau<-.5+sc$gamma*(X[,1]<=0)*(1+sc$a*(X[,2]>0)+(X[,2]*X[,3]>0))
  mu<-2*X[,3]+X[,4]-X[,5]+A*tau
  training<-sample.int(sc$n,sc$n%/%2L)
  list(X=X,A=A,B=B,eps=eps,tau=tau,mu=mu,training=training,
       holdout=setdiff(seq_len(sc$n),training))
}
potential_nodes <- function(D) {
  x<-"root"
  if(D>1L)for(j in seq_len(D-1L))
    x<-c(x,sort(apply(expand.grid(rep(list(c("L","R")),j),
                                stringsAsFactors=FALSE),1,paste0,collapse="")))
  x
}
prefix_tree <- function(tree,D) tree[vapply(tree,function(z)z$depth<D,logical(1))]
node_level <- function(name) ifelse(name=="root",0L,nchar(name))
control_scale <- function(Y,A,B) {
  ii<-which(A==0);qq<-qr(B[ii,,drop=FALSE]);df<-length(ii)-qq$rank
  if(df<=0L)stop("No residual degrees of freedom for control scale")
  s<-sqrt(sum(qr.resid(qq,Y[ii])^2)/df)
  if(!is.finite(s)||s<=0)stop("Nonpositive estimated scale")
  s
}

# All methods use their own actual finite-design projection target.
infer_from_support <- function(node,name,Y,mu,support,sigma,invert) {
  t<-sum(node$u*Y);truth_u<-sum(node$u*mu);truth<-truth_u/node$norm
  lo<-tn_logodds(t,0,sigma,support);lt<-tn_logodds(t,truth_u,sigma,support)
  ci<-c(NA_real_,NA_real_);err<-"";status<-"not_requested"
  if(invert) {
    ci<-tryCatch({
      if(nrow(support)==1L&&all(support[1,]==c(-Inf,Inf)))
        t+c(-1,1)*qnorm(.975)*sigma else tn_interval(t,sigma,support)
    },error=function(e){err<<-conditionMessage(e);c(NA_real_,NA_real_)})
    status<-if(all(is.finite(ci)))"ok" else "failed"
    if(status=="ok") {
      ee<-abs(c(tn_logodds(t,ci[1],sigma,support)-qlogis(.975),
                tn_logodds(t,ci[2],sigma,support)-qlogis(.025)))
      if(any(!is.finite(ee))||max(ee)>2e-5) {
        status<-"failed";err<-paste("Endpoint logodds residual",max(ee))
      }
    }
    ci<-ci/node$norm
  }
  endpoints<-support[is.finite(support)]
  data.frame(method="",node=name,available=TRUE,split=node$spec$id,
    variable=node$spec$j,threshold=node$spec$cut,estimate=t/node$norm,
    truth=truth,p_value=2*plogis(-abs(lo)),pivot_truth=plogis(lt),
    covered=abs(lt)<=qlogis(.975),lower=ci[1],upper=ci[2],width=ci[2]-ci[1],
    endpoint_covered=if(all(is.finite(ci)))ci[1]<=truth&&truth<=ci[2] else NA,
    ci_status=status,ci_error=err,components=nrow(support),
    boundary_distance=if(length(endpoints))min(abs(t-endpoints)) else Inf,
    inference_status="ok",inference_error="",stringsAsFactors=FALSE)
}
base_columns <- names(unavailable_row("naive","root","not grown"))
node_indices <- function(tree,name,n) {
  if(name=="root")return(seq_len(n))
  parent<-if(nchar(name)==1L)"root" else substr(name,1,nchar(name)-1L)
  if(is.null(tree[[parent]]))return(NULL)
  tree[[parent]]$spec[[if(substr(name,nchar(name),nchar(name))=="L")"left" else "right"]]
}
node_metadata <- function(tree,positions,n) {
  out<-data.frame(node=positions,reached=FALSE,candidate_count=NA_integer_,
                  grown_parent_n=NA_integer_,grown_split=NA_character_,
                  grown_variable=NA_integer_,grown_threshold=NA_real_,
                  grown_contrast_norm=NA_real_)
  for(i in seq_along(positions)) {
    nm<-positions[i];idx<-node_indices(tree,nm,n)
    if(!is.null(idx)) {
      out$reached[i]<-TRUE;out$candidate_count[i]<-0L
      out$grown_parent_n[i]<-length(idx)
    }
    nd<-tree[[nm]]
    if(!is.null(nd)) {
      out$candidate_count[i]<-length(nd$candidates$specs)
      out$grown_split[i]<-nd$spec$id;out$grown_variable[i]<-nd$spec$j
      out$grown_threshold[i]<-nd$spec$cut;out$grown_contrast_norm[i]<-nd$norm
    }
  }
  out
}
holdout_details <- function(tree,name,X,A,B,shape) {
  idx<-seq_len(nrow(X))
  for(ancestor in head(tree[[name]]$path,-1L)) {
    s<-tree[[ancestor]]$spec
    d<-substr(name,tree[[ancestor]]$depth+1L,tree[[ancestor]]$depth+1L)
    yes<-X[idx,s$j]<=s$cut;idx<-idx[if(d=="L")yes else !yes]
  }
  reason<-"";h<-X[idx,tree[[name]]$spec$j]<=tree[[name]]$spec$cut
  if(length(idx)<ncol(B)+4L)reason<-"parent_too_small"
  else if(min(table(factor(A[idx],0:1),factor(h,c(FALSE,TRUE))))<1L)
    reason<-"empty_treatment_child"
  nd<-if(reason=="")holdout_node(tree,name,X,A,B,variance_shape=shape) else NULL
  if(is.null(nd)&&reason=="")reason<-"unidentifiable"
  list(node=nd,parent_n=length(idx),reason=reason)
}

# Independent fixed evaluation covariates for procedure-comparable recovery.
.evaluation_X <- local({set.seed(6091031L);matrix(runif(4096L*5L,-1,1),4096L,5L)})
.true_split_labels <- function() {
  X<-.evaluation_X
  parents<-list(root=rep(TRUE,nrow(X)),L=X[,1]<=0,
                LL=X[,1]<=0&X[,2]<=0,LR=X[,1]<=0&X[,2]>0)
  variables<-c(root=1L,L=2L,LL=3L,LR=3L)
  lapply(names(parents),function(nm)list(name=nm,parent=parents[[nm]],
    labels=as.integer(X[parents[[nm]],variables[nm]]>0)))
}
.true_labels <- .true_split_labels()
.ari <- function(x,y) {
  tt<-table(x,y);N<-sum(tt);if(N<2L)return(NA_real_)
  c2<-function(z)sum(z*(z-1)/2)
  ab<-c2(tt);aa<-c2(rowSums(tt));bb<-c2(colSums(tt));tot<-N*(N-1)/2
  ex<-aa*bb/tot;den<-(aa+bb)/2-ex
  if(den==0)return(if(identical(as.integer(x),as.integer(y)))1 else 0)
  (ab-ex)/den
}
recovery_geometry <- function(tree) {
  if(!length(tree))return(matrix(numeric(),nrow=4L,ncol=0L,
                                dimnames=list(c("root","L","LL","LR"),NULL)))
  labels<-list();X<-.evaluation_X
  for(nm in names(tree)) {
    parent<-rep(TRUE,nrow(X))
    for(an in head(tree[[nm]]$path,-1L)) {
      s<-tree[[an]]$spec;dir<-substr(nm,tree[[an]]$depth+1L,tree[[an]]$depth+1L)
      yes<-X[,s$j]<=s$cut;parent<-parent&if(dir=="L")yes else !yes
    }
    s<-tree[[nm]]$spec
    labels[[nm]]<-ifelse(!parent,2L,as.integer(X[,s$j]>s$cut))
  }
  ans<-vapply(labels,function(y)vapply(.true_labels,function(z)
    .ari(z$labels,y[z$parent]),numeric(1)),numeric(4))
  rownames(ans)<-c("root","L","LL","LR");ans
}

one_version <- function(dat,sc,r,stage,lambdas,invert,plugin_lambda=4) {
  t0<-proc.time()[3];X<-dat$X;A<-dat$A;B<-dat$B;mu<-dat$mu
  shape<-ifelse(A==1,sc$variance_ratio,1);Y<-mu+sqrt(shape)*dat$eps
  train<-dat$training;hold<-dat$holdout
  grown<-fit_tree(X,A,B,Y,max_depth=3L,variance_shape=shape)
  trained<-fit_tree(X[train,,drop=FALSE],A[train],B[train,,drop=FALSE],Y[train],
                    max_depth=3L,variance_shape=shape[train])
  grow_seconds<-unname(proc.time()[3]-t0)
  gs<-setNames(lapply(names(grown),function(nm)selection_support(grown,nm,Y,"path")),names(grown))
  hs<-setNames(lapply(names(trained),function(nm)holdout_details(trained,nm,
    X[hold,,drop=FALSE],A[hold],B[hold,,drop=FALSE],shape[hold])),names(trained))
  sigma_hat<-control_scale(Y,A,B)
  honest_sigma_hat<-control_scale(Y[train],A[train],B[train,,drop=FALSE])
  plugin<-sc$shape%in%c("null","g10")
  geom<-if(sc$gamma>0)recovery_geometry(grown) else NULL
  hgeom<-if(sc$gamma>0)recovery_geometry(trained) else NULL
  cache<-new.env(parent=emptyenv());blocks<-list();recoveries<-list()
  inference<-function(node,nm,method,support,sigma) {
    key<-paste(if(method=="honest")"hold" else "full",nm,
       paste(format(support,digits=17),collapse=":"),format(sigma,digits=17),sep="|")
    if(exists(key,cache,inherits=FALSE))return(get(key,cache,inherits=FALSE))
    z<-tryCatch(infer_from_support(node,nm,if(method=="honest")Y[hold] else Y,
          if(method=="honest")mu[hold] else mu,support,sigma,invert),
      error=function(e)unavailable_row(method,nm,conditionMessage(e),"failed"))
    z$available<-TRUE;assign(key,z,cache);z
  }
  for(D in dense_depths) {
    tree<-prefix_tree(grown,D);tr<-prefix_tree(trained,D);pos<-potential_nodes(D)
    full_meta<-node_metadata(tree,pos,nrow(X));honest_meta<-node_metadata(tr,pos,length(train))
    for(lambda in lambdas) {
      pr<-prune_tree(tree,Y,lambda);hp<-prune_tree(tr,Y[train],lambda)
      modes<-if(plugin&&lambda==plugin_lambda)c("known","plugin") else "known"
      for(vm in modes)for(method in dense_methods) {
        ishonest<-method=="honest";tt<-if(ishonest)tr else tree
        kept<-if(ishonest)hp$kept else pr$kept
        z<-unavailable_row(method,"root","not grown")[rep(1L,length(pos)),,drop=FALSE]
        rownames(z)<-NULL;z$node<-pos
        z$grown<-pos%in%names(tt);z$retained<-pos%in%kept
        z$holdout_eligible<-NA;z$eligible<-z$retained;z$delta<-NA_real_
        z$holdout_parent_n<-NA_integer_;z$holdout_contrast_norm<-NA_real_
        z$holdout_failure_reason<-NA_character_
        z$inference_error[z$grown&!z$retained]<-"pruned"
        for(i in which(z$retained)) {
          nm<-pos[i];nd<-tt[[nm]]
          if(ishonest) {
            hd<-hs[[nm]];nd<-hd$node
            z$holdout_parent_n[i]<-hd$parent_n
            z$holdout_failure_reason[i]<-hd$reason
            z$holdout_eligible[i]<-!is.null(nd);z$eligible[i]<-!is.null(nd)
            if(is.null(nd)) {z$inference_error[i]<-hd$reason;next}
            z$holdout_contrast_norm[i]<-nd$norm
          }
          sigma<-if(vm=="known")1 else if(ishonest)honest_sigma_hat else sigma_hat
          support<-if(method%in%c("naive","honest"))whole_line() else gs[[nm]]$support
          if(method=="proposed") {
            ss<-selection_support_pruned(tree,nm,Y,lambda,prune=pr,grow_support=gs[[nm]])
            support<-ss$support;z$delta[i]<-ss$delta
          }
          rr<-inference(nd,nm,method,support,sigma)
          z[i,base_columns]<-rr[,base_columns,drop=FALSE]
          if(rr$inference_status=="failed"||rr$ci_status=="failed") {
            dir.create(file.path(simulation_run_dir,"failures"),recursive=TRUE,showWarnings=FALSE)
            ff<-file.path(simulation_run_dir,"failures",sprintf("s%02d_r%05d_D%d_l%g_%s_%s_%s.rds",
                            sc$id,r,D,lambda,method,vm,nm))
            if(!file.exists(ff))saveRDS(list(scenario=sc,replication=r,stage=stage,
              method=method,variance_mode=vm,max_depth=D,lambda=lambda,node=nm,
              result=rr,support=support,sigma=sigma,t=sum(nd$u*(if(ishonest)Y[hold] else Y)),
              norm=nd$norm),ff)
          }
        }
        z$method<-method;z$available<-z$eligible
        md<-if(ishonest)honest_meta else full_meta
        z<-cbind(z,md[,setdiff(names(md),"node"),drop=FALSE])
        z$node_depth<-node_level(pos);z$max_depth<-D;z$lambda<-lambda;z$variance_mode<-vm
        z$ci_planned<-invert
        blocks[[length(blocks)+1L]]<-z
        if(sc$gamma>0) {
          gg<-if(ishonest)hgeom else geom;kk<-intersect(kept,colnames(gg))
          rg<-data.frame(true_node=c("root","L","LL","LR"),true_depth=c(0L,1L,2L,2L),
             match_node=NA_character_,recovery_ari=NA_real_,recovered=FALSE,
             eligible=FALSE,rejected=FALSE)
          if(length(kk))for(i in 1:4) {
            ar<-gg[i,kk];j<-which.max(ar);nm<-kk[j];zz<-z[z$node==nm,,drop=FALSE]
            rg$match_node[i]<-nm;rg$recovery_ari[i]<-ar[j]
            rg$recovered[i]<-ar[j]>.75
            rg$eligible[i]<-zz$inference_status=="ok"
            rg$rejected[i]<-rg$recovered[i]&&rg$eligible[i]&&zz$p_value<.05
          }
          rg$method<-method;rg$variance_mode<-vm;rg$max_depth<-D;rg$lambda<-lambda
          recoveries[[length(recoveries)+1L]]<-rg
        }
      }
    }
  }
  add_meta<-function(x) {
    if(is.null(x)||!nrow(x))return(data.frame())
    x$scenario<-sc$id;x$base_group<-sc$base_group;x$n<-sc$n;x$shape<-sc$shape
    x$gamma<-sc$gamma;x$a<-sc$a;x$variance_ratio<-sc$variance_ratio
    x$allocation<-.5;x$replication<-r;x$seed<-study_seed(sc$base_group,r,stage)
    x$training_n<-length(train);x$holdout_n<-length(hold);rownames(x)<-NULL;x
  }
  out<-add_meta(do.call(rbind,blocks));rec<-add_meta(do.call(rbind,recoveries))
  trials<-add_meta(data.frame(trial_seconds=unname(proc.time()[3]-t0),
    grow_seconds=grow_seconds,infer_seconds=unname(proc.time()[3]-t0)-grow_seconds,
    full_sigma_hat=sigma_hat,honest_sigma_hat=honest_sigma_hat,
    full_root_candidates=if(length(grown))length(grown$root$candidates$specs) else 0L,
    honest_root_candidates=if(length(trained))length(trained$root$candidates$specs) else 0L,
    full_tree_bytes=as.numeric(object.size(grown)),training_tree_bytes=as.numeric(object.size(trained))))
  stopifnot(all(out$available==out$eligible),all(!out$eligible|out$retained),
            all(!out$retained|out$grown),!anyDuplicated(out[c("max_depth","lambda","method","variance_mode","node")]))
  list(records=out,trials=trials,recovery=rec)
}
