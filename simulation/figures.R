# Generate manuscript figures and tables from this working copy's simulation summaries.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
arg<-grep("^--file=",commandArgs(),value=TRUE)
if(length(arg))setwd(dirname(normalizePath(sub("^--file=","",arg[1L]))))
root<-normalizePath("..")
args<-commandArgs(trailingOnly=TRUE);stage<-"main"
if(length(args)) {
  stopifnot(length(args)==2L,args[1L]=="--stage");stage<-args[2L]
}
stopifnot(stage%in%c("main","quick"))
summary_dir<-file.path(root,"results","simulation",stage,"summaries")
revision<-file.path(root,"figures","simulation",stage)
tdir<-file.path(revision,"tables");fdir<-revision
for(d in c(tdir,file.path(fdir,"data")))dir.create(d,recursive=TRUE,showWarnings=FALSE)
tree<-fread(file.path(summary_dir,"tree_summary.csv"))
depth<-fread(file.path(summary_dir,"depth_summary.csv"))
recovery<-fread(file.path(summary_dir,"recovery_summary.csv"))
paired<-fread(file.path(summary_dir,"sigma_paired_summary.csv"))
one<-function(z) {stopifnot(nrow(z)==1L);z}
pp<-function(x)ifelse(is.na(x),"---",sprintf("%.2f",100*x))
decimal<-function(x)ifelse(is.na(x),"---",sprintf("%.2f",x))
int<-function(x)ifelse(is.na(x),"---",formatC(x,format="d",big.mark=","))
rate_se<-function(x,s)ifelse(is.na(x),"---",paste0(pp(x)," (",pp(s),")"))
signed_se<-function(x,s)ifelse(is.na(x),"---",paste0(sprintf("%+.2f",100*x)," (",pp(s),")"))
rowline<-function(x)paste0(paste(x,collapse=" & ")," \\\\")
write_table<-function(file,caption,label,columns,header,rows,notes=NULL) {
  if(stage=="quick")caption<-paste("Quick execution check (not manuscript results).",caption)
  x<-c("\\begin{table}[tbp]","\\centering",paste0("\\caption{",caption,"}\\label{",label,"}"),
       "\\begingroup\\footnotesize\\setlength{\\tabcolsep}{3pt}",
       paste0("\\begin{tabular}{",columns,"}"),"\\toprule",rowline(header),"\\midrule",
       rows,"\\bottomrule","\\end{tabular}","\\endgroup")
  if(!is.null(notes))x<-c(x,"\\par\\smallskip",paste0("\\parbox{\\linewidth}{\\footnotesize ",notes,"}"))
  writeLines(c(x,"\\end{table}"),file.path(tdir,file))
}

caltab<-tree[shape=="null"&variance_ratio==1&max_depth==3&
                variance_mode=="known"&lambda%in%c(4,8)]
rows<-character()
for(lam in c(4,8))for(nn in c(100,200,400,800)) {
  z<-caltab[lambda==lam&n==nn]
  s<-one(z[method=="proposed"]);g<-one(z[method=="growth_only"])
  nv<-one(z[method=="naive"]);h<-one(z[method=="honest"])
  stopifnot(s$available==g$available,s$available==nv$available)
  rows<-c(rows,rowline(c(nn,lam,int(s$available),int(h$available),
    rate_se(nv$reject_rate,nv$reject_rate_mcse),
    rate_se(g$reject_rate,g$reject_rate_mcse),
    rate_se(s$reject_rate,s$reject_rate_mcse),
    rate_se(h$reject_rate,h$reject_rate_mcse))))
}
write_table("calibration.tex",
  "Splitwise false-positive rates under a constant treatment effect with known variances, $D=3$ and $\\rho=1$. Rates and Monte Carlo standard errors in parentheses are in percent and percentage points, respectively.",
  "tab:dense-calibration","rrrrrrrr",
  c("$n$","$\\lambda$","$N_{\\mathrm{full}}$","$N_{\\mathrm H}$","Naive","Growth-only","Selective","Honest"),
  rows,
  "$N_{\\mathrm{full}}$ is the total number of successfully inferred retained splits shared by Naive, Growth-only and Selective; $N_{\\mathrm H}$ is the corresponding count for Honest's own trees and holdout targets. Each row uses the trial count recorded in tree_summary.csv. Pooled-node standard errors use trials as clusters, including trials without a retained target; the counts are not numbers of independent replications.")
fwrite(caltab,file.path(tdir,"calibration_data.csv"))
if(stage=="main")stopifnot(one(caltab[n==100&lambda==4&method=="proposed"])$available==7854,
          one(caltab[n==100&lambda==4&method=="honest"])$available==2366,
          round(one(caltab[n==100&lambda==4&method=="proposed"])$reject_rate*7854)==389,
          round(one(caltab[n==100&lambda==4&method=="honest"])$reject_rate*2366)==99)

wt<-depth[shape=="g10"&variance_ratio==1&max_depth==3&lambda==4&
              variance_mode=="known"&method%in%c("proposed","honest")&node_depth%in%c(0,2)]
wt[,stage:=node_depth+1L]
rows<-character()
for(nn in c(100,200,400,800))for(st in c(1L,3L))for(mm in c("proposed","honest")) {
  z<-one(wt[n==nn&stage==st&method==mm])
  nci<-int(z$ci_n)
  if(z$ci_n>0&&z$ci_n<30)nci<-paste0(nci,"$^{\\dagger}$")
  rows<-c(rows,rowline(c(nn,st,if(mm=="proposed")"Selective"else"Honest",
                         int(z$available),nci,decimal(z$width_median),decimal(z$width_p90))))
}
write_table("widths.tex",
  "Interval widths by actual split stage with known variances, for $\\gamma=a=\\rho=1$, $D=3$ and $\\lambda=4$. Stage 1 is the root and stage 3 corresponds to node depth two.",
  "tab:dense-widths","rrlrrrr",
  c("$n$","Stage","Method","$N$","$N_{\\mathrm{CI}}$","Median","90th percentile"),
  rows,
  "$N$ counts successfully inferred splits over all generated trials; $N_{\\mathrm{CI}}$ counts intervals at the prespecified trial-index subset. A dagger flags fewer than 30 computed intervals. Widths use all successful intervals in that subset without capping or trimming. The two procedures generally have different projection targets and survivor populations, including within the same stage. The Honest training sample at $n=100$ cannot form a third-stage split under the common minimum-cell rule; its unavailable quantiles are shown as dashes.")
fwrite(wt,file.path(tdir,"widths_data.csv"))
if(stage=="main")stopifnot(one(wt[n==100&method=="proposed"&stage==1])$ci_n==500,
          one(wt[n==100&method=="honest"&stage==1])$ci_n==478,
          one(wt[n==100&method=="honest"&stage==3])$available==0)

pt<-tree[shape=="g10"&max_depth==3&lambda==4&n%in%c(100,800)&
              method%in%c("proposed","honest")]
pd<-paired[shape=="g10"&max_depth==3&lambda==4&n%in%c(100,800)&
              method%in%c("proposed","honest")]
rows<-character()
for(nn in c(100,800))for(rho in c(.25,1,4))for(mm in c("proposed","honest")) {
  k<-one(pt[n==nn&variance_ratio==rho&method==mm&variance_mode=="known"])
  h<-one(pt[n==nn&variance_ratio==rho&method==mm&variance_mode=="plugin"])
  d<-one(pd[n==nn&variance_ratio==rho&method==mm])
  stopifnot(k$available==h$available,d$n_pairs==k$available)
  rows<-c(rows,rowline(c(nn,if(rho==.25)"$1/4$"else rho,
    if(mm=="proposed")"Selective"else"Honest",int(k$available),
    rate_se(k$coverage,k$coverage_mcse),rate_se(h$coverage,h$coverage_mcse),
    signed_se(d$coverage_difference,d$coverage_difference_mcse))))
}
write_table("plugin.tex",
  "Coverage after substituting the estimated common scale, for $\\gamma=a=1$, $D=3$ and $\\lambda=4$.",
  "tab:dense-plugin","rrlrrrr",
  c("$n$","$\\rho$","Method","$N$","Known scale","Estimated scale","Difference"),
  rows,
  "Coverage is in percent; parentheses give trial-cluster Monte Carlo standard errors in percentage points. Difference is estimated minus known coverage, with its paired trial-cluster standard error. Within each method the two scale versions use the same data, trees, targets and $N$ inferred splits. Honest is evaluated on its own holdout targets and estimates the scale using training controls only. These plug-in calculations retain the known variance ratio and do not carry a finite-sample exactness guarantee.")
fwrite(pt,file.path(tdir,"plugin_data.csv"))
fwrite(pd,file.path(tdir,"plugin_paired_data.csv"))

# Figures use readable labels at intended manuscript print scale.
colors<-c(naive="#B64850",growth_only="#B58323",proposed="#087E8B",honest="#5267A8")
labels<-c(naive="Naive",growth_only="Growth-only",proposed="Selective",honest="Honest")
decor<-function(z) {
  z<-copy(z);z[,method:=factor(method,names(colors))]
  z[,n_label:=factor(paste0("n = ",n),paste0("n = ",c(100,200,400,800)))]
  z[,rho_label:=factor(paste0("rho = ",variance_ratio),paste0("rho = ",c(.25,1,4)))]
  if("node_depth"%in%names(z))z[,stage:=node_depth+1L]
  z
}
theme_print<-theme_bw(base_size=9,base_family="Helvetica")+
 theme(panel.grid.minor=element_blank(),legend.position="bottom",legend.title=element_blank(),
       strip.background=element_rect(fill="#F1F3F5"),strip.text=element_text(size=8.5),
       legend.text=element_text(size=8),panel.spacing=grid::unit(.55,"lines"),
       axis.text=element_text(size=8),axis.title=element_text(size=9))
savefig<-function(name,p,data,width=7,height=5.5) {
  fwrite(data,file.path(fdir,"data",paste0(name,".csv")))
  pdf<-file.path(fdir,paste0(name,".pdf"))
  if(stage=="quick")p<-p+labs(subtitle="Quick execution check; not manuscript results")
  ggsave(pdf,p,width=width,height=height,units="in",device=grDevices::pdf,
         family="Helvetica",useDingbats=FALSE)
  cat("Wrote ",name,"\n",sep="")
}

# A: all maximum growth limits D1/2/3; full y-range retains Naive over-rejection.
nullfig<-decor(tree[shape=="null"&variance_mode=="known"&lambda==4])
stopifnot(nrow(nullfig)==144L,setequal(nullfig$max_depth,1:3))
p<-ggplot(nullfig,aes(max_depth,reject_rate,color=method,shape=method,group=method))+
 geom_hline(yintercept=.05,linetype=2,linewidth=.35,color="#45494D")+
 geom_line(linewidth=.4,na.rm=TRUE)+geom_point(size=1.6,na.rm=TRUE)+
 geom_errorbar(aes(ymin=reject_rate_lo,ymax=reject_rate_hi),width=.07,linewidth=.3,na.rm=TRUE)+
 scale_color_manual(values=colors,labels=labels)+
 scale_shape_manual(values=c(naive=1,growth_only=15,proposed=19,honest=17),labels=labels)+
 facet_grid(rho_label~n_label)+scale_x_continuous(breaks=1:3)+
 scale_y_continuous(breaks=c(0,.25,.5,.75,1),limits=c(0,1))+
 labs(x="Maximum growth depth",y="Pooled split false-positive rate")+theme_print
savefig("null_calibration",p,nullfig,7,5.5)

# B: maximum ARI match to the reference-root partition, not unique-tree recovery.
z<-decor(recovery[true_node=="root"&a==1&variance_mode=="known"&max_depth==3&
                  lambda==4&variance_ratio==1&method%in%c("proposed","honest")])
first<-copy(z);first[,event:="Reference agreement"]
first[,probability:=recovered_rate];first[,lo:=recovered_rate_lo];first[,hi:=recovered_rate_hi]
second<-copy(z);second[,event:="Agreement and rejection"]
second[,probability:=recovered_and_rejected_rate]
second[,lo:=recovered_and_rejected_rate_lo];second[,hi:=recovered_and_rejected_rate_hi]
rootfig<-rbindlist(list(first,second))
rootfig[,event:=factor(event,c("Reference agreement","Agreement and rejection"))]
pdodge<-position_dodge(width=.055)
p<-ggplot(rootfig,aes(gamma,probability,color=method,shape=method,linetype=event,
                group=interaction(method,event)))+
 geom_line(linewidth=.45,position=pdodge)+geom_point(size=1.8,position=pdodge)+
 geom_errorbar(aes(ymin=lo,ymax=hi),width=.05,linewidth=.3,position=pdodge)+
 scale_color_manual(values=colors,labels=labels)+
 scale_shape_manual(values=c(proposed=19,honest=17),labels=labels)+
 scale_linetype_manual(values=c(1,2))+
 facet_wrap(~n_label,nrow=2)+scale_x_continuous(breaks=c(.5,1,2))+
 scale_y_continuous(breaks=seq(0,1,.25),limits=c(0,1))+
 labs(x="Effect amplitude gamma (a = 1)",y="Probability per trial")+theme_print+
 theme(legend.box="vertical")
savefig("root_recovery",p,rootfig,7,4.9)

# C: all actual stages, retaining the observed upper-tail summaries.
z<-decor(depth[shape=="g10"&variance_mode=="known"&max_depth==3&lambda==4&
                 method%in%c("proposed","honest")])
widthfig<-melt(z,id.vars=setdiff(names(z),c("width_median","width_p90")),
  measure.vars=c("width_median","width_p90"),variable.name="quantile",value.name="interval_width")
widthfig[,quantile:=factor(quantile,c("width_median","width_p90"),c("Median","90th percentile"))]
p<-ggplot(widthfig,aes(stage,interval_width,color=method,linetype=quantile,
                group=interaction(method,quantile)))+
 geom_line(linewidth=.45,na.rm=TRUE)+geom_point(size=1.5,na.rm=TRUE)+
 geom_point(data=widthfig[ci_n<30&is.finite(interval_width)],shape=1,color="#25262A",size=2.5)+
 scale_color_manual(values=colors,labels=labels)+scale_linetype_manual(values=c(1,2))+
 scale_y_log10()+scale_x_continuous(breaks=1:3)+facet_grid(rho_label~n_label)+
 labs(x="Actual split stage (root = 1)",y="95% interval width (log scale)")+theme_print
savefig("interval_width",p,widthfig,7,5.5)

# D: rho1 reference panel; no data are capped, and the vertical axis is zoomed.
pluginfig<-decor(depth[shape=="g10"&max_depth==3&lambda==4&variance_ratio==1&
                      method%in%c("proposed","honest")])
pluginfig[,variance_mode:=factor(variance_mode,c("known","plugin"),c("Known scale","Plug-in scale"))]
if(stage=="main")stopifnot(all(pluginfig$coverage[is.finite(pluginfig$coverage)]>=.85))
p<-ggplot(pluginfig,aes(stage,coverage,color=method,linetype=variance_mode,
                group=interaction(method,variance_mode)))+
 geom_hline(yintercept=.95,linetype=2,linewidth=.35,color="#45494D")+
 geom_line(linewidth=.45,na.rm=TRUE)+geom_point(size=1.5,na.rm=TRUE)+
 geom_errorbar(aes(ymin=coverage_lo,ymax=coverage_hi),width=.07,linewidth=.3,na.rm=TRUE)+
 geom_point(data=pluginfig[small_cell==TRUE&is.finite(coverage)],shape=1,
            color="#25262A",size=2.8,na.rm=TRUE)+
 scale_color_manual(values=colors,labels=labels)+scale_linetype_manual(values=c(1,2))+
 scale_x_continuous(breaks=1:3)+scale_y_continuous(breaks=c(.85,.9,.95,1))+
 coord_cartesian(ylim=if(stage=="main")c(.85,1)else c(0,1))+facet_wrap(~n_label,nrow=2)+
 labs(x="Actual split stage (root = 1)",y="Coverage")+theme_print
savefig("plugin_coverage",p,pluginfig,7,4.9)

cat("Created three tables and four PDF figures from ",stage," summaries.\n",sep="")
