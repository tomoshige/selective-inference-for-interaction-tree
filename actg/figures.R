# Rebuild the ACTG table and two figures. Invoke: Rscript actg/figures.R
root<-normalizePath(getwd())
out<-file.path(root,"results/actg")
fd<-file.path(root,"figures/actg");td<-out
stopifnot(file.exists(file.path(out,"contrasts.csv")))
dir.create(fd,recursive=TRUE,showWarnings=FALSE)
results<-read.csv(file.path(out,"contrasts.csv"),stringsAsFactors=FALSE)
nodes<-read.csv(file.path(out,"nodes.csv"),stringsAsFactors=FALSE)
settings<-read.csv(file.path(out,"settings.csv"),stringsAsFactors=FALSE)
if(!"lambda"%in%names(nodes))nodes$lambda<-settings$lambda[match(nodes$setting,settings$setting)]
main<-results[results$lambda==0,,drop=FALSE]
stopifnot(nrow(main)>0,!anyNA(main$lambda),all(main$inference_status%in%c("ok","unavailable")))
main<-main[main$method%in%c("naive","proposed","honest"),,drop=FALSE]
full<-main[main$method=="proposed",,drop=FALSE];honest<-main[main$method=="honest",,drop=FALSE]
stopifnot(nrow(full)>0,nrow(honest)>0)
variable<-c(cd40="Baseline CD4",cd80="Baseline CD8",age="Age",preanti="Prior days",
            str2="Prior experience",symptom="Symptoms")
shortvariable<-c(cd40="CD4(0)",cd80="CD8(0)",age="Age",preanti="Prior days",
                str2="Prior experience",symptom="Symptoms")
texvariable<-c(cd40="\\mathrm{CD4}_0",cd80="\\mathrm{CD8}_0",age="\\mathrm{Age}",
               preanti="\\mathrm{PriorDays}",str2="\\mathrm{PriorExp}",symptom="\\mathrm{Symptoms}")
methods<-c("naive","proposed","honest")
colors<-c(naive="#787878",proposed="#1D628D",honest="#33816A")
method_labels<-c(naive="Naive",proposed="Selective",honest="Honest")
cutprint<-function(x)format(x,digits=12,scientific=FALSE,trim=TRUE)
label<-function(r,compact=FALSE)paste0((if(compact)shortvariable else variable)[r$variable]," <= ",cutprint(r$threshold))
texrule<-function(r)paste0("$",texvariable[r$variable],"\\leq ",cutprint(r$threshold),"$")
fmt<-function(x)ifelse(is.finite(x),sprintf("%.2f",x),"---")
pfmt<-function(x,math=TRUE) {
 if(!is.finite(x))return("---")
 if(x<.001)return(if(math)"$<0.001$" else "<0.001")
 if(math)sprintf("$%.3f$",x) else sprintf("%.3f",x)
}
rowline<-function(x)paste0(paste(x,collapse=" & ")," \\\\")
textp<-function(x)if(x<.001)"$p<0.001$" else sprintf("$p=%.3f$",x)
one<-function(z){stopifnot(nrow(z)==1L);z}

# Publish every main interval, while avoiding duplicated zero-penalty curves.
rows<-character()
append_row<-function(r,method,repeat_rule=TRUE) {
  nd<-if(r$node=="root")"Root" else r$node
  if(isTRUE(r$available)&&r$inference_status=="ok")
   rowline(c(if(repeat_rule)nd else "",if(repeat_rule)texrule(r)else "",method_labels[method],fmt(r$estimate),
     paste0("[",fmt(r$lower),", ",fmt(r$upper),"]"),pfmt(r$p_value),fmt(r$width)))
  else rowline(c(nd,texrule(r),method_labels[method],"---","Unavailable","---","---"))
}
for(i in seq_len(nrow(full))) {
 nm<-full$node[i]
 rows<-c(rows,append_row(one(main[main$method=="naive"&main$node==nm,]),"naive"),
          append_row(full[i,],"proposed",FALSE),"\\addlinespace[2pt]")
}
rows<-c(rows,"\\midrule")
for(i in seq_len(nrow(honest))) {
 rows<-c(rows,append_row(honest[i,],"honest"))
}
caption<-sprintf("ACTG 175 full-midpoint unpruned display at maximum depth three and $\\lambda=0$: %d full-data splits and %d Honest training splits. All intervals and $p$-values are individual nominal working-model calculations. Growth-only and Selective coincide at zero penalty and are displayed once as Selective. CD4 units are cells per cubic millimeter.",nrow(full),nrow(honest))
writeLines(c("\\begin{table}[tbp]\\centering\\footnotesize",paste0("\\caption{",caption,"}\\label{tab:actg-main}"),
 "\\setlength{\\tabcolsep}{3pt}","\\begin{tabular}{lllrlrr}\\toprule",
 rowline(c("Node","Rule","Method","Estimate","95\\% interval","$p$","Width")),"\\midrule",rows,
 "\\bottomrule\\end{tabular}",
 "\\par\\smallskip\\parbox{\\linewidth}{\\footnotesize Addresses identify each tree's own nodes. The full-data Naive and Selective results share targets, whereas Honest uses its own training-selected rules and holdout-design targets. Unavailable holdout contrasts are distinguished from numerical failure. Intervals are not simultaneous, and estimated variances do not establish exact coverage.}",
 "\\end{table}"),file.path(td,"main_results.tex"))

# Standard PDF fonts are used; no Ghostscript or PNG-generation dependency.
figure<-function(name,width,height,draw) {
 pdf(file.path(fd,paste0(name,".pdf")),width=width,height=height,
     pointsize=11,family="Helvetica",useDingbats=FALSE)
 tryCatch(draw(),error=function(e){dev.off();stop(e)});dev.off()
}
position<-function(name) {
 if(name=="root")return(c(.5,0));ss<-strsplit(name,"")[[1]]
 index<-sum(as.numeric(ss=="R")*2^rev(seq_along(ss)-1L));c((index+.5)/2^length(ss),-length(ss))
}
draw_tree<-function(which_tree,heading) {
 frame<-nodes[nodes$lambda==0&nodes$tree==which_tree,,drop=FALSE];kept<-frame[frame$retained,,drop=FALSE]
 stopifnot(nrow(kept)>0);col<-colors[if(which_tree=="full")"proposed" else "honest"]
 maxdepth<-max(kept$depth)+1L
 plot(0,0,type="n",xlim=c(-.025,1.025),ylim=c(-maxdepth-.4,.38),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
 title(sprintf("%s: %d internal splits",heading,nrow(kept)),adj=0,line=.65,cex.main=1.10)
 terminals<-list()
 for(i in seq_len(nrow(kept))) {
  nd<-kept[i,];p<-position(nd$node)
  for(side in c("L","R")) {
   child<-paste0(if(nd$node=="root")"" else nd$node,side);q<-position(child)
   segments(p[1],p[2]-.25,q[1],q[2]+.22,col="#8B939B",lwd=.8)
   text((p[1]+q[1])/2,(p[2]+q[2])/2+.015,if(side=="L")"Yes" else "No",cex=.92,col="#555B62")
   if(!child%in%kept$node)terminals[[child]]<-c(if(side=="L")nd$n0left else nd$n0right,if(side=="L")nd$n1left else nd$n1right)
  }
 }
 for(i in seq_len(nrow(kept))) {
  nd<-kept[i,];p<-position(nd$node)
  txt<-sprintf("%s\n%s\nn0=%d, n1=%d",if(nd$node=="root")"Root" else nd$node,label(nd,TRUE),nd$n0left+nd$n0right,nd$n1left+nd$n1right)
  # At the manuscript's approximately 6.2-inch display width, 0.92 x 11 pt
  # remains above 8 pt after scaling the 7.5-inch vector source.
  maxwidth<-if(nd$depth>=2).238 else .39;cex<-.98
  while(max(strwidth(strsplit(txt,"\n")[[1]],cex=cex))>maxwidth-.02&&cex>.92)cex<-cex-.02
  w<-max(strwidth(strsplit(txt,"\n")[[1]],cex=cex))+.025
  rect(p[1]-w/2,p[2]-.38,p[1]+w/2,p[2]+.38,border=col,col="white",lwd=1.1);text(p[1],p[2],txt,cex=cex)
 }
 for(nm in names(terminals)) {
  p<-position(nm);cc<-terminals[[nm]];txt<-sprintf("n0=%d\nn1=%d",cc[1],cc[2]);cex<-.94
  w<-max(strwidth(strsplit(txt,"\n")[[1]],cex=cex))+.025
  rect(p[1]-w/2,p[2]-.25,p[1]+w/2,p[2]+.25,col="#F0F2F4",border="#AEB6BE",lwd=.8);text(p[1],p[2],txt,cex=cex)
 }
}
figure("main_tree",7.5,7.1,function(){
 layout(matrix(1:3,3,1),heights=c(3.35,3.35,.4));par(cex=1,mar=c(.25,.3,2,.3),oma=c(.1,.1,.1,.1),xpd=NA)
 draw_tree("full","(a) Full sample; no pruning");draw_tree("Honest_training","(b) Honest training; no pruning")
 par(mar=rep(0,4));plot.new();text(.5,.55,"n0: ZDV+zal; n1: ZDV+ddI. Honest counts are training counts.",cex=.94)
})
shown_all<-main[main$available&main$inference_status=="ok",,drop=FALSE]
shared_lim<-range(c(0,shown_all$lower,shown_all$upper,shown_all$estimate),finite=TRUE)
shared_lim<-shared_lim+c(-1,1)*.045*diff(shared_lim)
shared_ticks<-pretty(shared_lim)
figure("main_intervals",7.8,7.4,function(){
 layout(matrix(1:2,2,1),heights=c(max(3.5,nrow(full)*.55),max(3.5,nrow(honest)*.55)))
 par(cex=1,mar=c(3.7,8.6,2.2,.7),oma=c(.2,0,0,0),mgp=c(2.15,.65,0),las=1,tcl=-.22)
 for(isfull in c(TRUE,FALSE)) {
  dat<-if(isfull)full else honest;ms<-if(isfull)c("naive","proposed") else "honest"
  lim<-shared_lim
  yy<-rev(seq_len(nrow(dat)))
  plot(0,0,type="n",xlim=lim,ylim=c(.4,nrow(dat)+.6),axes=FALSE,xlab=expression("GLS interaction contrast (CD4 cells/mm"^3*")"),ylab="",yaxs="i")
  abline(v=shared_ticks,col="#ECEFF1",lwd=.6);abline(v=0,col="#969EA5",lty=2,lwd=.8);axis(1,at=shared_ticks,cex.axis=.96);box(bty="l",col="#8F969D")
  axis(2,at=yy,labels=vapply(seq_len(nrow(dat)),function(i)paste0(if(dat$node[i]=="root")"Root"else dat$node[i],": ",label(dat[i,],TRUE)),""),tick=FALSE,cex.axis=.96)
  title(if(isfull)sprintf("(a) Full sample: %d shared targets",nrow(full))else sprintf("(b) Honest: %d training-selected rules",nrow(honest)),adj=0,line=.65,cex.main=1.10)
  for(i in seq_len(nrow(dat)))for(j in seq_along(ms)) {
   m<-ms[j];r<-one(main[main$node==dat$node[i]&main$method==m,,drop=FALSE]);y<-yy[i]+if(isfull)c(.14,-.14)[j] else 0
   if(isTRUE(r$available)&&r$inference_status=="ok") {
    segments(r$lower,y,r$upper,y,col=colors[m],lwd=1.5);segments(c(r$lower,r$upper),y-.04,c(r$lower,r$upper),y+.04,col=colors[m]);points(r$estimate,y,col=colors[m],pch=19,cex=.55)
   } else text(mean(lim),y,"Holdout unavailable",col="#A33E3E",cex=.96)
  }
  legend("bottomleft",method_labels[ms],col=colors[ms],lty=1,pch=19,pt.cex=.6,bty="n",horiz=TRUE,cex=.98,inset=c(.005,-.03),xpd=NA)
 }
})
cat("ACTG assets: results/actg/main_results.tex and figures/actg/main_tree.pdf, main_intervals.pdf.\n")
