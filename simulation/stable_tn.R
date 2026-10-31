# Stable TN log odds when the reference mean lies far beyond finite support.
# These identities change numerical evaluation only, not selection or inference.
if (!exists(".stable_tn_original_logodds", inherits=FALSE)) {
  .stable_tn_original_logodds <- tn_logodds
}
.stable_mills_coefficients <- c(1,-1,3,-15,105,-945,10395,-135135,
                               2027025,-34459425,654729075,-13749310575,
                               316234143225)
.stable_mills_series <- function(r) {
  p <- rep(tail(.stable_mills_coefficients,1L),length(r))
  for(j in rev(seq_len(length(.stable_mills_coefficients)-1L)))
    p <- .stable_mills_coefficients[j]+r*p
  p
}
.stable_mills_log <- function(z) {
  stopifnot(all(is.finite(z)),all(z>20))
  -log(z)+log(.stable_mills_series(1/z^2))
}
.stable_mills_log_ratio <- function(z, gap) {
  # log R(z+gap)-log R(z), R(z)=Phi_bar(z)/phi(z).
  # The stable divided difference avoids cancellation for extremely small gaps.
  rb <- 1/z^2
  lr <- log1p(gap/z)
  ra <- rb*exp(-2*lr)
  dr <- rb*expm1(-2*lr)
  pa <- rep(tail(.stable_mills_coefficients,1L),length(z))
  divided <- numeric(length(z))
  for(j in rev(seq_len(length(.stable_mills_coefficients)-1L))) {
    divided <- pa+rb*divided
    pa <- .stable_mills_coefficients[j]+ra*pa
  }
  -lr+log1p(dr*divided/.stable_mills_series(rb))
}
.stable_tilted_mass <- function(a,b,k,H) {
  ans <- rep(-Inf,length(a))
  ok <- a<b
  if(!any(ok))return(ans)
  aa <- a[ok];bb <- b[ok]
  zb <- k-bb
  zH <- k-H
  distance <- H-bb
  base <- -distance*(zH+distance/2)+.stable_mills_log(zb)
  finite_left <- is.finite(aa)
  adjustment <- numeric(length(aa))
  if(any(finite_left)) {
    gap <- bb[finite_left]-aa[finite_left]
    log_ratio <- -gap*(zb[finite_left]+gap/2)+
      .stable_mills_log_ratio(zb[finite_left],gap)
    adjustment[finite_left] <- log(-expm1(pmin(0,log_ratio)))
  }
  ans[ok] <- base+adjustment
  ans
}
.stable_positive_odds <- function(k,a,b) {
  H <- max(b)
  lower <- .stable_tilted_mass(a,pmin(0,b),k,H)
  upper <- .stable_tilted_mass(pmax(0,a),b,k,H)
  log_sum(lower)-log_sum(upper)
}
.stable_normal_log_mass <- function(lo,hi,mean,sd) {
  ans <- normal_log_mass(lo,hi,mean,sd)
  a <- (lo-mean)/sd;b <- (hi-mean)/sd
  right <- lo<hi & a>20 & is.finite(a)
  left <- lo<hi & b< -20 & is.finite(b)
  for(is_right in c(TRUE,FALSE)) {
    keep <- if(is_right)right else left
    if(!any(keep))next
    z <- if(is_right)a[keep] else -b[keep]
    gap <- (hi[keep]-lo[keep])/sd
    adjust <- numeric(length(z))
    finite <- is.finite(gap)
    if(any(finite)) {
      lr <- -gap[finite]*(z[finite]+gap[finite]/2)+
        .stable_mills_log_ratio(z[finite],gap[finite])
      adjust[finite] <- log(-expm1(pmin(0,lr)))
    }
    ans[keep] <- -.5*z^2-.5*log(2*pi)+.stable_mills_log(z)+adjust
  }
  # CDF subtraction also loses relative precision for very narrow intervals.
  # Centered five-point Gauss-Legendre integration has negligible remainder
  # under this width/density-variation bound and avoids subtracting CDFs.
  width <- (hi-lo)/sd
  narrow <- lo<hi & is.finite(width) & width<1e-3 &
    is.finite(a) & is.finite(b) & width*pmax(1,abs(a),abs(b))<1e-3
  nodes <- c(-.9061798459386639928,-.5384693101056830910,0,
              .5384693101056830910,.9061798459386639928)
  weights <- c(.2369268850561890875,.4786286704993664680,.5688888888888888889,
                .4786286704993664680,.2369268850561890875)
  for(i in which(narrow)) {
    half <- width[i]/2
    midpoint <- a[i]+half
    displacement <- half*nodes
    ans[i] <- -.5*midpoint^2-.5*log(2*pi)+log(width[i])+
      log_sum(log(weights/2)-midpoint*displacement-.5*displacement^2)
  }
  ans
}
tn_logodds <- function(t,mean,sd,support) {
  if(length(t)!=1L || !is.finite(t) || length(mean)!=1L || !is.finite(mean) ||
     length(sd)!=1L || !is.finite(sd) || sd<=0 || !is.matrix(support) ||
     ncol(support)!=2L || !nrow(support) || anyNA(support) ||
     any(support[,1L]>=support[,2L]))stop("Invalid truncated-normal arguments")
  k <- (mean-t)/sd
  a <- (support[,1L]-t)/sd
  b <- (support[,2L]-t)/sd
  high <- max(b);low <- min(a)
  if(is.finite(high) && is.finite(k) && k-high>20) {
    ans <- .stable_positive_odds(k,a,b)
  } else if(is.finite(low) && is.finite(k) && low-k>20) {
    ans <- -.stable_positive_odds(-k,-b,-a)
  } else {
    lower <- .stable_normal_log_mass(support[,1L],pmin(t,support[,2L]),mean,sd)
    upper <- .stable_normal_log_mass(pmax(t,support[,1L]),support[,2L],mean,sd)
    ans <- log_sum(lower)-log_sum(upper)
  }
  if(is.nan(ans))stop("Undefined truncated-normal probability")
  ans
}
tn_interval <- function(t,sd,support,alpha=.05) {
  if(length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("alpha must lie strictly between zero and one")
  solve_mean <- function(prob) {
    target <- qlogis(prob)
    f <- function(m)tn_logodds(t,m,sd,support)-target
    found <- FALSE
    for(k in 0:90) {
      radius <- sd*2^k
      lo <- t-radius;hi <- t+radius
      if(!all(is.finite(c(lo,hi))))break
      flo <- f(lo);fhi <- f(hi)
      if(is.finite(flo) && is.finite(fhi) && flo>=0 && fhi<=0) {
        found <- TRUE
        break
      }
    }
    if(!found)stop("Confidence endpoint could not be bracketed")
    uniroot(f,c(lo,hi),tol=1e-7*sd,maxiter=250L)$root
  }
  c(lower=solve_mean(1-alpha/2),upper=solve_mean(alpha/2))
}
