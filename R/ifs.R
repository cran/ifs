### file ifs/R/ifs.R
### copyright (C) 2001-2 S. M. Iacus



# The IFS distribution function estimator
#
# x : where to estimate the distribution function
# y : observation vector
# k : number of iteration, default = 5
# q : proportion of quantiles to use in the
#     definition of the functional
# f : the starting point in the space of
#     of continuous distribution
#     functions on [0,1]. If not specified
#     the uniform distribution is used.
# n : number of points of F to calculate
#
# This functions calls ifs() or ifs.flex() where it 
# is really defined the IFS estimator


IFS <- function (y, k = 5, q = 0.5, f = NULL, n = 512, maps = c("quantile", "wl1", "wl2")) 
{
    if ((q <= 0) && (q >= 1)) {
        q = 0.5
        warning("proportion q set to 0.5")
    }
    parms <- ifs.setup(y, maps)
    nm <- parms$n
    if (maps != "quantile") 
        nm <- 10

    s <- parms$s[1:nm]
    a <- parms$a[1:nm]
    p <- rep(1/length(a), length(a))

    if (maps != "quantile") {
     QF <- setQF(parms$m, parms$s, parms$a, nm)
     p <- ifs.cf(QF$Q, QF$b)
    }

    if (is.function(f)) {
        g <- function(x) ifs.flex(x, p, s, a, k, f)
    }
    else {
        g <- function(x) ifs(x, p, s, a, k)
    }
    xx <- seq(0, 1, length = n)
    yy <- apply(matrix(xx), 1, g)
    return(list(x = xx, y = yy))
}

# The IFS estimator engine
# x  : where to estimate the distribution function
# p  : the coefficients vector
# s  : the vector of coefficients s_i in w_i = s_i * x + a_i
# a  : the vector of coefficients a_i in w_i = s_i * x + a_i
# k  : number of iterations, default = 5
#

ifs <- function(x, p, s, a, k = 5)
{
   val <- .Call("ifs_df",x,p,s,a,as.integer(k))
   val
}

# The IFS estimator engine (flexible version)
# x  : where to estimate the distribution function
# p  : the coefficients vector
# s  : the vector of coefficients s_i in w_i = s_i * x + a_i
# a  : the vector of coefficients a_i in w_i = s_i * x + a_i
# k  : number of iterations, default = 5
# f  : the starting function

ifs.flex <- function(x, p, s, a, k = 5, f)
{
   val <- .Call("ifs_df_flex",x,p,s,a,as.integer(k),f,new.env())
   val
}

# Used to build the quadratic form in terms of the moments
#  and the coefficients of the affine maps.
# m  : the vector fo the moments
# s  : the vector of coefficients s_i in w_i = s_i * x + a_i
# a  : the vector of coefficients a_i in w_i = s_i * x + a_i
# n  : number of coefficients to be used in the IFS

setQF <- function(m, s, a, n = 10)
{

 if( length(m) - 1 < n)
  stop("`n' too big")

 if( length(s) != length(a) )
  stop("`s' and `a' must have the same length")
  
 if( length(s) != (length(m)-1) )
  stop("`m', `s' and `a' not of conformable length")

 return( .Call("ifs_setQF",m,s,a,as.integer(n)) )

}

# Calculate the coefficients of the IFS given the quadratic form

ifs.cf <- function(Q, b){

    ln <- dim(Q)[1]
    w <- rep(0,ln)

    fr <- function(x) {
        w[1] <- 1 - sum(x)
        w[2:ln] <- x[1:(ln-1)]        
        return(t(w) %*% Q %*% w + b %*% w)
    }
    sol <- optim(rep(1/(2 * ln), ln-1), fr, gr = NULL, 
        method = "L-BFGS-B", lower = rep(0, ln-1), upper = rep(1,ln-1))

    p <- c(1-sum(sol$par),sol$par)

    return(p)
}

# y : the vector of observations
# returns the coefficients a and s of w = s*x+a and the moments' vector m

# y : the vector of observations
# returns the coefficients a and s of w = s*x+a and the moments' vector m

ifs.setup <-function (y,maps=c("quantile","wl1","wl2"),qtl) 
{
   s <- NULL
   a <- NULL 

   if( maps == "quantile" ) {
    if(missing(qtl))
     qt <- as.numeric(c(0, quantile(y, probs = seq(0, 1, 1/(0.5 * 
        length(y)))), 1))
    else
     qt <- qtl
    np <- length(qt) - 1
    for (i in 1:np) 
     s <- c(s, qt[i + 1] - qt[i])
    a <- qt[1:np]
    }
    
    if( maps == "wl1") {
     M <- 4
     np <- sum(2^(1:M))
     for(i in 1:M)
      for(j in 1:(2^i)){
        s <- c(s,1/(2^i))
        a <- c(a,(j-1)/(2^i))
      }
    }
    
    if( maps == "wl2") {
     M <- 4
     np <- M*(m-1)/2
     for(i in 2:M)
      for(j in 2:i){
        s <- c(s,1/i)
        a <- c(a,(j-1)/i)
      }
    }
  
     
  m <- 1
  for (i in 1:np) 
   m <- c(m, mean(y^i))

  return(list(m = m, s = s, a = a, n = np))
}

# Fourier distribution function estimation
#
ifs.pf.FT <- function(x,b,nterms){
 val <- x/(2*pi)
 if(missing(nterms) | (nterms > length(b)) ) 
  nterms <- length(b)

 for(k in 1:nterms){
  val <- val + (Re(b[k])*sin(k*x)/k + Im(b[k])*(cos(k*x)-1)/k)/pi
 }
 ifelse(val<0,0,ifelse(val>1,1,val))
}

# Fourier density function estimation
#
ifs.df.FT <- function(x,b,nterms){
 val <- 1/(2*pi)
 if(missing(nterms) | (nterms > length(b)) ) 
 nterms <- length(b)

 for(k in 1:nterms){
  val <- val + (Re(b[k])*cos(k*x) -Im(b[k])*sin(k*x))/pi
 }
 ifelse(val<0,0,val)
}

# x : where to calculate the FT
ifs.FT <- function(x, p, s, a, k = 2){
 .Call("ifs_ft",x,p,s,a,as.integer(k))
}

# m : number of Fourier coefficients to estimate
# k : number of iterations
# p, s, a as in ifs
# cutoff = sqrt(2/(n+1)) where n is the sample size

ifs.setup.FT <- function(m, p, s, a, k = 2, cutoff){

 b <- NULL
 for(i in 1:m)
  b <- c(b, .Call("ifs_ft",i,p,s,a,as.integer(k)))

 if(missing(cutoff))
  nterms <- length(b)
 else{
  nterms <-0
  nfail <- 0
  for( k in 1:(length(b)) ) {
   if( abs(b[k]) > cutoff ) {
    nterms <- k;  nfail <- 0 }
   else  nfail <- nfail + 1
 
  if(nfail > 2)
   break;
  }
 }
 
 return(list(b = b, nterms = nterms))

}



IFS.pf.FT <- function (y, k = 2, n = 512, maps=c("quantile","wl1","wl2")) 
{
    
    parms <- ifs.setup(y, maps)

    a <- parms$a
    s <- parms$s
    p <- rep(1/length(a),length(a))
 
    nobs <- length(y)
    cutoff <- sqrt(2/(nobs+1))

    coeff <- ifs.setup.FT(20,p,s,a,k,cutoff) 

    b <- coeff$b
    nterms <- coeff$nterms

    
    xx <- seq(0, 1, length = n)
    yy <- ifs.pf.FT(xx,b,nterms)
     return(list(x = xx, y = yy, b=b, nterms=nterms))
}

IFS.df.FT <- function (y, k = 2, n = 512, maps=c("quantile","wl1","wl2")) 
{
    
    parms <- ifs.setup(y, maps)

    a <- parms$a
    s <- parms$s
    p <- rep(1/length(a),length(a))
 
    nobs <- length(y)
    cutoff <- sqrt(2/(nobs+1))

    coeff <- ifs.setup.FT(20,p,s,a,k,cutoff) 

    b <- coeff$b
    nterms <- coeff$nterms

    
    xx <- seq(0, 1, length = n)
    yy <- ifs.df.FT(xx,b,nterms)
    return(list(x = xx, y = yy, b=b, nterms=nterms))
}
