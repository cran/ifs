library(ifs)
load("gayser.rda")
y<-faithful$eruptions

 
 y <- (y-min(y))/(max(y)-min(y))

 qtl <- as.numeric(quantile(y,probs=seq(0,1,length=40)))
 parms <- ifs.setup(y,qtl=qtl)
    a <- parms$a
    s <- parms$s
    p <- rep(1/length(a), length(a))
    nobs <- length(y)
    cutoff <- sqrt(2/(nobs + 1))
    coeff <- ifs.setup.FT(30, p, s, a, 2, cutoff)
    b <- coeff$b
    nterms <- coeff$nterms
    den <- density(y,from=0,to=1,bw=0.03)
    xx <- den$x
    yy <- ifs.df.FT(xx, b, nterms)
    pdf("gayser.pdf")
   # 
   plot(den,xlim=c(0,1),ylim=c(0,2.5),main="Kernel vs IFS",sub="",xlab="x",ylab="",lty=3)
    lines(xx,yy,col="red")

  dev.off()


