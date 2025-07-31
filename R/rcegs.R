rcegs<-function(n=100,l=1,g=2)	{
	rgeom(n,1 / (rexp(n)^g / l + 1))
}

dcegs.val <- function(x, l, g, log=FALSE) {
  cegs.pdf <- function(E, x, l, g) {
    if(l<0) stop("l must be positive")
    if(g<=0) stop("g must be positive")
    p <- 1/(E^g/l + 1)
    dgeom(x,p)*dexp(E)
  }
  dens <- integrate(cegs.pdf, x=x, l=l,g=g,lower=1e-20,upper=100, 
                    stop.on.error=FALSE)$value
  if(log) dens <- log(dens)
  dens
}

dcegs <- function(x, l, g, log=FALSE) {
  if(length(x)==1) {
    res <- dcegs.val(x=x, l=l, g=g, log=log)
  } else {
    res <- sapply(x, dcegs.val, l=l, g=g, log=log)
  }
  res
}

# dcegs(x=1:4, l=l,g=g, log=TRUE)
