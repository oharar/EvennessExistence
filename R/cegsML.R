cegsML<-function(n, fitted=TRUE)	{
	if (length(table(n)) < 3)
		return(list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA))

	S <- length(n)
	s <- array(dim=max(n),data=0)
	t <- table(n)
	s[as.numeric(names(t))] <- t
	u <- which(s > 0)
	px<-function(U,l,g,i)	{
		p <- 1 / ((-log(U))^g / l + 1)
		p[p == 0] <- 1e-8
		dgeom(i,p)
	}
	p0<-function(U,l,g)	{
		1 / ((-log(U))^g / l + 1)
	}
	like<-function(l,g)	{
		g <- exp(g)
		if (l <= 0 || g <= -10)
			return(1e10)
		p <- array()
		for (i in 1:length(u))
			p[i] <- integrate(px,l=l,g=g,i=u[i],lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		if (is.nan(p[1]) || min(p) == 0)
			return(1e10)
		p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
		p <- p / (1 - p0)
		ll <- -sum(s[u] * log(p))
		if (is.infinite(ll) || is.nan(ll))
			return(1e10)
		ll
	}
#	cf <- try(stats4::mle(like,lower=list(l=0,g=-10),upper=list(l=1e8,g=10),start=list(l=1,g=2)), silent=T)
	cf <- tryCatch(stats4::mle(like,lower=list(l=0,g=-10),upper=list(l=1e8,g=10),start=list(l=1,g=2)), 
	               error=function(e) NULL, silent=TRUE)
	if(!is.null(cf)) {
	  l <- coef(cf)[1]
	  g <- coef(cf)[2]
	} else {
	  l <- NA
	  g <- NA
	  
	}
	if (is.na(l) || l == 0 || l == 1e8 || g == -10 || g == 10) {
	  res <- list('richness' = NA, 'scale' = NA, 'shape' = NA, 'AICc' = NA, 'fitted.RAD' = NA, 'fitted.SAD' = NA)
	} else {
	  aicc <- 2 * like(l,g) + 4 + 12 / (S - 3)
	  g <- exp(g)
	  p <- array()
	  mx <- max(2^12,2^ceiling(log2(max(n))))
	  for (i in 1:mx)
	    p[i] <- integrate(px,l=l,g=g,i=i,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	  p0 <- integrate(p0,l=l,g=g,lower=1e-20,upper=1 - 1e-20,stop.on.error=F)$value
	  p <- p / (1 - p0)
	  
	  res <- list('richness' = as.numeric(S / (1 - p0)), 'scale' = as.numeric(l), 
	              'shape' = as.numeric(g), 'AICc' = aicc)
	  if(fitted) {
	    res$fitted.RAD <-  sadrad(length(n),p)
	    res$fitted.SAD <- p[1:2^12]
	  }
	  
	}
	return(res)
}
