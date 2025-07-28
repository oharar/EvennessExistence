shannon<-function(n,exponentiate=TRUE)	{
	f <- n / sum(n)
	if (exponentiate == TRUE)
		return(exp(-sum(f * log(f))))
	return(-sum(f * log(f)))
}
