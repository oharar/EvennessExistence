fisher<-function(n)	{
	if (max(n) == 1)
		return(NA)
	S <- length(n)
	N <- sum(n)
	a <- 1
	lasta <- 0
	z <- 0
	while (abs(a - lasta) > 0.0000001 && z < 1000)	{
		z <- z + 1
		lasta <- a
		a <- S / log(1 + N / a)
	}
	return(a)
}
