simpson<-function(n)	{
	N <- sum(n)
	S <- sum(n * (n - 1) / (N * (N - 1)))
	if (is.na(S) || S == 0)
		return(NA)
	return(1 / S)
}
