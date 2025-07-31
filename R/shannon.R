shannon<-function(n,exponentiate=TRUE)	{
	f <- n / sum(n)
	ifelse(exponentiate, exp(-sum(f * log(f))), -sum(f * log(f)) )
}
