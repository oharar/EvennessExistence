sadrad<-function(S,p)	{
	p <- p / sum(p)
	r <- array()
	q <- 1:S / (S + 1)
	cs <- cumsum(p)
	w <- 1
	for (i in 1:S)	{
		for (j in w:length(p))
			if (cs[j] > q[i])
				break
		w <- j
		if (w == 1)
			r[i] <- 1
		else
			r[i] <- w
	}
	return(r)
}
