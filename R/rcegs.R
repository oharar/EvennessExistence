rcegs<-function(n=100,l=1,g=2)	{
	rgeom(n,1 / (rexp(n)^g / l + 1))
}
