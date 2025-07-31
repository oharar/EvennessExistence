CalcParSR <- function(n, trim=TRUE, removeDom=FALSE) {
  if(removeDom & length(unique(n))>1) n <- n[n!=max(n)]
  if(trim) n <- trim(n)
  
  # These are in case anything crashes
  ln <- list(p=as.numeric(NA), par=as.numeric(c(NA,NA)), logLval=as.numeric(NA))
  eg <- list(richness=as.numeric(NA), shape=as.numeric(NA))
  nb <- list(richness=as.numeric(NA), size=as.numeric(NA), 
             prob=as.numeric(NA), loglik=as.numeric(NA))
  
  try(ln <- poilog::poilogMLE(n), silent = TRUE)
  try(nb <- FitNB(n), silent = TRUE)
  try(eg <- cegsML(n, fitted = FALSE), silent = TRUE)
  cegs_loglik <- ifelse(is.na(eg$scale), NA, sum(dcegs(x=n, l=eg$scale,g=eg$shape, log=TRUE)))
  
  res <- c(
    pln_richness = length(n)/ln$p, pln_mu = ln$par[1], pln_sigma = ln$par[2], pln_loglik = ln$logLval,
    cegs_richness = eg$richness, cegs_gamma = eg$shape, cegs_lambda = eg$scale, cegs_loglik = cegs_loglik, 
    nb_richness = nb$richness, nb_size=nb$size, nb_loglik = nb$loglik,
    nb_prob=nb$prob
  )
  res
}
