# Neg bin
# would like to try to use the VGAM package, but the package has a bug, and the update
# doesn't work on a Mac. In theory VGAM::vglm(n ~ 1, posnegbinomial) should work

FitNB <- function(n) {
  fit <- MASS::fitdistr(trim(n), densfun=actuar::dztnbinom,
                        start=list(size=0.2, prob=0.5), lower=c(5e-3,1e-5))  
  p0 <- dnbinom(0, size=fit$estimate["size"], prob=fit$estimate["prob"])
  
  k <- length(fit$estimate)
  return(list('richness' = length(n) / (1 - p0), 
              'size' = fit$estimate["size"], 
              'prob' = fit$estimate["prob"], 
              'loglik' = fit$loglik,
              'AICc' = -2*fit$loglik + (2*k^2 + 2*k)/(length(n)-k-1)))
              # 'fitted.RAD' = sadrad(length(n),p), 
              # 'fitted.SAD' = p[1:2^12]))
}


# this isn't working! If it did, we could use it to get CIs for the species richness estimates
FitNB2 <- function(n) {
  GetNB <- function(nzero, counts, lhood=TRUE) {
    tab <- table(counts)
    nnn <- c(0,as.numeric(names(tab)))
    wt <- c(nzero, tab)
    
    fit <- MASS::glm.nb(nnn~1, weights = wt)
    if(lhood) {
      # Likelihood for n
      res <- sum(wt[-1]*dnbinom(nnn[-1], size=fit$theta, mu=exp(coef(fit)), log=TRUE))
    } else {
      res <- fit
    }
    res
  }
  
  opt <- optimize(GetNB, interval=c(1, 1e6), counts=n, maximum = TRUE)
  
  tab <- table(n)
  nnn <- c(0,as.numeric(names(tab)))
  wt <- c(round(opt$maximum), tab)
  
  fit <- MASS::glm.nb(nnn~1, weights = wt)
  lhood <- sum(wt[-1]*dnbinom(nnn[-1], size=fit$theta, mu=exp(coef(fit)), log=TRUE))
  k <- length(fit$estimate) + 1
  return(list('richness' = round(opt$maximum), 
              'mu' = coef(fit)["(Intercept)"], 
              'theta' = fit$theta, 
              'loglik' = lhood, 
              'AICc' = -2*lhood + (2*k^2 + 2*k)/(length(n)-k-1)))
  # 'fitted.RAD' = sadrad(length(n),p), 
  # 'fitted.SAD' = p[1:2^12]))
}
