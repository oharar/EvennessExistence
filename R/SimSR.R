SimPLN <- function(S, mu, sigma) {
  loglambda <- rnorm(S, mu, sigma)
  N <- rpois(length(loglambda), exp(loglambda))
  N
}

SimCEGS <- function(S, gamma, lambda) {
  p <- 1/((rexp(S, 1)^gamma)/lambda+1)
  N <- rgeom(length(p), p)
  N
}
SimNB <- function(S, size, prob) {
  lambda <- rgamma(S, shape=size, prob/(1-prob))
  N <- rpois(length(lambda), lambda)
  N
}

SimSR <- function(pars, model=c("PLN", "CEGS", "NB")) {
  # Check for NAs in pars
  AnyNAs  <- switch(model,
                    PLN = any(is.na(pars[c("pln_richness", "pln_mu.mu", "pln_sigma.sig")])),
                    CEGS = any(is.na(pars[c("cegs_richness", "cegs_gamma", "cegs_lambda")])),
                    NB = any(is.na(pars[c("nb_richness", "nb_size.size", "nb_prob.prob")])))
  if(AnyNAs) {
    snames <- switch(model,
                     PLN = c("pln_mu.mu", "pln_sigma.sig"),
                     CEGS = c("cegs_gamma", "cegs_lambda"),
                     NB = c("nb_size.size", "nb_prob.prob"))
    
    Names <- c(paste0("Diff.", snames), "SimSR", "pln_richness", "cegs_richness", "nb_richness", snames)
    res <- rep(NA, times=length(Names)); names(res) <- Names
  } else {
    # Fit models
    # simulate from model
    n.sim  <- switch(model,
                     PLN = SimPLN(S=round(pars["pln_richness"]), mu=pars["pln_mu.mu"], 
                                  sigma=pars["pln_sigma.sig"]),
                     CEGS = SimCEGS(S=round(pars["cegs_richness"]), gamma=pars["cegs_gamma"], 
                                    lambda=pars["cegs_lambda"]),
                     NB = SimNB(S=round(pars["nb_richness"]), size=pars["nb_size.size"], 
                                prob=pars["nb_prob.prob"])
    )
    simpars <- CalcParSR(n.sim[n.sim>0], trim=FALSE)
    params <- switch(model,
                     PLN = pars[c("pln_mu.mu", "pln_sigma.sig")],
                     CEGS = pars[c("cegs_gamma", "cegs_lambda")],
                     NB = pars[c("nb_size.size", "nb_prob.prob")])
    simparams <- switch(model,
                        PLN = simpars[c("pln_mu.mu", "pln_sigma.sig")],
                        CEGS = simpars[c("cegs_gamma", "cegs_lambda")],
                        NB = simpars[c("nb_size.size", "nb_prob.prob")])
    #  c(Sim=length(n.sim), pars[c("pln_richness", "cegs_richness", "nb_richness")], params)
    Diff <- simparams-params; names(Diff) <- paste0("Diff.", names(Diff))
    res <- c(Diff, SimSR=length(n.sim), simpars[c("pln_richness", "cegs_richness", "nb_richness")], simparams)
  }
  res
}


SimulateSpRich <- function(n, model="CEGS", nsims=10) {
  pars <- CalcParSR(n[n>0], trim=FALSE)
  # AnyNAs  <- switch(model,
  #                   PLN = any(is.na(pars[grep("pln_", names(pars))])),
  #                   CEGS = any(is.na(pars[grep("cegs_", names(pars))])),
  #                   NB = any(is.na(pars[grep("nb_", names(pars))])))
  # if(AnyNAs) {
  #   res <- matrix(rep(NA, times=8*nsims), ncol=nsims)
  # } else {
    res <- replicate(nsims, SimSR(pars, model=model))
#  }
  res
}



