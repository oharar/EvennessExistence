# This fits each distribution to the data, and then 
ParametricBootstraps <- function(n, nsim=5) {
  Sims <- list(
    CEGS = SimulateSpRich(n, model="CEGS", nsim=nsim),
    NB = SimulateSpRich(n, model="NB", nsim=nsim),
    PLN = SimulateSpRich(n, model="PLN", nsim=nsim)
  )
  
  Res.l <- lapply(Sims, function(mat) {
    Bias <- apply(mat[grep("Diff", rownames(mat)),], 1, function(x) mean(x, na.rm=TRUE))
    logSR <- apply(mat[grep("richness", rownames(mat)),], 1, function(x) mean(x[x<10e5], na.rm=TRUE))
    SimSR <- mat["SimSR",1]
    BiaslogSR <- logSR - log(SimSR)
    c(Bias, logSR, SimSR)
  })
  
  Res <- unlist(Res.l)
  Res
}
