CalcBasicStats <- function(n, removeDom=FALSE, vegan=FALSE) {
  if(removeDom & length(unique(n))>1) n <- n[n!=max(n)]
  
  if(vegan) {
    require(vegan)
    res <- c(N = sum(n), S = length(n), alpha = vegan::fisher.alpha(n), 
             H = vegan::diversity(n, index="shannon"), 
             D = 1/(1-vegan::simpson.unb(n)), 
             UniqueCounts = length(unique(n)))
  } else {
    res <- c(N = sum(n), S = length(n), alpha = fisher(n), H = shannon(n), 
             D = simpson(n), UniqueCounts = length(unique(n)))
  }
  # compute Pielou's J & the log-transformed Hill ratio of H and D
  res <- c(res, J=log(res['H'])/log(res['S']), HD = log(res['H']/res['D']))
  names(res) <- gsub("\\.H", "", names(res))
  
  res
}
