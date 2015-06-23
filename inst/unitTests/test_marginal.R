test_selectK_easy <- function(){
  library(GenomicRanges)
  set.seed(1000)
  truth <- simulateData(N=2500, p=c(1/4, 1/2, 1-1/2-1/4), theta=means, sds=sds)
  x2 <- computeMarginalLik(y(truth), nchains=3, K=1:4, T=5000, T2=1000, burnin=1000)
  m2 <- orderModels(x2)
  checkTrue(k(m2)[1] >= 3)
}
