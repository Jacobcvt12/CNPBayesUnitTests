test_senseless_batch <- function(){
  set.seed(2000)
  library(oligoClasses)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.1, 0.1),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  m1 <- computeMarginalLik(y(truth), nchains=5,
                           burnin=300,
                           T2=500, T=1000)
  models <- orderModels(m1)
  ## might select k=2 due to variability of marginal lik estimates
  checkTrue(k(models)[1] == 3)
  if(FALSE){
    ## Make up a bogus batch
    m2 <- computeMarginalLik(y(truth), batch=rep(1:3, length.out=2500),
                             nchains=3, burnin=300, T2=300, T=500)
    ## we get the right model even when batch is independent of the data
    checkTrue(is(m2$models, "BatchModelList"))
    checkTrue(k(orderModels(m2))[1]==3)
  }
}
