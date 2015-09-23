test_that("senseless_batch"{
  set.seed(2000)
  library(oligoClasses)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.1, 0.1),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  mp <- McmcParams(iter=1000, burnin=300, nStarts=5)
  model3 <- MarginalModel(y(truth), mcmc.params=mp, k=3)
  model3 <- posteriorSimulation(model3)
  model2 <- MarginalModel(y(truth), mcmc.params=mp, k=2)
  model2 <- posteriorSimulation(model2)
  mlist <- list(model2, model3)
  x <- marginalLikelihood(mlist)
  models <- mlist[order(x, decreasing=TRUE)]
  ## might select k=2 due to variability of marginal lik estimates
  expect_true(k(models[[1]]) == 3)
  if(FALSE){
    ## Make up a bogus batch
    m2 <- computeMarginalLik(y(truth), batch=rep(1:3, length.out=2500),
                             nchains=3, burnin=300, T2=300, T=500)
    ## we get the right model even when batch is independent of the data
    checkTrue(is(m2$models, "BatchModelList"))
    checkTrue(k(orderModels(m2))[1]==3)
  }
}
