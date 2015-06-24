test_kbatch <- function(){
  library(oligoClasses)
  set.seed(123)
  k <- 3
  means <- matrix(c(rnorm(5, -1, 0.1),
                    rnorm(5, 0, 0.1),
                    rnorm(5, 1, 0.1)), 5, k, byrow=FALSE)
  sds <- matrix(0.1, 5, k)
  ##N <- 1500
  N <- 9000
  probs <- c(1/3, 1/3, 3/10, 0.02, 0.013)
  probs <- probs/sum(probs)
  batch <- sample(1:5, size=N, prob=probs, replace=TRUE)

  p <- c(1/5, 1/3)
  p <- c(p, 1-sum(p))
  truth <- simulateBatchData(N=N,
                             batch=batch,
                             theta=means,
                             sds=sds,
                             p=p)
  mp <- McmcParams(iter=1000, burnin=250, nStarts=5)
  kmod <- BatchModel(y(truth), batch(truth), k=3, mcmc.params=mp)
  kmod <- posteriorSimulation(kmod)
  cn <- map(kmod)

  set.seed(1000)
  index <- sample(seq_len(N), 3000)
  kmod2 <- BatchModel(y(truth)[index], batch(truth)[index], k=3, mcmc.params=mp)
  kmod2 <- posteriorSimulation(kmod2)
  yy <- setNames(y(truth), seq_along(y(truth)))
  df <- CNPBayes:::imputeFromSampledData(kmod2, yy, index)
  cn2 <- df$cn
  mean(cn != cn2)

  cn2 <- map(kmod2)
  pz <- probz(kmod2)
  pz <- mapCnProbability(kmod2)





  if(FALSE){
    ##
    ## With 10k iterations and 1k iterations burnin, all models
    ## provide consistent estimates of the marginal likelihood and the
    ## K=4 model is best.  This unit test is too time consuming to be
    ## run during each package update.
    ##
    fit <- computeMarginalLik(y(truth), batch(truth), K=1:4,
                              burnin=1000,
                              T2=1000, T=10000,
                              nchains=3)
    prz <- probz(fit$models[[4]])
    cn <- map(fit$models[[4]])
    plot(r, cn, pch=20, cex=0.3)
    ## The 4th state is the map estimate in 2 individuals, so even if
    ## K = 4 has a higher marginal likelihood the copy number
    ## inference is not effected.
    trace(cnProbability, browser)
    prz <- cnProbability(prz, 4)
    plot(jitter(prz, amount=0.05), jitter(cn, amount=0.05), pch=20, cex=0.3)
    table(cn)

    pz <- cnProbability(probz(fit$models[[4]]), 4)
    r <- y(fit$models[[4]])
    plot(r, pz, pch=".")
    checkTrue(k(orderModels(fit))[1] == 3)
  }
}
