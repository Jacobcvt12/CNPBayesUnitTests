test_that("marginal_galaxy", {
  library(MASS)
  data(galaxies)
  set.seed(42)
  galaxies[78] <- 26960
  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100)
  mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=20)
  model <- MarginalModel(data=galaxies / 1000, k=3,
                         hypp=hypp,
                         mcmc.params=mp)
  model <- posteriorSimulation(model)
  if(FALSE){
    plot(model, breaks=50)
  }
  lik.marginal <- marginalLikelihood(model, 1000L)
  expect_equal(lik.marginal, -227, tolerance=1)
})

test_that("full_theta_pooled", {
  library(MASS)
  data(galaxies)
  set.seed(42)
  galaxies[78] <- 26960
  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100, mu.0=20)
  mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)
  model <- CNPBayes:::SingleBatchPooledVar(data=galaxies / 1000, k=3,
                                           hypp=hypp,
                                           mcmc.params=mp)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1000,
                                                burnin=1000, thin=10)
  model <- useModes(model)
  p.thetas.cpp <- CNPBayes:::full_theta_pooled(model)

  thetastar <- modes(model)[["theta"]]
  probs <- full_theta_pooled(model)
  Z <- z(chains(model))
  tau2c <- tau2(chains(model))
  muc <- mu(chains(model))
  sigma2 <- sigma2(chains(model))
  p_theta <- rep(NA, 500)
  for(i in 1:500){
    zz <- Z[i, ]
    z(model) <- zz ;
    nn <- as.integer(tablez(model))
    data_mean <- compute_means(model)
    tau2_tilde = 1/tau2c[i]
    sigma2_tilde = 1/sigma2[i, 1]
    prod <- 1
    for(k in 1:3){
      post_prec = tau2_tilde + sigma2_tilde * nn[k];
      tau_n = sqrt(1/post_prec);
      w1 <- tau2_tilde/post_prec
      w2 <- nn[k]*sigma2_tilde/post_prec
      mu_n <- w1*muc[i] + w2*data_mean[k];
      tmp = dnorm(thetastar, mu_n, tau_n) ;
      prod <- prod * tmp[k] ;
    }
    p_theta[i] <- prod
  }
  ## small values occur because of label switching (good mixing)
  expect_equal(p_theta, p.thetas.cpp[1:500])
})

test_that("reduced_sigma_pooled", {
  library(MASS)
  data(galaxies)
  set.seed(42)
  galaxies[78] <- 26960
  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100, mu.0=20)
  mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)
  model <- CNPBayes:::SingleBatchPooledVar(data=galaxies / 1000, k=3,
                                           hypp=hypp,
                                           mcmc.params=mp)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1000,
                                                burnin=1000, thin=10)
  model <- useModes(model)
  mp <- McmcParams(iter=1000, burnin=0)
  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp
  model.psigma2 <- reduced_sigma_pooled(model.reduced)
  psigma.star <- p_sigma_reduced_pooled(model.psigma2)

  Z <- z(chains(model.psigma2))
  s20chain <- sigma2.0(chains(model.psigma2))
  nu0chain <- nu.0(chains(model.psigma2))
  sigma2star <- modes(model)[["sigma2"]]
  prec = 1.0/sigma2star;
  p.prec <- rep(NA, 1000)
  thetastar <- modes(model.psigma2)[["theta"]]
  ##for (int s = 0; s < S; ++s) {
  for(s in 1:1000){
    zz = Z[s, ]
    nn = as.integer(table(zz));
    if(length(nn) !=3 ) stop()
    s20 = s20chain[s];
    nu0 = nu0chain[s];

    squared <- (y(model.psigma2)-thetastar[zz])^2
    ss <- sum(squared)
    nu.n <- 0.5*(nu0 + length(zz))
    sigma2_n <- 0.5*(nu0 * s20 + ss)
    p.prec[s] = dgamma(prec, nu.n, rate=sigma2_n)
  }
  expect_equal(psigma.star, p.prec)
})


test_that("marginalLikelihoodPooledVar", {
  library(MASS)
  data(galaxies)
  set.seed(42)
  galaxies[78] <- 26960
  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100, mu.0=20)
  mp <- McmcParams(iter=1000, burnin=2000, thin=10, nStarts=10)
  model <- CNPBayes:::SingleBatchPooledVar(data=galaxies / 1000, k=3,
                                           hypp=hypp,
                                           mcmc.params=mp)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1000,
                                                burnin=5000, thin=10)
  set.seed(123)
  (ml.cpp <- marginalLikelihood(model))

  mp <- McmcParams(iter=1000L)
  logLik <- modes(model)[["loglik"]] ## includes 2nd stage
  stage2.loglik <- stageTwoLogLik_pooled(model2)
  logPrior <- modes(model)[["logprior"]]
  pstar <- blockUpdatesPooledVar(model2, mp)
  set.seed(123)
  m.y <- logLik + stage2.loglik + logPrior - sum(pstar) + log(factorial(k(model)))
  ## there is a lot of label-switching with pooled variance model
  ## -- no need to incorporate biase correction
  ##log(factorial(k(model)))
  expect_equal(ml.cpp, m.y, tolerance=0.01)
  published.ml <- -226.803
  expect_equal(ml.cpp, published.ml, tolerance=1)


  ## 2 components, equal variance
  set.seed(123)
  hypp <- Hyperparameters(type="marginal", k=2, m2.0=100, mu.0=20)
  mp <- McmcParams(iter=1000, burnin=2000, thin=10, nStarts=10)
  model <- CNPBayes:::SingleBatchPooledVar(data=galaxies / 1000, k=2,
                                           hypp=hypp,
                                           mcmc.params=mp)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1000,
                                                burnin=5000, thin=10)
  set.seed(123)
  (ml.cpp <- marginalLikelihood(model))
  expect_equal(ml.cpp, expect=-239.764, tolerance=0.1)
})

test_that("marginalLikelihoodPooledVar", {
  library(MASS)
  data(galaxies)
  set.seed(42)
  galaxies[78] <- 26960
  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100, mu.0=20)
  mp <- McmcParams(iter=1000, burnin=1000, thin=5, nStarts=10)
  model <- CNPBayes:::SingleBatchPooledVar(data=galaxies / 1000, k=3,
                                           hypp=hypp,
                                           mcmc.params=mp)
  model <- CNPBayes:::posteriorSimulationPooled(model, iter=1000,
                                                burnin=1000, thin=10)

  model <- useModes(model)
  model.reduced <- model
  mcmcParams(model.reduced, force=TRUE) <- mp
  model.s20star <- reduced_s20(model.reduced)
  (p.cpp <- p_s20_reduced_pooled(model.s20star))

  K <- k(model.s20star)
  a <- a(hyperParams(model.s20star))
  b <- b(hyperParams(model.s20star))
  nu0star <- modes(model.s20star)[["nu0"]]
  s20star <- modes(model.s20star)[["sigma2.0"]]
  sigma2star <- modes(model.s20star)[["sigma2"]]
  a_k <- a + 0.5*K*nu0star
  b_k  <- b + 0.5*nu0star*1/sigma2star
  (p_s20 = dgamma(s20star, a_k, rate=b_k) )
  expect_equal(p.cpp, p_s20)
})



test_that("compare_galaxy_batch", {
  library(MASS)
  data(galaxies)
  set.seed(42)

  galaxies[78] <- 26960
  galaxies <- c(galaxies, galaxies + 5000)

  mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=20)
  hypp <- Hyperparameters(type="batch", k=2, m2.0=100)

  model <- BatchModel(data=galaxies / 1000, k=2,
                      hypp=hypp,
                      batch=rep(1:2, each=length(galaxies) / 2),
                      mcmc.params=mp)

  mlist.batch <- list(posteriorSimulation(model, k=1),
                      posteriorSimulation(model, k=2),
                      posteriorSimulation(model, k=3),
                      posteriorSimulation(model, k=4),
                      posteriorSimulation(model, k=5))

  lik.batch <- marginalLikelihood(mlist.batch)

  hypp <- Hyperparameters(type="marginal", k=2, m2.0=100)

  model <- MarginalModel(data=galaxies / 1000, k=2,
                         hypp=hypp,
                         mcmc.params=mp)

  mlist.marginal <- list(posteriorSimulation(model, k=1),
                         posteriorSimulation(model, k=2),
                         posteriorSimulation(model, k=3),
                         posteriorSimulation(model, k=4),
                         posteriorSimulation(model, k=5))

  lik.marginal <- marginalLikelihood(mlist.marginal)
  expect_true(lik.batch[3] >= max(lik.batch))
  expect_true(lik.batch[3] >= lik.marginal[3])
})

test_that("marginal_is_preferred", {
  library(MASS)
  data(galaxies)
  set.seed(42)

  galaxies[78] <- 26960
  galaxies <- c(galaxies, galaxies + 10)

  mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=20)
  hypp <- Hyperparameters(type="batch", k=3, m2.0=100)

  model <- BatchModel(data=galaxies / 1000, k=3,
                      hypp=hypp,
                      batch=rep(1:2, each=length(galaxies) / 2),
                      mcmc.params=mp)

  mlist.batch <- posteriorSimulation(model)
  lik.batch <- batchLikelihood(mlist.batch)

  hypp <- Hyperparameters(type="marginal", k=3, m2.0=100)

  model <- MarginalModel(data=galaxies / 1000, k=3,
                         hypp=hypp,
                         mcmc.params=mp)

  mlist.marginal <- posteriorSimulation(model)
  lik.marginal <- marginalLikelihood(mlist.marginal)
  expect_true(lik.batch <= lik.marginal)
})
