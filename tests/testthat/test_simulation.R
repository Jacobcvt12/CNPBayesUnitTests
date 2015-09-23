test_that("simulation", {
  ##library(devtools)
  ##load_all()
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  dat <- CNPBayes:::simulateProbeLevel(cnvs=1, K=4, probes=10,
                                       arguments=arguments,
                                       qual="easy")

  ## dimensions are samples x probes x cnp x components
  x <- dat[[1]]
  ## data for 15th CNP under 3-component mixture
  if(FALSE)
    hist(rowMeans(x[, ,1, 3]), col="gray", breaks=80)
  K <- 3
  xx <- x[, , 1, K]
  mns <- rowMeans(xx)
  pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
  if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
  if(FALSE)
    hist(pc, breaks=100, col="gray", border="gray")

  mp <- McmcParams(iter=1000, nStarts=10, burnin=1000, thin=5)
  model <- MarginalModel(data=pc, k=2,
                         mcmc.params=mp)
  singlebatch.models <- list(posteriorSimulation(model, k=1),
                             posteriorSimulation(model, k=2),
                             posteriorSimulation(model, k=3),
                             posteriorSimulation(model, k=4))
  ml <- marginalLikelihood(singlebatch.models)
  expect_true(which.max(ml) >= 3)
}
