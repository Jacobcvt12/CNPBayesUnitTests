
library(RUnit)
##library(CNPBayes)
library(MASS)
library(devtools)
load_all("~/Software/CNPBayes")
data(galaxies)
galaxies2 <- c(galaxies, galaxies + 5000)
batch <- rep(1:2, each=length(galaxies))
mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=20)
##
## Must make the priors much more uninformative
##
hypp <- Hyperparameters(type="batch", k=3, m2.0=6, eta.0=1.8)
model <- BatchModel(data=galaxies2/1000,
                    batch=batch,
                    k=3,
                    hypp=hypp,
                    mcmc.params=mp)
bmodel <- posteriorSimulation(model)
plot(bmodel, breaks=80)

loglik <- modes(bmodel)[["loglik"]]
bmodel2 <- useModes(bmodel)
s2.loglik <- .Call("stageTwoLogLikBatch", bmodel2)
complete.loglik <- loglik + s2.loglik
logprior <- modes(bmodel2)[["logprior"]]
loglikAndPrior <- complete.loglik + logprior

## check:  values at or near zero
ptheta.star <- .Call("marginal_theta_batch", bmodel2)
(p.theta.rb <- log(mean(ptheta.star)))

fit.pi.star <- fit
mcmcParams(fit.pi.star, force=TRUE) <- mp.reduced
fit.pi.star <- .Call("reduced_pi", fit.pi.star)
identical(modes(fit.pi.star), modes(fit))
p.pi.star <- .Call("p_pmix_reduced", fit.pi.star)
(p.pi.rb <- log(mean(p.pi.star)))
## check
zz <- z(chains(fit.pi.star))
gtools::ddirichlet(modes(fit)[["mixprob"]], alpha(hypp) + table(zz[2,]))
mp <- modes(fit)[["mixprob"]]
ztab <- tableZ(3, z(fit))
##ddirichlet(mp, 1+ztab)
