library(gtools)
library(MASS)

### data
data(galaxies)
galaxies[78] <-26960

#  nsim <- 5e7
 nsim <- 5e6

## hyperparams
mu.0 <- 20
tau2.0 <- 100
eta.0 <- 1
m2.0 <- 0.1
alpha <- c(1,1,1)
beta <- 0.1
a <- 1.8
b <- 6
eta.0 <- 1
m2.0 <- 0.1

set.seed(123)
## simulate hyperpriors
mu <- rnorm(nsim, mu.0, tau2.0)
tau2 <- 1/rgamma(nsim, 0.5*eta.0, 0.5*eta.0*m2.0)
sigma2.0 <- rgamma(nsim, a, b)
nu.0 <- rgeom(nsim, beta)
nu.0[nu.0 == 0] <- 1

## simulate theta, sigma2, p
theta <- cbind(rnorm(nsim, mu, sqrt(tau2)),
               rnorm(nsim, mu, sqrt(tau2)),
               rnorm(nsim, mu, sqrt(tau2)))
sigma2 <- cbind(rgamma(nsim, 0.5*nu.0, 0.5*nu.0*sigma2.0),
                rgamma(nsim, 0.5*nu.0, 0.5*nu.0*sigma2.0),
                rgamma(nsim, 0.5*nu.0, 0.5*nu.0*sigma2.0))
p <- rdirichlet(nsim, alpha)


l <- rep(1, nsim) # vector of likelihoods
for(y in galaxies) {
    l <- l *  (p[,1]*sqrt(1/(2*sigma2[,1]*pi)) * exp(-(y - theta[,1])^2/(2*sigma2[,1]))
              +p[,2]*sqrt(1/(2*sigma2[,2]*pi)) * exp(-(y - theta[,2])^2/(2*sigma2[,2])) 
              +p[,3]*sqrt(1/(2*sigma2[,3]*pi)) * exp(-(y - theta[,3])^2/(2*sigma2[,3])) )
}

## dnorm possibly less stable than writing out normal density
# for(y in yg) {
#     l <- l* ( p[,1]*dnorm(y, mean=theta[,1], sd=sqrt(sigma2[,1])) + 
#              p[,2]*dnorm(y, mean=theta[,2], sd=sqrt(sigma2[,2])) + 
#              p[,3]*dnorm(y, mean=theta[,3], sd=sqrt(sigma2[,3])))
# }

mlik <- log(mean(l))
e <- sqrt(var(l)/nsim)

print(mlik)
print(e)
print(log(mean(l) + c(-1,1) * 1.96 * e))

attr(sim.data, "mloglik") <- mlik
attr(sim.data, "mlikerror") <- e
saveRDS(sim.data, "sim_marglik.RDS")
