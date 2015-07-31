library(MASS)
data(galaxies)
galaxies[78] <- 26960

N <- 1000000

p <- runif(N)
v <- rgamma(N,1/3)/20
u1 <- rnorm(N,20,10)
u2 <- rnorm(N,20,10)

l <- rep(1,N)
for(x in galaxies) {
    l <- l*( p * sqrt(v/(2*pi)) * exp(-v*(x-u1)^2/2)
            + (1-p)* sqrt(v/(2*pi)) * exp(-v*(x-u2)^2/2))
}

m <- mean(l)
e <- sqrt(var(l)/N)
