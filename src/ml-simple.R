#!/usr/bin/env R
##install.package("plot3D")               # if necessary
library(plot3D)
y <- c(11, 7, 12)                       # the observations
m <- mean(y)
s <- sd(y)
x <- seq(m-1, m+1, .01)                 # grid centers about mean(y)
y <- seq(s-1, s+1, .01)                 # and about sd(y)
g <- mesh(x, y)                         # the grid
## probability density
pd <- matrix(vector(mode = "numeric", length = length(g$x)), nrow=length(x))
for (i in 1:length(x))
    for (j in 1:length(x))
        pd[i, j] = prod(pnorm(y, mean=g$x[i, j], sd=g$y[i, j]))
surf3D(g$x, g$y, pd)
