


setwd("../../fresh/glpath/code/")
system("shl")
dyn.load("gamlr.so")
for(f in Sys.glob("~/Projects/packages/gamlr/R/*.R")) source(f)
setwd("../../../packages/textir")
source("R/collapse.R")
source("R/mnlm.R")

library(parallel)
load("data/we8there.rda")

B <- mnlm(we8thereCounts,we8thereRatings,bins=5,cores=4, k=2)

C <- collapse(we8thereCounts,we8thereRatings,5)
v <- C$v
x <- C$x
mu <- rowSums(x)/ncol(x)
fit <- gamlr(v, x[,"serv best"], family="poisson", fix=mu, gam=1, k=2)
