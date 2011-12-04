
source("code/sample.R")

require(RColorBrewer)
palette(brewer.pal(7, "OrRd"))

h <- 6; w <- 7; mu <- w * h - w - h; n <- w * h
psamp <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
s <- 2^seq.int(0, 6)

#m <- mu * psamp
m <- 10 * psamp

alpha <- sqrt(s) * 2

#nrep <- 100
#tau <- sapply(1:nrep, function(seed) sapply(s, function(si) demo(nswaps=si, seed=seed)$tau))
#tau <- matrix(tau, nrep, length(s), byrow=TRUE)
#tau.mean <- colMeans(tau)

pdf("plots/optimal-theory.pdf", width = 5, height = 5)
par(mar=c(5,6,4,2) + .1)
plot(range(m), c(0,1), log="x", t="n",
     xlab=expression(m),
     ylab=expression(frac(m, m + 2 * (alpha - 1))))
#alpha <- c(1, 2, 5, 10, 20, 50, 100)
for (i in seq_along(alpha)) {
    lines(m, (1 + 2 * (alpha[i] - 1)/m)^(-1), col=i, lty=i)
}

ix <- rev(seq_along(alpha))
legend("topleft", inset = 0.05, title = expression(alpha),
       legend = format(alpha, digits = 1)[ix],
       col = seq_along(alpha)[ix], lty = seq_along(alpha)[ix])
dev.off()

