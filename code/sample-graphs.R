source("code/graph.R")

set.seed(0)

width <- 7
height <- 6
wrap <- FALSE
mu <- 2 * width * height - width - height
s <- 4

psamp <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0) # must be ascending
m <- ceiling(mu * psamp)




adj0 <- grid.2d(width, height, wrap = wrap)
adj0 <- shuffle.edges(adj0, s)


adj <- as.list(rep(NA, length(s)))
m1 <- 0
adj1 <- sample.edges(0, adj0)
for (i in seq_along(psamp)) {
	delta <- sample.edges(m[i] - m1, adj0)
	adj1 <- adj1 + delta
	m1 <- m[i]
	adj[[i]] <- adj1
}

pdf("plots/sample-graphs-long.pdf", width=15, height=2.5)
par(mfrow=c(1, length(psamp)))
#pdf("plots/sample-graphs.pdf", width=7.5, height=5)
#par(mfrow=c(2, (length(psamp) + 1) %/% 2))
for (i in seq_along(psamp)) {
	plot(adj[[i]], xlab="", ylab="", axes=FALSE,
	     main=bquote(m/mu == .(psamp[i])))
}
plot(20 * adj0, xlab="", ylab="", axes=FALSE,
     main=expression(m/mu == infinity))
dev.off()


