source("code/graph.R")

set.seed(0)

width <- 7
height <- 6
wrap <- FALSE
s <- (0:7)^2 # must be sorted in ascending order


adj0 <- grid.2d(width, height, wrap = wrap)


adj <- as.list(rep(NA, length(s)))
adj1 <- adj0
s1 <- 0
for (i in seq_along(s)) {
	adj1 <- shuffle.edges(adj1, s[i] - s1)
	s1 <- s[i]
	adj[[i]] <- adj1
}

pdf("plots/interpolate-graphs-long.pdf", width=15, height=2.5)
par(mfrow=c(1, length(s)))
#pdf("plots/interpolate-graphs.pdf", width=7.5, height=5)
#par(mfrow=c(2, (length(s) + 1) %/% 2))
for (i in seq_along(s)) {
	plot(adj[[i]], xlab="", ylab="", axes=FALSE,
	     main=paste("s =", s[i]))
}
dev.off()


