# code/guattery.R
# ---------------

source("code/graph.R")
source("code/laplacian.R")

guattery <- function(k) {
    adj <- matrix(0, 4 * k, 4 * k)

    # upper path
    for (i in 1:(2 * k - 1)) {
        adj[i, i + 1] <- 1
    }

    # lower path
    for (i in (2 * k + 1):(4*k - 1)) {
        adj[i, i + 1] <- 1
    }

    # body
    for (i in (k + 1):(2 * k)) {
        adj[i, i + 2 * k] <- 1
    }

    # symmetrize
    adj <- adj + t(adj)

    class(adj) <- "adj"
    attr(adj, "x") <- c(1:(2*k), 1:(2*k))
    attr(adj, "y") <- c(rep(1, 2*k), rep(0, 2*k))
    adj
}



k <- 10
adj <- guattery(k)
legs <-c(1:k, 2*k + 1:k)

lap <- laplacian(adj)
ilap <- pseudoinv(lap)

set.seed(0)
nreps <- 500
nreps <- 4*k
s1 <- matrix(NA, nreps, 4*k)
for (r in seq_len(nreps)) {
    #s0 <- rnorm(nrow(ilap))
    s0 <- rep(0, nrow(lap)); s0[r] <- 1
    s <- solve(lap + 1/nrow(lap), s0) - mean(s0)
    #s <- ilap %*% s0
    sign <- ifelse(sum(s[legs] > median(s)) < k, -1, 1)
    s1[r,] <- sign * s
}


pdf("plots/guattery.pdf", width=8, height=1.25)
par(mfrow=c(1,2))
par(oma=rep(0,4), mar=rep(1,4))
plot(adj, axes=FALSE, xlab="", ylab="", main="Global cut", asp=1.5)
plotcl(specclust(lap), adj)

plot(adj, axes=FALSE, xlab="", ylab="", main="Average local cut", asp=1.5)
plotcl(colMeans(s1) > median(colMeans(s1)), adj)
dev.off()


