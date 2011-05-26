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



k <- 4
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


#pdf("plots/guattery.pdf", width=8, height=1.25)



ramp <- (function() {
            gamma <- 1
            function(u) { gray((1-u)^(1/gamma)) }
        })()

 pdf("plots/guattery-global.pdf", width=2, height=1)
par(mfrow=c(1,1), oma=rep(0,4), mar=rep(1, 4))
u <- eigen(lap)$vectors[,4*k-1]
u <- u - min(u)
u <- u / max(u)
plot(adj, axes=FALSE, xlab="", ylab="", main="", asp=1)
points(attr(adj, "x"), attr(adj, "y"), col = ramp(u), pch=16, cex = 1.2)
points(attr(adj, "x"), attr(adj, "y"), col = 1, pch=1, cex = 1.2)
dev.off()


pdf("plots/guattery-local.pdf", width=8, height=2)
v <-  seq(1, 2 *k, by = 1)
par(mfrow=c(2, length(v)/2), oma=rep(0, 4), mar=rep(1, 4))

for (i in seq_along(v)) {
    s0 <- rep(0, nrow(lap)); s0[v[i]] <- 1
    u <- solve(lap + 1/nrow(lap), s0) - mean(s0)
    u <- u - min(u)
    u <- u / max(abs(u))
    plot(adj, axes=FALSE, xlab="", ylab="", main="", asp=1)
    points(attr(adj, "x"), attr(adj, "y"), col = ramp(u), pch=16, cex=2)
    points(attr(adj, "x"), attr(adj, "y"), col = 1, pch=1, cex=2)
}
dev.off()



pr <- (function() {
    deg <- colSums(adj)
    mat <- adj %*% diag(1 / deg)
    eye <- diag(nrow(adj))
    
    function(s, gamma) {
        gamma * solve(eye - (1 - gamma) * mat, s)
    }
})()

pdf("plots/guattery-pagerank.pdf", width=8, height=2)
par(mfrow=c(2, length(v)/2), oma=rep(0, 4), mar=rep(1, 4))

for (i in seq_along(v)) {
    s0 <- rep(0, nrow(lap)); s0[v[i]] <- 1
    u <- pr(s0, 0.1)
    u <- u - min(u)
    u <- u / max(abs(u))
    plot(adj, axes=FALSE, xlab="", ylab="", main="", asp=1)
    points(attr(adj, "x"), attr(adj, "y"), col = ramp(u), pch=16, cex=2)
    points(attr(adj, "x"), attr(adj, "y"), col = 1, pch=1, cex=2)
}
dev.off()

par(mfrow=c(1,2))
par(oma=rep(0,4), mar=rep(1,4))
plot(adj, axes=FALSE, xlab="", ylab="", main="Global cut", asp=1.5)
plotcl(specclust(lap), adj)

plot(adj, axes=FALSE, xlab="", ylab="", main="Average local cut", asp=1.5)
plotcl(colMeans(s1) > median(colMeans(s1)), adj)
#dev.off()


