# code/interpolate.R
# ------------------

source("code/graph.R")

demo2.swaps <- function(width = 7, height = 7, wrap = TRUE,
                        nswaps = unique(c(0, 1, 2, 4,
                                          round(seq(0, sqrt(2 * width * height),
                                                    length.out = 16)^2)))[-1],
                        seed = 0,
                        nreps = 1000)
{
    adj <- grid.2d(width, height, wrap = wrap)
    n <- nrow(adj)
    deg <- colSums(adj)
    tot <- sum(deg) / 2
    lap <- diag(n) - diag(sqrt(1/deg)) %*% adj %*% diag(sqrt(1/deg))
    lambda0 <- eigen(lap, TRUE, TRUE)$values[-n]
    ilambda0 <- (1/lambda0) / sum(1/lambda0)
    
    lambda <- array(NA, c(nreps, length(nswaps), nrow(adj) - 1))
    for (i in seq_along(nswaps)) {
        set.seed(seed)
        for (r in seq_len(nreps)) {
            adj1 <- shuffle.edges(adj, nswaps[i])
            deg1 <- colSums(adj1)
            lap1 <- diag(n) - diag(sqrt(1/deg1)) %*% adj1 %*% diag(sqrt(1/deg1))
            lambda[r,i,] <- eigen(lap1, TRUE, TRUE)$values[-n]
        }
    }
    ilambda <- 1/lambda
    scale <- apply(ilambda, c(1,2), sum)
    ilambda.scale <- ilambda / rep(scale, n - 1)
    ilambda.mean <- apply(ilambda.scale, c(2, 3), mean)
    lambda.mean <- apply(lambda, c(2,3), mean)

    #par(mfrow=c(1,2))
    #plot(sqrt(c(0, max(nswaps))), c(0, 8), t='n',
    #     ylab="Laplacian Eigenvalues",
    #     xlab=expression(sqrt(Swaps)))
    #points(rep(0, nrow(adj) - 1), lambda0)
    #for (i in seq_along(nswaps)) {
    #    points(rep(sqrt(nswaps[i]), nrow(adj) - 1), lambda.mean[i,])
    #}

    plot(sqrt(c(0, max(nswaps)) / tot), c(0, max(ilambda0, max(ilambda.mean))),t='n',
           ylab="Inverse Laplacian Eigenvalues",
           xlab=expression(sqrt(Swaps / Edges)))
    points(rep(0, nrow(adj) - 1), ilambda0)
    for (i in seq_along(nswaps)) {
        points(rep(sqrt(nswaps[i] / tot), nrow(adj) - 1), ilambda.mean[i,])
    }
}


pdf("plots/interpolate.pdf", 4, 4)
demo2.swaps(7,6,FALSE)
dev.off()

