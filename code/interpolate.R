# code/interpolate.R
# ------------------

source("code/graph.R")

demo2.swaps <- function(width = 7, height = 7,
                        nswaps = unique(c(0, 1, 2, 4,
                                          round(seq(0, sqrt(2 * width * height),
                                                    length.out = 16)^2)))[-1],
                        seed = 0,
                        nreps = 200)
{
    adj <- grid.2d(width, height, wrap = TRUE)
    lambda0 <- (4 - eigen(adj, TRUE, TRUE)$values[-1])
    
    lambda <- array(NA, c(nreps, length(nswaps), nrow(adj) - 1))
    for (i in seq_along(nswaps)) {
        set.seed(seed)
        for (r in seq_len(nreps)) {
            adj1 <- shuffle.edges(adj, nswaps[i])
            lambda[r,i,] <- (4 - eigen(adj1, TRUE, TRUE)$values[-1])
        }
    }
    lambda.mean <- apply(lambda, c(2,3), mean)

    par(mfrow=c(1,2))
    plot(c(0, 8), sqrt(c(0, max(nswaps))), t='n',
         xlab="Laplacian Eigenvalue",
         ylab=expression(sqrt(Swaps)))
    points(lambda0, rep(0, nrow(adj) - 1))
    for (i in seq_along(nswaps)) {
        points(lambda.mean[i,], rep(sqrt(nswaps[i]), nrow(adj) - 1))
    }

    plot(c(1/8, 1/min(lambda.mean)), sqrt(c(0, max(nswaps))), t='n',
           xlab="Inverse Laplacian Eigenvalue",
           ylab=expression(sqrt(Swaps)))
    points(1/lambda0, rep(0, nrow(adj) - 1))
    for (i in seq_along(nswaps)) {
        points(1/lambda.mean[i,], rep(sqrt(nswaps[i]), nrow(adj) - 1))
    }
}


pdf("plots/interpolate.pdf", 8, 4)
demo2.swaps()
dev.off()

