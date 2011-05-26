# code/interpolate.R
# ------------------

source("code/graph.R")

demo2.swaps <- function(width = 7, height = 7, wrap = TRUE,
                        nswaps = unique(c(0, 1, 2, 4,
                                          round(seq(0, sqrt(2 * width * height),
                                                    length.out = 16)^2)))[-1],
                        seed = 0,
                        nreps = 500)
{
    adj <- grid.2d(width, height, wrap = wrap)
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
    ilambda.mean <- apply(1/lambda, c(2,3), mean)

    par(mfrow=c(1,2))
    plot(sqrt(c(0, max(nswaps))), c(0, 8), t='n',
         ylab="Laplacian Eigenvalues",
         xlab=expression(sqrt(Swaps)))
    points(rep(0, nrow(adj) - 1), lambda0)
    for (i in seq_along(nswaps)) {
        points(rep(sqrt(nswaps[i]), nrow(adj) - 1), lambda.mean[i,])
    }

    plot(sqrt(c(0, max(nswaps))), c(0, max(1/lambda0, max(1/lambda.mean))),t='n',
           ylab="Inverse Laplacian Eigenvalues",
           xlab=expression(sqrt(Swaps)))
    points(rep(0, nrow(adj) - 1), 1/lambda0)
    for (i in seq_along(nswaps)) {
        points(rep(sqrt(nswaps[i]), nrow(adj) - 1), 1/lambda.mean[i,])
    }
}


pdf("plots/interpolate.pdf", 8, 4)
#demo2.swaps()
demo2.swaps(8,7,FALSE)
dev.off()

