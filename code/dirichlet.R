# code/dirichlet.R
# ----------------

pdf("plots/dirichlet.pdf", 4, 4)

n <- 35
nreps <- 500
shape <- seq(1.01, 10, len=16)
seed <- 0

lambda <- array(NA, c(nreps, length(shape), n))

for (i in seq_along(shape)) {
    set.seed(seed)
    for (r in seq_len(nreps)) {
        y <- rgamma(n, shape = shape[i])
        lambda[r,i,] <- sort(y) / sum(y)
    }
}

lambda.mean <- apply(lambda, c(2,3), mean)


par(mfrow=c(1,1))

plot(range(shape), range(lambda.mean), t='n',
     ylab = "Order statistics",
     xlab = expression(Shape))
for (i in seq_along(shape)) {
    points( rep(shape[i], n), lambda.mean[i,])
}

#plot(range(0, 1), range(shape), t='n',
#     xlab = "Order statistics",
#     ylab = "Dirichlet parameter")
#for (i in seq_along(shape)) {
#    x <- lambda.mean[i,]
#    x <- x - min(x)
#    x <- x / max(x)
#    points(x, rep(shape[i], n))
#}


#plot(range(1 / lambda.mean), range(shape), t='n',
#     xlab = "Inverse Order statistics",
#     ylab = "Dirichlet parameter")
#for (i in seq_along(shape)) {
#    points(1 / lambda.mean[i,], rep(shape[i], n))
#}
#ilambda.mean <- lambda.mean
#for (i in seq_along(shape)) {
#    scale <- max(ilambda.mean[i,])
#    ilambda.mean[i,] <- ilambda.mean[i,] / scale
#}
#plot(range(ilambda.mean), range(shape), t='n',
#     xlab = "Scaled Inverse Order statistics",
#     ylab = "Dirichlet parameter")
#for (i in seq_along(shape)) {
#    points(ilambda.mean[i,], rep(shape[i], n))
#}

dev.off()
