# sample.R
# --------

source("code/graph.R")
source("code/laplacian.R")

rank1 <- function(x) {
    ev <- eigen(x)
    u <- ev$vectors[,1]
    u %*% t(u)
}

loss.lap <- function(lap0, lap1) {
    il0 <- pseudoinv(lap0); il0 <- il0 / sum(diag(il0))
    il1 <- pseudoinv(lap1); il1 <- il1 / sum(diag(il1))
    loss(il0, il1)
}

loss <- function(ilap0, ilap1) {
    ev <- eigen(ilap0 - ilap1)$values
    l2 <- max(abs(ev))
    lF <- sqrt(sum(ev^2))
    list(spectral = l2, frobenius = lF)
}

constant <- function(lap) {
    n <- nrow(lap)
    diag(1/n, n)
}

rank1 <- function(lap) {
    iln <- pseudoinv(lap)
    iln.eig <- eigen(iln)
    iln.1 <- iln.eig$vectors[,1] %*% t(iln.eig$vectors[,1,drop=FALSE])
}

mle <- function(lap) {
    ilap <- pseudoinv(lap)
    ilap / sum(diag(ilap))
}


demo <- function(width = 10, height = 5, nswaps = 0, psamp = 0.99,
                 penalty = seq(0.1, 50, length.out = 100), seed = 0, ev.tol = 1e-8) {
    set.seed(seed)

    # generate population graph and laplacian
    adj0 <- grid.2d(width, height)
    adj <- shuffle.edges(adj0, nswaps)
    nedge <- sum(adj > 0) / 2
    lap <- laplacian(adj)
    
    # generate sample graph and laplacian
    n <- floor(nedge * psamp)
    adj.samp <- sample.edges(n, adj)
    lap.samp <- laplacian(adj.samp)

    # compute the true normalized inverse laplacian
    iln <- pseudoinv(lap, scale = TRUE)

    # compute the MLE of normalized inverse laplacian
    iln.samp <- pseudoinv(lap.samp, scale = TRUE)

    # compute the pagerank estimate of normalized inverse laplacian
    npenalty <- length(penalty)
    loss.frobenius <- rep(NA, npenalty)
    loss.spectral <- rep(NA, npenalty)
    eigen.lap.samp <- eigen(lap.samp, symmetric = TRUE)
    ev.ok <- abs(eigen.lap.samp$values) > ev.tol
    if (sum(!ev.ok) > 1)
        warning("sample laplacian is not connected")

    for (i in seq_len(npenalty)) {
        iln.pr <- pagerank(lap.samp, eigen.lap.samp, penalty = penalty[[i]], ev.tol = ev.tol)
        l <- loss(iln, iln.pr)
        loss.frobenius[[i]] <- l$frobenius
        loss.spectral[[i]] <- l$spectral
    }

    penalty.mle = sum(1/eigen.lap.samp$values[ev.ok])
    loss.mle <- loss(iln, iln.samp)

    
    res <- list(penalty = penalty, loss.frobenius = loss.frobenius, loss.spectral = loss.spectral,
                penalty.mle = penalty.mle, loss.frobenius.mle = loss.mle$frobenius,
                loss.spectral.mle = loss.mle$spectral)

    par(mfrow=c(1,2))
    plot(penalty, loss.spectral, main="Spectral",
         ylim=c(min(c(loss.spectral, loss.mle$spectral)),
                max(c(loss.spectral, loss.mle$spectral))))
    abline(v=penalty[which.min(loss.spectral)], col=2)
    abline(h=loss.mle$spectral, col=3, lty=2)
    plot(penalty, loss.frobenius, main="Frobenius",
         ylim=c(min(c(loss.frobenius, loss.mle$frobenius)),
                max(c(loss.frobenius, loss.mle$frobenius))))
    abline(v=penalty[which.min(loss.frobenius)], col=2)
    abline(h=loss.mle$frobenius, col=3, lty=2)

    cat("penalty.mle: ", penalty.mle, "\n")
    cat("sample l2: ", loss(iln, iln.samp)$spectral, "\n")
    cat("sample lF: ", loss(iln, iln.samp)$frobenius, "\n")
    cat("best regularized l2: ", min(loss.spectral), "\n")
    cat("best regularized lF: ", min(loss.frobenius), "\n")

    invisible(res)
}

repdemo <- function(nreps = 200,  penalty = seq(0.1, 100, length.out = 100), ...) {
    npenalty <- length(penalty)

    seed <- seq_len(nreps)
    loss.spectral <- array(NA, c(nreps, npenalty))
    loss.frobenius <- array(NA, c(nreps, npenalty))
    loss.spectral.mle <- array(NA, c(nreps))
    loss.frobenius.mle <- array(NA, c(nreps))
    penalty.mle <- array(NA, c(nreps))

    for (r in seq_len(nreps)) {
        d <- demo(..., penalty = penalty, seed = seed[[r]])
        loss.spectral[r,] <- d$loss.spectral
        loss.frobenius[r,] <- d$loss.frobenius
        loss.spectral.mle[r] <- d$loss.spectral.mle
        loss.frobenius.mle[r] <- d$loss.frobenius.mle
        penalty.mle[r] <- d$penalty.mle
    }

    loss.spectral.mean <- apply(loss.spectral, 2, mean)
    loss.frobenius.mean <- apply(loss.frobenius, 2, mean)
    loss.spectral.mle.mean <- mean(loss.spectral.mle)
    loss.frobenius.mle.mean <- mean(loss.frobenius.mle)
    penalty.mle.mean <- mean(penalty.mle)

    loss.spectral.sd <- apply(loss.spectral, 2, sd)
    loss.frobenius.sd <- apply(loss.frobenius, 2, sd)
    loss.spectral.mle.sd <- sd(loss.spectral.mle)
    loss.frobenius.mle.sd <- sd(loss.frobenius.mle)
    penalty.mle.sd <- sd(penalty.mle)
    
    par(mfrow=c(1,2))
    plot(penalty, loss.spectral.mean, main="Spectral",
         ylim=c(min(c(loss.spectral.mean, loss.spectral.mle.mean)),
                max(c(loss.spectral.mean, loss.spectral.mle.mean))))
    lines(penalty, loss.spectral.mean + loss.spectral.sd)
    lines(penalty, loss.spectral.mean - loss.spectral.sd)    
    abline(h=loss.spectral.mle.mean, col=3, lty=2)
    abline(h=loss.spectral.mle.mean + loss.spectral.mle.sd, col=3, lty=1)
    abline(h=loss.spectral.mle.mean - loss.spectral.mle.sd, col=3, lty=1)
 
    plot(penalty, loss.frobenius.mean, main="Frobenius",
         ylim=c(min(c(loss.frobenius.mean, loss.frobenius.mle.mean)),
                max(c(loss.frobenius.mean, loss.frobenius.mle.mean))))
    lines(penalty, loss.frobenius.mean + loss.frobenius.sd)
    lines(penalty, loss.frobenius.mean - loss.frobenius.sd)    
    abline(h=loss.frobenius.mle.mean, col=3, lty=2)
    abline(h=loss.frobenius.mle.mean + loss.frobenius.mle.sd, col=3, lty=1)
    abline(h=loss.frobenius.mle.mean - loss.frobenius.mle.sd, col=3, lty=1)
}
