# sample.R
# --------

source("code/kregular.R")


sample.adj <- function(n, adj) {
    e.ix <- which(adj > 0) # includes duplicate edges (both i ~ j and j ~ i)
    ix <- suppressWarnings(sample(e.ix, n, replace = TRUE, prob = adj[e.ix]))
    ix.tab <- tabulate(ix, nbins = nrow(adj) * ncol(adj))
    adj.samp <- matrix(ix.tab, nrow(adj), ncol(adj))
    adj.samp <- adj.samp + t(adj.samp)
    adj.samp
}

laplacian <- function(adj) {
    d <- rowSums(adj)
    l <- -adj
    diag(l) <- diag(l) + d
    l
}


pseudoinv <- function(a, tol = 1e-8, scale = FALSE) {
    udv <- svd(a)
    active <- udv$d > tol
    ainv <- udv$u[,active] %*% ((1/udv$d[active]) * t(udv$v[,active]))
    if (scale)
        ainv <- ainv / sum(diag(ainv))
    ainv
}

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

pagerank <- function(lap, penalty = 1.0) {
    if (penalty <= 0)
        stop("'penalty' is negative")
    
    p <- nrow(lap)
    ev <- eigen(lap)
    ev$values <- ev$values[-p]
    ev$vectors <- ev$vectors[,-p]

    f <- function(lambda) { sum(-1/(penalty*(lambda - ev$values))) - 1 }
    lower <- -1; while(is.finite(lower) && f(lower) >= 0) { lower <- lower * 10 }
    upper <- min(ev$values) - 1e-10
    lambda <- uniroot(f, c(lower, upper))$root

    est <- ev$vectors %*% (-1/(penalty*(lambda - ev$values)) * t(ev$vectors))
    est
}


demo <- function(width = 4, height = 4, nswaps = 1024, psamp = 0.99, seed = 0) {
    set.seed(seed)

    # generate population graph and laplacian
    adj0 <- grid.2d(width, height)
    adj <- shuffle.edges(adj0, nswaps)
    nedge <- sum(adj > 0) / 2
    lap <- laplacian(adj)
    
    # generate sample graph and laplacian
    n <- floor(nedge * psamp)
    adj.samp <- sample.adj(n, adj)
    lap.samp <- laplacian(adj.samp)

    # compute the true normalized inverse laplacian
    iln <- pseudoinv(lap, scale = TRUE)

    # compute the MLE of normalized inverse laplacian
    iln.samp <- pseudoinv(lap.samp, scale = TRUE)

    # compute the pagerank estimate of normalized inverse laplacian
    penalty <- seq(0.1, 20, len = 100)
    npenalty <- length(penalty)    
    loss.frobenius <- rep(NA, npenalty)
    loss.spectral <- rep(NA, npenalty)
    for (i in seq_len(npenalty)) {
        iln.pr <- pagerank(lap.samp, penalty = penalty[[i]])
        l <- loss(iln, iln.pr)
        loss.frobenius[[i]] <- l$frobenius
        loss.spectral[[i]] <- l$spectral
    }

    par(mfrow=c(1,2))
    plot(penalty, loss.spectral, main="Spectral",
         ylim=c(min(c(loss.spectral, loss(iln, iln.samp)$spectral)),
                max(c(loss.spectral, loss(iln, iln.samp)$spectral))))
    abline(v=penalty[which.min(loss.spectral)], col=2)
    abline(h=loss(iln, iln.samp)$spectral, col=3, lty=2)
    plot(penalty, loss.frobenius, main="Frobenius",
         ylim=c(min(c(loss.frobenius, loss(iln, iln.samp)$frobenius)),
                max(c(loss.frobenius, loss(iln, iln.samp)$frobenius))))
    abline(v=penalty[which.min(loss.frobenius)], col=2)
    abline(h=loss(iln, iln.samp)$frobenius, col=3, lty=2)    

    cat("sample l2: ", loss(iln, iln.samp)$spectral, "\n")
    cat("sample lF: ", loss(iln, iln.samp)$frobenius, "\n")
    cat("best regularized l2: ", min(loss.spectral), "\n")
    cat("best regularized lF: ", min(loss.frobenius), "\n")

}
