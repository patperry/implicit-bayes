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


demo <- function(width = 6, height = 7, wrap = FALSE, nswaps = 0, psamp = 0.99,
                 penalty = seq(0.1, 200, length.out = 100), seed = 0, ev.tol = 1e-8,
                 warn.connected = FALSE, quiet = TRUE)
{
    set.seed(seed)

    # generate population graph and laplacian
    adj0 <- grid.2d(width, height, wrap = wrap)
    adj <- shuffle.edges(adj0, nswaps)
    nedge <- sum(adj > 0) / 2
    lap0 <- laplacian(adj, normalized = FALSE)
    lap <- laplacian(adj, normalized = TRUE)
    
    # generate sample graph and laplacian
    n <- floor(nedge * psamp)
    adj.samp <- sample.edges(n, adj)
    lap0.samp <- laplacian(adj.samp, normalized = FALSE)
    lap.samp <- laplacian(adj.samp, normalized = TRUE)

    # compute the true normalized inverse laplacian
    iln0 <- pseudoinv(lap0, scale = TRUE)
    iln <- pseudoinv(lap, scale = TRUE)
    tau <- sum(diag(pseudoinv(lap, scale = FALSE)))

    # compute the MLE of normalized inverse laplacian
    iln0.samp <- pseudoinv(lap0.samp, scale = TRUE)
    iln.samp <- pseudoinv(lap.samp, scale = TRUE)    

    # compute the pagerank estimate of normalized inverse laplacian
    npenalty <- length(penalty)
    loss.frobenius <- rep(NA, npenalty)
    loss.spectral <- rep(NA, npenalty)
    eigen.lap.samp <- eigen(lap.samp, symmetric = TRUE)
    ev.ok <- abs(eigen.lap.samp$values) > ev.tol
    if (sum(!ev.ok) > 1 && warn.connected)
        warning("sample laplacian is not connected")

    for (i in seq_len(npenalty)) {
        iln.pr <- pagerank(lap.samp, eigen.lap.samp, penalty = penalty[[i]], ev.tol = ev.tol)
        l <- loss(iln, iln.pr)
        loss.frobenius[[i]] <- l$frobenius
        loss.spectral[[i]] <- l$spectral
    }

    penalty.mle = sum(1/eigen.lap.samp$values[ev.ok])
    loss.mle <- loss(iln, iln.samp)

    
    res <- list(tau = tau,
                penalty = penalty, loss.frobenius = loss.frobenius, loss.spectral = loss.spectral,
                penalty.mle = penalty.mle, loss.frobenius.mle = loss.mle$frobenius,
                loss.spectral.mle = loss.mle$spectral)

    if (quiet)
        return(res)

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

repdemo <- function(nreps = 200,  penalty = seq(0.1, 100, length.out = 100),
                    type=c("frobenius", "spectral"), quiet = FALSE, ...)
{
    type <- match.arg(type)
    npenalty <- length(penalty)

    seed <- seq_len(nreps)
    loss.spectral <- array(NA, c(nreps, npenalty))
    loss.frobenius <- array(NA, c(nreps, npenalty))
    loss.spectral.mle <- array(NA, c(nreps))
    loss.frobenius.mle <- array(NA, c(nreps))
    penalty.mle <- array(NA, c(nreps))
    tau <- rep(NA, npenalty)

    d <- NULL    
    for (r in seq_len(nreps)) {
        d <- demo(..., penalty = penalty, seed = seed[[r]])
        loss.spectral[r,] <- d$loss.spectral
        loss.frobenius[r,] <- d$loss.frobenius
        loss.spectral.mle[r] <- d$loss.spectral.mle
        loss.frobenius.mle[r] <- d$loss.frobenius.mle
        penalty.mle[r] <- d$penalty.mle
        tau[r] <- d$tau
    }

    penalty <- penalty / mean(tau)
    
    loss.spectral.mean <- apply(loss.spectral / rep(loss.spectral.mle, npenalty), 2, mean)
    loss.frobenius.mean <- apply(loss.frobenius / rep(loss.frobenius.mle, npenalty), 2, mean)
    #loss.spectral.mle.mean <- mean(loss.spectral.mle)
    #loss.frobenius.mle.mean <- mean(loss.frobenius.mle)
    #penalty.mle.mean <- mean(penalty.mle)

    loss.spectral.sd <- apply(loss.spectral / rep(loss.spectral.mle, npenalty), 2, sd)
    loss.frobenius.sd <- apply(loss.frobenius / rep(loss.frobenius.mle, npenalty), 2, sd)
    #loss.spectral.mle.sd <- sd(loss.spectral.mle)
    #loss.frobenius.mle.sd <- sd(loss.frobenius.mle)
    #penalty.mle.sd <- sd(penalty.mle)

    res <- NULL
    if (quiet)
        return(res)
    
    par(mfrow=c(1,length(type)))
    if ("spectral" %in% type) {
        plot(range(penalty), range(loss.spectral.mean), t = 'n',
             main = "",
             xlab = "Regularization",
             ylab = "Rel. Spectral Error",
             ylim=range(c(0,
                          loss.spectral.mean - loss.spectral.sd,
                          loss.spectral.mean + loss.spectral.sd)))
        abline(h = 1, lty = 2)        
        lines(penalty, loss.spectral.mean + loss.spectral.sd)
        lines(penalty, loss.spectral.mean - loss.spectral.sd)
        points(penalty, loss.spectral.mean)        
    }

    if ("frobenius" %in% type) {
        plot(range(penalty), range(loss.frobenius.mean), t = 'n',
             main = "",
             xlab = "Regularization",
             ylab = "Rel. Frobenius Error",
             ylim=range(c(0,
                          loss.frobenius.mean - loss.frobenius.sd,
                          loss.frobenius.mean + loss.frobenius.sd)))
        abline(h = 1, lty = 2)        
        lines(penalty, loss.frobenius.mean + loss.frobenius.sd)
        lines(penalty, loss.frobenius.mean - loss.frobenius.sd)
        points(penalty, loss.frobenius.mean)
    }

    invisible(res)
}
