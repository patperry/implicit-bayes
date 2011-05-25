# laplacian.R
# -----------

source("code/graph.R")

laplacian <- function(adj, ...) {
  d <- degrees(adj)
  lap <- -adj
  diag(lap) <- diag(lap) + d
  lap
}


pseudoinv <- function(a, tol = 1e-8, scale = FALSE) {
    udv <- svd(a)
    active <- udv$d > tol
    ainv <- udv$u[,active] %*% ((1/udv$d[active]) * t(udv$v[,active]))
    if (scale)
        ainv <- ainv / sum(diag(ainv))
    ainv
}


pagerank <- function(lap, lap.eigen = eigen(lap), penalty = 1.0, ev.tol = 1e-8) {
    if (penalty <= 0)
        stop("'penalty' is negative")
    
    p <- nrow(lap)
    ev <- lap.eigen
    ok <- abs(ev$values) > 1e-8
    ev$values <- ev$values[ok]
    ev$vectors <- ev$vectors[,ok]

    f <- function(lambda) { sum(-1/(penalty*(lambda - ev$values))) - 1 }
    lower <- -1; while(is.finite(lower) && f(lower) >= 0) { lower <- lower * 10 }
    upper <- min(ev$values) - 1e-10
    lambda <- uniroot(f, c(lower, upper))$root

    est <- ev$vectors %*% (-1/(penalty*(lambda - ev$values)) * t(ev$vectors))
    est
}


specclust <- function(lap, ev.tol = 1e-9) {
    eig <- eigen(lap, symm = TRUE)
    val <- eig$values
    ok <- abs(eig$values) > 1e-8
    vec <- eig$vectors
    u <- vec[,max(which(ok))]
    threshold <- median(u)
    cl <- u > threshold
    if (!cl[1])
        cl <- !cl
    cl
}

powclust <- function(ilap, steps = 10, init = rnorm(nrow(ilap))) {
    x <- init
    x <- x / as.numeric(sqrt(t(x) %*% x))
    for (i in seq_len(steps)) {
        x <- ilap %*% x
        x <- x / as.numeric(sqrt(t(x) %*% x))
    }
    threshold <- median(x)
    cl <- x > threshold
    if (!cl[1])
        cl <- !cl
    cl
}

plotcl <- function(cl, adj, col1 = "black", col2 = "white") {
    cols <- c(col1, col2)
    points(attr(adj, "x"), attr(adj, "y"), col = cols[cl + 1], pch=16)
}
