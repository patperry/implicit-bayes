# code/guattery.R
# ---------------


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
