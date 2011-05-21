# graph.R
# -------

grid.2d <- function(width, height, wrap = FALSE) {
  if (width < 0 || height < 0)
    warning("width and height should be positive")
  
  n <- width * height
  adj <- matrix(0, n, n)

  ix <- matrix(seq_len(n), height, width)

  # handle the interior
  for (j in 1:width) {
      for (i in 1:height) {
          if (i < height)
              adj[ix[i,j], ix[i+1,j  ]] <- 1
          if (j < width)
              adj[ix[i,j], ix[i  ,j+1]] <- 1
          if (i > 1)
              adj[ix[i,j], ix[i-1,j  ]] <- 1
          if (j > 1)
              adj[ix[i,j], ix[i  ,j-1]] <- 1
      }
  }

  if (wrap) {
      if (width > 1) {
          for (i in 1:height) {
              adj[ix[i,1    ], ix[i,width]] <- 1
              adj[ix[i,width], ix[i,1    ]] <- 1
          }
      }
      if (height > 1) {
          for (j in 1:width) {
              adj[ix[1     ,j], ix[height,j]] <- 1
              adj[ix[height,j], ix[1     ,j]] <- 1
          }
      }
  }

  class(adj) <- "adj"
  attr(adj, "x") <- rep(1:width, each = height)
  attr(adj, "y") <- rep(1:height, times = width)
  adj
}

plot.adj <- function(object, ...) {
    adj <- object
    plot(attr(adj, "x"), attr(adj, "y"), ...)
    for (i in seq_len(nrow(adj))) {
        for (j in seq_len(ncol(adj))) {
            if (adj[i,j] > 0)
                segments(attr(adj, "x")[i], attr(adj, "y")[i],
                         attr(adj, "x")[j], attr(adj, "y")[j], ...)
        }
    }
}

plot.adj.compact <- function(object, ...) {
    adj.compact <- object
    plot(attr(adj.compact, "x"), attr(adj.compact, "y"), ...)
    for (i in seq_len(nrow(adj.compact))) {
        for (k in seq_len(ncol(adj.compact))) {
            j <- adj.compact[i,k]
            if (!is.na(j)) {
                segments(attr(adj, "x")[i], attr(adj, "y")[i],
                         attr(adj, "x")[j], attr(adj, "y")[j], ...)
            }
        }
    }
}

degrees <- function(adj, ...) UseMethod("degrees")
degrees.adj <- function(adj, ...) colSums(adj)
degrees.adj.compact <- function(adj, ...) rowSums(!is.na(adj))

to.compact <- function(adj) {
  n <- nrow(adj)
  d <- rowSums(adj)
  dmax <- max(d)
  adj.compact <- matrix(NA, n, dmax)
  for (i in seq_len(n)) {
    nbhd <- which(adj[i,] > 0)
    adj.compact[i,seq_along(nbhd)] <- nbhd
  }
  mostattributes(adj.compact) <- attributes(adj)
  class(adj.compact) <- "adj.compact"
  attr(adj.compact, "dim") <- c(n, dmax)
  adj.compact
}

from.compact <- function(adj.compact) {
  n <- nrow(adj.compact)
  adj <- matrix(0, n, n)
  for (i in seq_len(n)) {
    adj[i,adj.compact[i,]] <- 1
  }
  mostattributes(adj) <- attributes(adj.compact)
  class(adj) <- "adj"
  attr(adj, "dim") <- c(n,n)
  adj
}

dual.edge <- function(adj.compact, node, edge) {
  j <- adj.compact[node,edge]
  f <- which(adj.compact[j,] == node)[1]
  list(node = j, edge = f)
}

swap.edge <- function(adj.compact, node1, edge1, node2, edge2) {
  if (node1 != node2) {
    dual1 <- dual.edge(adj.compact, node1, edge1)
    dual2 <- dual.edge(adj.compact, node2, edge2)

    if (dual1$node != dual2$node) {
      adj.compact[node1,edge1] <- dual2$node
      adj.compact[dual2$node, dual2$edge] <- node1
      adj.compact[node2,edge2] <- dual1$node
      adj.compact[dual1$node, dual1$edge] <- node2
    }
  }

  adj.compact
}


shuffle.edges <- function(adj, nswaps) {
  if (nswaps == 0)
    return(adj)
  
  adj.compact <- to.compact(adj)
  n <- nrow(adj.compact)
  k <- ncol(adj.compact)

  d <- degrees(adj.compact)
  i1 <- sample.int(n, nswaps, replace = TRUE, prob = d)
  e1 <- 1 + floor(runif(nswaps) * d[i1])
  i2 <- sample.int(n, nswaps, replace = TRUE, prob = d)
  e2 <- 1 + floor(runif(nswaps) * d[i2])

  invalid <- 0
  for (s in seq_len(nswaps)) {
      if(adj.compact[i1[s],e1[s]] %in% c(i2[s], adj.compact[i2[s],])
         || adj.compact[i2[s],e2[s]] %in% c(i1[s], adj.compact[i1[s],])) {
          invalid <- invalid + 1
      } else {
          adj.compact <- swap.edge(adj.compact, i1[s], e1[s], i2[s], e2[s])
      }
  }

  if (!all(degrees(adj.compact) == d))
      browser()
  
  adj1 <- from.compact(adj.compact)
  attributes(adj1) <- attributes(adj)
  shuffle.edges(adj1, invalid)
}


add.edges <- function(adj, nadd) {
    if (nadd == 0)
        return(adj)

    n <- nrow(adj)
    i1 <- sample.int(n, nadd, replace = TRUE)
    i2 <- 1 + (sample.int(n - 1, nadd, replace = TRUE) + i1 - 1) %% n

    for (i in seq_len(nadd)) {
        adj[i1[i], i2[i]] <- adj[i1[i], i2[i]] + 1
        adj[i2[i], i1[i]] <- adj[i2[i], i1[i]] + 1
    }

    adj
}

sample.edges <- function(n, adj) {
    e.ix <- which(adj > 0) # includes duplicate edges (both i ~ j and j ~ i)
    ix <- suppressWarnings(sample(e.ix, n, replace = TRUE, prob = adj[e.ix]))
    ix.tab <- tabulate(ix, nbins = nrow(adj) * ncol(adj))
    adj.samp <- matrix(ix.tab, nrow(adj), ncol(adj))
    adj.samp <- adj.samp + t(adj.samp)
    attributes(adj.samp) <- attributes(adj)
    adj.samp
}


demo.swaps <- function(width = 25, height = 25, nswaps = 2^seq(0,10), seed = 0,
                 hist = TRUE) {
  require("RColorBrewer")
  pal <- brewer.pal(4, "OrRd")
  col <- colorRampPalette(pal, space="Lab")(length(nswaps) + 1)
  
  adj <- grid.2d(width, height, wrap = TRUE)
  lambda <- (4 - eigen(adj, TRUE, TRUE)$values[-1])
  plot(lambda, col=col[1], t="l",
       xlab="Index", ylab="(Unnormalized) Laplacian Eigenvalue")
  
  for (i in seq_along(nswaps)) {
    set.seed(seed)
    adj1 <- shuffle.edges(adj, nswaps[i])
    lambda1 <- (4 - eigen(adj1, TRUE, TRUE)$values[-1])
    lines(lambda1, col=col[1+i])
  }

  axis(3, labels=FALSE)
  axis(4, labels=FALSE)

  legend("bottomright", legend=c(0,nswaps), col=col, lty=1, lwd=3,
         title="Edge Swaps")
}
