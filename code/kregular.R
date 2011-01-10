# kregular.R
# ----------

nnodes <- 10

grid.2d <- function(width, height) {
  if (width < 3 || height < 3)
    warning("width and height should be at least 3")
  
  n <- width * height
  x <- matrix(0, n, n)

  ix <- matrix(seq_len(n), height, width)
    
  for (j in seq_len(width)) {
    jprev <- ifelse(j == 1, width, j - 1)
    jnext <- ifelse(j == width, 1, j + 1)

    for (i in seq_len(height)) {
      iprev <- ifelse(i == 1, height, i - 1)
      inext <- ifelse(i == height, 1, i + 1)

      x[ix[i,j], ix[inext,j    ]] <- 1
      x[ix[i,j], ix[i    ,jnext]] <- 1
      x[ix[i,j], ix[iprev,j    ]] <- 1
      x[ix[i,j], ix[i    ,jprev]] <- 1
    }
  }

  x
}

laplacian <- function(adj) {
  d <- colSums(adj)
  l <- -adj
  diag(l) <- diag(l) + d
  l
}

# only works for regular graphs
to.compact <- function(adj) {
  n <- nrow(adj)
  d <- rowSums(adj)
  dmax <- max(d)
  adj.compact <- matrix(NA, n, dmax)
  for (i in seq_len(n)) {
    nbhd <- which(adj[i,] > 0)
    adj.compact[i,seq_along(nbhd)] <- nbhd
  }
  adj.compact
}

from.compact <- function(adj.compact) {
  n <- nrow(adj.compact)
  adj <- matrix(0, n, n)
  for (i in seq_len(n)) {
    adj[i,adj.compact[i,]] <- 1
  }
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

  i1 <- sample.int(n, nswaps, replace = TRUE)
  e1 <- sample.int(k, nswaps, replace = TRUE)
  i2 <- sample.int(n, nswaps, replace = TRUE)
  e2 <- sample.int(k, nswaps, replace = TRUE)

  invalid <- 0
  for (s in seq_len(nswaps)) {
    if(adj.compact[i1[s],e1[s]] %in% c(i2[s], adj.compact[i2[s],])
       || adj.compact[i2[s],e2[s]] %in% c(i1[s], adj.compact[i1[s],])) {
      invalid <- invalid + 1
    } else {
      adj.compact <- swap.edge(adj.compact, i1[s], e1[s], i2[s], e2[s])
    }
  }

  adj <- from.compact(adj.compact)
  shuffle.edges(adj, invalid)
}


demo <- function(width = 25, height = 25, nswaps = 2^seq(0,10), seed = 0) {
  require("RColorBrewer")
  pal <- brewer.pal(4, "OrRd")
  col <- colorRampPalette(pal, space="Lab")(length(nswaps) + 1)
  
  adj <- grid.2d(width, height)
  plot((4 - eigen(adj, TRUE, TRUE)$values[-1]), col=col[1], t="l",
       xlab="Index", ylab="(Unnormalized) Laplacian Eigenvalue")
  
  for (i in seq_along(nswaps)) {
    set.seed(seed)
    adj1 <- shuffle.edges(adj, nswaps[i])
    lines((4 - eigen(adj1, TRUE, TRUE)$values[-1]), col=col[1+i])
  }

  axis(3, labels=FALSE)
  axis(4, labels=FALSE)

  legend("bottomright", legend=c(0,nswaps), col=col, lty=1, lwd=3,
         title="Edge Swaps")
}
