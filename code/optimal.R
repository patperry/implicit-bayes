source("code/sample.R")

set.seed(0)

require(RColorBrewer)
palette(brewer.pal(7, "OrRd"))


optdemo <- function(psamp = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                    swaps = 2^seq.int(0, 6),
                    type=c("frobenius", "spectral"), ...)
{
    type <- match.arg(type)
    npsamp <- length(psamp)
    nswaps <- length(swaps)

    eta.frob <- matrix(NA, nswaps, npsamp)
    eta.spec <- matrix(NA, nswaps, npsamp)    
    
    for (i in seq_along(swaps)) {
        for (j in seq_along(psamp)) {
            d <- repdemo(nswaps = swaps[i], psamp = psamp[j],
                         type = c("frobenius", "spectral"), quiet = TRUE, ...)
            eta.frob[i,j] <- d[1]
            eta.spec[i,j] <- d[2]
        }
    }

    adj <- grid.2d(7, 6, wrap=FALSE)
    tot <- sum(adj) / 2

    plot(range(psamp), range(eta.frob), t='n', log = "x",
         xlab = "Sample Proportion",
         ylab = "Optimal Penalty")
    for (i in seq_along(swaps)) {
        points(psamp, eta.frob[i,], col = i)
        lines(psamp, eta.frob[i,], col = i, lty = i)        
    }
    ix <- nswaps:1
    legend("topleft", inset = 0.05, title = "Swaps/Edges",
           legend = format(swaps / tot, digits = 1)[ix],
           col = seq_along(swaps)[ix], lty = seq_along(swaps)[ix])
}
    

pdf("plots/optimal.pdf", width = 5, height = 5)
optdemo(width = 6, height = 7, wrap = FALSE,
        nreps = 100,
        type="frobenius")
dev.off()

