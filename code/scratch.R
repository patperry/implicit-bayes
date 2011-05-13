
source("code/graph.R")
source("code/laplacian.R")

PAUSE <- TRUE

nrep <- 50
nstart <- 5
npenalty <- 50

psamp = 0.5

penalty <- seq(0.1, 5000, len = npenalty)
pmid <- floor(npenalty / 2) + 1


width <- 10
height <- 5

steps <- 1

acc.samp <- array(NA, c(nrep, nstart))
acc.pr <- array(NA, c(nrep, nstart, npenalty))

accuracy <- function(cl.true, cl) {
    err <- min(mean(abs(cl.true - cl)), mean(abs(!cl.true - cl)))
    1 - err
}

for (r in seq_len(nrep)) {
    set.seed(r)

    adj <- grid.2d(width, height)
    adj.samp <- sample.edges(round(psamp * width * height * 4), adj)
    init <- matrix(rnorm(nstart * nrow(adj)), nstart, nrow(adj))
    
    lap <- laplacian(adj)
    cl.true <- specclust(lap)
    lap.samp <- laplacian(adj.samp)
    ilap.samp <- pseudoinv(lap.samp)
    cl.samp <- array(NA, c(nstart, nrow(lap)))
    cl.pr <- array(NA, c(nstart, nrow(lap), npenalty))
    
    for (s in seq_len(nstart)) {
        x0 <- init[s,,drop=TRUE]
        cl.samp[s,] <- powclust(ilap.samp, init = x0, steps = steps)
        #cl.samp[s,] <- specclust(lap.samp)
        acc.samp[r,s] <- accuracy(cl.samp[s,], cl.true)
        
        for (p in seq_len(npenalty)) {
            ilap.pr <- pagerank(lap.samp, penalty = penalty[p])
            cl.pr[s,,p] <- powclust(ilap.pr, init = x0, steps = steps)
            acc.pr[r,s,p] <- accuracy(cl.pr[s,,p], cl.true)
        }
    }

    acc.pr.avg <- apply(acc.pr[1:r,,,drop=FALSE], 3, mean)
    pbest <- which.max(acc.pr.avg)
    sbest.samp <- which.max(acc.samp[r,])
    sbest.pr <- which.max(acc.pr[r,,pbest])

    if (r > 1 && PAUSE)
        browser()
    
    par(mfrow=c(2,2))
    plot(adj.samp, main = paste("sample (", round(acc.samp[r,sbest.samp], 3), ")"))
    plotcl(cl.samp[sbest.samp,], adj.samp)
    
    plot(adj.samp, main = paste("pagerank", round(penalty[pbest]), " (",
                   round(acc.pr[r,sbest.samp,pbest], 3), ")"))
    plotcl(cl.pr[sbest.samp,,pbest], adj.samp)

    plot(adj.samp, main = paste("sample (", round(acc.samp[r,sbest.pr], 3), ")"))
    plotcl(cl.samp[sbest.pr,], adj.samp)
    
    plot(adj.samp, main = paste("pagerank", round(penalty[pbest]), " (",
                   round(acc.pr[r,sbest.pr,pbest], 3), ")"))
    plotcl(cl.pr[sbest.pr,,pbest], adj.samp)

}

cat("Sample:\n")
print(overall.samp <- summary(as.numeric(acc.samp)))
cat("Pagerank:\n")
print(overall.pr <- summary(as.numeric(acc.pr[,,pbest])))
