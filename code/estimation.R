source("code/sample.R")

pdf("plots/estimation-frob-p010.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p020.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p050.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.50,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p100.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 1,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p200.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 2,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p500.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 5,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p1000.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-spec-p010.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p020.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p050.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.50,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p100.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 1,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-spec-p200.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 2,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p500.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 5,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p1000.pdf", 4, 4)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()
