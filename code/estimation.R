source("code/sample.R")

require(RColorBrewer)
palette(brewer.pal(8, "Dark2"))

pdf("plots/estimation-frob-s0-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 0,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 0.2, s == 0)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-s0-p100.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 0,
        psamp = 1.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 1, s == 0)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-s0-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 0,
        psamp = 2.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 2, s == 0)),
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 0.2, s == 4)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-p100.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 1,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 1, s == 4)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 2,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 2, s == 4)),
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-s32-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 32,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 0.2, s == 32)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-s32-p100.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 32,
        psamp = 1.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 1, s == 32)),
        type="frobenius")
dev.off()
pdf("plots/estimation-frob-s32-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 32,
        psamp = 2.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
	label=expression(list(m/mu == 2, s == 32)),
        type="frobenius")
dev.off()


pdf("plots/estimation-frob-p010.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p050.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.50,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()


pdf("plots/estimation-frob-p500.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 5,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-frob-p1000.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-spec-p010.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p050.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 0.50,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p100.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 1,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="frobenius")
dev.off()

pdf("plots/estimation-spec-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 2,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p500.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 5,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-p1000.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 4,
        psamp = 10,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()



pdf("plots/estimation-spec-s32-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 32,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()
pdf("plots/estimation-spec-s32-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 32,
        psamp = 2.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

pdf("plots/estimation-spec-s0-p020.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 0,
        psamp = 0.20,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()
pdf("plots/estimation-spec-s0-p200.pdf", 5, 5)
repdemo(width = 6, height = 7, wrap = FALSE,
        nswaps = 0,
        psamp = 2.00,
        penalty = seq(0.1, 200, length.out = 100),
        nreps = 100,
        type="spectral")
dev.off()

