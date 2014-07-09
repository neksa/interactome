
library(data.table)

f <- read.table("real_scores_bin47.tab")
png("plot_47_nat.png")
hist(f$V6, main="0.5A bin, 2-30A\n(natual interfaces)", xlab="score", breaks = 50)
dev.off()

f <- read.table("shuffl_scores_bin47.tab")
png("plot_47_back.png")
hist(f$V6, main="0.5A bin, 2-30A\n (shuffl. interfaces)", xlab="score",  breaks = 50)
dev.off()

f <-read.table("real_scores_bin20.tab")
png("plot_20_nat.png")
hist(f$V6, main="bin47 (natural interfaces)", xlab="score",  breaks = 50)
dev.off()

f <- read.table("shuffl_scores_bin20.tab")
png("plot_20_back.png")
hist(f$V6, main="bin47 (shuffl. interfaces)", xlab="score",  breaks = 50)
dev.off()
