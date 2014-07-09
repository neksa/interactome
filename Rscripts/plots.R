
library(ggplot2)
library(data.table)
f <- fread("matches_human.processed.tab", header=T)
setkey(f, "complex_type")

# image(as.matrix(table(f$id, f$interface_id)))
# image(as.matrix(log(table(f$id, f$interface_id))))
# cor(f$id, f$interface_id)

png("bs_score1_homo.png")
hist(f["Homo"]$bs_score1, main="Score1, Homo")
dev.off()

png("bs_score1_hetero.png")
hist(f["Hetero"]$bs_score1, main="Score1, Hetero")
dev.off()


# p <- qplot(identity, data=f, geom="histogram")
# ggsave("identity.png")

# p <- qplot(bs_len, data=f, geom="histogram")
# ggsave("bs_len.png")

# p <- qplot(bs_err, data=f, geom="histogram")
# ggsave("bs_err.png")

# p <- qplot(bs_covered, data=f, geom="histogram")
# ggsave("bs_covered.png")

# p <- qplot(bs_aligned, data=f, geom="histogram")
# ggsave("bs_aligned.png")

# p <- qplot(bs_identical, data=f, geom="histogram")
# ggsave("bs_identical.png")

# # ntempl <- fread("ntempl.tab", header=T)