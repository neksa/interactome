library(ggplot2)
library(GGally)
library(data.table)

d <- fread("/Users/agoncear/projects/Interactome/Workflow/Alignments/matches_human_25_analysis.tab", header=TRUE)
setkey(d, "complex_type")

# 
# y <- cbind(
#         # complex=d$complex_type,
#         id=d$identity,
#         bid=d$bs_identical/d$bs_len,
#         baln=d$bs_aligned/d$bs_len,
#         score1=d$bs_score1,
#         blosum=d$bs_BLOSUM)

y <- data.table(
        # complex=d$complex_type,
        id=d$identity,
        bid=d$bs_identical/d$bs_len,
        baln=d$bs_aligned/d$bs_len,
        score1=d$bs_score1,
        blosum=d$bs_BLOSUM)






ntempl = fread("/Users/agoncear/projects/Interactome/Workflow/Alignments/ntempl.tab", header=TRUE)

png("ntempl.png")
hist(ntempl$ntempl[ntempl$ntempl < 600], main="Number of templates per query", breaks=50)
dev.off()

# png("bs_score1_homo.png")
# hist(f["Homo"]$bs_score1, main="Score1, Homo")
# dev.off()


png("bs_score1_homo.png")
hist(d["Homo"]$bs_score1, main="Score1, Homo", breaks=50)
dev.off()

png("bs_score1_hetero.png")
hist(d["Hetero"]$bs_score1, main="Score1, Hetero", breaks=50)
dev.off()




png("bs_score1_zoom_homo.png")
hist(d["Homo"]$bs_score1[abs(d["Homo"]$bs_score1)<10000], main="Score1, zoom, Homo", breaks=50)
dev.off()

png("bs_score1_zoom_hetero.png")
hist(d["Hetero"]$bs_score1[abs(d["Hetero"]$bs_score1)<10000], main="Score1, zoom, Hetero", breaks=50)
dev.off()




png("bs_blosum_homo.png")
hist(d["Homo"]$bs_BLOSUM, main="BLOSUM, Homo", breaks=50)
dev.off()

png("bs_blosum_hetero.png")
hist(d["Hetero"]$bs_BLOSUM, main="BLOSUM, Hetero", breaks=50)
dev.off()


png("bs_aln_homo.png")
hist(d["Homo"]$bs_aligned / d["Homo"]$bs_len, main="Aligned interface, Homo", breaks=50)
dev.off()


png("bs_aln_hetero.png")
hist(d["Hetero"]$bs_aligned / d["Hetero"]$bs_len, main="Aligned interface, Hetero", breaks=50)
dev.off()



png("bs_ident_homo.png")
hist(d["Homo"]$bs_identical / d["Homo"]$bs_len, main="Identity on the interface, Homo", breaks=50)
dev.off()


png("bs_ident_hetero.png")
hist(d["Hetero"]$bs_identical / d["Hetero"]$bs_len, main="Identity on the interface, Hetero", breaks=50)
dev.off()




png("ident_homo.png")
hist(d["Homo"]$identity, main="Identity of the hit, Homo", breaks=50)
dev.off()


png("ident_hetero.png")
hist(d["Hetero"]$identity, main="Identity of the hit, Hetero", breaks=50)
dev.off()



png("bs_alignment_homo.png")
hist(d["Homo"]$bs_aligned, main="Binding site alignment length, Homo")
dev.off()


png("bs_alignment_hetero.png")
hist(d["Hetero"]$bs_aligned, main="Binding site alignment length, Hetero")
dev.off()


tpl <- fread("../template_analysis.tab", header=T)

summary(as.factor(tpl$complex_type))



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


