#!/usr/bin/env Rscript --vanilla
#
library(MASS)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)

f <- fread("/Users/agoncear/projects/Interactome/PREDRES_DUMP.tab", header=TRUE)
setkey(f, label)

f$complex_type <- as.factor(f$complex_type)
levels(f$complex_type) = c("Hetero", "Homo")

# pairs(~ model_minus_avg + bs_similarity + bs_aln_similarity + full_seq_id + complex_type + protein_size + bs_size, data=f)

r <- lda(formula = label ~ model_minus_avg + bs_similarity + bs_aln_similarity + full_seq_id + protein_size + bs_size, data = f)


lda_predict <- predict(r, f)
lda_label <- lda_predict$class
round(table(f$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f)*100, 2)

r.reduced <- lda(formula = label ~ model_minus_avg + bs_similarity, data = f)


lda_predict.reduced <- predict(r.reduced, f)
lda_label.reduced <- lda_predict.reduced$class
round(table(f$label, lda_label.reduced, dnn=c("Structure", "LDA classifier"))/nrow(f)*100, 2)


############################
f.subset <- subset(f, bs_similarity >= 0.4 & model_minus_avg >= 5.0)

table(f.subset$label)
table(f.subset$complex_type)

r <- lda(formula = label ~ model_minus_avg + bs_similarity + bs_aln_similarity + full_seq_id + protein_size + bs_size, data = f.subset)
lda_predict <- predict(r, f.subset)
lda_label <- lda_predict$class
# round(table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f.subset)*100, 2)
table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))

r <- lda(formula = label ~ model_minus_avg + bs_similarity + bs_aln_similarity + full_seq_id + protein_size + bs_size, data = f.subset)
lda_predict <- predict(r, f.subset)
lda_label <- lda_predict$class
round(table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f.subset)*100, 2)


r <- lda(formula = label ~ model_minus_avg + bs_similarity + complex_type, data = f.subset)
lda_predict <- predict(r, f.subset)
lda_label <- lda_predict$class
round(table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f)*100, 2)


f.subset$bs_size_norm <- (f.subset$bs_size - mean(f.subset$bs_size)) / sd(f.subset$bs_size)
f.subset$protein_size_norm <- (f.subset$protein_size - f.subset$protein_size) / sd(f.subset$protein_size)
r <- lda(formula = label ~ model_minus_avg + bs_similarity , data = f.subset)
lda_predict <- predict(r, f.subset)
lda_label <- lda_predict$class
round(table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f.subset)*100, 2)

r <- qda(formula = label ~ model_minus_avg + bs_similarity , data = f.subset)
lda_predict <- predict(r, f.subset)
lda_label <- lda_predict$class
# round(table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))/nrow(f.subset)*100, 2)
table(f.subset$label, lda_label, dnn=c("Structure", "LDA classifier"))



#     0     1
# 65290 22382

#  Hetero    Homo Unknown
#   51497   36175       0

#   Hetero  Homo Unknown
# 0  38218 27072       0
# 1  13279  9103       0

