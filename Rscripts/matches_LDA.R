#!/usr/bin/env Rscript --vanilla
#
# LDA analysis for Ecoli matches
#
# Finding the best combination of factors

# source("../../Rscripts/matches_LDA.R")

library(MASS)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)

f <- fread("/Users/agoncear/projects/Interactome/Workflow/Alignments/labeled_matches_ecoli.tab", header=TRUE)
setkey(f, label)

f$template_type <- as.factor(f$template_type)
f$query_type <- as.factor(f$query_type)
levels(f$query_type) = c("Hetero", "Homo", "Unknown")

f <- mutate(f, identityA=as.double(identicalA) / as.double(aln_lenA))
f <- mutate(f, identityB=as.double(identicalB) / as.double(aln_lenB))

f <- mutate(f, similarityA=as.double(positiveA) / as.double(aln_lenA))
f <- mutate(f, similarityB=as.double(positiveB) / as.double(aln_lenB))

f <- mutate(f, bs_identityA=as.double(bs_identicalA) / as.double(bs_lenA))
f <- mutate(f, bs_identityB=as.double(bs_identicalB) / as.double(bs_lenB))

f <- mutate(f, bs_similarityA=as.double(bs_positiveA) / as.double(bs_lenA))
f <- mutate(f, bs_similarityB=as.double(bs_positiveB) / as.double(bs_lenB))

f <- mutate(f, bs_coverageA=as.double(bs_alignedA) / as.double(bs_lenA))
f <- mutate(f, bs_coverageB=as.double(bs_alignedB) / as.double(bs_lenB))

f <- mutate(f, query_as_template= (template_type == query_type))

# see also base::pmin
f <- f %>% 
  rowwise() %>% 
  mutate(
    min_identity=min(identityA, identityB),
    min_similarity=min(similarityA, similarityB),
    min_bs_identity=min(bs_identityA, bs_identityB),
    min_bs_similarity=min(bs_similarityA, bs_similarityB),
    min_bs_coverage=min(bs_coverageA, bs_coverageB),
    min_bs_BLOSUM=min(bs_BLOSUMA, bs_BLOSUMB),
    min_bs_score1=min(bs_score1A, bs_score1B))

f <- f %>% 
  rowwise() %>% 
  mutate(
    min_len=min(aln_lenA, aln_lenB),
    min_bs_len=min(bs_lenA, bs_lenB)
    )


# ggplot(L, aes(score, fill = source)) + geom_density(alpha = 0.2)
# ggplot(L, aes(score, fill = source)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
# ggplot(d, aes(V1, V2, color = filter)) + geom_point()

# r <- lda(formula = label ~ ., data = f )#, prior = c(1,1,1)/3)
s <- sample_n(f, 100000)
table(s$label)

r <- lda(
    formula = label ~ score + zscore + score_minus_avg + min_len + min_bs_len + min_identity + min_similarity + min_bs_identity + min_bs_similarity + min_bs_coverage + min_bs_BLOSUM + min_bs_score1 + query_as_template,
    data = f)


r.reduced1 <- lda(
    formula = label ~ score + zscore + score_minus_avg + min_bs_len +  min_bs_identity + min_bs_similarity + min_bs_coverage + query_as_template,
    data = f)

r.reduced2 <- lda(
    formula = label ~ score_minus_avg + min_bs_similarity,
    data = f)


m <- lm(label ~ score + zscore + score_minus_avg + min_len + min_bs_len + min_identity + min_similarity + min_bs_identity + min_bs_similarity + min_bs_coverage + min_bs_BLOSUM + min_bs_score1 + query_as_template,
    data = f)

m.reduced1 <- lm(label ~ score + zscore + score_minus_avg + min_bs_len +  min_bs_identity + min_bs_similarity + min_bs_coverage + query_as_template,
    data = f)
m.reduced2 <- lm(label ~ score_minus_avg + min_bs_similarity, data = f)


s <- sample_n(f[f$label==FALSE,], 10000)
s <- rbind(s, f[f$label==TRUE,])
rs <- lda(
    formula = label ~ score + zscore + score_minus_avg + min_len + min_bs_len + min_identity + min_similarity + min_bs_identity + min_bs_similarity + min_bs_coverage + min_bs_BLOSUM + min_bs_score1 + query_as_template,
    data = s,CV=TRUE)
x <- predict(rs, s)
table(s$label, x$class)




rs.score <- lda(label ~ min_bs_similarity + score_minus_avg, data=s, )
x.score <- predict(rs.score, s)
table(s$label, x.score$class)
plot(rs.score, asp=1)


m <- glm(label ~ score + zscore + score_minus_avg + min_len + min_bs_len + min_identity + min_similarity + min_bs_identity + min_bs_similarity + min_bs_coverage + min_bs_BLOSUM + min_bs_score1 + query_as_template,
family=binomial(logit),
data = s)


m <- glm(label ~ score_minus_avg * min_bs_len + min_identity + min_bs_identity + min_bs_similarity ,
family=binomial(logit),
data = s)

