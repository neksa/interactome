#!/usr/bin/env Rscript --vanilla

library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)

f <- fread("matches_human_033bs_cov.tab", header=TRUE, drop=c(12,13))
f$template_type <- as.factor(f$template_type)
f$query_type <- as.factor(f$query_type)

# with(f, by(template_type, CS1, summary))
# with(f, by(template_type, CS2, summary))
# with(f, by(query_type, CS1, summary))
# with(f, by(query_type, CS2, summary))


bg <- data.frame(score = f$CS1)
fg <- data.frame(score = f$CS2)

#Now, combine your two dataframes into one.  First make a new column in each.
bg$source <- 'Template'
fg$source <- 'Model'

#and combine into your new data frame vegLengths
L <- rbind(bg, fg)

ggplot(L, aes(score, fill = source)) + geom_density(alpha = 0.2)
ggplot(L, aes(score, fill = source)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')


ggplot(d, aes(V1, V2, color = filter)) + geom_point()
