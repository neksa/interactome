#!/usr/bin/env Rscript --vanilla

library(ggplot2)
library(data.table)
library(reshape2)

f <- fread("/Users/agoncear/projects/Interactome/Workflow/Potential/potential_1.index-warp", header=F, na.strings=NULL, drop=17:21)
setnames(f, c("aa", "1.5-2.5", "2.5-3.5", "3.5-4.5", "4.5-5.5", "5.5-6.5", "6.5-7.5", "7.5-8.5", "8.5-9.5", "9.5-10.5", "10.5-11.5",
              "11.5-12.5", "12.5-13.5", "13.5-14.5", "14.5-15.5", "15.5-16.5"))
setkey(f, "aa")

mf <- melt(f, variable.name = "dCalpha", value.name="value", id.vars="aa")
mf$dCalpha <- as.factor(mf$dCalpha)
setkey(mf, "aa")

pl <- function(a) {
    png(paste("potential_1/pairs/", a, ".png", sep=""))
    barplot(mf[a]$value, ylim=c(-2, +2))
    dev.off()
}

dist <- function(d) {
    png(paste("potential_1/distances/", d, ".png", sep=""))
    hist(mf[dCalpha == d]$value, xlim=c(-2, +2), main=d)
    dev.off()
}

apply(as.array(mf$aa), 1, pl)
apply(as.array(levels(mf$dCalpha)), 1, dist)


# d <- fread("/Users/agoncear/projects/Interactome/Workflow/Potential/distance_stats.tab")
# # setnames(d, c("pdb", "chain", "resn", "resi", "sCa", "sCb", "svCb", "squaternion"))
# setnames(d, c("res", "dCa", "dCb", "dvCb", "cos_theta", "theta", "dtheta", "cos_omega", "omega", "domega"))

# q <- fread("/Users/agoncear/projects/Interactome/Workflow/Potential/distance_stats_shuffled.tab")
# setnames(q, c("res", "dCa", "dCb", "dvCb", "cos_theta", "theta", "dtheta", "cos_omega", "omega", "domega"))

