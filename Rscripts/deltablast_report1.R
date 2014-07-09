#!/usr/local/bin/Rscript
# 
# Plot reports on Deltablast for human dimers
#

require(data.table)
require(ggplot2)
library(plyr)

blast <- fread("deltablast_report.tab")
# setnames(distances, c("V1", "V2", "V3"), c("aa", "dHeavy", "dCalpha"))
setkey(blast, "query")


# query               hit                score          bit_score           evalue             identity          positive            gaps
#  Length:3320786     Length:3320786     Min.   :  48.0   Min.   :  22.78   Min.   : 0.000000   Min.   :   0.00   Min.   :   2.00   Min.   :  0.00
#  Class :character   Class :character   1st Qu.:  95.0   1st Qu.:  41.09   1st Qu.: 0.000000   1st Qu.:  17.00   1st Qu.:  31.00   1st Qu.:  3.00
#  Mode  :character   Mode  :character   Median : 230.0   Median :  93.22   Median : 0.000000   Median :  37.00   Median :  65.00   Median : 10.00
#                                        Mean   : 311.6   Mean   : 124.40   Mean   : 0.387498   Mean   :  54.77   Mean   :  80.66   Mean   : 14.86
#                                        3rd Qu.: 494.0   3rd Qu.: 194.72   3rd Qu.: 0.000221   3rd Qu.:  71.00   3rd Qu.: 118.00   3rd Qu.: 20.00
#                                        Max.   :4023.0   Max.   :1553.91   Max.   : 9.999880   Max.   :1068.00   Max.   :1068.00   Max.   :191.00
#    align_len
#  Min.   :   8.0
#  1st Qu.:  92.0
#  Median : 165.0
#  Mean   : 169.2
#  3rd Qu.: 230.0
#  Max.   :1096.0

g <- ggplot(blast, aes(x=align_len)) +
    geom_histogram() +
    theme_bw()
ggsave("align_len_hist.png")

g <- ggplot(blast, aes(x=identity, y=evalue)) +
    geom_point() +
    scale_y_log10() +
    theme_bw()
ggsave("identity_evalue.png")


g <- ggplot(blast, aes(x=identity/align_len, y=evalue)) +
    geom_point() +
    scale_y_log10() +
    theme_bw()
ggsave("percent_identity_evalue.png")



blast_summary <- ddply(blast, .(query), summarize, 
    hits = length(unique(hit)))
    # evalue = min(evalue),
    # identity = max(identity)


g <- ggplot(blast_summary, aes(x=hits)) +
    geom_histogram() +
    theme_bw()
ggsave("num_hits_hist.png")

