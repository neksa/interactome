#!/usr/bin/env Rscript --vanilla
#
# Analyze 
# Queries per template
# Templates per query
# 

require(ggplot2)
require(data.table)
require(dplyr)

theme_set(theme_classic(base_size = 18)) 

fname  <- "/Users/agoncear/projects/Interactome/Workflow/BLAST-results/Hsapiens_merged.tab"

f <- fread(fname, header=T)
setkey(f, "query")
templates_per_query <- f %>% 
    group_by(query) %>%
    summarise(N_templates = n_distinct(template))

setkey(f, "template")
queries_per_template <- f %>% 
    group_by(template) %>%
    summarise(N_queries = n_distinct(query))
# summarise_each(funs(mean, max, sum, n_distinct), aa)

f <- mutate(f, pdb = substr(template, 1, 4))
setkey(f, "pdb")
queries_per_pdb <- f %>% 
    group_by(pdb) %>%
    summarise(N_queries = n_distinct(query))

p <- ggplot(templates_per_query, aes(x = N_templates)) +
    geom_histogram(colour = "black", fill = "white", binwidth = 50) +
    ylab('Queries') + xlab('Templates / Query')
ggsave("human_hits/templates_per_query_human.png", dpi=150)

p <- ggplot(queries_per_template, aes(x = N_queries)) +
    geom_histogram(colour = "black", fill = "white", binwidth = 50) +
    ylab('Templates') + xlab('Queries / Template')
ggsave("human_hits/queries_per_template_human.png", dpi=150)


p <- ggplot(queries_per_pdb, aes(x = N_queries)) +
    geom_histogram(colour = "black", fill = "white", binwidth = 50) +
    ylab('PDBs') + xlab('Queries / PDB')
ggsave("human_hits/queries_per_pdb_human.png", dpi=150)

write.csv(select(filter(queries_per_pdb, N_queries > 200), pdb, N_queries), sep="\t", file="PDB_list.txt")

# image(as.matrix(table(f$id, f$interface_id)))
# image(as.matrix(log(table(f$id, f$interface_id))))
# cor(f$id, f$interface_id)

# png("bs_score1_homo.png")
# hist(f["Homo"]$bs_score1, main="Score1, Homo")
# dev.off()

# png("bs_score1_hetero.png")
# hist(f["Hetero"]$bs_score1, main="Score1, Hetero")
# dev.off()


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