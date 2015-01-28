#!/usr/local/bin/Rscript
# 
# Template statistics

require(data.table)
require(ggplot2)
require("gridExtra")

tpl <- fread("/Users/agoncear/projects/Interactome/Workflow/Structures/template_analysis.tab", header=TRUE)
# setnames(distances, c("V1", "V2", "V3"), c("aa", "dHeavy", "dCalpha"))
setkey(tpl, "template")
tpl$complex_type <- as.factor(tpl$complex_type)
tpl$redundant <- as.factor(tpl$redundant)

table(tpl$redundant)
table(tpl$complex_type)

# template\tpdb\tA\tB\tprot_A\tprot_B\tcomplex_type\tnsubunits\tndirect_int\tredundant\tbs_lenA\tncontactsA\tbs_lenB\tncontactsB
p1 <- ggplot(aes(x = nsubunits, y = ndirect_int, color=complex_type), data = tpl) +
  geom_point(alpha=0.5, size=0.8) +
  # xlim(0, 20) + ylim(0, 5) +
  labs(title="Subunits in dimer template library", x="total subunits", y="directly interacting subunits") +
  theme_classic()

p2 <- ggplot(tpl, aes(nsubunits)) + 
    # guides(colour=FALSE) +
    # geom_density(size=0.2, alpha=0.5) +
    geom_histogram(binwidth=1, colour="black", fill="white") +
    labs(title="Distribution of total number of subunits") +
    theme_classic()

p3 <- ggplot(tpl, aes(ndirect_int)) + 
    # guides(colour=FALSE) +
    # geom_density(size=0.2, alpha=0.5) +
    geom_histogram(binwidth=1, colour="black", fill="white") +
    labs(title="Distribution of directly interacting subunits") +
    theme_classic()

p4 <- ggplot(tpl, aes(x=bs_lenA, y=bs_lenB, color=complex_type)) + 
      geom_point(alpha=0.5, size=0.8) +
      labs(title="Lengths of binding sites in residues", x="site A", y="site B") +
    theme_classic()


p5 <- ggplot(tpl, aes(x=ncontactsA, y=ncontactsB, color=complex_type)) + 
      geom_point(alpha=0.5, size=0.8) +
      labs(title="Total number of contacts in binding sites", x="site A", y="site B") +
    theme_classic()


multiplot <- arrangeGrob(p1, p2, p3, p4, ncol=2, nrow=2)
ggsave(filename="templates.png", plot=multiplot, width=12, height=10)


# # Multiple plot function
# #
# # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# # - cols:   Number of columns in layout
# # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# #
# # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# # then plot 1 will go in the upper left, 2 will go in the upper right, and
# # 3 will go all the way across the bottom.
# #
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   require(grid)

#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)

#   numPlots = length(plots)

#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                     ncol = cols, nrow = ceiling(numPlots/cols))
#   }

#  if (numPlots==1) {
#     print(plots[[1]])

#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

# # multiplot(p1, p2, p3, p4)

# # ggsave("template_analysis.png") #, width=7, height=7)

