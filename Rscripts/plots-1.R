#!/usr/local/bin/Rscript
# 
# Plot panels of distance-dependent statistics
#

require(data.table)
require(ggplot2)

distances <- fread("aa-distances.tab")
setnames(distances, c("V1", "V2", "V3"), c("aa", "dHeavy", "dCalpha"))
setkey(distances, "aa")

# aa dHeavy dCalpha
# 
aa_list <- c("ALA", "LEU", "PRO", "GLY", "ASP", "ASN", "TYR", "HIS", "GLU", "CYS", "PHE", "VAL", "ILE", "ARG", "THR", "LYS", "SER", "GLN", "MET", "TRP")

if (FALSE) {
    build_bg_plot <- function(amino_acid, data) {
        p <- ggplot(aes(x = dCalpha, y = dHeavy), data = data[grep(amino_acid, aa)]) +
            geom_point(alpha=0.3, size=0.8) +
            xlim(0, 20) + ylim(0, 5) +
            labs(title=paste(amino_acid,"with ANYTHING")) +
            theme_classic()
        ggsave(paste("plots/distance_AA_with_anything/", amino_acid, ".png", sep=""), width=5, height=5)
        print(amino_acid)
    }

    lapply(aa_list, build_bg_plot, distances[dCalpha<20])

    build_plot <- function(amino_acids, data) {
        p <- ggplot(aes(x = dCalpha, y = dHeavy), data = data[amino_acids]) +
          geom_point(alpha=0.3, size=0.8) +
          xlim(0, 20) + ylim(0, 5) +
          labs(title=amino_acids) +
          theme_classic()
        ggsave(paste("plots/distance_AA_pairs/", amino_acids, ".png", sep=""), width=5, height=5)
        print(amino_acids)
    }

    combinations <- levels(as.factor(distances$aa))
    lapply(combinations, build_plot, distances[dCalpha<20])
}



# ggplot(aes(x = dCalpha, y = dHeavy), data = distances[dCalpha<20][grep("ASP", aa)]) +
#     geom_point(alpha=0.3, size=0.8) +
#     facet_wrap(~ aa) +
#     theme_classic()


# ggplot(aes(x = dCalpha, y = dHeavy), data = distances[dCalpha<20][grep("ASP", aa)]) +
#     geom_point(alpha=0.3, size=0.8) +
#     facet_wrap(~ aa) +
#     theme_classic()

# ggplot(aes(x = dCalpha, y = dHeavy), data = distances[dCalpha<20]["CYS_CYS"]) +
#     geom_point(alpha=0.3, size=0.8) +
#     facet_wrap(~ aa) +
#     theme_classic()

# # PHE with any
# ggplot(aes(x = dCalpha, y = dHeavy), data = distances[dCalpha<20][grep("PHE", aa)]) + geom_point(alpha=0.4, size=0.8) + theme_classic()

#     facet_wrap(~ aa) + geom_point(alpha=0.2)

# ggplot(aes(x = dCalpha, y = dHeavy), data = distances) +
#     facet_wrap(~ aa) + geom_point(alpha=0.2)


contacts <- fread("aa-contacts.tab")
setnames(contacts, c("V1", "V2", "V3"), c("aa", "Ncontacts", "dCalpha"))
setkey(contacts, "aa")

# contacts.agr <- contacts[, mean(Ncontacts, na.rm = TRUE), by=aa]
contacts.agr <- contacts[, list(
    Contacts_mean = mean(Ncontacts, na.rm = TRUE),
    Contacts_sd=sd(Ncontacts, na.rm = TRUE)),
    by=aa]

p <- ggplot(contacts[Ncontacts<25], aes(Ncontacts, colour=aa)) + 
    guides(colour=FALSE) +
    geom_density(size=0.2, alpha=0.5) +
    labs(title="Distribution of #contacts for all amino acid pairs") +
    theme_classic()
ggsave("plots/Num_contacts_histogram.png", width=7, height=7)


######## 
if (FALSE) {
    build_bg_plot_contacts <- function(amino_acid, data) {
        p <- ggplot(aes(x = dCalpha, y = Ncontacts), data = data[grep(amino_acid, aa)]) +
            geom_point(alpha=0.5, size=1) +
            labs(title=paste(amino_acid,"with ANYTHING")) +
            xlim(0, 20) + ylim(0, 70) +
            theme_classic()
        ggsave(paste("plots/Ncontacts_AA_with_anything/", amino_acid, ".png", sep=""), width=5, height=5)
        print(amino_acid)
    }

    lapply(aa_list, build_bg_plot_contacts, contacts[dCalpha<20])

    build_plot_contacts <- function(amino_acids, data) {
        p <- ggplot(aes(x = dCalpha, y = Ncontacts), data = data[amino_acids]) +
          geom_point(alpha=0.5, size=1) +
          labs(title=amino_acids) +
          xlim(0, 20) + ylim(0, 70) +
          theme_classic()
        ggsave(paste("plots/Ncontacts_AA_pairs/", amino_acids, ".png", sep=""), width=5, height=5)
        print(amino_acids)
    }

    combinations_contacts <- levels(as.factor(contacts$aa))
    lapply(combinations_contacts, build_plot_contacts, contacts[dCalpha<20])
}

###################

if (FALSE) {
    dist <- fread("asp_lys-distances.tab")
    setnames(dist, c("V1", "V2", "V3", "V4"), c("aa", "atoms", "dHeavy", "dCalpha"))
    setkey(dist, "atoms")


    build_plot2 <- function(atoms, data) {
        p <- ggplot(aes(x = dCalpha, y = dHeavy), data = data[atoms]) +
          geom_point(alpha=0.3, size=0.8) +
          xlim(0, 20) + ylim(0, 5) +
          labs(title=atoms) +
          theme_classic()
        ggsave(paste("plots/", atoms, ".png", sep=""), width=5, height=5)
        print(atoms)
    }

    build_plot2("O_N", dist)
    build_plot2("CB_CB", dist)
}
####################

binary <- fread("binary_interactions.tab")
setnames(binary, c("V1", "V2", "V3"), c("local", "aa", "dCalpha"))
setkey(binary, "aa")

build_plot3 <- function(aa, data, aloc) {
    filename <- paste("plots/contact_probabilities/", aloc, "/", aa, ".png", sep="")
    if (is.na(data[data$loc == aloc][aa][1]$dCalpha)) {
        print(paste("No data for", aa))
        file.remove(filename)
        return()
    }

    print(aa)
    p <- ggplot(aes(x = dCalpha), data = data[data$loc == aloc][aa]) +
      geom_density() +
      xlim(0, 20) +
      labs(title=paste("Probability density of a contact (HA, 4A)", aa)) +
      theme_classic()
    ggsave(filename, width=5, height=5)
}

combinations <- levels(as.factor(binary$aa))
dev.null <- lapply(combinations, build_plot3, binary, "backbone")
dev.null <- lapply(combinations, build_plot3, binary, "sidechain")
dev.null <- lapply(combinations, build_plot3, binary, "mixed")

# build_plot3("GLY_GLY", binary[loc == "backbone"])


