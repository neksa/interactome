from itertools import islice
from collections import namedtuple, defaultdict
from benchmark_common import get_scored_pairs, roc, multiroc


# def load_gi_uniprot_mapping(fname):
#     mapping = {}
#     with open(fname) as f:
#         for line in f:
#             gis, uniprot = line.strip().split()
#             for gi in gis.split(","):
#                 mapping[gi] = uniprot
#     return mapping


def load_gene_protein_mapping(mapping_fname):
    # gene_proteins = defaultdict(set)
    # protein_genes = defaultdict(set)
    protein_genes = {}
    with open(mapping_fname) as f:
        for line in f:
            fields = line.strip().split()
            # print fields
            if len(fields) != 2:
                continue
            gene = int(fields[0])
            protein = fields[1]
            # gene_proteins[gene].add(protein)
            # protein_genes[protein].add(gene)
            protein_genes[protein] = gene
    # return gene_proteins, protein_genes
    return protein_genes


def load_interactions(fname):
    interactions = set()
    with open(fname) as f:
        for line in islice(f, 1, None):
            if len(line) == 0: continue
            print "LINE", line
            gA, nameA, gB, nameB = line.strip().split("\t")
            # print gA, nameA, gB, nameB
            gA = int(gA)
            gB = int(gB)
            pair = (gA, gB)
            if gA > gB:
                pair = (gB, gA)
            interactions.add(pair)
    return interactions


def scored(predicted_pairs, observed_pairs, mapping):
    true_labels = []
    scores = []

    for unmapped_pair, score in predicted_pairs:
        A, B = unmapped_pair
        # print A, B
        mappedA = mapping.get(A)
        if mappedA is None:
            print "skipping unmapped", A
            continue
        mappedB = mapping.get(B)
        if mappedB is None:
            print "skipping unmapped", B
            continue
        # print A, B, mappedA, mappedB

        label = 0
        if (mappedA, mappedB) in observed_pairs:
            label = 1
        if (mappedB, mappedA) in observed_pairs:
            label = 1
            # print pair, score, label
        # print  label, score
        true_labels.append(label)
        scores.append(score)
    return true_labels, scores


##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    roc_d = d + "/Benchmarks/Hsapiens"
    matches_fname = d + "/Alignments/matches_human.tab"
    # validation_fname = "/Users/agoncear/data/Hpylori/409211A0_S1.txt"
    # mapping_fname = "/Users/agoncear/data/Hpylori/helpy.txt"
    # mapping_fname = "/Users/agoncear/data/Hpylori/M2014110493WE3322RA.tab"

    validation_fname= "/Users/agoncear/data/Vidal/Human/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/Human/unique_gene_ids_mapped.tab"

    # Gene pairs from EBI (Uetz)
    print "Loading Interactome Y2H Human"
    pairs_obs = load_interactions(validation_fname)
    # print pairs_obs
    # mapping = load_swissprot(mapping_fname)
    # mapping = load_mapping(mapping_fname)
    # gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    mapping = load_gene_protein_mapping(mapping_fname)
    # mapping = protein_genes
    # print mapping

    # Uniprot pairs
    print "Processing the interactome"
    pairs, a_pairs = get_scored_pairs(matches_fname)

    print "Comparing observed with predicted"
    r_identity = scored(pairs["identity"], pairs_obs, mapping)
    r_bs_identity = scored(pairs["bs_identity"], pairs_obs, mapping)
    r_positive = scored(pairs["positive"], pairs_obs, mapping)
    r_bs_positive = scored(pairs["bs_positive"], pairs_obs, mapping)
    r_bs_coverage = scored(pairs["bs_coverage"], pairs_obs, mapping)
    r_score = scored(pairs["model_score"], pairs_obs, mapping)
    r_zscore = scored(pairs["model_zscore"], pairs_obs, mapping)
    r_model_minus_avg = scored(pairs["model_minus_avg"], pairs_obs, mapping)
    r_random = scored(pairs["random"], pairs_obs, mapping)

    # print r_identity

    print "Making plots"
    roc(r_identity, roc_d + "/identity.png", "H. sapiens: Alignment identity")
    roc(r_bs_identity, roc_d + "/bs_identity.png", "H. sapiens: BS alignment identity")
    roc(r_positive, roc_d + "/positive.png", "H. sapiens: Alignment BLOSUM62 positive")
    roc(r_bs_positive, roc_d + "/bs_positive.png", "H. sapiens: BS alignment BLOSUM62 positive")
    roc(r_bs_coverage, roc_d + "/bs_coverage.png", "H. sapiens: BS coverage")
    roc(r_score, roc_d + "/bs_score.png", "H. sapiens: Compatibility score 1")
    roc(r_zscore, roc_d + "/bs_zscore.png", "H. sapiens: Compatibility Z-score")
    roc(r_model_minus_avg, roc_d + "/bs_model_minus_avg.png", "H. sapiens: Compatibility score-avg(decoy)")
    roc(r_random, roc_d + "/bs_random.png", "H. sapiens: Random control")

    multiroc(
        (r_identity, r_bs_identity, r_positive, r_bs_positive),
        ('id', 'bs_id', 'pos', 'bs_pos'),
        roc_d + "/id_pos.png",
        "Basic ROC in H. sapiens, Y2H Vidal")

    multiroc(
        (r_bs_positive, r_model_minus_avg),
        ('Binding site sequence similarity', 'Interface compatibily score'),
        roc_d + "/Hsapiens_combined.eps",
        "H. sapiens, Y2H Vidal")
