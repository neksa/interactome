from itertools import islice
from collections import defaultdict
from copy import deepcopy

def best_pairs(matches_fname, results_fname, scoring=1):
    pairs = defaultdict(list)
    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            queryA, queryB, tpl, query_type, template_type, SC1, SC2, SC3, SC4, SC5, SC6, siteA, siteB = line.strip().split("\t")
            SC1, SC2, SC3, SC4, SC5, SC6 = map(float, (SC1, SC2, SC3, SC4, SC5, SC6))
            
            # scoring == 1
            score = SC2
            if scoring == 2:
                score = SC5*SC2
            if scoring == 3:
                score = SC5*(2*SC2/20.0 + (SC1 - SC2)/20.0 + SC3/700.0 + SC4/100000.0) / 4.0
            if scoring == 4:
                score = (SC1 - SC2)/20.0

            pair = (queryA, queryB)

            # if scoring == 0:
            #     if query_type != "Hetero" or template_type != "Hetero":
            #         continue
            #     if float(SC5) < 0.4:
            #         continue
            #     if SC4 < 0:
            #         continue
            #     if SC1 - SC2 < 0:
            #         continue

            pairs[pair].append(score)

    with open(results_fname, 'w') as o:
        for pair, scores in pairs.iteritems():
            score = max(scores)
            o.write("{}\t{}\t{}\n".format(pair[0], pair[1], score))

    return pairs


def load_HPA_filter_results(fname, strict=False):
    pairs = {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            A, B, protein, RNA, protein_locs, RNA_locs = line.strip().split("\t")
            pair = (A, B)
            outcome = False
            if strict:
                if protein == "true" or RNA == "true":
                    outcome = True
            else:
                outcome = True
                if protein == "false" and RNA == "false":
                    outcome = False
            pairs[pair] = outcome
    return pairs


def load_gene_protein_mapping(mapping_fname):
    gene_proteins = defaultdict(set)
    protein_genes = defaultdict(set)
    with open(mapping_fname) as f:
        for line in f:
            fields = line.strip().split()
            # print fields
            if len(fields) != 2:
                continue
            gene = int(fields[0])
            protein = fields[1]
            gene_proteins[gene].add(protein)
            protein_genes[protein].add(gene)
    return gene_proteins, protein_genes


def load_interactions(fname):
    interactions = []
    with open(fname) as f:
        for line in islice(f, 1, None):
            gA, nameA, gB, nameB = line.strip().split("\t")
            # print gA, nameA, gB, nameB
            gA = int(gA)
            gB = int(gB)
            pair = (gA, gB)
            interactions.append(pair)
    return interactions


def filtered_pairs(pairs, HPA_results, strict=False):
    new_pairs = deepcopy(pairs)
    for pair in pairs:
        outcome = True
        if strict:
            outcome = HPA_results.get(pair, False)
        else:
            outcome = HPA_results.get(pair, True)
        if not outcome:
            del new_pairs[pair]
    return new_pairs


def inverse_scored_labels(predicted, observed, gene_proteins, output):
    with open(output, 'w') as o:
        for (A, B) in observed.iteritems():
            proteinsA = gene_proteins[A]
            proteinsB = gene_proteins[B]
            # print (pA, pB), genesA, genesB
            scores = []
            for pA in proteinsA:
                for pB in proteinsB:
                    scores.extend(predicted.get((pA, pB), []))
                    scores.extend(predicted.get((pB, pA), []))
            score = max(scores)
            label = 0 if score is None else 1
            if score is None:
                score = -1e6
            o.write("{}\t{}\n".format(score, label))


def scored_labels(predicted, observed, protein_genes, output):
    with open(output, 'w') as o:
        for (pA, pB), scores in predicted.iteritems():
            score = max(scores)
            genesA = protein_genes[pA]
            genesB = protein_genes[pB]
            # print (pA, pB), genesA, genesB
            outcome = None
            for gA in genesA:
                for gB in genesB:
                    if (gA, gB) in observed or (gB, gA) in observed:
                        outcome = True
                    if outcome is not True:
                        outcome = False
            if outcome is None:
                continue
            label = 1 if outcome else 0
            o.write("{}\t{}\n".format(score, label))


if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    matches_fname = d + "/Alignments/matches_human_033bs_cov.tab"
    pairs_fname = d + "/Alignments/human_033bs_cov.tab"
    HPA_fname = d + "/Coexpression/human_033bs_cov_HPA.tab"
    vidal_fname = "/Users/agoncear/data/Vidal/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"
    # vs Vidal
    benchmark1_fname = d + "/Benchmarks/benchmark1.scored-label"
    # filtered vs Vidal
    benchmark2_fname = d + "/Benchmarks/benchmark2.scored-label"
    # alternative scoring vs Vidal
    benchmark3_fname = d + "/Benchmarks/benchmark3.scored-label"
    # alternative scoring, filtered vs Vidal
    benchmark4_fname = d + "/Benchmarks/benchmark4.scored-label"
    benchmark5_fname = d + "/Benchmarks/benchmark5.scored-label"

    benchmark_inv = d + "/Benchmarks/vidal_coverage.txt"

    gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    # print protein_genes
    
    HPA_strict_results = load_HPA_filter_results(HPA_fname, strict=True)
    vidal_interactome = load_interactions(vidal_fname)
    print "V:", len(vidal_interactome)

    # pairs_scoring_0 = best_pairs(matches_fname, pairs_fname, scoring=0)
    # scored_labels(pairs_scoring_0, vidal_interactome, protein_genes, d + "/Benchmarks/bench_score0.txt")


    pairs_scoring_II = best_pairs(matches_fname, pairs_fname, scoring=2)
    """
    pairs_scoring_I = best_pairs(matches_fname, pairs_fname, scoring=1)
    pairs_scoring_III = best_pairs(matches_fname, pairs_fname, scoring=3)
    # run HPAfilter java
    HPA_results = load_HPA_filter_results(HPA_fname, strict=False)
    HPA_strict_results = load_HPA_filter_results(HPA_fname, strict=True)

    # scored_labels(pairs_scoring_I, vidal_interactome, protein_genes, benchmark1_fname)

    # pairs_filtered = filtered_pairs(pairs_scoring_I, HPA_results)
    # scored_labels(pairs_filtered, vidal_interactome, protein_genes, benchmark2_fname)

    # pairs_strict_filtered = filtered_pairs(pairs_scoring_I, HPA_strict_results)
    # scored_labels(pairs_strict_filtered, vidal_interactome, protein_genes, benchmark3_fname)
    """

    pairs_strict_filtered_sII = filtered_pairs(pairs_scoring_II, HPA_strict_results)
    scored_labels(pairs_scoring_II, vidal_interactome, protein_genes, benchmark1_fname)
    scored_labels(pairs_strict_filtered_sII, vidal_interactome, protein_genes, benchmark2_fname)
    # inverse_scored_labels(pairs_scoring_II, vidal_interactome, gene_proteins, benchmark_inv)

    """"
    pairs_strict_filtered_sIII = filtered_pairs(pairs_scoring_III, HPA_strict_results)
    scored_labels(pairs_strict_filtered_sIII, vidal_interactome, protein_genes, benchmark5_fname)
    """
