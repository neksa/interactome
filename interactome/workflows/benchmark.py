from itertools import islice
from collections import defaultdict
from copy import deepcopy

def best_pairs(matches_fname, scoring=1):
    pairs = defaultdict(list)
    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            queryA, queryB, tpl, query_type, template_type, SC1, SC2, SC3, SC4, SC5, SC6, siteA, siteB = line.strip().split("\t")
            SC1, SC2, SC3, SC4, SC5, SC6 = map(float, (SC1, SC2, SC3, SC4, SC5, SC6))
            
            score = SC2  # compatibility score
            if scoring == 4:
                score = SC6

            if scoring == 5:
                score = SC5*SC2

            # if scoring == 3:
            #     score = SC5*(2*SC2/20.0 + (SC1 - SC2)/20.0 + SC3/700.0 + SC4/100000.0) / 4.0
            if scoring == 6:
                score = (SC1 - SC2)

            if scoring == 7:
                score = (SC1 - SC2) / (SC5*(len(siteA.split(";") * len(siteB.split(";")))))

            pair = (queryA, queryB)
            if queryA > queryB:
                pair = (queryB, queryA)

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

    # with open(results_fname, 'w') as o:
    #     for pair, scores in pairs.iteritems():
    #         score = max(scores)
    #         o.write("{}\t{}\t{}\n".format(pair[0], pair[1], score))
    return {pair: max(scores) for pair, scores in pairs.iteritems()}


def load_HPA_filter_results(fname, strict=False):
    pairs = {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            A, B, protein, RNA, protein_locs, RNA_locs = line.strip().split("\t")
            pair = (A, B)
            if A > B:
                pair = (B, A)
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


def load_uniprot_filter(fname):
    ids = set()
    with open(fname) as f:
        for line in islice(f, 1, None):
            uniprot = line.strip().split()[0]
            ids.add(uniprot)
    # for i in ids:
    #     for j in ids:
    #         pair = i, j
    #         if i > j:
    #             pair = j, i
    #         pairs[pair] = False
    return ids


def load_interactions(fname):
    interactions = []
    with open(fname) as f:
        for line in islice(f, 1, None):
            gA, nameA, gB, nameB = line.strip().split("\t")
            # print gA, nameA, gB, nameB
            gA = int(gA)
            gB = int(gB)
            pair = (gA, gB)
            if gA > gB:
                pair = (gB, gA)
            interactions.append(pair)
    return interactions


def filter_pairs(pairs, HPA_results):
    new_pairs = deepcopy(pairs)
    for pair in pairs:
        if not HPA_results.get(pair, False):
            del new_pairs[pair]
    return new_pairs


def filter_pairs_blacklist(pairs, blacklist):
    new_pairs = deepcopy(pairs)
    for x, y in pairs:
        if x in blacklist or y in blacklist:
            del new_pairs[x, y]
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


def scored_labels(predicted, observed, gene_proteins, output_scores, output_network):
    with open(output_scores, 'w') as o, open(output_network, 'w') as o_net:
        gene_pairs = defaultdict(list)

        for (pA, pB), score in predicted.iteritems():
            genesA = protein_genes[pA]
            genesB = protein_genes[pB]
            for gA in genesA:
                for gB in genesB:
                    pair = gA, gB
                    if gA > gB:
                        pair = gB, gA
                    gene_pairs[pair].append(score)

        for (gA, gB), scores in gene_pairs.iteritems():
            score = max(scores)
            label = 0
            if (gA, gB) in observed:
                label = 1
            # print score, label
            o.write("{}\t{}\n".format(score, label))
            o_net.write("{}\t{}\t{}\t{}\n".format(gA, gB, score, label))
    return gene_pairs

##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    matches_fname = d + "/Alignments/matches_human_033bs_cov.tab"
    HPA_fname = d + "/Coexpression/human_033bs_cov_HPA.tab"
    Exclude_Uniprot_fname = d + "/Benchmarks/Immunoglobulin_uniprot.tab"
    vidal_fname = "/Users/agoncear/data/Vidal/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"

    # # Our vs Vidal
    b1_fname = d + "/Benchmarks/b1.scored-label"
    b1_netname = d + "/Benchmarks/b1.net"
    # Our Filtered vs Vidal
    b2_fname = d + "/Benchmarks/b2.scored-label"
    b2_netname = d + "/Benchmarks/b2.net"

    b3_fname = d + "/Benchmarks/b3.scored-label"
    b3_netname = d + "/Benchmarks/b3.net"

    b4_fname = d + "/Benchmarks/b4.scored-label"
    b4_netname = d + "/Benchmarks/b4.net"

    b5_fname = d + "/Benchmarks/b5.scored-label"
    b5_netname = d + "/Benchmarks/b5.net"

    b6_fname = d + "/Benchmarks/b6.scored-label"
    b6_netname = d + "/Benchmarks/b6.net"

    b7_fname = d + "/Benchmarks/b7.scored-label"
    b7_netname = d + "/Benchmarks/b7.net"

    # Mapping gene->prot
    gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    
    # HPA results for uniprot pairs
    HPA_strict = load_HPA_filter_results(HPA_fname, strict=True)

    # Exclude uniprot proteins: Ig, antigens, antibodies,...
    uniprot_filter = load_uniprot_filter(Exclude_Uniprot_fname)

    # Gene pairs from Vidal
    print "Loading Vidal"
    vidal_pairs = load_interactions(vidal_fname)
    # print "V:", len(vidal_pairs)

    # Uniprot pairs
    print "Loading our"
    # pairs1 = best_pairs(matches_fname, scoring=1)
    # print pairs1
    # pairs4 = best_pairs(matches_fname, scoring=4)
    # pairs5 = best_pairs(matches_fname, scoring=5)
    pairs6 = best_pairs(matches_fname, scoring=6)
    pairs7 = best_pairs(matches_fname, scoring=7)
    # Uniprot pairs
    print "Filtering our set with HPA"
    # pairs_filtered2 = filter_pairs(pairs1, HPA_strict)
    # print pairs_filtered2
    # pairs_filtered4 = filter_pairs(pairs4, HPA_strict)
    # pairs_filtered5 = filter_pairs(pairs5, HPA_strict)
    pairs_filtered6 = filter_pairs(pairs6, HPA_strict)
    pairs_filtered7 = filter_pairs(pairs7, HPA_strict)

    # Filter with uniprot
    print "Filtering with a Uniprot blacklist"
    # pairs_filtered_uniprot3 = filter_pairs_blacklist(pairs_filtered2, uniprot_filter)
    # print pairs_filtered_uniprot3
    # pairs_filtered_uniprot4 = filter_pairs_blacklist(pairs_filtered4, uniprot_filter)
    # pairs_filtered_uniprot5 = filter_pairs_blacklist(pairs_filtered5, uniprot_filter)
    pairs_filtered_uniprot6 = filter_pairs_blacklist(pairs_filtered6, uniprot_filter)
    pairs_filtered_uniprot7 = filter_pairs_blacklist(pairs_filtered7, uniprot_filter)

    # print "Benchmarking..."
    # gene_pairs = scored_labels(pairs, vidal_pairs, gene_proteins, b1_fname, b1_netname)
    # print "Benchmarking filtered..."
    # gene_pairs_filtered = scored_labels(pairs_filtered, vidal_pairs, gene_proteins, b2_fname, b2_netname)

    print "Benchmarking filtered uniprot..."
    # scored_labels(pairs1, vidal_pairs, gene_proteins, b1_fname, b1_netname)
    # scored_labels(pairs_filtered2, vidal_pairs, gene_proteins, b2_fname, b2_netname)
    # scored_labels(pairs_filtered_uniprot3, vidal_pairs, gene_proteins, b3_fname, b3_netname)
    # scored_labels(pairs_filtered_uniprot4, vidal_pairs, gene_proteins, b4_fname, b4_netname)
    # scored_labels(pairs_filtered_uniprot5, vidal_pairs, gene_proteins, b5_fname, b5_netname)
    scored_labels(pairs_filtered_uniprot6, vidal_pairs, gene_proteins, b6_fname, b6_netname)
    scored_labels(pairs_filtered_uniprot7, vidal_pairs, gene_proteins, b7_fname, b7_netname)

    # # gene_pairs = set(gene_pairs)
    # vidal_pairs = set(vidal_pairs)
    # uniprot_filtered_pairs = set(gene_pairs_filtered_uniprot)

    # print "Vidal:", len(vidal_pairs)
    # # print "Our:", len(gene_pairs)
    # print "Our filtered uniprot:", len(uniprot_filtered_pairs)
    # # print "Overlap:", len(vidal_pairs & gene_pairs)
    # print "Overlap:", len(vidal_pairs & uniprot_filtered_pairs)
