
from benchmark_common import get_scored_pairs, load_MITAB, scored_labels, roc, multiroc


def load_gi_uniprot_mapping(fname):
    mapping = {}
    with open(fname) as f:
        for line in f:
            gis, uniprot = line.strip().split()
            for gi in gis.split(","):
                mapping[gi] = uniprot
    return mapping


##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    roc_d = d + "/Benchmarks/Ecoli"
    matches_fname = d + "/Alignments/matches_ecoli.tab"
    validation_fname = "/Users/agoncear/data/Ecoli/MITAB-ecoli.txt"
    mapping_fname = "/Users/agoncear/data/Ecoli/gi_uniprot_mapping.tab"

    # Gene pairs from EBI (Uetz)
    print "Loading MIAB interactome ECOLI"
    pairs_uetz = load_MITAB(validation_fname)
    # mapping = load_gi_uniprot_mapping(mapping_fname)
    mapping = None
    # print mapping

    # Uniprot pairs
    print "Processing the interactome"
    pairs, a_pairs = get_scored_pairs(matches_fname)
    # pairs_homo, a_pairs = get_scored_pairs(matches_fname, mode="homo")
    # pairs_hetero, a_pairs = get_scored_pairs(matches_fname, mode="hetero")
    # pairs_short, a_pairs = get_scored_pairs(matches_fname, mode="short")
    # pairs_long, a_pairs = get_scored_pairs(matches_fname, mode="long")

    print "Comparing observed with predicted"
    # r_identity = scored_labels(pairs["identity"], pairs_uetz, mapping)
    # r_bs_identity = scored_labels(pairs["bs_identity"], pairs_uetz, mapping)
    # r_positive = scored_labels(pairs["positive"], pairs_uetz, mapping)
    # r_bs_positive = scored_labels(pairs["bs_positive"], pairs_uetz, mapping)
    # r_bs_coverage = scored_labels(pairs["bs_coverage"], pairs_uetz, mapping)
    # r_score = scored_labels(pairs["model_score"], pairs_uetz, mapping)
    # r_zscore = scored_labels(pairs["model_zscore"], pairs_uetz, mapping)
    # r_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_uetz, mapping)
    # r_random = scored_labels(pairs["random"], pairs_uetz, mapping)

    r_score = scored_labels(pairs["score"], pairs_uetz, mapping)
    
    # r_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_uetz, mapping)
    # r_m_homo = scored_labels(pairs_homo["model_minus_avg"], pairs_uetz, mapping)
    # r_m_hetero = scored_labels(pairs_hetero["model_minus_avg"], pairs_uetz, mapping)
    # r_m_short = scored_labels(pairs_short["model_minus_avg"], pairs_uetz, mapping)
    # r_m_long = scored_labels(pairs_long["model_minus_avg"], pairs_uetz, mapping)

    print "Making plots"
    # roc(r_identity, roc_d + "/identity.png", "Ecoli: Alignment identity")
    # roc(r_bs_identity, roc_d + "/bs_identity.png", "Ecoli: BS alignment identity")
    # roc(r_positive, roc_d + "/positive.png", "Ecoli: Alignment BLOSUM62 positive")
    # roc(r_bs_positive, roc_d + "/bs_positive.png", "Ecoli: BS alignment BLOSUM62 positive")
    # roc(r_bs_coverage, roc_d + "/bs_coverage.png", "Ecoli: BS coverage")
    # roc(r_score, roc_d + "/bs_score.png", "Ecoli: Compatibility score 1")
    # roc(r_zscore, roc_d + "/bs_zscore.png", "Ecoli: Compatibility Z-score")
    # roc(r_model_minus_avg, roc_d + "/bs_model_minus_avg.png", "Ecoli: Compatibility score-avg(decoy)")
    # roc(r_random, roc_d + "/bs_random.png", "Ecoli: Random control")

    # multiroc(
    #     (r_identity, r_bs_identity, r_positive, r_bs_positive),
    #     ('id', 'bs_id', 'pos', 'bs_pos'),
    #     roc_d + "/id_pos.png",
    #     "Basic ROC in Ecoli, MITAB Uetz 2014")

    # multiroc(
    #     (r_bs_positive, r_model_minus_avg),
    #     ('Binding site sequence similarity', 'Interface compatibily score'),
    #     roc_d + "/ecoli_combined.eps",
    #     "E.coli, Uetz 2014")

    # multiroc((r_identity, r_bs_positive),
    #     ("Total sequence identity", "Binding site similarity"),
    #     roc_d + "/ecoli_sequence.png", "E.coli, Uetz 2014")

    roc(r_score, roc_d + "/combined_score.png", "E.coli: combined score vs Uetz 2014")

    # multiroc((r_model_minus_avg, r_m_homo, r_m_hetero, r_m_short, r_m_long),
    #     ("Score (all)", "Score (homo)", "Score (hetero)", "Score (short)", "Score (long)"),
    #     roc_d + "/ecoli_templates.eps", "E.coli, Uetz 2014")

    # multiroc((r_model_minus_avg, r_m_homo, r_m_hetero, r_m_short, r_m_long),
    #     ("Score (all)", "Score (homo)", "Score (hetero)", "Score (short)", "Score (long)"),
    #     roc_d + "/ecoli_templates.eps", "E.coli, Uetz 2014")

