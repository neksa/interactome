from itertools import islice
from benchmark_common import get_scored_pairs, scored_labels, multiroc, roc


def load_vidal_interactions(fname):
    interactions = []
    with open(fname) as f:
        for line in islice(f, 0, None):
            A, B = line.strip().split("\t")
            pair = (A, B)
            if A > B:
                pair = (B, A)
            interactions.append(pair)
    return set(interactions)


##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    roc_d = d + "/Benchmarks/Yeast2"
    matches_fname = d + "/Alignments/matches_yeast.tab"
    validation_fname = "/Users/agoncear/data/Vidal/Yeast/Y2H_union.txt"
    golden_fname = "/Users/agoncear/data/Vidal/Yeast/Binary-GS.txt"

    # Gene pairs from Vidal
    print "Loading Vidal interactome Y2H"
    pairs_vidal = load_vidal_interactions(validation_fname)
    pairs_golden = load_vidal_interactions(golden_fname)

    # Uniprot pairs
    print "Processing the interactome"
    # max pairs, all pairs
    pairs, a_pairs = get_scored_pairs(matches_fname)

    print "Comparing observed with predicted: Union"
    # r_identity = scored_labels(pairs["identity"], pairs_vidal)
    # r_bs_identity = scored_labels(pairs["bs_identity"], pairs_vidal)
    # r_positive = scored_labels(pairs["positive"], pairs_vidal)
    r_bs_positive = scored_labels(pairs["bs_positive"], pairs_vidal)
    # r_bs_coverage = scored_labels(pairs["bs_coverage"], pairs_vidal)
    # r_score = scored_labels(pairs["model_score"], pairs_vidal)
    # r_zscore = scored_labels(pairs["model_zscore"], pairs_vidal)
    r_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_vidal)
    # r_random = scored_labels(pairs["random"], pairs_vidal)

    print "Comparing observed with predicted: Golden"
    # g_identity = scored_labels(pairs["identity"], pairs_golden)
    # g_bs_identity = scored_labels(pairs["bs_identity"], pairs_golden)
    # g_positive = scored_labels(pairs["positive"], pairs_golden)
    g_bs_positive = scored_labels(pairs["bs_positive"], pairs_golden)
    # g_bs_coverage = scored_labels(pairs["bs_coverage"], pairs_golden)
    # g_score = scored_labels(pairs["model_score"], pairs_golden)
    # g_zscore = scored_labels(pairs["model_zscore"], pairs_golden)
    g_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_golden)
    # g_random = scored_labels(pairs["random"], pairs_golden)

    # print "Comparing observed with predicted: Union"
    # ra_identity = scored_labels(a_pairs["identity"], pairs_vidal)
    # ra_bs_identity = scored_labels(a_pairs["bs_identity"], pairs_vidal)
    # ra_positive = scored_labels(a_pairs["positive"], pairs_vidal)
    # ra_bs_positive = scored_labels(a_pairs["bs_positive"], pairs_vidal)
    # ra_bs_coverage = scored_labels(a_pairs["bs_coverage"], pairs_vidal)
    # ra_score = scored_labels(a_pairs["model_score"], pairs_vidal)
    # ra_zscore = scored_labels(a_pairs["model_zscore"], pairs_vidal)
    # ra_model_minus_avg = scored_labels(a_pairs["model_minus_avg"], pairs_vidal)
    # ra_random = scored_labels(a_pairs["random"], pairs_vidal)

    # print "Comparing observed with predicted: Golden"
    # ga_identity = scored_labels(a_pairs["identity"], pairs_golden)
    # ga_bs_identity = scored_labels(a_pairs["bs_identity"], pairs_golden)
    # ga_positive = scored_labels(a_pairs["positive"], pairs_golden)
    # ga_bs_positive = scored_labels(a_pairs["bs_positive"], pairs_golden)
    # ga_bs_coverage = scored_labels(a_pairs["bs_coverage"], pairs_golden)
    # ga_score = scored_labels(a_pairs["model_score"], pairs_golden)
    # ga_zscore = scored_labels(a_pairs["model_zscore"], pairs_golden)
    # ga_model_minus_avg = scored_labels(a_pairs["model_minus_avg"], pairs_golden)
    # ga_random = scored_labels(a_pairs["random"], pairs_golden)

    # roc(g_bs_positive, roc_d + "/roc_bs_positive.png", "Yeast, golden set (best model), BS positive (BLOSUM62)")
    # roc(g_bs_coverage, roc_d + "/roc_bs_coverage.png", "Yeast, golden set (best model), BS coverage")
    # roc(g_score, roc_d + "/roc_score.png", "Yeast, golden set (best model), Score")
    # roc(g_zscore, roc_d + "/roc_zscore.png", "Yeast, golden set (best model), Z-score")
    # roc(g_model_minus_avg, roc_d + "/roc_model_minus_avg.png", "Yeast, golden set (best model), Score-avg(decoys)")

    print "Making plots"
    # legend = ('Y2H union (best model)', 'Gold set from literature (best model)', 'Y2H union (all models)', 'Gold set from literature (all models)')
    legend = ('Y2H union (Yu et.al.)', 'Gold set from literature (Yu et.al.)')
    # multiroc((r_identity, g_identity, ra_identity, ga_identity), legend, roc_d + "/identity.png", "Yeast: Alignment identity")
    # multiroc((r_bs_identity, g_bs_identity, ra_bs_identity, ga_bs_identity), legend, roc_d + "/bs_identity.png", "Yeast: BS alignment identity")
    # multiroc((r_positive, g_positive, ra_positive, ga_positive), legend, roc_d + "/positive.png", "Yeast: Alignment BLOSUM62 positive")
    # multiroc((r_bs_positive, g_bs_positive, ra_bs_positive, ga_bs_positive), legend, roc_d + "/bs_positive.png", "Yeast: BS alignment BLOSUM62 positive")
    multiroc((r_bs_positive, g_bs_positive), legend, roc_d + "/bs_positive.png", "Yeast: Binding site similarity score")
    # # multiroc((r_bs_positive_homo, g_bs_positive_homo), legend, roc_d + "/bs_positive_homo.png", "Yeast: Preserve homo/hetero type. BS BLOSUM62 positive")
    # multiroc((r_bs_coverage, g_bs_coverage, ra_bs_coverage, ga_bs_coverage), legend, roc_d + "/bs_coverage.png", "Yeast: BS coverage")
    # multiroc((r_score, g_score, ra_score, ga_score), legend, roc_d + "/score.png", "Yeast: Compatibility score 1")
    # multiroc((r_zscore, g_zscore, ra_zscore, ga_zscore), legend, roc_d + "/zscore.png", "Yeast: Compatibility Z-score")
    # multiroc((r_model_minus_avg, g_model_minus_avg, ra_model_minus_avg, ga_model_minus_avg), legend, roc_d + "/model_minus_avg.png", "Yeast: Compatibility score-avg(decoy)")
    # multiroc((r_random, g_random, ra_random, ga_random), legend, roc_d + "/random.png", "Yeast: Random control")
    # # roc(r_random, roc_d + "/bs_random.png", "Random control")

    # multiroc(
    #     (r_identity, r_bs_identity, r_positive, r_bs_positive),
    #     ('id', 'bs_id', 'pos', 'bs_pos'),
    #     roc_d + "/id_pos.png",
    #     "Basic ROC in Yeast, Y2H_union")
    multiroc(
        (r_bs_positive, g_bs_positive, r_model_minus_avg, g_model_minus_avg),
        ('Bind. site similarity (Y2H)', 'Bind. site similarity (Gold set)', 'Interface compatibily (Y2H)', 'Interface compatibility (Gold set)'),
        roc_d + "/yeast_combined.eps",
        "Yeast, Y2H and Gold set (Literature)")
