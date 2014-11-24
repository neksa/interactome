from itertools import islice
from collections import namedtuple
from benchmark_common import get_scored_pairs, roc, multiroc


# def load_gi_uniprot_mapping(fname):
#     mapping = {}
#     with open(fname) as f:
#         for line in f:
#             gis, uniprot = line.strip().split()
#             for gi in gis.split(","):
#                 mapping[gi] = uniprot
#     return mapping


def load_mapping(fname):
    mapping = {}
    with open(fname) as f:
        for line in f:
            a, b = line.strip().split()
            if not b.startswith("hpy:"):
                continue
            HP = b.split(":")[1]
            mapping[a] = HP
    return mapping


def load_interactions(fname):
    interactions = set()
    with open(fname) as f:
        for line in islice(f, 8, None):
            line = line.strip()
            if not line.startswith("HP"):
                continue
            pair = line.split()
            pair = pair[0], pair[1]
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


"""
C694_01785; HP_0352      FLIG_HELPY  O25119      343     fliG
"""
# Swissprot = namedtuple('Swissprot', 'swissprot, uniprot, length, genes')
def load_swissprot(fname):
    loci = {}
    with open(fname) as f:
        for line in islice(f, 30, None):
            line = line.strip()
            if not line.startswith('HP') and not line.startswith('C6'):
                continue
            locus = ''
            for l in line[:25].strip().split(';'):
                l = l.strip()
                if l.startswith('HP_'):
                    locus = l.split('.')[0]
                    locus = locus.replace('_', '')
                    break
            if locus == '':
                continue
            # swissprot = line[25:36].strip()
            uniprot = line[37:48].strip()
            # length = line[47:56].strip()
            # genes = line[57:].strip()
            # loci[locus] = Swissprot(swissprot=swissprot, uniprot=uniprot, length=length, genes=genes)
            loci[uniprot] = locus
            # print uniprot, locus
    return loci


##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    roc_d = d + "/Benchmarks/Hpylori"
    matches_fname = d + "/Alignments/matches_HPY_26695_PPI.tab"
    validation_fname = "/Users/agoncear/data/Hpylori/409211A0_S1.txt"
    # mapping_fname = "/Users/agoncear/data/Hpylori/helpy.txt"
    mapping_fname = "/Users/agoncear/data/Hpylori/M2014110493WE3322RA.tab"

    # Gene pairs from EBI (Uetz)
    print "Loading Interactome Hpylori from Rain et al., Nature 409 p211, 2000"
    pairs_obs = load_interactions(validation_fname)
    # print pairs_obs
    # mapping = load_swissprot(mapping_fname)
    mapping = load_mapping(mapping_fname)
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
    roc(r_identity, roc_d + "/identity.png", "H. pylori: Alignment identity")
    roc(r_bs_identity, roc_d + "/bs_identity.png", "H. pylori: BS alignment identity")
    roc(r_positive, roc_d + "/positive.png", "H. pylori: Alignment BLOSUM62 positive")
    roc(r_bs_positive, roc_d + "/bs_positive.png", "H. pylori: BS alignment BLOSUM62 positive")
    roc(r_bs_coverage, roc_d + "/bs_coverage.png", "H. pylori: BS coverage")
    roc(r_score, roc_d + "/bs_score.png", "H. pylori: Compatibility score 1")
    roc(r_zscore, roc_d + "/bs_zscore.png", "H. pylori: Compatibility Z-score")
    roc(r_model_minus_avg, roc_d + "/bs_model_minus_avg.png", "H. pylori: Compatibility score-avg(decoy)")
    roc(r_random, roc_d + "/bs_random.png", "H. pylori: Random control")

    multiroc(
        (r_identity, r_bs_identity, r_positive, r_bs_positive),
        ('id', 'bs_id', 'pos', 'bs_pos'),
        roc_d + "/id_pos.png",
        "Basic ROC in H.pylori, Rain et al. 2000")

    multiroc(
        (r_bs_positive, r_model_minus_avg),
        ('Binding site sequence similarity', 'Interface compatibily score'),
        roc_d + "/ecoli_combined.eps",
        "H.pylori, Rain et al. 2000")
