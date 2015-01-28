from itertools import islice
from collections import defaultdict
import random

import pylab as pl
from sklearn.metrics import roc_curve, auc


def get_scored_pairs(matches_fname, mode="max"):
    pairs = defaultdict(lambda: defaultdict(list))
    random.seed(800)
    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            # print "L", len(line.strip().split("\t"))
            # for i, a in enumerate(line.strip().split("\t")):
            #     print i, "==", a, "=="
            # 11 + 2 + 11 + 11 = 35
            try:
                # queryA, queryB, tpl, query_type, template_type,  #template_nsubunits,\
                queryA, queryB, tpl, query_type, template_type, \
                    SC1, SC2, SC3, SC4, SC5, SC6, \
                    identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                    identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                    siteA, siteB = line.strip().split("\t")
                    # _, _, _, _, \

                score_template_full, score_template, score_model, scaled_score_template_full, scaled_score_template, scaled_score = map(
                    float, (SC1, SC2, SC3, SC4, SC5, SC6))

                zscore, model_minus_avg = score_template_full, score_template

                identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A = map(
                    float, (identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A))

                identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B = map(
                    float, (identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B))

                if mode == "homo":
                    if template_type != "Homo":
                        continue

                if mode == "hetero":
                    if template_type != "Hetero":
                        continue

                if mode == "short":
                    if bs_alignedA > 10 or bs_alignedB > 10:
                        continue

                if mode == "long":
                    if bs_alignedA <= 10 or bs_alignedB <= 10:
                        continue

                if template_type != query_type:
                    continue

                # score = score_template_full
                score = score_model
                bs_positive = min(
                    bs_positiveA / bs_lenA,
                    bs_positiveB / bs_lenB)

                full_identity = min(
                    identicalA / aln_lenA,
                    identicalB / aln_lenB)

                score = (1-full_identity) * 0.067 * score + full_identity * 5.07 * bs_positive

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                pair = (queryA, queryB)
                if queryA > queryB:
                    pair = (queryB, queryA)

                pairs["identity"][pair].append(full_identity)
                pairs["bs_identity"][pair].append(min(
                    bs_identicalA / bs_lenA,
                    bs_identicalB / bs_lenB))
                pairs["positive"][pair].append(min(
                    positiveA / aln_lenA,
                    positiveB / aln_lenB))
                pairs["bs_positive"][pair].append(bs_positive)
                pairs["bs_coverage"][pair].append(min(
                    bs_coveredA / bs_lenA,
                    bs_coveredB / bs_lenB))
                pairs["random"][pair].append(random.random())
                pairs["score"][pair].append(score) # COMBINED
                pairs["model_score"][pair].append(score_model)
                pairs["model_zscore"][pair].append(zscore)
                pairs["model_minus_avg"][pair].append(model_minus_avg)

                # print pair, pairs["model_score"][pair]

                # pairs["bs_positive_homo"][pair].append()
                #     if query_type == "Unknown" or template_type == "Unknown" or query_type == template_type:
                #         score = min(
                #             bs_positiveA / bs_lenA,
                #             bs_positiveB / bs_lenB)
                #     else:
                #         score = 0.0
            except:
                print "Error reading the whole matches file. Skipping the line with an error"

    max_pairs = defaultdict(list)
    for scoring, pair_scores in pairs.iteritems():
        for pair, scores in pair_scores.iteritems():
            max_pairs[scoring].append((pair, max(scores)))

    all_pairs = defaultdict(list)
    for scoring, pair_scores in pairs.iteritems():
        for pair, scores in pair_scores.iteritems():
            for score in set(scores):
                all_pairs[scoring].append((pair, score))

    return max_pairs, all_pairs
    # return {pair: max(scores) for pair, scores in pairs.iteritems()}


def load_MITAB(fname):
    """
    uniprotkb:P64503    uniprotkb:P0A6R3    intact:EBI-9126792|uniprotkb:Q2MB15|uniprotkb:P76274    intact:EBI-550170|uniprotkb:Q2M8V4|uniprotkb:P37404|uniprotkb:P11028    psi-mi:yebv_ecoli(display_long)|uniprotkb:JW5302(locus name)|uniprotkb:b1836(locus name)|uniprotkb:yebV(gene name)|psi-mi:yebV(display_short)   psi-mi:fis_ecoli(display_long)|uniprotkb:fis(gene name)|psi-mi:fis(display_short)|uniprotkb:b3261(locus name)|uniprotkb:JW3229(locus name)|uniprotkb:Factor-for-inversion stimulation protein(gene name synonym)|uniprotkb:Hin recombinational enhancer-binding protein(gene name synonym)  psi-mi:"MI:0096"(pull down) Rajagopala et al. (2014)    pubmed:24561554|imex:IM-22059   taxid:83333(ecoli)|taxid:83333("Escherichia coli (strain K12)") taxid:83333(ecoli)|taxid:83333("Escherichia coli (strain K12)") psi-mi:"MI:0915"(physical association)  psi-mi:"MI:0469"(IntAct)    intact:EBI-9200188|imex:IM-22059-2470   intact-miscore:0.51 -   psi-mi:"MI:0499"(unspecified role)  psi-mi:"MI:0499"(unspecified role)  psi-mi:"MI:0498"(prey)  psi-mi:"MI:0496"(bait)  psi-mi:"MI:0326"(protein)   psi-mi:"MI:0326"(protein)   interpro:IPR009950(Protein of unknown function DUF1480)|go:"GO:0005515"(protein binding)|refseq:YP_490098.1|refseq:NP_416350.4|uniprotkb:P64503(identity)   interpro:IPR002197(Helix-turn-helix, Fis-type)|interpro:IPR009057(Homeodomain-like)|interpro:IPR005412(DNA-binding protein Fis)|rcsb pdb:1ETQ|rcsb pdb:1ETO|rcsb pdb:1ETK|rcsb pdb:4FIS|rcsb pdb:3JRI|rcsb pdb:3JRH|rcsb pdb:3JRG|rcsb pdb:3JRF|rcsb pdb:3JRE|rcsb pdb:3JRD|rcsb pdb:3JRC|rcsb pdb:3JRB|rcsb p
    """
    interactions = []
    with open(fname) as f:
        for line in islice(f, 1, None):
            A, B, _ = line.strip().split("\t", 2)
            A = A.split(":")[1]
            B = B.split(":")[1]
            pair = (A, B)
            if A > B:
                pair = (B, A)
            interactions.append(pair)
    return set(interactions)


def scored_labels(predicted_pairs, observed_pairs, mapping=None):
    true_labels = []
    scores = []
    # mapped_pairs = []
    # mapping
    for unmapped_pair, score in predicted_pairs:
        A, B = pair = unmapped_pair
        if mapping is not None:
            # giA = A.split("|")[1]
            mappedA = mapping.get(A)
            if mappedA is None:
                print "skipping unmapped", A
                continue
            # giB = B.split("|")[1]
            mappedB = mapping.get(B)
            if mappedB is None:
                print "skipping unmapped", B
                continue
            pair = mappedA, mappedB
            if mappedA > mappedB:
                pair = (mappedB, mappedA)
            # if pair not in mapped_pairs:
            #     mapped_pairs[pair] = score
            # else:
            #     mapped_pairs[pair] = max(mapped_pairs[pair], score)
        # mapped_pairs.append((pair, score))
        # # benchmarking
        # for pair, score in mapped_pairs:
        label = 0
        if pair in observed_pairs:
            label = 1
            # print pair, score, label
        true_labels.append(label)
        scores.append(score)
    return true_labels, scores


def roc(r, fname, title=None):
    test_labels, scores = r
    fpr, tpr, thresholds = roc_curve(test_labels, scores)

    with open(fname.split(".")[0] + ".txt", 'w') as o:
        o.write("tpr\tfpr\tthreshold\n")
        for i in xrange(len(fpr)):
            o.write("{}\t{}\t{}\n".format(tpr[i], fpr[i], thresholds[i]))
    roc_auc = auc(fpr, tpr)
    print "Area under the ROC curve : %f" % roc_auc

    # Plot ROC curve
    pl.clf()
    pl.plot(fpr, tpr, label='AUC = %0.3f' % roc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title(title)
    pl.legend(loc="lower right")
    pl.savefig(fname)


def multiroc(rs, curves, fname, title=None):
    pl.clf()
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.plot([0, 1], [0, 1], 'k--')
    for i, r in enumerate(rs):
        curve = curves[i]
        test_labels, scores = r
        fpr, tpr, thresholds = roc_curve(test_labels, scores)
        roc_auc = auc(fpr, tpr)
        print "ROC %s AUC : %f" % (curve, roc_auc)
        # Plot ROC curve
        pl.plot(fpr, tpr, label='%s (AUC = %0.2f)' % (curve, roc_auc))
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title(title)
    pl.legend(loc="lower right")
    pl.savefig(fname)
