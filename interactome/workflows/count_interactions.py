from itertools import islice
from collections import defaultdict
import random


def count_scored_pairs(matches_fname, t_bs=-1000, t_score=-1000):
    pairs = defaultdict(lambda: defaultdict(list))
    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            try:
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

                # if mode == "homo":
                #     if template_type != "Homo":
                #         continue

                # if mode == "hetero":
                #     if template_type != "Hetero":
                #         continue

                # if mode == "short":
                #     if bs_alignedA > 10 or bs_alignedB > 10:
                #         continue

                # if mode == "long":
                #     if bs_alignedA <= 10 or bs_alignedB <= 10:
                #         continue

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                pair = (queryA, queryB)
                if queryA > queryB:
                    pair = (queryB, queryA)


                bs = min(
                    bs_positiveA / bs_lenA,
                    bs_positiveB / bs_lenB)

                if bs < t_bs:
                    continue

                if model_minus_avg < t_score:
                    continue

                # pairs["identity"][pair].append(min(
                #     identicalA / aln_lenA,
                #     identicalB / aln_lenB))
                # pairs["bs_identity"][pair].append(min(
                #     bs_identicalA / bs_lenA,
                #     bs_identicalB / bs_lenB))
                # pairs["positive"][pair].append(min(
                #     positiveA / aln_lenA,
                #     positiveB / aln_lenB))
                pairs["bs_positive"][pair].append(min(
                    bs_positiveA / bs_lenA,
                    bs_positiveB / bs_lenB))
                # pairs["bs_coverage"][pair].append(min(
                #     bs_coveredA / bs_lenA,
                #     bs_coveredB / bs_lenB))
                # pairs["random"][pair].append(random.random())
                # pairs["model_score"][pair].append(score_model)
                # pairs["model_zscore"][pair].append(zscore)
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

    # max_pairs = defaultdict(list)
    pairs = set()
    for pair, scores in pairs["model_minus_avg"]:
        pairs.add(pair)
    return pairs


def main():
    d = "/Users/agoncear/projects/Interactome/Workflow"
    # roc_d = d + "/Benchmarks/Ecoli"
    matches_ecoli = d + "/Alignments/matches_ecoli.tab"
    pairs = count_scored_pairs(matches_ecoli, 0.25)  # , 3.0)
    p = set()
    for pair in pairs:
        p.add(pair[0])
        p.add(pair[1])
    print "N proteins = ", len(p)
    print "N interactions = ", len(pairs)


if __name__ == '__main__':
    main()

