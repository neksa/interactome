from itertools import islice
from collections import defaultdict


def get_scored_pairs(matches_fname):
    # pairs = defaultdict(lambda: defaultdict(list))
    max_pairs = defaultdict(dict)
    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            # print "L", len(line.strip().split("\t"))
            # for i, a in enumerate(line.strip().split("\t")):
            #     print i, "==", a, "=="
            # 11 + 2 + 11 + 11 = 35
            try:
                queryA, queryB, tpl, query_type, template_type, \
                    SC1, SC2, SC3, SC4, SC5, SC6, \
                    identicalA, poitiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                    identicalB, poitiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                    siteA, siteB = line.strip().split("\t")
                    # _, _, _, _, \

                score_template_full, score_template, score_model, scaled_score_template_full, scaled_score_template, scaled_score = map(
                    float, (SC1, SC2, SC3, SC4, SC5, SC6))

                zscore, model_minus_avg = score_template_full, score_template

                identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A = map(
                    float, (identicalA, poitiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A))

                identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B = map(
                    float, (identicalB, poitiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B))

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                pair = (queryA, queryB)
                if queryA > queryB:
                    pair = (queryB, queryA)

                m = max_pairs[pair].get("model_minus_avg")
                if m is None or model_minus_avg > m:
                    max_pairs[pair]["model_minus_avg"] = model_minus_avg
                    max_pairs[pair]["bs_positive"] = min(
                        bs_positiveA / bs_lenA,
                        bs_positiveB / bs_lenB)
                    max_pairs[pair]["bs_coverage"] = min(
                        bs_coveredA / bs_lenA,
                        bs_coveredB / bs_lenB)
                    max_pairs[pair]["sites"] = (siteA, siteB)

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

    # for scoring, pair_scores in pairs.iteritems():
    #     for pair, scores in pair_scores.iteritems():
    #         max_pairs[pair][scoring] = max(scores)
    return max_pairs


def savefile(infile, outfile):
    pairs = get_scored_pairs(infile)
    with open(outfile, 'w') as o:
        o.write("A\tB\tSequence similarity\tInterface compatibility\tsiteA\tsiteB\n")
        for pair in pairs:
            S1 = pairs[pair]["bs_positive"]
            S2 = pairs[pair]["model_minus_avg"]
            siteA, siteB = pairs[pair]["sites"]
            # S3 = pairs[pair]["bs_coverage"]
            o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pair[0], pair[1], S1, S2, siteA, siteB))


if __name__ == '__main__':
    DIR = "/Users/agoncear/projects/Interactome/Workflow/"

    # infile = DIR + "Alignments/matches_HPY_26695_PPI.tab"
    # outfile = DIR + "Interactomes/bacteria/HPY.tab"
    # savefile(infile, outfile)

    # infile = DIR + "Alignments/matches_SPN_TIGR4_PPI.tab"
    # outfile = DIR + "Interactomes/bacteria/SPN.tab"
    # savefile(infile, outfile)

    infile = DIR + "Alignments/matches_ecoli.tab"
    outfile = DIR + "Interactomes/bacteria/ecoli.tab"
    savefile(infile, outfile)



