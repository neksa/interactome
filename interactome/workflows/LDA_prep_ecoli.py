
from benchmark_common import get_scored_pairs, load_MITAB, scored_labels, roc, multiroc
from collections import defaultdict
from itertools import islice


def main():
    d = "/Users/agoncear/projects/Interactome/Workflow"
    matches_fname = d + "/Alignments/matches_ecoli.tab"
    validation_fname = "/Users/agoncear/data/Ecoli/MITAB-ecoli.txt"
    mapping_fname = "/Users/agoncear/data/Ecoli/gi_uniprot_mapping.tab"
    results_fname = d + "/Alignments/labeled_matches_ecoli.tab"

    print "Load GI -> Uniprot mapping"
    mapping = load_gi_uniprot_mapping(mapping_fname)

    # # Gene pairs from EBI (Uetz)
    print "Loading MITAB interactome ECOLI"
    pairs_uetz = load_MITAB(validation_fname)

    # M = namedtuple('M',
    #         "queryA queryB template query_type template_type " +
    #         "n1 m1 n2 m2 " +
    #         "score_template_full score_template score " +
    #         "scaled_score_template_full scaled_score_template scaled_score " +
    #         "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +
    #         "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +
    #         "siteA siteB")

    print "Writing results to", results_fname
    with open(matches_fname) as f, open(results_fname, 'w') as o:  # , open(aln_stringent, 'w') as o_aln:
        # header = "query_A query_B template query_type template_type " +\
        #     "score_template_full score_template score " +\
        #     "scaled_score_template_full scaled_score_template scaled_score " +\
        #     "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +\
        #     "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +\
        #     "siteA siteB"

        header = "label query_A query_B template query_type template_type " +\
            "zscore score_minus_avg score " +\
            "score4 score5 score6 " +\
            "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +\
            "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B"

        o.write("\t".join(header.split(" ")) + "\n")

        for i, line in enumerate(f):
            fields = line.strip().split("\t")
            # if len(fields) != 39:
            #     continue
            # I = M._make(fields)
            uniprot_A = mapping.get(fields[0])  # query_A)
            uniprot_B = mapping.get(fields[1])  # query_B)
            if uniprot_A is None or uniprot_B is None:
                continue
            pair = uniprot_A, uniprot_B
            if uniprot_A > uniprot_B:
                pair = uniprot_B, uniprot_A

            label = False
            if pair in pairs_uetz:
                label = True

            output = [label]
            output.extend(fields[:-2])
            format_str = "\t".join(['{}']*len(output)) + "\n"
            o.write(format_str.format(*output))

            """
            ### STRINGENT FILTER ####
            seq_id = min(float(I.identicalA) / float(I.aln_lenA), float(I.identicalB) / float(I.aln_lenB))
            bs_seq_id = min(float(I.bs_identicalA) / float(I.bs_alignedA), float(I.bs_identicalB) / float(I.bs_alignedB))
            bs_seq_aligned = min(float(I.bs_alignedA)/float(I.bs_lenA), float(I.bs_alignedB)/float(I.bs_lenB))

            # print seq_id, bs_seq_id, bs_seq_aligned, I.bs_BLOSUMA, I.bs_BLOSUMB, I.bs_score1A, I.bs_score1B, I.template_type, I.query_type

            if seq_id < 0.5:
                continue
            if bs_seq_id < 0.5:
                continue
            if bs_seq_aligned < 1.0:
                continue

            if I.template_type == 'Homo' and I.query_type != 'Homo':
                continue
            if I.template_type == 'Hetereo' and I.query_type != 'Hetero':
                continue

            if I.bs_BLOSUMA < 0 or I.bs_BLOSUMB < 0:
                continue
            if I.bs_score1A < 0 or I.bs_score1B < 0:
                continue

            #########################

            siteA = set()
            for site in I.siteA.split(';'):
                s = site.split(',')
                siteA.add(s[3] + s[4])

            siteB = set()
            for site in I.siteB.split(';'):
                s = site.split(',')
                siteB.add(s[3] + s[4])

            #########################
            """


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
    main()
 