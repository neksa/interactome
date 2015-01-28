from itertools import islice
from benchmark_common import get_scored_pairs, scored_labels, multiroc, roc
from itertools import islice
from collections import defaultdict
import random

"""
Residue-level ROC calculations
Score (s) is tunable

Main data structures for residue-level ROC:

P = goodres[uniprot] = set(), goodresnum[uniprot] = num_residues:int
N = badresnum[uniprot] = num_residues:int

Hits -> predres[uniprot][resi:int] = max_score:float

Mapping structures:

* From yeast.txt (Uniprot database):
locus_uniprot[locus] = uniprot
uniprot_nresidues[uniprot] nresidues:int

* From source.idx (PDB database):
yeast_pdbs = list(4 letter pdb in lowercase)

* From pdb_templates_5A.tab
  From pdb_proteins.tab

  (pdb,chain,resi) => goodres[uniprot] = list(resi:int)
  goodresnum[uniprot] = count(list)

* From  goodres:
  badresnum[uniprot in goodres] = uniprot_nresidues[uniprot] - goodresnum[uniprot]

* From matches.tab:
  uniprot = locus_uniprot[locus]
  for resi in siteA and in siteB
  predres[uniprot][resi:int] = max(score, score, score, ...)

"""


"""
AAC1                                                                       YMR056C             P04710     ADT1_YEAST  S000004660   309     13
AAC3                                                                       YBR085W             P18238     ADT3_YEAST  S000000289   307 (3)  2

We need OLN [76-83), Swissprot [96-102) and size [130,135)
Starting with line #60
"""
def load_uniprot():
    locus_uniprot = {}
    uniprot_nresidues = {}
    fname = "/Users/agoncear/data/Yeast/yeast.txt"
    with open(fname) as f:
        for i, line in enumerate(islice(f, 60, 6781)):
            locus = line[75:82].strip()
            uniprot = line[95:101].strip()
            if uniprot == '':
                continue
            # print i, "OK", line
            length = int(line[129:134].strip())
            locus_uniprot[locus] = uniprot
            uniprot_nresidues[uniprot] = length
    return locus_uniprot, uniprot_nresidues


def load_structures():
    pdbs = []
    fname = "/Users/agoncear/data/PDB/derived_data/source.idx"
    with open(fname) as f:
        for line in islice(f, 4, None):
            fields = line.strip().split("\t")
            if len(fields) == 2:
                if fields[1].startswith("SACCHAROMYCES CEREVISIAE"):
                    pdbs.append(fields[0].lower())
    return pdbs

"""
PDB:
  resi = 10
  uniprot_offset = 1
  resi_uniprot = 10 - 1 + 1

  resi = 10
  uniprot_offset = 5
  resi_uniprot = 10 - 5 + 1

"""
def load_templates(pdbs):
    d = "/Users/agoncear/projects/Interactome/Workflow/Structures/"
    fname_pdb_proteins = d + "pdb_proteins.tab"
    fname_templates = d + "pdb_templates_5A.tab"
    # pdb_protein = {}
    pdb_protein = defaultdict(list)
    goodres = defaultdict(set)

    with open(fname_pdb_proteins) as f_pdb, open(fname_templates) as f_template:
        for line in islice(f_pdb, 1, None):
            (pdb, chain, chain_author, uniprot,
                seq_aln_begin, seq_aln_begin_ins, seq_aln_end, seq_aln_end_ins,
                db_aln_begin, db_aln_end,
                auth_aln_begin, auth_aln_end) = line.strip().split("\t")
            # if uniprot == '':
            #     continue
            pdb_protein[(pdb, chain)].append((uniprot, int(seq_aln_begin), int(seq_aln_end), int(db_aln_begin), int(db_aln_end)))

        for line in f_template:
            pdb, chainA, chainB, siteA, siteB = line.strip().split("\t")
            if pdb not in pdbs:
                continue
            A = chainA.split("_")[0]
            B = chainB.split("_")[0]

            resmapA = {}
            for uniprotA, seq_offsetA, seq_endA, db_offsetA, db_endA in pdb_protein[(pdb, A)]:
                resmapA[(seq_offsetA, seq_endA)] = uniprotA, db_offsetA

            resmapB = {}
            for uniprotB, seq_offsetB, seq_endB, db_offsetB, db_endB in pdb_protein[(pdb, B)]:
                resmapB[(seq_offsetB, seq_endB)] = uniprotB, db_offsetB

            # examples of Uniprot-PDB discrepancies: 1a6r (Nterm), 4c0b (Cterm) and 4c46 (Completely artificial sequence)
            # if pdb == "1a6r":
            #     print "MAP A", A, resmapA
            #     print "MAP B", B, resmapB

            for s in siteA.split(";"):
                resn, resi, ncontacts = s.split(",")
                resi = int(resi)
                if int(ncontacts) == 1:
                    continue
                # uniprot_resi = int(resi) - seq_offsetA + db_offsetA
                # goodres[uniprotA].add(uniprot_resi)

                ok = False
                for key, val in resmapA.iteritems():
                    start, end = key
                    if start <= resi <= end:
                        uniprotA, db_offsetA = val
                        # skip unknown proteins
                        if uniprotA == '':
                            ok = True
                            continue
                        uniprot_resi = int(resi) - start + db_offsetA
                        goodres[uniprotA].add(uniprot_resi)
                        ok = True
                if not ok:
                    # print "UNKNOWN", pdb, chainA, resn, resi
                    pass
                # if uniprotA == "P25604":
                #     print pdb, chainA, A, uniprotA, resn, resi, "offset", offsetA, uniprot_resi

            for s in siteB.split(";"):
                resn, resi, ncontacts = s.split(",")
                resi = int(resi)
                if int(ncontacts) == 1:
                    continue
                # uniprot_resi = int(resi) - seq_offsetB + db_offsetB
                # goodres[uniprotB].add(uniprot_resi)
                ok = False
                for key, val in resmapB.iteritems():
                    start, end = key
                    if start <= resi <= end:
                        uniprotB, db_offsetB = val
                        # skip unknown proteins
                        if uniprotB == '':
                            ok = True
                            continue
                        uniprot_resi = int(resi) - start + db_offsetB
                        goodres[uniprotB].add(uniprot_resi)
                        ok = True
                if not ok:
                    # print "UNKNOWN", pdb, chainB, resn, resi
                    pass

    return goodres


def get_predicted_residues(
        matches_fname, locus_uniprot, goodres, uniprot_nresidues,
        pdbs=None, scoring="score1", calibration=False, calib_fname=None, dump_fname=None,
        threshold_identity=1.0, threshold_nsubunits=None, same_complex_type=None):
    fname_pdb_proteins = "/Users/agoncear/projects/Interactome/Workflow/Structures/pdb_proteins.tab"
    pdb_protein = defaultdict(list)
    with open(fname_pdb_proteins) as f_pdb:
        for line in islice(f_pdb, 1, None):
            (pdb, chain, chain_author, uniprot,
                seq_aln_begin, seq_aln_begin_ins, seq_aln_end, seq_aln_end_ins,
                db_aln_begin, db_aln_end,
                auth_aln_begin, auth_aln_end) = line.strip().split("\t")
            if uniprot == '':
                continue
            pdb_protein[(pdb, chain)].append((uniprot, int(seq_aln_begin), int(db_aln_begin)))

    predres = defaultdict(lambda: defaultdict(float))

    pred_dump = {
        'score': defaultdict(lambda: defaultdict(float)),
        'template_type': defaultdict(lambda: defaultdict(float)),
        'template_nsubunits': defaultdict(lambda: defaultdict(float)),
        'model_size': defaultdict(lambda: defaultdict(float)),
        'model_bs_size': defaultdict(lambda: defaultdict(float)),
        'bs_sim': defaultdict(lambda: defaultdict(float)),
        'bs_aln_sim': defaultdict(lambda: defaultdict(float)),
        'full_id': defaultdict(lambda: defaultdict(float))
    }

    with open(matches_fname) as f:
        for c, line in enumerate(islice(f, 1, None)):
            if c % 100000 == 0:
                print c
            try:
                # # file_format = 1  # single scoring
                file_format = 2  # multiscoring with number of subunits

                # if file_format == 1:
                #     queryA, queryB, tpl, query_type, template_type, \
                #         SC1, SC2, SC3, SC4, SC5, SC6, \
                #         identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                #         identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                #         siteA, siteB = line.strip().split("\t")
                #     SC1, SC2, SC3, SC4, SC5, SC6 = map(float, (SC1, SC2, SC3, SC4, SC5, SC6))
                #     # zscore, model_minus_avg, score_model = SC1, SC2, SC3
                #     template_nsubunits = 0

                queryA, queryB, tpl, query_type, template_type, template_nsubunits, \
                    SC1, SC2, SC3, SC4, SC5, SC6, \
                    identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                    identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                    siteA, siteB = line.strip().split("\t")
                SC1, SC2, SC3, SC4, SC5, SC6 = map(float, (SC1, SC2, SC3, SC4, SC5, SC6))
                # Convert Energies to scores (the higher - the better)
                SC2 = 0.0 - SC2
                SC3 = 0.0 - SC3
                template_nsubunits = int(template_nsubunits)

                identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A = map(
                    float, (identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A))

                identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B = map(
                    float, (identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B))

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                bs_similarityA = bs_positiveA / bs_lenA  # BS SIMILARITY
                bs_similarityB = bs_positiveB / bs_lenB  # BS SIMILARITY
                bs_similarity = min(bs_similarityA, bs_similarityB)  # BS SIMILARITY

                aligned_bs_similarityA = bs_positiveA / bs_alignedA
                aligned_bs_similarityB = bs_positiveB / bs_alignedB
                aligned_bs_similarity = min(aligned_bs_similarityA, aligned_bs_similarityB)  # ALIGNED BS SIMILARITY

                aligned_bs_identityA = bs_identicalA / bs_alignedA
                aligned_bs_identityB = bs_identicalB / bs_alignedB
                aligned_bs_idenity = min(aligned_bs_identityA, aligned_bs_identityB)  # ALIGNED BS IDENTITY

                full_identityA = identicalA / aln_lenA
                full_identityB = identicalB / aln_lenB
                full_identity = min(full_identityA, full_identityB)  # FULL SEQUENCE ALN

                # print SC1, SC3, SC3, SC4, SC5, SC6, bs_similarity
                scoring_map = {
                    'score1': SC1,  # was: zscore. now: model_score
                    'score2': SC2,  # was: model_minus_decoy_avg, now: 0-SKOLNICK potential
                    'score3': SC3,  # was: model_score, now: 0-DFIRE
                    'score4': SC4,
                    'score5': SC5,
                    'score6': SC6,

                    'random': random.random(),
                    'bs_similarity': bs_similarity,
                    'full_identity': full_identity,
                    'score_combined1': 0.067 * SC1 * 5.07 * bs_similarity,
                    'score_combined2': (1-full_identity) * 0.067 * SC1 + full_identity * 5.07 * bs_similarity,
                }

                score = scoring_map[scoring]
                # print "SCORE", score

                # scoring system
                # score = model_minus_avg  # SCORE(MODEL) - AVG(DECOYS)

                if same_complex_type is True:
                    if query_type != template_type:
                        continue

                if threshold_nsubunits is not None:
                    if not (threshold_nsubunits[0] <= template_nsubunits <= threshold_nsubunits[1]):
                        continue

                # print "OK"

                # Unless in calibration mode: make sure the template is not same as query (close homologs - by identity)
                if calibration is False:
                    # threshold_identity = 1.0
                    # threshold_identity = 0.90
                    # threshold_identity = 0.30
                    if identicalA / aln_lenA > threshold_identity:
                        continue
                    if identicalB / aln_lenB > threshold_identity:
                        continue
                    # we can further remove all yeast paralogous structures and only use orthologs
                    if pdbs is not None:
                        if tpl.split("|")[0] in pdbs:
                            continue
                else:
                    if full_identity < 0.98:
                        continue

                # make sure that Uniprot IDs for the predicted residues have good residues
                # otherwise residue predictions will not be comparable

                ####### A ########
                uniprot = queryB
                if file_format == 1:
                    if queryA not in locus_uniprot:
                        continue
                    uniprot = locus_uniprot[queryA]

                if uniprot not in goodres:
                    continue
                # print "A", siteA
                # z = (tpl.split("|")[0], tpl.split("|")[1].split("_")[0])
                # uniprotA, seq_offsetA, db_offsetA = pdb_protein[z]

                sA = siteA.split(";")
                sB = siteB.split(";")
                # bs_total_number_of_residues = len(sA) + len(sB)

                for s in sA:
                    x = s.split(",")  # ['E', '13', '56', 'E', '1364']
                    if x[3] in ('*', '-'):
                        continue
                    if int(x[2]) == 1:
                        continue
                    # DEBUGGING complex cases with INDELS in sequences of templates
                    # if int(x[1]) - seq_offsetA + db_offsetA != int(x[4]):
                    #     print "Score", score
                    #     print "Alignment for residue", x
                    #     print "Template:", z, uniprotA, "offset seq", seq_offsetA, "offset db", db_offsetA
                    #     print "Query:", queryA
                    #     print "side A:", "template", x[0], int(x[1]), "in Uniprot coord:", int(x[1]) - seq_offsetA + db_offsetA, "Query coord:", int(x[4])
                    #     print "-------"
                    resi = int(x[4])  # query residue index
                    if resi not in predres[uniprot] or predres[uniprot][resi] < score:
                        predres[uniprot][resi] = score
                        # print "PREDRES", uniprot, resi, predres[uniprot][resi]


                    # DUMP:
                    if resi not in pred_dump['score'][uniprot] or pred_dump['score'][uniprot][resi] < score:
                        pred_dump['score'][uniprot][resi] = score
                        pred_dump['template_type'][uniprot][resi] = template_type
                        pred_dump['template_nsubunits'][uniprot][resi] = template_nsubunits
                        pred_dump['model_bs_size'][uniprot][resi] = len(sB)
                        pred_dump['model_size'][uniprot][resi] = uniprot_nresidues[uniprot]
                        pred_dump['bs_sim'][uniprot][resi] = bs_similarityA
                        pred_dump['bs_aln_sim'][uniprot][resi] = aligned_bs_similarityA
                        pred_dump['full_id'][uniprot][resi] = full_identityA

                ####### same for B ########
                uniprot = queryB
                if file_format == 1:
                    if queryB not in locus_uniprot:
                        continue
                    uniprot = locus_uniprot[queryB]
                if uniprot not in goodres:
                    continue
                # print "A", siteA

                for s in sB:
                    x = s.split(",")
                    if x[3] in ('*', '-'):
                        continue
                    if int(x[2]) == 1:
                        continue
                    resi = int(x[4])  # query residue index
                    if resi not in predres[uniprot] or predres[uniprot][resi] < score:
                        predres[uniprot][resi] = score

                    # DUMP:
                    if resi not in pred_dump['score'][uniprot] or pred_dump['score'][uniprot][resi] < score:
                        pred_dump['score'][uniprot][resi] = score
                        pred_dump['template_type'][uniprot][resi] = template_type
                        pred_dump['template_nsubunits'][uniprot][resi] = template_nsubunits
                        pred_dump['model_size'][uniprot][resi] = uniprot_nresidues[uniprot]
                        pred_dump['model_bs_size'][uniprot][resi] = len(sB)
                        pred_dump['bs_sim'][uniprot][resi] = bs_similarityB
                        pred_dump['bs_aln_sim'][uniprot][resi] = aligned_bs_similarityB
                        pred_dump['full_id'][uniprot][resi] = full_identityB

                # print queryA, queryB, locus_uniprot[queryA], locus_uniprot[queryB]
            except:
                print "Error reading the whole matches file. Skipping the line with an error"
                raise

    if calibration is True:
        with open(calib_fname, 'w') as o:
            o.write('uniprot\tresi\n')
            for uniprot, resi_list in predres.iteritems():
                s = [str(r) for r in resi_list if r in goodres[uniprot]]
                if len(s) == 0:
                    continue
                o.write("{}\t{}\n".format(uniprot, ";".join(s)))
    else:
        if dump_fname is not None:
            with open(dump_fname, 'w') as o:
                o.write("label\tuniprot\tresi\tscore\tbs_similarity\tbs_aln_similarity\tfull_seq_id\tcomplex_type\tcomplex_nsubunits\tprotein_size\tbs_size\n")
                for uniprot, resi_list in predres.iteritems():
                    for resi in resi_list:
                        label = 0
                        if resi in goodres[uniprot]:
                            label = 1
                        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            label, uniprot, resi,
                            pred_dump['score'][uniprot][resi],
                            pred_dump['bs_sim'][uniprot][resi],
                            pred_dump['bs_aln_sim'][uniprot][resi],
                            pred_dump['full_id'][uniprot][resi],
                            pred_dump['template_type'][uniprot][resi],
                            pred_dump['template_nsubunits'][uniprot][resi],
                            pred_dump['model_size'][uniprot][resi],
                            pred_dump['model_bs_size'][uniprot][resi]))
    return predres


def load_calibration(fname):
    calibration = defaultdict(set)
    with open(fname) as f:
        for line in islice(f, 1, None):
            x, y = line.split("\t")
            y = y.strip()
            if len(y) == 0:
                continue
            calibration[x] = set(map(int, y.split(";")))
    return calibration

# def load_vidal_interactions(fname):
#     interactions = []
#     with open(fname) as f:
#         for line in islice(f, 0, None):
#             A, B = line.strip().split("\t")
#             pair = (A, B)
#             if A > B:
#                 pair = (B, A)
#             interactions.append(pair)
#     return set(interactions)

def scored_residue_labels(predres, goodres, badresnum=None):
    true_labels = []
    scores = []
    for uniprot in predres.iterkeys():
        for res, score in predres[uniprot].iteritems():
            label = 0
            if res in goodres[uniprot]:
                label = 1
            true_labels.append(label)
            scores.append(score)
    return true_labels, scores


# def add_pseudo_scores(labels, predres, goodres, badresnum):
#     """
#     in goodres, not in predres (FN) => min score, label 1
#     badresnum - (goodres | predres) => min score, label 0
#     in predres TP or FP => real score
#     """
#     lab, scores = labels
#     min_score = min(scores)
#     max_score = max(scores)
#     new_lab = []
#     new_scores = []
#     # galse negatives
#     for uniprot in goodres:
#         for r in goodres[uniprot]:
#             if r not in predres[uniprot]:
#                 new_lab.append(1)
#                 new_scores.append(min_score)
#         for r in goodres[uniprot] | predres[uniprot]:


def plots_sens_spec(r, fname, title=None):
    from sklearn.metrics import roc_curve, auc
    import pylab as pl

    test_labels, scores = r
    fpr, tpr, thresholds = roc_curve(test_labels, scores)
    sens = tpr
    spec = [1-x for x in fpr]

    pl.clf()
    # remove first and last entries (these were artificially added by roc_curve() calculations)
    pl.plot(thresholds[1:-1], sens[1:-1], color="green")
    pl.plot(thresholds[1:-1], spec[1:-1], color="red")
    pl.ylim([0.0, 1.0])
    pl.xlabel('Score')
    pl.ylabel('Sensitivity (green) and Specificity (red)')
    pl.title(title)
    pl.savefig(fname)


def plots_score(predres, goodres, fname, title=None):
    from sklearn.metrics import roc_curve, auc
    import pylab as pl
    vals = []
    for resscore in predres.itervalues():
        for score in resscore.itervalues():
            vals.append(score)
    count = len(vals)
    min_score = min(vals)
    max_score = max(vals)
    print "Count", count, min_score, max_score
    thresholds = []
    tpr = []
    tnr = []
    threshold = min_score
    step = 0.05
    while threshold < max_score:
        tp = 0
        tn = 0
        print "T", threshold,
        for uniprot, res_score in predres.iteritems():
            for res, score in res_score.iteritems():
                if score >= threshold:
                    if res in goodres[uniprot]:
                        tp += 1
                else:
                    if res not in goodres[uniprot]:
                        tn += 1
        print tp, tn, tp / float(count), tn / float(count)
        tpr.append(tp / float(count))
        tnr.append(tn / float(count))
        thresholds.append(threshold)
        threshold += step

    # test_labels, scores = r
    # fpr, tpr, thresholds = roc_curve(test_labels, scores)

    # Plot ROC curve
    pl.clf()
    pl.plot(thresholds, tnr, color="red")
    # pl.plot([0, 1], [0, 1], 'k--')
    pl.plot(thresholds, tpr, color="green")
    # pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('Score')
    pl.ylabel('Sensitivity (green) and Specificity (red)')
    pl.title(title)
    # pl.legend(loc="lower right")
    pl.savefig(fname)

##########################################
if __name__ == '__main__':
    d = "/Users/agoncear/projects/Interactome/Workflow"
    roc_d = d + "/Benchmarks/Yeast_residues"
    # matches_fname = d + "/Alignments/matches_yeast.tab"
    matches_fname = d + "/Alignments/matches_Scerevisiae.tab"
    calib_fname = roc_d + "/CALIB.tab"

    locus_uniprot, uniprot_nresidues = load_uniprot()
    pdbs = load_structures()

    try:
        with open(calib_fname) as f:
            pass
    except:
        print "Calibrating..."
        goodres_ = load_templates(pdbs)  # Contains error
        get_predicted_residues(
            matches_fname, locus_uniprot, goodres_, uniprot_nresidues,
            pdbs=None, calibration=True, calib_fname=calib_fname,
            scoring='score3',
            threshold_identity=1.0, threshold_nsubunits=None, same_complex_type=True)

    print "Loading calibration..."
    goodres = load_calibration(calib_fname)

    print "Evaluating..."
    goodresnum = {}
    badresnum = {}
    for uniprot in goodres.iterkeys():
        goodresnum[uniprot] = len(goodres[uniprot])
        if uniprot not in uniprot_nresidues:
            # print "Protein not in Swissprot (possibly not Yeast or not the right strain, or not reviewed) -> gene unknown", uniprot
            continue
        badresnum[uniprot] = uniprot_nresidues[uniprot] - goodresnum[uniprot]

    print "GOOD", sum(goodresnum.itervalues())
    print "BAD", sum(badresnum.itervalues())

    # for prefix in ('', 'dimers', 'twilight', 'twilight_dimers'):
    for prefix in ('', 'twilight', 'dark', ):
        scoring_labels = []
        # scorings = ("score1", "score2", "score3", "random", "bs_similarity", "full_identity")
        # scoring_legends = 'Our score', 'Skolnick', 'DFIRE', 'random', 'bs similarity', 'tot. seq. identity'
        # scorings = ("score_combined1", "score_combined2", "score1", "score2", "score3", "bs_similarity")
        scorings = ("score_combined2", "score2", "score3", "bs_similarity")
        # scoring_legends = ('comb1', 'comb2', 'Our score', 'Skolnick', 'DFIRE', 'bs similarity')
        scoring_legends = ('Our combined score', 'Skolnick potential', 'DFIRE potential', 'BS similarity')
        for scoring in scorings:
            # scoring = "score1"
            print "Predicting with (scoring = {}, prefix = {})".format(scoring, prefix)

            TID = 0.98
            if prefix.startswith("twilight"):
                TID = 0.25

            if prefix.startswith("dark"):
                TID = 0.15

            TSUB = None
            if prefix.endswith("dimers"):
                TSUB = (2, 2)

            predres = get_predicted_residues(
                matches_fname, locus_uniprot, goodres, uniprot_nresidues,
                pdbs=None, calibration=False, calib_fname=None,
                scoring=scoring,
                threshold_identity=TID, threshold_nsubunits=TSUB, same_complex_type=True)

            predresnum = {uniprot: len(resdict.keys()) for uniprot, resdict in predres.iteritems()}

            print "GOOD", sum(goodresnum.itervalues())
            print "BAD", sum(badresnum.itervalues())
            print "PRED", sum(predresnum.itervalues())

            labels = scored_residue_labels(predres, goodres, badresnum)
            scoring_labels.append(labels)
            roc(labels, roc_d + "/{}_{}_ROC.png".format(prefix, scoring), "Yeast benchmark, {} {}".format(prefix, scoring))
            plots_sens_spec(labels, roc_d + "/{}_{}_SS.png".format(prefix, scoring), "Yeast benchmark, {} {}".format(prefix, scoring))

        multiroc(
            scoring_labels,
            scoring_legends,
            roc_d + "/{}_yeast_residues_combined.png".format(prefix),  #.eps
            "Yeast residue-level benchmark")

    # plots_score(predres, goodres, roc_d + "/scores.png", "Scores")

    # # Uniprot pairs
    # print "Processing the interactome"
    # # max pairs, all pairs
    # pairs, a_pairs = get_scored_pairs(matches_fname)

    # print "Comparing observed with predicted: Union"
    # # r_identity = scored_labels(pairs["identity"], pairs_vidal)
    # # r_bs_identity = scored_labels(pairs["bs_identity"], pairs_vidal)
    # # r_positive = scored_labels(pairs["positive"], pairs_vidal)
    # r_bs_positive = scored_labels(pairs["bs_positive"], pairs_vidal)
    # # r_bs_coverage = scored_labels(pairs["bs_coverage"], pairs_vidal)
    # # r_score = scored_labels(pairs["model_score"], pairs_vidal)
    # # r_zscore = scored_labels(pairs["model_zscore"], pairs_vidal)
    # r_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_vidal)
    # # r_random = scored_labels(pairs["random"], pairs_vidal)

    # print "Comparing observed with predicted: Golden"
    # # g_identity = scored_labels(pairs["identity"], pairs_golden)
    # # g_bs_identity = scored_labels(pairs["bs_identity"], pairs_golden)
    # # g_positive = scored_labels(pairs["positive"], pairs_golden)
    # g_bs_positive = scored_labels(pairs["bs_positive"], pairs_golden)
    # # g_bs_coverage = scored_labels(pairs["bs_coverage"], pairs_golden)
    # # g_score = scored_labels(pairs["model_score"], pairs_golden)
    # # g_zscore = scored_labels(pairs["model_zscore"], pairs_golden)
    # g_model_minus_avg = scored_labels(pairs["model_minus_avg"], pairs_golden)
    # # g_random = scored_labels(pairs["random"], pairs_golden)

    # # print "Comparing observed with predicted: Union"
    # # ra_identity = scored_labels(a_pairs["identity"], pairs_vidal)
    # # ra_bs_identity = scored_labels(a_pairs["bs_identity"], pairs_vidal)
    # # ra_positive = scored_labels(a_pairs["positive"], pairs_vidal)
    # # ra_bs_positive = scored_labels(a_pairs["bs_positive"], pairs_vidal)
    # # ra_bs_coverage = scored_labels(a_pairs["bs_coverage"], pairs_vidal)
    # # ra_score = scored_labels(a_pairs["model_score"], pairs_vidal)
    # # ra_zscore = scored_labels(a_pairs["model_zscore"], pairs_vidal)
    # # ra_model_minus_avg = scored_labels(a_pairs["model_minus_avg"], pairs_vidal)
    # # ra_random = scored_labels(a_pairs["random"], pairs_vidal)

    # # print "Comparing observed with predicted: Golden"
    # # ga_identity = scored_labels(a_pairs["identity"], pairs_golden)
    # # ga_bs_identity = scored_labels(a_pairs["bs_identity"], pairs_golden)
    # # ga_positive = scored_labels(a_pairs["positive"], pairs_golden)
    # # ga_bs_positive = scored_labels(a_pairs["bs_positive"], pairs_golden)
    # # ga_bs_coverage = scored_labels(a_pairs["bs_coverage"], pairs_golden)
    # # ga_score = scored_labels(a_pairs["model_score"], pairs_golden)
    # # ga_zscore = scored_labels(a_pairs["model_zscore"], pairs_golden)
    # # ga_model_minus_avg = scored_labels(a_pairs["model_minus_avg"], pairs_golden)
    # # ga_random = scored_labels(a_pairs["random"], pairs_golden)

    # # roc(g_bs_positive, roc_d + "/roc_bs_positive.png", "Yeast, golden set (best model), BS positive (BLOSUM62)")
    # # roc(g_bs_coverage, roc_d + "/roc_bs_coverage.png", "Yeast, golden set (best model), BS coverage")
    # # roc(g_score, roc_d + "/roc_score.png", "Yeast, golden set (best model), Score")
    # # roc(g_zscore, roc_d + "/roc_zscore.png", "Yeast, golden set (best model), Z-score")
    # # roc(g_model_minus_avg, roc_d + "/roc_model_minus_avg.png", "Yeast, golden set (best model), Score-avg(decoys)")

    # print "Making plots"
    # # legend = ('Y2H union (best model)', 'Gold set from literature (best model)', 'Y2H union (all models)', 'Gold set from literature (all models)')
    # legend = ('Y2H union (Yu et.al.)', 'Gold set from literature (Yu et.al.)')
    # # multiroc((r_identity, g_identity, ra_identity, ga_identity), legend, roc_d + "/identity.png", "Yeast: Alignment identity")
    # # multiroc((r_bs_identity, g_bs_identity, ra_bs_identity, ga_bs_identity), legend, roc_d + "/bs_identity.png", "Yeast: BS alignment identity")
    # # multiroc((r_positive, g_positive, ra_positive, ga_positive), legend, roc_d + "/positive.png", "Yeast: Alignment BLOSUM62 positive")
    # # multiroc((r_bs_positive, g_bs_positive, ra_bs_positive, ga_bs_positive), legend, roc_d + "/bs_positive.png", "Yeast: BS alignment BLOSUM62 positive")
    # multiroc((r_bs_positive, g_bs_positive), legend, roc_d + "/bs_positive.png", "Yeast: Binding site similarity score")
    # # # multiroc((r_bs_positive_homo, g_bs_positive_homo), legend, roc_d + "/bs_positive_homo.png", "Yeast: Preserve homo/hetero type. BS BLOSUM62 positive")
    # # multiroc((r_bs_coverage, g_bs_coverage, ra_bs_coverage, ga_bs_coverage), legend, roc_d + "/bs_coverage.png", "Yeast: BS coverage")
    # # multiroc((r_score, g_score, ra_score, ga_score), legend, roc_d + "/score.png", "Yeast: Compatibility score 1")
    # # multiroc((r_zscore, g_zscore, ra_zscore, ga_zscore), legend, roc_d + "/zscore.png", "Yeast: Compatibility Z-score")
    # # multiroc((r_model_minus_avg, g_model_minus_avg, ra_model_minus_avg, ga_model_minus_avg), legend, roc_d + "/model_minus_avg.png", "Yeast: Compatibility score-avg(decoy)")
    # # multiroc((r_random, g_random, ra_random, ga_random), legend, roc_d + "/random.png", "Yeast: Random control")
    # # # roc(r_random, roc_d + "/bs_random.png", "Random control")

    # # multiroc(
    # #     (r_identity, r_bs_identity, r_positive, r_bs_positive),
    # #     ('id', 'bs_id', 'pos', 'bs_pos'),
    # #     roc_d + "/id_pos.png",
    # #     "Basic ROC in Yeast, Y2H_union")
    # multiroc(
    #     (r_bs_positive, g_bs_positive, r_model_minus_avg, g_model_minus_avg),
    #     ('Bind. site similarity (Y2H)', 'Bind. site similarity (Gold set)', 'Interface compatibily (Y2H)', 'Interface compatibility (Gold set)'),
    #     roc_d + "/yeast_combined.eps",
    #     "Yeast, Y2H and Gold set (Literature)")
