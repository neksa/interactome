# import copy
# import pprint
import math
import signal
from itertools import islice
from collections import defaultdict, namedtuple
from multiprocessing import Pool

import numpy as np
import random
# MIN: 20%, 0, 0.5: Total alignments: 1115249

debug = False
# debug = True
# if debug:
#     import cProfile

# from interactome.sequences.blast import BLASTReport
from interactome.structures.complexes import Complexes  # , SiteResidue
# from interactome.sequences.BLOSUM import read_BLOSUM_matrix


def get_root():
    return "/Users/agoncear/projects/Interactome/Workflow"


def site_to_list(site):
    # print  [s.split(",") for s in site.split(";")]
    return map(lambda s: (s[0], int(s[1]), int(s[2]), s[3], int(s[4])), [s.split(",") for s in site.split(";")])


def site_to_dict(site_str):
    # print "A", {s for site in site_str.split(';') for s in site.split(',')}
    # return {(s[0], int(s[1])): (s[3], int(s[4]), int(s[2])) for site in site_str.split(';') for s in site.split(',')}
    # d = {}
    # gen = (site for site in site_str.split(';'))
    # site.split(';') for site in gen]

    # 4240836 of alignment consume 32 Gb RAM
    # we convert strings to ints with ord()

    d = {}
    for site in site_str.split(";"):
        s = site.split(",")
        d[(s[0], int(s[1]))] = (s[3], int(s[4]), int(s[2]))  # substitution: from key to value: resn, resi, ncontacts (temp. removed)
    return d


def dict_to_site(d):
    l = []
    for key, value in d.iteritems():
        resn, resi = key
        subst_resn, subst_resi, ncontacts = value
        l.append("{},{},{},{},{}".format(resn, resi, ncontacts, subst_resn, subst_resi))
    return ";".join(l)


TPLPROP = namedtuple('TPLPROP', 'ctype nsubunits redundant bs_lenA bs_lenB ncontacts')
def load_complex_props(fname):
    complex_props= {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            tpl, pdb, chA, chB, protA, protB, ctype, nsubunits, ndirect_int, redundant, bs_lenA, ncontactsA, bs_lenB, ncontactsB = line.strip().split("\t")
            redundant = bool(int(redundant))
            nsubunits = int(nsubunits)
            bs_lenA = int(bs_lenA)
            bs_lenB = int(bs_lenB)
            ncontacts = int(ncontactsA)
            complex_props[tpl] = TPLPROP(
                ctype=ctype, nsubunits=nsubunits, redundant=redundant, bs_lenA=bs_lenA, bs_lenB=bs_lenB, ncontacts=ncontacts)
    return complex_props


ALN = namedtuple('ALN', 'query template identical positive aln_len bs_len bs_covered bs_aligned bs_identical bs_positive bs_contacts bs_BLOSUM bs_score1 site')
def load_alignments(fname, templates=None, bs_cov_threshold=0.7, aln_identity_threshold=0.2, skip_problematic=None, sel_proteins=None):
    global debug
    alignments = defaultdict(list)
    with open(fname) as f:
        limit = None
        if debug:
            limit = 100000

        for c, line in enumerate(islice(f, 1, limit)):
            if c % 100000 == 0:
                print c
            fields = line.strip().split("\t")
            fields[2:11] = map(int, fields[2:11])
            fields[11:13] = map(float, fields[11:13])
            # fields[-1] = site_to_list(fields[-1])
            fields[-1] = site_to_dict(fields[-1])
            # print fields
            aln = ALN._make(fields)

            # check that the protein is in the set of selected proteins
            if sel_proteins and aln.query not in sel_proteins:
                continue
            # check the coverage
            if (aln.bs_covered / float(aln.bs_len)) < bs_cov_threshold:
                continue

            # TEMP:
            if aln.bs_identical == 0:
                continue
            if aln.identical / float(aln.aln_len) < aln_identity_threshold:
                continue

            if templates is not None and aln.template not in templates:
                # print "skip tpl", aln.template
                continue

            # skip alignments where SEQRES and ATOMRES where not the same (or problems in gi SEQRES at NCBI in BLAST pdbaa database)
            if skip_problematic is not None:
                skip = False
                for s in aln.site.itervalues():
                    if s[0] == '@':
                        skip = True
                        break
                if skip:
                    continue

            alignments[aln.template].append(aln)
    return alignments


def ca_to_bin(ca):
    b = int(round(ca)) - 2
    if b < 0:
        return 0
    elif b > 15-2:
        return 15-2
    return b


def load_potential(fname):
    potential = {}
    with open(fname) as f:
        for line in f:
            fields = line.strip().split("\t")
            aa = fields[0]
            scores = fields[1:]
            scores_float = map(float, scores)
            potential[aa] = scores_float
    return potential


def load_mutipot_potential(fname):
    # aaL = "GLY ALA SER CYS VAL THR ILE PRO MET ASP ASN LEU LYS GLU GLN ARG HIS PHE TYR TRP".split()
    # aaS = "G A S C V T I P M D N L K E Q R H F Y W".split()
    # aa_map = dict(zip(aaL, aaS))
    aa_map = "GASCVTIPMDNLKEQRHFYW"
    potential = {}
    with open(fname) as f:
        for i, line in enumerate(islice(f, 1, None)):
            for j, v in enumerate(islice(line.split(), 1, None)):
                potential[aa_map[i] + aa_map[j]] = [float(v)] * 15  # repeat 15 times the same value
    return potential


def load_dfire_potential(fname, atom):
    aaL = "GLY ALA SER CYS VAL THR ILE PRO MET ASP ASN LEU LYS GLU GLN ARG HIS PHE TYR TRP".split()
    aaS = "G A S C V T I P M D N L K E Q R H F Y W".split()
    aa_map = dict(zip(aaL, aaS))
    potential = {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            s = line.split()
            rA, aA, rB, aB = s[0:4]
            rA = aa_map[rA]
            rB = aa_map[rB]
            if aA != atom or aB != atom:
                continue
            potential[rA + rB] = map(float, s[4:])
            potential[rB + rA] = potential[rA + rB]
    return potential


def get_atoms(pdb, chainA, chainB):
    pdb_file = get_root() + '/Interfaces/' + pdb[1:3].lower() + '/' + pdb.lower() + '_interface_residue_atoms_5.0A.tab'
    atoms_A = defaultdict(dict)
    atoms_B = defaultdict(dict)
    with open(pdb_file) as f:
        for line in islice(f, 1, None):
            chain, resn, resi, atomn, x, y, z = line.strip().split()
            resi = int(resi)
            coord = np.array(map(float, (x, y, z)))
            # if atomn == "CA":
            if chain == chainA:
                atoms_A[(resn, resi)][atomn] = coord
            if chain == chainB:
                atoms_B[(resn, resi)][atomn] = coord
    return atoms_A, atoms_B


def calc_multiple_scores(siteA, siteB, atomsA, atomsB, alnA, alnB):
    """ We have to calculate 6 compatibility scores """
    global potentials

    max_distance = 15.0
    distances = [x for x in ca_distances(siteA, siteB, atomsA, atomsB, max_distance)]

    score1 = our_score(potentials['our'], siteA, siteB, atomsA, atomsB, alnA, alnB, distances)
    score2 = mutipot_score(potentials['mutipot'], siteA, siteB, atomsA, atomsB, alnA, alnB, distances)
    score3 = dfire_score(potentials['dfire_ca'], siteA, siteB, atomsA, atomsB, alnA, alnB, distances)
    score4 = .0
    score5 = .0
    score6 = .0
    # print "S", score1, score2

    # if error
    if score1 is None:
        return None
    return score1, score2, score3, score4, score5, score6


def ca_distances(siteA, siteB, atomsA, atomsB, max_distance):
    for k, rA in enumerate(siteA):
        Ca_A = atomsA.get((rA.resn, rA.resi), {}).get("CA")
        if Ca_A is None:
            continue
        for l, rB in enumerate(siteB):
            Ca_B = atomsB.get((rB.resn, rB.resi), {}).get("CA")
            if Ca_B is None:
                continue
            d_Ca = math.sqrt(
                (Ca_A[0] - Ca_B[0]) ** 2 +
                (Ca_A[1] - Ca_B[1]) ** 2 +
                (Ca_A[2] - Ca_B[2]) ** 2)
            if d_Ca > max_distance:
                continue
            yield (k, l), d_Ca


aa_short = set("ALPGDNYHECFVIRTKSQMW")
blacklist = set("-*@")
def our_score(potential, siteA, siteB, atomsA, atomsB, alnA, alnB, distances=None):
    model_score = .0
    if distances is None:
        max_distance = 15.0
        distances = ca_distances(siteA, siteB, atomsA, atomsB, max_distance)
    bins = {k: ca_to_bin(d) for k, d in distances}

    def calc(resnA, resnB):
        total_score = .0
        summations = 0
        for (k, l), d in bins.iteritems():
            a = resnA[k]
            b = resnB[l]
            # print a, b, d
            if a not in aa_short or b not in aa_short:
                continue
            # if a in blacklist or b in blacklist:
            #     continue
            # if d > limit for distance between a and b: continue
            total_score += potential[a+b][d]  # resnA resnB distance_bin
            summations += 1
        return total_score, summations

    resnA = [alnA.site[(s.resn, s.resi)][0] for s in siteA]
    resnB = [alnB.site[(s.resn, s.resi)][0] for s in siteB]
    model_score, pairs_summed_up = calc(resnA, resnB)

    # DECOYS:
    # resnA = rA[:]
    # resnB = rB[:]
    # if pairs_summed_up > 0:
    #     decoys = []
    #     for i in xrange(100):
    #         if i % 2 == 0:
    #             random.shuffle(resnA)
    #         else:
    #             random.shuffle(resnB)
    #         decoy_score, _ = calc(resnA, resnB)
    #         decoys.append(decoy_score)
    #     model_avg_score = model_score - np.mean(decoys)
    #     model_zscore = model_avg_score / np.std(decoys)

    if pairs_summed_up == 0:
        return None
    else:
        return model_score


def mutipot_score(potential, siteA, siteB, atomsA, atomsB, alnA, alnB, distances=None):
    model_score = .0
    if distances is None:
        max_distance = 15.0
        distances = ca_distances(siteA, siteB, atomsA, atomsB, max_distance)
    bins = {k: 0 for k, d in distances}

    def calc(resnA, resnB):
        total_score = .0
        summations = 0
        for (k, l), d in bins.iteritems():
            a = resnA[k]
            b = resnB[l]
            if a not in aa_short or b not in aa_short:
                continue
            # if d > limit for distance between a and b: continue
            total_score += potential[a+b][d]  # resnA resnB distance_bin
            summations += 1
        return total_score, summations

    resnA = [alnA.site[(s.resn, s.resi)][0] for s in siteA]
    resnB = [alnB.site[(s.resn, s.resi)][0] for s in siteB]
    model_score, pairs_summed_up = calc(resnA, resnB)

    if pairs_summed_up == 0:
        return None
    else:
        return model_score


def dfire_score(potential, siteA, siteB, atomsA, atomsB, alnA, alnB, distances=None):
    model_score = .0
    if distances is None:
        max_distance = 15.0
        distances = ca_distances(siteA, siteB, atomsA, atomsB, max_distance)
    bins = {k: int(d * 2) for k, d in distances}

    def calc(resnA, resnB):
        total_score = .0
        summations = 0
        for (k, l), d in bins.iteritems():
            a = resnA[k]
            b = resnB[l]
            if a not in aa_short or b not in aa_short:
                continue
            # if d > limit for distance between a and b: continue
            total_score += potential[a+b][d]  # resnA resnB distance_bin
            summations += 1
        return total_score, summations

    resnA = [alnA.site[(s.resn, s.resi)][0] for s in siteA]
    resnB = [alnB.site[(s.resn, s.resi)][0] for s in siteB]
    model_score, pairs_summed_up = calc(resnA, resnB)

    if pairs_summed_up == 0:
        return None
    else:
        return model_score


########################################################################################################################
def runner_parallel(template):
    """ This is a collector function. It joins the results of runner generator in a packet of results for each template """
    return [x for x in runner(template)]


def runner(template):
    global alignments
    # pdb, chainA, chainB = template.split("|")
    (pdb, chainA, chainB), (site1, site2) = template

    tplA = pdb.upper() + '|' + chainA
    tplB = pdb.upper() + '|' + chainB
    # print "T:", tplA, tplB

    all_alnA = alignments[tplA]
    all_alnB = alignments[tplB]
    if len(all_alnA) == 0 or len(all_alnB) == 0:
        return

    atomsA = None
    atomsB = None

    site1_resi = set([s.resi for s in site1])
    site2_resi = set([s.resi for s in site2])

    for alnA in all_alnA:
        # if it is an alignment for a different site in the same template -> skip it
        alnA_resi = set([a[1] for a in alnA.site.keys()])
        if alnA_resi != site1_resi:
            continue

        for alnB in all_alnB:
            # if it is an alignment for a different site in the same template -> skip it
            alnB_resi = set([a[1] for a in alnB.site.keys()])
            if alnB_resi != site2_resi:
                continue

            # delayed loading of atoms, load only once
            if atomsA is None or atomsB is None:
                atomsA, atomsB = get_atoms(pdb, chainA, chainB)

            scores = calc_multiple_scores(site1, site2, atomsA, atomsB, alnA, alnB)
            if scores is not None:
                yield '|'.join((pdb, chainA, chainB)), alnA.query, alnB.query, scores, alnA, alnB


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


alignments = None
potentials = None
def predict(templates_fname, alignment_fname, matches_fname):
    global alignments
    global potentials

    pdb_templates = get_root() + "/Structures/" + templates_fname
    alignments_fname = get_root() + "/BLAST-results/" + alignment_fname
    matches_fname = get_root() + "/Alignments/" + matches_fname
    template_analysis_fname = get_root() + "/Structures/template_analysis.tab"
    template_pdb_blacklist = get_root() + "/Structures/pdb_blacklist.tab"

    print "Loading potentials..."
    potentials = {}
    potential_fname = get_root() + "/Potential/potential_1.index-warp"
    potential_mutipot_fname = get_root() + "/Potential/mutipot_lu.index.txt"
    potential_dfire_fname = get_root() + "/Potential/dfire_pair.lib"
    potentials['our'] = load_potential(potential_fname)
    potentials['mutipot'] = load_mutipot_potential(potential_mutipot_fname)
    potentials['dfire_ca'] = load_dfire_potential(potential_dfire_fname, 'CA')
    # print potential["EE"][ca_to_bin(6.5)]

    print "Loading templates..."
    complexes = Complexes()
    # templates = complexes.loadTemplates(pdb_templates)
    iter_templates = complexes.loadTemplates(pdb_templates)

    # Blacklist contains template structure with huge number of sequence alignments (>200 in human)
    pdb_blacklist = set()
    with open(template_pdb_blacklist) as f:
        for line in islice(f, 1, None):
            # print line.split()
            pdb, n = line.split()
            n = int(n)
            if n > 200:
                pdb_blacklist.add(pdb)

    chains = set()
    templates = []
    for t in iter_templates:
        (pdb, chainA, chainB), (site1, site2) = t
        pdb = pdb.upper()
        if pdb in pdb_blacklist:
            # print "skip blacklisted template", (pdb, chainA, chainB)
            continue
        tplA = pdb + '|' + chainA
        tplB = pdb + '|' + chainB
        chains.add(tplA)
        chains.add(tplB)
        templates.append(t)
    print "Total templates (interfaces):", len(templates)

    print "Loading complex properties..."
    complex_props = load_complex_props(template_analysis_fname)
    # print complex_props

    print "Loading alignments..."
    # alignments = load_alignments(alignments_fname, bs_cov_threshold=0.5)  # , sel_proteins=selected_proteins)
    # THRESHOLDS REMOVED:
    alignments = load_alignments(alignments_fname, templates=chains, bs_cov_threshold=0.5, aln_identity_threshold=0.0, skip_problematic=True)  # , sel_proteins=selected_proteins)
    print "Total alignments:", sum((len(x) for x in alignments.itervalues()))

    pool = None
    per_template_results = None
    if debug:
        per_template_results = (runner_parallel(template) for template in islice(templates, 0, 100))  # sequential (good for debugging)
    else:
        N_CPU = 7
        pool = Pool(N_CPU, init_worker)
        per_template_results = pool.imap_unordered(runner_parallel, templates)  # parallel

    c = 0
    print "Writing results to", matches_fname
    with open(matches_fname, 'w') as o:
        header = "query_A query_B template query_type template_type template_nsubunits " +\
            "score1 score2 score3 score4 score5 score6 " +\
            "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +\
            "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +\
            "siteA siteB"
        o.write("\t".join(header.split(" ")) + "\n")
        try:
            for results_block in per_template_results:
                for result in results_block:
                    if c % 1000 == 0:
                        print c
                    c += 1
                    tpl, queryA, queryB, scores, alnA, alnB = result

                    if scores is None:
                        continue

                    score1, score2, score3, score4, score5, score6 = map(lambda x: round(x, 4), scores)  # round 4

                    props = complex_props.get(tpl.upper(), None)
                    if props is None:
                        continue
                    template_nsubunits = props.nsubunits
                    if props.redundant:
                        continue

                    template_type = "Unknown"
                    if props.ctype is not None:
                        template_type = props.ctype

                    query_type = "Hetero"
                    if queryA == queryB:
                        query_type = "Homo"

                    if queryA.startswith("gi"):
                        queryA = queryA.split("|")[1]
                        queryB = queryB.split("|")[1]
                    else:
                        queryA = queryA.split()[0]
                        queryB = queryB.split()[0]

                    # print queryA, queryB, score1, score2, score3, score4, score5, score6
                    output = (
                        queryA, queryB, tpl,
                        query_type, template_type, template_nsubunits,
                        score1, score2, score3, score4, score5, score6,
                        alnA.identical, alnA.positive, alnA.aln_len, alnA.bs_len, alnA.bs_covered, alnA.bs_aligned, alnA.bs_identical, alnA.bs_positive, alnA.bs_contacts, alnA.bs_BLOSUM, alnA.bs_score1,
                        alnB.identical, alnB.positive, alnB.aln_len, alnB.bs_len, alnB.bs_covered, alnB.bs_aligned, alnB.bs_identical, alnB.bs_positive, alnB.bs_contacts, alnB.bs_BLOSUM, alnB.bs_score1,
                        dict_to_site(alnA.site), dict_to_site(alnB.site))
                    # print len(output)  # 35 ?
                    format_str = "\t".join(['{}']*len(output)) + "\n"
                    o.write(format_str.format(*output))

        except KeyboardInterrupt:
            print "Caught KeyboardInterrupt, terminating workers"
            if pool is not None:
                pool.terminate()
                pool.join()
                pool = None
            else:
                raise
        if pool is not None:
            pool.close()
            pool.join()

