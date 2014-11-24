"""
1. Take the human structures with self-hits (100% identity hit)
2. For each structure A,B:
    3. score the native structure. report score S_native
    4. get the list of templates X:Z Y:W (each binding site not longer than in A, B)
       and (temlates not related to A or B)
    5. remap A,B to X,Y template
    6. score the decoy X:A, Y:B, report the scores
7. calculate z-scores relative to the native structure  z_native = s_native - mean(s) / sigma(s)

"""
import math
import signal
import random
from multiprocessing import Pool
from itertools import islice

from interactome.structures.complexes import Complexes  # , SiteResidue
from interactome.workflows.predict_interactions import get_root, load_potential, get_atoms, aa_short, ca_to_bin
#  load_complex_types, load_alignments
from interactome.structures.nr_filter import NRFilter


def calc_score(scoring_function, siteA, siteB, resnA, resnB, atomsA, atomsB, distance_cutoff=15.0):
    max_distance = distance_cutoff
    total_score = .0
    # total_template_score = .0

    for k, rA in enumerate(siteA):
        Ca_A = atomsA.get((rA.resn, rA.resi))
        if Ca_A is None:
            # raise Exception("Atoms for residues not found in scoring {} {}".format(rA.resn, rA.resi))
            print "Atoms for residues not found in scoring {} {}".format(rA.resn, rA.resi)
            continue
        # print alnA
        # alignment A:
        resnA_substituted = resnA[k]
        # if resnA_substituted in ('*', '-', '@'):
        skip_pair = False  # skip the pair of residues from the pairwise calculations in the aligned site
        if resnA_substituted not in aa_short:
            skip_pair = True

        for l, rB in enumerate(siteB):
            Ca_B = atomsB.get((rB.resn, rB.resi))
            if Ca_B is None:
                # raise Exception("Atoms for residues not found in scoring {} {}".format(rB.resn, rB.resi))
                print "Atoms for residues not found in scoring {} {}".format(rB.resn, rB.resi)
                continue

            # d_Ca = np.linalg.norm(Ca_A - Ca_B)  # too slow!
            d_Ca = math.sqrt(sum([(Ca_A[z] - Ca_B[z])**2 for z in range(3)]))
            if d_Ca > max_distance:
                continue

            # alignment B:
            resnB_substituted = resnB[l]
            # if resnB_substituted in ('*', '-', '@'):
            if resnB_substituted not in aa_short:
                skip_pair = True

            if skip_pair:
                continue

            # res_template = rA.resn + rB.resn
            res = resnA_substituted + resnB_substituted
            # bin_num = ca_to_bin(d_Ca)

            # s_template = scoring_function(res_template, d_Ca)
            # total_template_score += s_template

            s = scoring_function(res, d_Ca)
            total_score += s
            # print "LB={}\tA={}\tB={}\tNcont_A={}\tNcont_B={}\td_CA={}\tScore f={}".format(label, (rA.resn, rA.resi), (rB.resn, rB.resi), resnA_ncontacts, resnB_ncontacts, d_Ca, s)

    return total_score  # , total_template_score


def score_distance_potential1(res_pair, d_Ca):
    global potential
    bin_num = ca_to_bin(d_Ca)
    s = potential[res_pair][bin_num]
    return s


def load_mutipot_potential(fname):
    # GLY ALA SER CYS VAL THR ILE PRO MET ASP ASN LEU LYS GLU GLN ARG HIS PHE TYR TRP
    aa_order = "GASCVTIPMDNLKEQRHFYW"
    p = {}
    with open(fname) as f:
        for i, line in enumerate(islice(f, 1, None)):
            fields = line.strip().split()[1:]
            print fields
            for j in range(len(aa_order)):
                pair = aa_order[i] + aa_order[j]
                invpair = aa_order[j] + aa_order[i]
                p[pair] = float(fields[j])
                p[invpair] = p[pair]
    return p


def score_potential2(res_pair, d_Ca=None):
    global potential
    return potential[res_pair]


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def run_parallel_inverse(native):


    structure, value = native
    # structure = '3C5W|A|B'
    # structure = '4BD2|A|B'
    params = structure.split("|")
    atomsA, atomsB = get_atoms(*params)
    siteA, siteB = value  # templates_by_pdb[structure]
    resnA = [s.resn for s in siteA]
    resnB = [s.resn for s in siteB]
    score = calc_score(scoring_function, siteA, siteB, resnA, resnB, atomsA, atomsB, distance_cutoff)
    print "Native", structure, "score", score
    resnA = [s.resn for s in siteA]
    resnB = [s.resn for s in siteB]

    outputs = []
    output = (structure, score, "native")
    outputs.append(output)
    for i in range(300):
        random.shuffle(resnA)
        random.shuffle(resnB)
        score_decoy = calc_score(scoring_function, siteA, siteB, resnA, resnB, atomsA, atomsB, distance_cutoff)
        output = (structure, score_decoy, "decoy")
        outputs.append(output)
    return outputs


# def run_parallel(native):
#     global templates_by_pdb
#     # setup:
#     scoring_function = score_distance_potential1
#     distance_cutoff = 15.0

#     native_structure, native_value = native
#     # native_structure = '3C5W|A|B'
#     # native_structure = '4BD2|A|B'
#     params = native_structure.split("|")
#     native_atomsA, native_atomsB = get_atoms(*params)
#     native_siteA, native_siteB = native_value  # templates_by_pdb[native_structure]
#     native_resnA = [s.resn for s in native_siteA]
#     native_resnB = [s.resn for s in native_siteB]
#     score_native = score(scoring_function, native_siteA, native_siteB, native_resnA, native_resnB, native_atomsA, native_atomsB, distance_cutoff)
#     print "Native", native_structure, "score", score_native
#     # o.write("{}\t{}\t{}\n".format(native_structure, native_structure, score_native))

#     outputs = []
#     for tpl_key, value in templates_by_pdb.iteritems():
#         pdb, chainA, chainB = tpl_key.split('|')
#         site1, site2 = value

#         atomsA, atomsB = get_atoms(pdb, chainA, chainB)

#         if len(site1) > len(native_resnA):
#             continue
#         if len(site2) > len(native_resnB):
#             continue

#         score_decoy = score(scoring_function, site1, site2, native_resnA, native_resnB, atomsA, atomsB, distance_cutoff)
#         # print "scoring decoy T:", tpl_key,
#         # print "score", score_decoy

#         model_type = "decoy"
#         if native_structure == tpl_key:
#             model_type = "native"

#         output = (tpl_key, native_structure, score_decoy, model_type)
#         outputs.append(output)
#     return outputs


def run():
    global alignments
    global potential
    global templates_by_pdb

    global scoring_function
    global distance_cutoff
    scoring_function = score_distance_potential1
    distance_cutoff = 15.0

    templates_fname = get_root() + "/Structures/pdb_templates_5A.tab"
    # template_analysis_fname = get_root() + "/Structures/template_analysis.tab"
    # alignments_fname = get_root() + "/BLAST-results/human_deltablast_report_merged-bs.tab"
    potential_fname = get_root() + "/Potential/potential_1.index-warp"
    results_fname = get_root() + "/ScoringBenchmark/potential1.tab"

    print "Loading potential..."
    potential = load_potential(potential_fname)
    print "Loading templates..."
    NR = NRFilter()

    complexes = Complexes()
    templates = complexes.loadTemplates(templates_fname)
    templates_by_pdb = {}
    # for tpl in islice(templates, 0, 30000):
    for tpl in templates:
        (pdb, chainA, chainB), (site1, site2) = tpl

        if not NR.isNR(pdb):
            continue

        key = pdb.upper() + '|' + chainA + '|' + chainB
        templates_by_pdb[key] = site1, site2
    # print "Loading complex properties..."
    # complex_types = load_complex_types(template_analysis_fname)
    print "Loading alignments..."

    print "Writing results to", results_fname
    with open(results_fname, 'w') as o:
        header = "AB score model_type"
        o.write("\t".join(header.split(" ")) + "\n")

        pool = Pool(8, init_worker)
        # per_native_result = pool.imap_unordered(run_parallel, templates_by_pdb.iteritems())
        per_native_result = pool.imap_unordered(run_parallel_inverse, templates_by_pdb.iteritems())
        try:
            for results in per_native_result:
                for output in results:
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


def run_potential2():
    global alignments
    global potential
    global templates_by_pdb

    global scoring_function
    global distance_cutoff
    scoring_function = score_potential2
    distance_cutoff = 6.0

    templates_fname = get_root() + "/Structures/pdb_templates_5A.tab"
    potential_fname = get_root() + "/Potential/mutipot_lu.index.txt"
    results_fname = get_root() + "/ScoringBenchmark/potential2_6A.tab"

    print "Loading potential..."
    potential = load_mutipot_potential(potential_fname)
    print "Loading templates..."
    NR = NRFilter()

    complexes = Complexes()
    templates = complexes.loadTemplates(templates_fname)
    templates_by_pdb = {}
    # for tpl in islice(templates, 0, 30000):
    for tpl in templates:
        (pdb, chainA, chainB), (site1, site2) = tpl

        if not NR.isNR(pdb):
            continue

        key = pdb.upper() + '|' + chainA + '|' + chainB
        templates_by_pdb[key] = site1, site2
    # print "Loading complex properties..."
    # complex_types = load_complex_types(template_analysis_fname)
    print "Loading alignments..."

    print "Writing results to", results_fname
    with open(results_fname, 'w') as o:
        header = "AB score model_type"
        o.write("\t".join(header.split(" ")) + "\n")

        pool = Pool(6, init_worker)
        # per_native_result = pool.imap_unordered(run_parallel, templates_by_pdb.iteritems())
        per_native_result = pool.imap_unordered(run_parallel_inverse, templates_by_pdb.iteritems())
        try:
            for results in per_native_result:
                for output in results:
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


if __name__ == '__main__':
    # run()
    run_potential2()
