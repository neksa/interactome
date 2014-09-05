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
from itertools import islice
# from collections import defaultdict, namedtuple
# import numpy as np
from interactome.structures.complexes import Complexes  # , SiteResidue

from interactome.workflows.predict_interactions import get_root, load_potential, get_atoms, aa_short, load_complex_types, load_alignments, ca_to_bin


def score(scoring_function, siteA, siteB, resnA, resnB, atomsA, atomsB, distance_cutoff=15.0):
    max_distance = distance_cutoff
    total_score = .0
    total_template_score = .0

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

            res_template = rA.resn + rB.resn
            res = resnA_substituted + resnB_substituted
            # bin_num = ca_to_bin(d_Ca)

            s_template = scoring_function(res_template, d_Ca)
            total_template_score += s_template
            
            s = scoring_function(res, d_Ca)
            total_score += s
            # print "LB={}\tA={}\tB={}\tNcont_A={}\tNcont_B={}\td_CA={}\tScore f={}".format(label, (rA.resn, rA.resi), (rB.resn, rB.resi), resnA_ncontacts, resnB_ncontacts, d_Ca, s)

    return total_score  #, total_template_score


def score_contact_potential1(res_pair, d_Ca):
    global potential
    bin_num = ca_to_bin(d_Ca)
    s = potential[res_pair][bin_num]
    return s


def run():
    global alignments
    global potential

    templates_fname = get_root() + "/Structures/pdb_templates_5A.tab"
    template_analysis_fname = get_root() + "/Structures/template_analysis.tab"
    alignments_fname = get_root() + "/BLAST-results/human_deltablast_report_merged-bs.tab"
    potential_fname = get_root() + "/Potential/potential_1.index-warp"
    results_fname = get_root() + "/ScoringBenchmark/scores1_10.tab"

    print "Loading potential..."
    potential = load_potential(potential_fname)
    print "Loading templates..."
    complexes = Complexes()
    templates = complexes.loadTemplates(templates_fname)
    templates_by_pdb = {}
    for tpl in islice(templates, 0, 30000):
        (pdb, chainA, chainB), (site1, site2) = tpl
        key = pdb.upper() + '|' + chainA + '|' + chainB
        templates_by_pdb[key] = site1, site2
    # print "Loading complex properties..."
    # complex_types = load_complex_types(template_analysis_fname)
    print "Loading alignments..."

    # a couple of human protein complexes (observed!)
    # selected_proteins = set('P30153', 'P67775', 'P55957', 'Q07812')
    # selected_pdbs = ('3C5W|A|C', '4BD2|A|C')
    # selected_pdbs = ('3C5W|A|B', '4BD2|A|B')  # mmcif names

    # alignments = load_alignments(alignments_fname, bs_cov_threshold=0.5, sel_proteins=selected_proteins)
    # print "Total alignments:", sum((len(x) for x in alignments.itervalues()))

    print "Writing results to", results_fname
    with open(results_fname, 'w') as o:
        header = "AB orig_AB score"
        o.write("\t".join(header.split(" ")) + "\n")

        scoring_function = score_contact_potential1
        distance_cutoff = 10.0

        # ALN_test = alignments[0]
        native_structure = '3C5W|A|B'
        # native_structure = '4BD2|A|B'
        params = native_structure.split("|")
        native_atomsA, native_atomsB = get_atoms(*params)
        native_siteA, native_siteB = templates_by_pdb[native_structure]
        native_resnA = [s.resn for s in native_siteA]
        native_resnB = [s.resn for s in native_siteB]
        score_native = score(scoring_function, native_siteA, native_siteB, native_resnA, native_resnB, native_atomsA, native_atomsB, distance_cutoff)
        print "score native", score_native
        o.write("{}\t{}\t{}\n".format(native_structure, native_structure, score_native))

        for tpl_key, value in templates_by_pdb.iteritems():
            pdb, chainA, chainB = tpl_key.split('|')
            site1, site2 = value

            print "scoring T:", tpl_key
            atomsA, atomsB = get_atoms(pdb, chainA, chainB)

            if len(site1) > len(native_resnA):
                continue
            if len(site2) > len(native_resnB):
                continue

            score_decoy = score(scoring_function, site1, site2, native_resnA, native_resnB, atomsA, atomsB, distance_cutoff)
            print "score decoy", score_decoy

            output = (tpl_key, native_structure, score_decoy)
            format_str = "\t".join(['{}']*len(output)) + "\n"
            o.write(format_str.format(*output))


if __name__ == '__main__':
    run()
