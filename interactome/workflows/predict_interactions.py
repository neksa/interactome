# import copy
# import pprint
import math
import signal
from itertools import islice
from collections import defaultdict, namedtuple
from multiprocessing import Pool

import numpy as np

# MIN: 20%, 0, 0.5: Total alignments: 1115249

debug = False
# debug = True

if debug:
    import cProfile

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
    d = {}
    for site in site_str.split(";"):
        s = site.split(",")
        d[(s[0], int(s[1]))] = (s[3], int(s[4]), int(s[2]))  # substitution: from key to value: resn, resi, ncontacts
    return d


def dict_to_site(d):
    l = []
    for key, value in d.iteritems():
        resn, resi = key
        subst_resn, subst_resi, ncontacts = value
        l.append("{},{},{},{},{}".format(resn, resi, ncontacts, subst_resn, subst_resi))
    return ";".join(l)


def load_complex_types(fname):
    complex_type = {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            tpl, pdb, chA, chB, protA, protB, c_type = line.strip().split("\t", 6)
            complex_type[tpl] = c_type
    return complex_type


# def gene_to_protein():
#     d = "/Users/agoncear/data/pdbsws/"
#     uniprot_gene_fname = d + "pir-id-mapping.tab"
#     gene_proteins = defaultdict(set)
#     with open(uniprot_gene_fname) as f:
#         for line in f:
#             fields = line.strip().split()
#             if len(fields) == 2:
#                 uniprot, gene = fields
#                 gene = int(gene)
#                 gene_proteins[gene].add(uniprot)
#     return gene_proteins

# def get_selected_proteins(gene_protein):
#     fname = "/Users/agoncear/data/Vidal/unique_gene_ids.tab"
#     genes = set()
#     proteins = set()
#     with open(fname) as f:
#         for line in f:
#             gene = int(line.strip())
#             genes.add(gene)
#             for protein in gene_protein.get(gene, []):
#                 proteins.add(protein)
#     return proteins, genes

def get_Vidal_proteins():
    fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"
    genes = set()
    proteins = set()
    with open(fname) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) != 2:
                continue
            gene = int(fields[0])
            uniprot = fields[1]
            genes.add(gene)
            proteins.add(uniprot)
    return proteins, genes


def get_PINA_proteins():
    # fname = "/Users/agoncear/data/PINA/Homo sapiens-20140521.sif"
    fname = "/Users/agoncear/data/PINA/unique_uniprot_ids.tab"
    proteins = set()
    with open(fname) as f:
        for line in f:
            uniprot = line.strip()
            if len(uniprot) == 0:
                continue
            proteins.add(uniprot)
    return proteins


ALN = namedtuple('ALN', 'query template identical positive aln_len bs_len bs_covered bs_aligned bs_identical bs_positive bs_contacts bs_BLOSUM bs_score1 site')
def load_alignments(fname, bs_cov_threshold=0.7, sel_proteins=None):
    alignments = defaultdict(list)
    with open(fname) as f:  
        c = 0
        limit = None
        if debug:
            limit = 100000

        for line in islice(f, 1, limit):
            if c % 100000 == 0:
                print c
            c += 1
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
            if aln.identical / float(aln.aln_len) < 0.2:
                continue

            alignments[aln.template].append(aln)
    return alignments


def ca_to_bin(ca):
    b = int(round(ca)) - 2
    if b < 0:
        b = 0
    if b > 15-2:
        b = 15-2
    return b


def load_potential(fname):
    potential = defaultdict(int)
    with open(fname) as f:
        for line in f:
            fields = line.strip().split("\t")
            aa = fields[0]
            scores = fields[1:]
            scores_float = map(float, scores)
            potential[aa] = scores_float
    return potential


def get_atoms(pdb, chainA, chainB):
    pdb_file = get_root() + '/Interfaces/' + pdb[1:3].lower() + '/' + pdb.lower() + '_interface_residue_atoms_5.0A.tab'
    atoms_A = {}
    atoms_B = {}
    with open(pdb_file) as f:
        for line in islice(f, 1, None):
            chain, resn, resi, atomn, x, y, z = line.strip().split()
            resi = int(resi)
            coord = np.array(map(float, (x, y, z)))
            if atomn == "CA":
                if chain == chainA:
                    atoms_A[(resn, resi)] = coord
                if chain == chainB:
                    atoms_B[(resn, resi)] = coord
    return atoms_A, atoms_B


aa_short = "ALPGDNYHECFVIRTKSQMW"
def score(potential, siteA, siteB, atomsA, atomsB, alnA, alnB):
    max_distance = 15.0
    total_score_template_full = .0
    total_scaled_score_template_full = .0
    total_scaled_score_template = .0
    total_scaled_score = .0
    total_score_template = .0
    total_score = .0
    pairs_summed_up = 0

    for k, rA in enumerate(siteA):
        Ca_A = atomsA.get((rA.resn, rA.resi))
        if Ca_A is None:
            # raise Exception("Atoms for residues not found in scoring {} {}".format(rA.resn, rA.resi))
            print "Atoms for residues not found in scoring {} {}".format(rA.resn, rA.resi)
            continue
        # print alnA
        # alignment A:
        resnA_substituted, _, resnA_ncontacts = alnA.site[(rA.resn, rA.resi)][0:3]
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
            resnB_substituted, _, resnB_ncontacts = alnB.site[(rB.resn, rB.resi)][0:3]  # get substituted resn, ncont
            # if resnB_substituted in ('*', '-', '@'):
            if resnB_substituted not in aa_short:
                skip_pair = True

            bin_num = ca_to_bin(d_Ca)
            res_template = rA.resn + rB.resn

            scaling_factor = math.sqrt(float(min(resnA_ncontacts, resnB_ncontacts)))

            # the score of the FULL template (not only the aligned portion)
            s = potential[res_template][bin_num]
            ss = s * scaling_factor
            total_score_template_full += s
            total_scaled_score_template_full += ss

            if skip_pair:
                continue
            # Process the aligned portion of the binding sites:
            # the score of the template
            s = potential[res_template][bin_num]
            ss = s * scaling_factor
            total_score_template += s
            total_scaled_score_template += ss

            # the score of the match
            res = resnA_substituted + resnB_substituted
            # print res, bin_num
            s = potential[res][bin_num]
            ss = s * scaling_factor
            total_score += s
            total_scaled_score += ss

            pairs_summed_up += 1

            # interesting, how the difference of match - template scores will work, will it be normally distributed?

            # print k, rA, Ca_A, "<====>", l, rB, Ca_B, " === ", res_template, d_Ca, bin_num, potential_score_template
            # print res_template, round(d_Ca, 2), potential_score_template

    # print "Summary", total_score_template, total_score
    if pairs_summed_up == 0:
        return
    else:
        return total_score_template_full, total_score_template, total_score, \
            total_scaled_score_template_full, total_scaled_score_template, total_scaled_score


def runner_parallel(template):
    global alignments
    global potential
    return [x for x in runner(alignments, potential, template)]


def runner(alignments, potential, template):
    # pdb, chainA, chainB = template.split("|")
    (pdb, chainA, chainB), (site1, site2) = template

    if chainA != "A" or chainB != "B":
        return

    tplA = pdb.upper() + '|' + chainA
    tplB = pdb.upper() + '|' + chainB
    print "T:", tplA, tplB

    # print atomsA, atomsB
    if len(alignments[tplA]) == 0 or len(alignments[tplB]) == 0:
        return

    atomsA, atomsB = get_atoms(pdb, chainA, chainB)

    # print alignments["4HC1|B"]
    site1_resi = set([s.resi for s in site1])
    site2_resi = set([s.resi for s in site2])

    for alnA in alignments[tplA]:
        # if alnA.query != 'Q9NSI6':
        #     return

        # if it is an alignment for a different site in the same template -> skip it
        alnA_resi = set([a[1] for a in alnA.site.keys()])
        if alnA_resi != site1_resi:
            continue

        for alnB in alignments[tplB]:
            # if alnB.query != 'O75916':
            #     return

            # if it is an alignment for a different site in the same template -> skip it
            alnB_resi = set([a[1] for a in alnB.site.keys()])
            if alnB_resi != site2_resi:
                continue

            scores = score(potential, site1, site2, atomsA, atomsB, alnA, alnB)
            if scores is None:
                continue
            result = '|'.join((pdb, chainA, chainB)), alnA.query, alnB.query, scores, alnA, alnB
            yield result


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def predict():
    global alignments
    global potential

    pdb_templates = get_root() + "/Structures/pdb_templates_5A.tab"
    template_analysis_fname = get_root() + "/Structures/template_analysis.tab"
    alignments_fname = get_root() + "/BLAST-results/human_deltablast_report_merged-bs.tab"
    potential_fname = get_root() + "/Potential/potential_1.index-warp"
    # matches_fname = get_root() + "/Alignments/matches_human_05.tab"
    matches_fname = get_root() + "/Alignments/matches_human_05_20.tab"
    if debug:
        matches_fname = get_root() + "/Alignments/matches_human_05_prof.tab"

    print "Loading potential..."
    potential = load_potential(potential_fname)
    # print potential["EE"][ca_to_bin(6.5)]

    print "Loading templates..."
    complexes = Complexes()
    templates = complexes.loadTemplates(pdb_templates)

    print "Loading complex properties..."
    complex_types = load_complex_types(template_analysis_fname)

    print "Loading selected proteins..."
    # gene_protein = gene_to_protein()
    # selected_proteins, selected_genes = get_selected_proteins(gene_protein)
    # Vidal_proteins, Vidal_genes = get_Vidal_proteins()
    # # print "Total number of selected genes:", len(selected_genes)
    # PINA_proteins = get_PINA_proteins()
    # selected_proteins = PINA_proteins | Vidal_proteins
    # print "Total number of selected proteins:", len(selected_proteins)

    print "Loading alignments..."
    alignments = load_alignments(alignments_fname, bs_cov_threshold=0.5)  # , sel_proteins=selected_proteins)
    print "Total alignments:", sum((len(x) for x in alignments.itervalues()))

    pool = None
    per_template_results = None
    if debug:
        per_template_results = (runner_parallel(template) for template in templates)  # sequential (good for debugging)
    else:
        pool = Pool(8, init_worker)
        per_template_results = pool.imap_unordered(runner_parallel, templates)  # parallel

    c = 0
    print "Writing results to", matches_fname
    with open(matches_fname, 'w') as o:
        header = "query_A query_B template query_type template_type " +\
            "score_template_full score_template score " +\
            "scaled_score_template_full scaled_score_template scaled_score " +\
            "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +\
            "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +\
            "siteA siteB"
        o.write("\t".join(header.split(" ")) + "\n")
        try:
            for results in per_template_results:
                for result in results:
                    if c % 1000 == 0:
                        print c
                    c += 1
                    # if result is None:
                    #     continue
                    # print result
                    tpl, queryA, queryB, scores, alnA, alnB = result
                    # continue

                    # ALN. query template identical positive aln_len bs_len bs_covered bs_aligned bs_identical bs_positive bs_contacts bs_BLOSUM bs_score1 site
                    # s1 = scores[0]  # tempalate compat score
                    # s2 = scores[1]  # model compat score
                    # s1 = scores[2]  # scaled tempalate compat score
                    # s2 = scores[3]  # scaled model compat score

                    # s3 = alnA.bs_BLOSUM + alnB.bs_BLOSUM
                    # s4 = alnA.bs_score1 + alnB.bs_score1

                    # ncont_alnA = sum([site[2] for site in alnA.site.itervalues() if site[1] not in ('-', '*', '@', 'X')])
                    # ncont_alnB = sum([site[2] for site in alnB.site.itervalues() if site[1] not in ('-', '*', '@', 'X')])
                    # s3 = alnA.bs_BLOSUM/float(ncont_alnA) + alnB.bs_BLOSUM/float(ncont_alnB)
                    # s4 = alnA.bs_score1/float(ncont_alnA) + alnB.bs_score1/float(ncont_alnB)

                    # s5 = (alnA.bs_covered / float(alnA.bs_len)) * (alnB.bs_covered / float(alnB.bs_len))
                    # s6 = (alnA.identical / float(alnA.aln_len)) * (alnB.identical / float(alnB.aln_len))
                    if scores is None:
                        continue

                    score_template_full, score_template, score, \
                        scaled_score_template_full, scaled_score_template, scaled_score = scores

                    template_type = complex_types.get(tpl.upper(), "Unknown")
                    query_type = "Hetero"
                    if queryA == queryB:
                        query_type = "Homo"

                    print queryA, queryB  # , tpl, query_type, template_type, s1, s2, s3, s4, s5, s6
                    output = (
                        queryA, queryB, tpl,
                        query_type, template_type,
                        score_template_full, score_template, score,
                        scaled_score_template_full, scaled_score_template, scaled_score,
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


if __name__ == '__main__':
    if debug:
        cProfile.run('predict()')
    else:
        predict()
