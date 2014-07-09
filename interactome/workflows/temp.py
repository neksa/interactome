
from collections import defaultdict
from itertools import islice

from BLOSUM import read_BLOSUM_matrix

import signal
import time
import math

import multiprocessing as mp


def protein_genes():
    d = "/Users/agoncear/data/pdbsws/"
    uniprot_gene_fname = d + "pir-id-mapping.tab"
    # pdb_uniprot_fname = d + "pdb_uniprot_chain_map.lst.2"
    uniprot_gene = defaultdict(list)
    # pdb_uniprot = defaultdict(list)
    with open(uniprot_gene_fname) as f1:#, open(pdb_uniprot_fname) as f2:
        for line in f1:
            # print line
            fields = line.strip().split()
            if len(fields) == 2:
                uniprot, gene = fields
                gene = int(gene)
                uniprot_gene[uniprot].append(gene)
        # for line in f2:
        #     # print line
        #     fields = line.strip().split()
        #     if len(fields) == 3:
        #         pdb, chain, uniprot = fields
        #         if uniprot.endswith("?"): continue
        #         pdb_chain = pdb.upper() + chain.upper()
        #         pdb_uniprot[pdb_chain].append(uniprot)
    return uniprot_gene


def best(in_fname):
    best_protein_template = defaultdict(dict)
    # best_query_template = defaultdict(dict)
    with open(in_fname) as f, open('best_matches_human.tab', 'w') as o:
        for i, line in enumerate(f):
            if i == 0:
                o.write(line)
                continue
            query, template, id_A, q_A, t_B, site_A, id_B, q_B, t_B, site_B = line.strip().split()
            score = (float(id_A) + float(id_B)) / 2.0
            if query not in best_protein_template:
                best_protein_template[query]['score'] = score
                best_protein_template[query]['line'] = line
            else:
                if best_protein_template[query]['score'] < score:
                    best_protein_template[query]['score'] = score
                    best_protein_template[query]['line'] = line
        
        for query, d in best_protein_template.iteritems():
            o.write(d['line'])

def template_analysis():
    complex_types = {}
    with open("template_analysis.tab", 'r') as f:
        for line in islice(f, 1, None):
            # template\tpdb\tA\tB\tprot_A\tprot_B\tcomplex_type
            template, pdb, A, B, prot_A, prot_B, complex_type = line.strip().split("\t")
            complex_types[template] = complex_type
    return complex_types



def calc_properties(line):
    global uniprot_gene
    global blosum
    global complex_types

    # if i == 0: return
    fields = ()
    try:
        fields = line.strip().split()
    except:
        return

    if len(fields) != 11:
        return

    query_A, query_B, template, id_A, q_A, t_B, site_A, id_B, q_B, t_B, site_B = fields

    id_A = int(id_A)
    id_B = int(id_B)

    # if id_A < 25 or id_B < 25: return

    complex_type = 'Unknown'
    try:
        complex_type = complex_types[template]
    except:
        pass

    pdb, _A, _B = template.split("|")
    template_A = pdb + _A
    template_B = pdb + _B

    # print template, complex_type

    # template_A = template[0:4] + template[4:5] #NOTE: error here with multichar chains
    # template_B = template[0:4] + template[5:6]

    ok_A = None
    ok_B = None

    # print id_A
    if True: #int(id_B) != 100:
        # print "A"
        residues = site_A.split(';')
        
        bs_len = len(residues)
        bs_err = 0 # @
        bs_covered = 0 # L,-
        bs_aligned = 0 # L
        bs_identical = 0
        bs_contacts = 0

        bs_score1 = 0.0
        bs_BLOSUM = 0.0
        
        for rec in residues:
            resn1, resi, ncont, resn2 = rec.split(",")

            ncont = int(ncont)
            bs_contacts += ncont

            blos = 0.0
            if resn2 == "*": # not aligned
                pass
                blos = -4.0
            elif resn2 == "@": # error
                bs_err += 1
                blos = -4.0
            elif resn2 == "-": # gap
                bs_covered += 1
                blos = -4.0
            else:
                bs_covered += 1
                bs_aligned += 1
                try:
                    blos = blosum[resn1][resn2]
                except:
                    blos = 0.0

            if resn1 == resn2:
                bs_identical += 1

            bs_BLOSUM += blos
            bs_score1 += ncont * blos

        if bs_covered < 3: return
        ok_A = "A", query_A, query_B, template, complex_type, id_A, bs_len, bs_err, bs_covered, bs_aligned, bs_identical, bs_contacts, bs_BLOSUM, bs_score1
        
        # o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        #     "A", query_A, query_B, template, complex_type, id_A, bs_len, bs_err, bs_covered, bs_aligned, bs_identical, bs_contacts, bs_BLOSUM, bs_score1))
    

    # print id_A
    if True: #int(id_B) != 100:
        # print "A"
        residues = site_B.split(';')
        
        bs_len = len(residues)
        bs_err = 0 # @
        bs_covered = 0 # L,-
        bs_aligned = 0 # L
        bs_identical = 0
        bs_contacts = 0

        bs_score1 = 0.0
        bs_BLOSUM = 0.0
        
        for rec in residues:
            resn1, resi, ncont, resn2 = rec.split(",")

            ncont = int(ncont)
            bs_contacts += ncont

            blos = 0.0
            if resn2 == "*": # not aligned
                pass
                blos = -4.0
            elif resn2 == "@": # error
                bs_err += 1
                blos = -4.0
            elif resn2 == "-": # gap
                bs_covered += 1
                blos = -4.0
            else:
                bs_covered += 1
                bs_aligned += 1
                try:
                    blos = blosum[resn1][resn2]
                except:
                    blos = 0.0

            if resn1 == resn2:
                bs_identical += 1

            bs_BLOSUM += blos
            bs_score1 += ncont * blos
        
        if bs_covered < 3: return
        ok_B = "B", query_A, query_B, template, complex_type, id_B, bs_len, bs_err, bs_covered, bs_aligned, bs_identical, bs_contacts, bs_BLOSUM, bs_score1

        # o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        #     "B", query_A, query_B, template, complex_type, id_B, bs_len, bs_err, bs_covered, bs_aligned, bs_identical, bs_contacts, bs_BLOSUM, bs_score1))

    if ok_A is not None and ok_B is not None:
        return (ok_A, ok_B)



uniprot_gene = None
blosum = None
complex_types = None
def init_worker():
    global uniprot_gene
    global blosum
    global complex_types
    uniprot_gene = protein_genes()
    blosum = read_BLOSUM_matrix("BLOSUM62")
    complex_types = template_analysis()

    signal.signal(signal.SIGINT, signal.SIG_IGN)
    print "Worker initialized."


def gen_lines(fname):
    with open(fname) as f:
        for line in islice(f, 1, None):
            yield line


def main(in_fname, dd):
    n_templates = defaultdict(int)
    template_bslens = {}

    lines = gen_lines(in_fname)
    pool = mp.Pool(8, init_worker)
    gen_results = pool.imap_unordered(calc_properties, lines, 100) # 5000

    c = 0
    with open(dd + "/matches_human.processed.tab-", "w") as o, open(dd + "/ntempl.tab", 'w') as o_ntempl, open(dd + "/bslen.tab", 'w') as o_bslen:
        o_ntempl.write("query\tntempl\n")
        o_bslen.write("tpl\t\complex_type\tbs_len\n")
        o.write("side\tquery_A\tquery_B\ttemplate\tcomplex_type\tidentity\tbs_len\tbs_err\tbs_covered\tbs_aligned\tbs_identical\tbs_contacts\tbs_BLOSUM\tbs_score1\n")

        try:
            for result in gen_results:

                c += 1
                if c % 10000 == 0: print c

                if result is None: continue
                A, B = result
                # print A
                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*A))
                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*B))
                template_bslens[A[3]] = (A[4], A[6])
                template_bslens[B[3]] = (B[4], B[6])

                # number of templates per query pair
                query_A = A[1]            
                query_B = A[2]            
                if query_A < query_B:
                    n_templates[query_A + '+' + query_B] += 1
                else:
                    n_templates[query_B + '+' + query_A] += 1

            for key, val in n_templates.iteritems():
               o_ntempl.write("{}\t{}\n".format(key, val))

            for key, val in template_bslens.iteritems():
                # for complex_type, bslen in val.iteritems():
                complex_type, length = val
                o_bslen.write("{}\t{}\t{}\n".format(key, complex_type, length))
        # except KeyboardInterrupt:
        except:
            print "Caught KeyboardInterrupt or other exception, terminating workers"
            pool.terminate()
            pool.join()
            raise

    pool.close()
    pool.join()



if __name__ == '__main__':
    # best("matches_human.tab")
    # main("matches_human.tab", 'all')
    main("matches_human_25.tab", '25')
    main("matches_human_20-25.tab", '20-25')
    main("matches_human_15-20.tab", '15-20')
    # main("best_matches_human.tab", 'bestp')

    