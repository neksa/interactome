#

import os
import fnmatch
from itertools import islice
from complexes import Complexes
from collections import defaultdict

def pdb_proteins(pdb_path, fname):
    # chain, m.chain_author, m.uniprot, m.begin, m.end
    pdb_uniprot = defaultdict(list)
    chain_author_to_mmcif = {}

    try:
        with open(fname, 'r') as f: pass
    except:
        try:
            with open(fname, 'w') as o:
                o.write("pdb\tchain\tchain_author\tuniprot\tbegin\tend\n")
                for root, dirnames, filenames in os.walk(pdb_path):
                    for filename in fnmatch.filter(filenames, '*_chain_protein_mapping.tab'):
                        pdb, _ = os.path.basename(filename).lower().split("_", 1)
                        fname_mapping = root + "/" + filename

                        print "Collecting the PDB-Uniprot mapping from", pdb
                        with open(fname_mapping, 'r') as f:
                            for line in islice(f, 1, None):
                                fields = line.strip().split("\t")
                                fields.insert(0, pdb)
                                o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*fields))
        except Exception, e:
            os.remove(fname)
            raise e

    with open(fname, 'r') as f:
        for line in islice(f, 1, None):
            pdb, chain, chain_author, uniprot, begin, end = line.strip().split("\t")
            pdb_chain = pdb.upper() + '|' + chain

            chain_author_to_mmcif[(pdb.upper(), chain_author)] = chain

            if len(uniprot) == 0:
                continue
            uniprot = uniprot.split("-")[0]
            uniprot = uniprot.split(".")[0]
            pdb_uniprot[pdb_chain].append(uniprot)

    return pdb_uniprot, chain_author_to_mmcif


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


def main():
    pdb_templates = "pdb_templates_5A.tab"
    pdb_path = "../Workflow/Interfaces/"

    print "Loading PDB-Uniprot mapping..."
    pdb_uniprot = pdb_proteins(pdb_path)[0]
    
    # print "Loading Uniprot-Entrez Gene ID mapping..."
    # uniprot_gene = protein_genes()

    complexes = Complexes()
    print "Loading Complexes-templates..."
    templates = complexes.loadTemplates(pdb_templates)#, mapping)
    # print "Info: Number of templates (PDB+chain+chain) = ", len(templates.keys())


    with open("template_analysis.tab", 'w') as o:
        o.write("template\tpdb\tA\tB\tprot_A\tprot_B\tcomplex_type\n")
        for (pdb, A, B), _ in templates:
            pdb_chain_A = pdb.upper() + '|' + A.split("_",1)[0]
            pdb_chain_B = pdb.upper() + '|' + B.split("_",1)[0]
            template = pdb.upper() + '|' + A + '|' + B
            # print pdb, A, B, "???"

            uniprots_A = set(pdb_uniprot.get(pdb_chain_A, list()))
            # genes_A = set()
            # for u in uniprots_A:
            #     genes_A |= set(uniprot_gene.get(u, list()))

            uniprots_B = set(pdb_uniprot.get(pdb_chain_B, list()))
            # genes_B = set()
            # for u in uniprots_B:
            #     genes_B |= set(uniprot_gene.get(u, list()))

            complex_type = "Unknown"
            if len(uniprots_A) > 0 and len(uniprots_B) > 0:
                complex_type = "Hetero"
                if len(uniprots_A.intersection(uniprots_B)) > 0:
                    complex_type = "Homo"

            prot_A = ""
            prot_B = ""
            if len(uniprots_A) > 0:
                prot_A = list(uniprots_A)[0]
            if len(uniprots_B) > 0:
                prot_B = list(uniprots_B)[0]

            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(template, pdb, A.split("_",1)[0], B.split("_",1)[0], prot_A, prot_B, complex_type))

# with open("matches.processed.tab", 'r') as f:
#     for line in islice(f, 1, None):
#         q, t, identity, bs_len, bs_err, bs_covered, bs_aligned, bs_identical = line.strip().split("\t")
#         q1 = q[0:4] + q[4:5]
#         q2 = q[0:4]


if __name__ == '__main__':
    main()
