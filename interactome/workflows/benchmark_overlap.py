# from itertools import islice
# from collections import defaultdict, namedtuple

from benchmark import load_interactions, load_gene_protein_mapping
from filter_interactome import get_PINA_proteins, load_uniprot_mapping


def main():
    fname_PINA = "/Users/agoncear/data/PINA/Homo sapiens-20140521.sif"
    network = "/Users/agoncear/projects/Interactome/Workflow/Interactomes/stringent.tab"
    vidal_fname = "/Users/agoncear/data/Vidal/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"

    # uniprot_pdbs, mmcif_pdb_chain = pdb_mapping(pdb_proteins_fname)
    pina = get_PINA_proteins(fname_PINA)

    gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    records, old_new_ids = load_uniprot_mapping()
    vidal_gene_gene = load_interactions(vidal_fname)
    vidal = set()
    for gA, gB in vidal_gene_gene:
        pA = gene_proteins[gA]
        pB = gene_proteins[gB]
        for p in pA:
            for q in pB:
                pair = p, q
                if p > q:
                    pair = q, p
                vidal.add(pair)

    our = set()
    with open(network) as f:
        for line in f:
            a, b, _ = line.strip().split("\t", 2)
            pair = a, b
            if a > b:
                pair = b, a
            our.add(pair)

    # print list(pina)[0]
    # print list(vidal)[0]
    # print list(our)[0]

    print "Total PINA:", len(pina)
    print "Total Vidal:", len(vidal)
    print "Total Our:", len(our)

    print "PINA & Vidal & Our:", len(pina & vidal & our)

    print "PINA & Vidal -Our:", len(pina & vidal - our)
    print "Our & Vidal -PINA:", len(our & vidal - pina)
    print "Our & PINA -Vidal:", len(our & pina - vidal)

    print "only in PINA:", len(pina - vidal - our)
    print "only in Vidal:", len(vidal - pina - our)
    print "only in Our:", len(our - pina - vidal)


if __name__ == '__main__':
    main()
