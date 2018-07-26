
from complexes import Complexes


NUCLEOSOMES_FILE = "tests/data/nucleosomes/nucleosomes_pdb.csv"
MAPPING = "tests/data/mapping/pdb_chain_uniprot_201807.tsv"


def main():
    d = "/net/pan1/interactomes/pipeline/Interactome/Workflow/Structures"


    # roc_d = d + "/Benchmarks/Ecoli"
    # matches_ecoli = d + "/Alignments/matches_ecoli.tab"
    # pairs = count_scored_pairs(matches_ecoli, 0.25)  # , 3.0)
    # p = set()
    # for pair in pairs:
    #     p.add(pair[0])
    #     p.add(pair[1])
    # print "N proteins = ", len(p)
    # print "N interactions = ", len(pairs)

    pdb_templates = d + "/pdb_templates_5A.tab"

    complexes = Complexes()
    templates = complexes.loadTemplates(pdb_templates, mapping)



if __name__ == '__main__':
    main()
