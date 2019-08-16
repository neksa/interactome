
from structures.complexes import Complexes
from structures.SIFTS import SIFTS, PDBChain


PDB_MAPPING = "tests/data/mapping/pdb_chain_uniprot_201807.tsv"

NUCLEOSOMES_FILE = "tests/data/nucleosomes/nucleosomes_pdb.csv"
UNIPROT_MAPPING = "tests/data/nucleosomes/uniprot.tab"


def read_chain_mapping(fname):
    chain_mapping = {}
    with open(fname) as f:
        for line in f:
            pdb, chain, chain_author = line.strip().split()[:3]
            chain_mapping[(pdb, chain)] = chain_author
    return chain_mapping


def read_uniprot_mapping(fname):
    mapping = {}
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            _, uniprot, uniprot_name, _, protein_name, gene, organism = line.strip().split("\t")[:7]
            mapping[uniprot] = protein_name
    return mapping


def main():
    d = "/net/pan1/interactomes/pipeline/"

    # roc_d = d + "/Benchmarks/Ecoli"
    # matches_ecoli = d + "/Alignments/matches_ecoli.tab"
    # pairs = count_scored_pairs(matches_ecoli, 0.25)  # , 3.0)
    # p = set()
    # for pair in pairs:
    #     p.add(pair[0])
    #     p.add(pair[1])
    # print "N proteins = ", len(p)
    # print "N interactions = ", len(pairs)

    pdb_templates = d + "Interactome/Workflow/Structures/pdb_templates_5A.tab"
    pdb_proteins = d + "Interactome/Workflow/Structures/pdb_proteins.tab"
    pdb_mapping_fname = d + "interactome/" + PDB_MAPPING
    uniprot_mapping_fname = d + "interactome/" + UNIPROT_MAPPING
    nucleosomes = d + "interactome/" + NUCLEOSOMES_FILE

    selected_pdbs = set()
    with open(nucleosomes) as f:
        for line in f:
            selected_pdbs.add(line.strip().split()[0].lower())
    # print(len(selected_pdbs))

    chain_mapping = read_chain_mapping(pdb_proteins)
    uniprot_mapping = read_uniprot_mapping(uniprot_mapping_fname)

    complexes = Complexes()
    templates = complexes.loadTemplates(pdb_templates) # , mapping)

    sifts = SIFTS()
    sifts_mapping = sifts.getPDBMapping(pdb_mapping_fname)

    for template in templates:
        (pdb, chainA, chainB), (site1, site2) = template
        if "_" in chainA:
            chainA = chainA.split("_")[0]
        if "_" in chainB:
            chainB = chainB.split("_")[0]
        # print(pdb)
        if pdb in selected_pdbs:
            chainA_author = chain_mapping[(pdb, chainA)]
            chainA_uniprot = sifts_mapping.get((pdb, chainA_author), pdb + "_" + chainA_author)
            if type(chainA_uniprot) is PDBChain:
                chainA_uniprot = chainA_uniprot.protein
                # chainA_uniprot = uniprot_mapping.get(chainA_uniprot, chainA_uniprot)

            chainB_author = chain_mapping[(pdb, chainB)]
            chainB_uniprot = sifts_mapping.get((pdb, chainB_author), pdb + "_" + chainB_author)
            if type(chainB_uniprot) is PDBChain:
                chainB_uniprot = chainB_uniprot.protein
                # chainB_uniprot = uniprot_mapping.get(chainB_uniprot, chainB_uniprot)

            # print(pdb, chainA_uniprot, chainB_uniprot)
            # print(chainA_uniprot)
            # print(chainB_uniprot)
            # print(template)
            print("{} pp {}". format(chainA_uniprot, chainB_uniprot))


if __name__ == '__main__':
    main()
