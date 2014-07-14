"""
filter chains according to PDB NR-blast dataset

Write a dataset for residue-residue contacts, where A-B and B-A represent the same thing
Calculate the number of contacts beteen heavy atoms given Calpha threshold.

"""

# DIR = "/Users/agoncear/data/PDB/clusters/"
# DIR_scoring = "/Users/agoncear/projects/Interactome/scoring/"

from interactome.structures.mapping import pdb_proteins


class NRFilter():
    def __init__(self):
        NR_file = "/Users/agoncear/projects/Interactome/Workflow/Potential/NR_pdb_chains_40.tab"
        mapping_file = "/Users/agoncear/projects/Interactome/Workflow/Structures/pdb_proteins.tab"
        pdb_root = "/Users/agoncear/projects/Interactome/Workflow/Interfaces/"

        pdb_uniprot, author_to_mmcif = pdb_proteins(pdb_root, mapping_file)

        self.nr_set = set()
        self.nr_pdb_set = set()
        with open(NR_file, 'r') as f:
            for line in f:
                pdb, chain_author = line.strip().split()
                mmcif_chain = author_to_mmcif.get((pdb.upper(), chain_author), None)
                if mmcif_chain is None:
                    continue
                pdb = pdb.lower()
                self.nr_set.add((pdb, mmcif_chain))
                self.nr_pdb_set.add(pdb)

    def isNR(self, pdb, chain=None):
        if chain is None:
            return pdb in self.nr_pdb_set
        return (pdb, chain) in self.nr_set


# def nr_filter(identity):
#     fname_in = "{}bc-{}.out".format(DIR, identity)
#     fname_out = "{}NR_pdb_chains_{}.tab".format(DIR_scoring, identity)
#     # cluster_id = 0
#     # clusters = {}
#     # pdb_codes = []
#     with open(fname_in, 'r') as f, open(fname_out, 'w') as o:
#         for line in f:
#             pdbs = line.strip().split()
#             for i, pdb_chain in enumerate(pdbs):
#                 pdb, chain = pdb_chain.split("_")
#                 if i == 0:
#                     o.write("{}\t{}\n".format(pdb, chain))
#                 # clusters[pdb_chain] = cluster_id
#                 # if pdb not in pdb_codes: pdb_codes.append(pdb)
#             # cluster_id += 1


if __name__ == '__main__':
    pass
    # nr_filter("40")
