"""
Analyze complexes in PDB structures

Assembly #1

Interface is defined as contact between heavy atoms in two different chains < 5A ???

d_HA = heavy atom (our cutoff parameter). Contact is determined according to 5A threshold
d_CA = R1_CA - R2_CA distance

I would store all HA contacts within 5A.

Store:
pdb
chain 1 + model
chain 2 + model
resi 1 atom A
resi 2 atom B
res name 1 atom A
res name 2 atom B
distance A <-> B
distance res1_CA <-> res2_CA


PROBLEMS with certain PDB files.

1. result of residue names with letters. Ca belongs to the last mention of the residue
Example:

Alexanders-MacBook-Pro:scoring agoncear$ find results -name "*.int" |xargs grep "118.41"
results/hq/1hqm.int:1hqm    C   VAL 1099    CG1 D   VAL    8    CG1 3.57    118.41


"""
import fnmatch
import os
import gzip
import math
import string
import signal
import time

from collections import defaultdict
from scipy.spatial import cKDTree
import multiprocessing as mp

from non_redundant_filter import NRFilter


class Contacts:

    def __init__(self, NR=None, threshold=5.0):
        self.threshold_distance = threshold
        self.NR = NR

    def findAll(self, pdb, in_fname, out_fname):
        # chain_numbers = []
        # residue_numbers = {}
        # atom_ids = []
        # atom_xyz = []

        chains = defaultdict(lambda: defaultdict(list))
        residues_with_CA = set()
        chain_suffix = 0

        with open(in_fname, 'r') as f:
            for line in f:
                if line.startswith("MODEL"):
                    model = int(line.strip().split()[1])
                    # print "Model {}".format(model)
                    chain_suffix = model - 1

                if line.startswith("ATOM"):
                    pdb_atomn = line[13:16].strip()  # 'CA' # questions about strip!!! Calcium short vs CA
                    # pdb_atomi = int(line[4:12].strip())
                    pdb_resn = line[17:20].strip()
                    pdb_chain = line[21:22]
                    pdb_resi = line[22:26]
                    pdb_element = line[76:78].strip()

                    if pdb_element == "H":
                        continue  # skip all hydrogen atoms
                    if pdb_resn == "UNK":
                        continue  # skip unknown residues
                    if line[16] in string.letters:
                        continue  # skip atoms with letters in atom numbers
                    if line[26] in string.letters:
                        continue  # skip atoms with letters in residue numbers

                    chain = pdb_chain
                    if chain_suffix > 0:
                        chain = "{}_{}".format(pdb_chain, chain_suffix)

                    if pdb_atomn == "CA":
                        residues_with_CA.add((chain, pdb_resi))

                    # atom = (pdb, chain, pdb_resi, pdb_resn, pdb_atomi, pdb_atomn)
                    atom = (pdb, chain, pdb_resi, pdb_resn, pdb_atomn)  # Do not change the order of attributes in the tuple
                    xyz = [float(line[30+8*i:38+8*i]) for i in range(3)]
                    # print atom
                    # print xyz
                    # [[l[12:26],[float(l[30+8*i:38+8*i]) for i in range(3)]]
                    chains[chain]["atom"].append(atom)
                    chains[chain]["xyz"].append(xyz)

        for chain in sorted(chains.keys()):
            atoms = chains[chain]["atom"]  # ref, not a copy
            coordinates = chains[chain]["xyz"]  # ref
            atoms_to_remove = [i for i, atom in enumerate(atoms) if (chain, atom[2]) not in residues_with_CA]
            atoms[:] = [atom for i, atom in enumerate(atoms) if i not in atoms_to_remove]
            coordinates[:] = [xyz for i, xyz in enumerate(coordinates) if i not in atoms_to_remove]

            if len(atoms) == 0 and len(coordinates) == 0:
                # print "Removing the whole chain {} because there are no CA atoms. Possibly a nucleic acid chain".format(chain)
                del chains[chain]

        for chain in chains.iterkeys():
            chains[chain]["CA_atom_indices"] = {}
            for i, atom in enumerate(chains[chain]["atom"]):
                if atom[4] == "CA":
                    chains[chain]["CA_atom_indices"][atom[:-2]] = i

        for chain in chains.iterkeys():
            chains[chain]["tree"] = cKDTree(chains[chain]["xyz"])

        n_contacts = 0
        processed_chains = 0
        with open(out_fname, "w") as o:
            o.write("pdb\tchain1\tresn1\tresi1\tatm1\tchain2\tresn2\tresi2\tatm2\td12\tdCA12\n")

            for i, chain1 in enumerate(chains.iterkeys()):

                if self.NR is None:
                    # print "NR is NONE!"
                    pass
                if self.NR is not None and not self.NR.isNR(pdb, chain1):
                    # print pdb, "skip chain", chain1
                    continue
                # print pdb, "analyze chain", chain1,

                processed_chains += 1
                for j, chain2 in enumerate(chains.iterkeys()):
                    if j > i:
                        # print chain2,
                        atom1 = chains[chain1]["atom"]
                        atom2 = chains[chain2]["atom"]
                        xyz1 = chains[chain1]["xyz"]
                        xyz2 = chains[chain2]["xyz"]
                        CA1 = chains[chain1]["CA_atom_indices"]
                        CA2 = chains[chain2]["CA_atom_indices"]
                        tree1 = chains[chain1]["tree"]
                        tree2 = chains[chain2]["tree"]
                        results = tree1.query_ball_tree(tree2, self.threshold_distance)
                        # print results
                        for k, contacts in enumerate(results):
                            for l in contacts:
                                d12 = self.dist(xyz1[k], xyz2[l])
                                CA_k = CA1[atom1[k][:-2]]
                                CA_l = CA2[atom2[l][:-2]]
                                dCA12 = self.dist(xyz1[CA_k], xyz2[CA_l])
                                # print "{}: {} <-> {} {} <-> {} [d={:.2f} A] [dCA={:.2f} A] {} <-> {}".format(
                                #     pdb, chain1, chain2, k, l, d12, dCA12, str(atom1[k]), str(atom2[l]))
                                resi1 = atom1[k][2]
                                resi2 = atom2[l][2]
                                resn1 = atom1[k][3]
                                resn2 = atom2[l][3]
                                atm1 = atom1[k][4]
                                atm2 = atom2[l][4]

                                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\n".format(
                                    pdb, chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12, dCA12))
                                n_contacts += 1

        print "Processed {}, {} chains, {} contacts".format(pdb, processed_chains, n_contacts)

    def dist(self, a, b):
        return math.sqrt(sum([(a[i] - b[i])**2 for i in range(3)]))


def extract_contacts(params):
    pdb_code, gz_fname, fname_pdb, fname_int = params
    # print pdb_code
    try:
        with open(fname_pdb):
            pass
    except:
        with gzip.open(gz_fname, 'r') as in_gz:
            with open(fname_pdb, 'w') as out_pdb:
                raw_pdb = in_gz.read()
                # pdb = PDB(raw_pdb)
                # clean_pdb = pdb.cleanUp()
                out_pdb.write(raw_pdb)
                # print "IN", gz_fname, "OUT", fname_pdb
    NR = NRFilter()
    c = Contacts(NR)
    c.findAll(pdb_code, fname_pdb, fname_int)


# class Structure:

#     def __init__(self):
#         root = "/Users/agoncear/projects/Interactome/scoring/"
#         self.results_dir = root+"results/"
#         self.structures_dir = "/Users/agoncear/data/PDB/biounit/"  # root+"PDB/"

#     def dirAll(self, NR=None, only=None):
#         """
#             if results file exists - skip it
#         """
#         for root, dirnames, filenames in os.walk(self.structures_dir):
#             for filename in fnmatch.filter(filenames, '*.pdb1.gz'):
#                 pdb_code = ""
#                 asm = 0
#                 # print filename
#                 if filename.endswith("pdb1.gz"):
#                     pdb_code, pdb_asm, _ = os.path.basename(filename).lower().split(".", 2)
#                     asm = int(pdb_asm[3:])
#                 else:
#                     continue

#                 if only is not None and pdb_code[1:3] != only:
#                     continue

#                 if NR is not None and not NR.isNR(pdb_code):
#                     # print "Skip {}".format(pdb_code)
#                     continue

#                 gz_fname = root + "/" + filename
#                 fname_base = self.results_dir + pdb_code[1:3] + "/" + pdb_code
#                 fname_pdb = fname_base + ".pdb"
#                 fname_int = fname_base + ".int"

#                 new_dir = self.results_dir + pdb_code[1:3]
#                 try:
#                     os.mkdir(new_dir)
#                 except:
#                     pass

#                 try:
#                     with open(fname_int):
#                         raise Exception()  # TEMP!!!
#                         # pass
#                 except:
#                     yield((pdb_code, gz_fname, fname_pdb, fname_int))

#     def runAll(self, only=None):
#         for (pdb_code, gz_fname, fname_pdb, fname_int) in self.dirAll(only):
#             # print pdb_code, gz_fname, fname_pdb, fname_int
#             extract_pdb_gz_contacts(pdb_code, gz_fname, fname_pdb, fname_int)


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


if __name__ == '__main__':
    struct = Structure()

    pool = mp.Pool(8, init_worker)
    NR = NRFilter()
    pool.imap_unordered(extract_contacts, struct.dirAll(NR))
    # pool.imap_unordered(extract_pdb_gz_contacts, it)

    try:
        while(True):
            print "Watchdog... every 60 seconds (Ctrl-C to interrupt)"
            time.sleep(60)

    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()


# def process_pdb():

#     with open("output.tab", 'w') as f:
#         f.write("pdb\tchain1\tchain2\tresi1\tresi2\tresn1\tresn2\tatm1\tatm2\td12\tdCA12\n")
#         for root, dirnames, filenames in os.walk("PDB"):
#             for filename in fnmatch.filter(filenames, '*.pdb'):
#                 pdb, asm, _ = os.path.basename(filename).lower().split(".", 2)
#                 subdir = pdb[1:3]
#                 f.write("%s\t%s\t%s\n" % (pdb, subdir, asm))

#                 path = OUTPUT_PATH + subdir + "/"
#                 try:
#                     os.makedirs(path)
#                 except:
#                     pass


# if __name__ == '__main__':
#     process_pdb()
