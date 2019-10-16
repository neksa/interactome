import math
import os
from collections import defaultdict
from scipy.spatial import cKDTree


class Interface:
    def __init__(self, PATH, contact_threshold=None):
        """ Path to results directory and contat threshold in Angstrem between heavy atoms (default: 5A) """
        self.path = PATH
        self.contact_threshold = 5.0
        if contact_threshold is not None:
            self.contact_threshold = contact_threshold

    # def get(self, code):
    #     """ Read precalculated data from files. If files do not exist, run findContacts. The results can be empty """
    #     return atomic_contacts, residues, atoms

    def getOutFilenames(self, code):
        """ Out filenames. Creates directory. Files should be deleted on failure """
        results_dir = "{}/{}".format(self.path, code[1:3])
        if not os.path.isdir(results_dir):
            try:
                os.makedirs(results_dir)
            except:
                pass
        atomic_contacts = "{}/{}_atomic_contacts_{}A.tab".format(results_dir, code, round(self.contact_threshold, 1))
        interface_residue_atoms = "{}/{}_interface_residue_atoms_{}A.tab".format(results_dir, code, round(self.contact_threshold, 1))

        # print "DIRS", atomic_contacts, interface_residue_atoms
        return atomic_contacts, interface_residue_atoms

    def dist(self, a, b):
        """Euclidean distance in 3D"""
        return math.sqrt(sum([(a[i] - b[i])**2 for i in range(3)]))

    def findContacts(self, code, iteratoms):
        """Write interface, write interface atoms"""
        # I. */atomic_contacts_5A.tab
        # II. */interface_residue_atoms_5A.tab
        fname_atomic_contacts, fname_interface_residue_atoms = self.getOutFilenames(code)

        atoms = defaultdict(list)
        coords = defaultdict(list)
        kdtrees = defaultdict(list)
        # chain_suffix = 0

        try:
            for atom in iteratoms:
                atoms[atom.chain].append(atom)
                coords[atom.chain].append(atom.xyz)
                # print(len(atoms[atom.chain]))
        except Exception as e:
            print(e)
            # traceback.print_exc()
            # raise
        print(atoms.keys())
        chains = sorted(atoms.keys())

        for chain in chains:
            kdtrees[chain] = cKDTree(coords[chain])

        n_contacts = 0
        processed_chains = 0
        contacting_residues = set()
        try:
            buf = ""
            for i, c1 in enumerate(chains):
                # print code, "analyze chain", c1, "interactions with:",
                processed_chains += 1
                for j, c2 in enumerate(chains):
                    if j > i:
                        # print c2,
                        results = kdtrees[c1].query_ball_tree(kdtrees[c2], self.contact_threshold)
                        for k, contacts in enumerate(results):
                            for l in contacts:
                                # d12 = self.dist(coords[c1][k], coords[c2][l])
                                d12 = math.sqrt(sum([(coords[c1][k][x] - coords[c2][l][x])**2 for x in range(3)]))
                                a1 = atoms[c1][k]
                                a2 = atoms[c2][l]

                                contacting_residues.add((c1, a1.resi))
                                contacting_residues.add((c2, a2.resi))

                                buf += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n".format(
                                    c1, a1.resn_short, a1.resi, a1.atomn, c2, a2.resn_short, a2.resi, a2.atomn, d12)
                                n_contacts += 1

            if len(buf) > 0:
                with open(fname_atomic_contacts, "w") as o:
                    o.write("chain1\tresn1\tresi1\tatomn1\tchain2\tresn2\tresi2\tatomn2\td12\n")
                    o.write(buf)
                    buf = ""

                with open(fname_interface_residue_atoms, "w") as o:
                    o.write("chain\tresn\tresi\tatomn\tx\ty\tz\n")
                    for chain in chains:
                        for atom in atoms[chain]:
                            if (chain, atom.resi) in contacting_residues:
                                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    chain, atom.resn_short, atom.resi, atom.atomn, atom.xyz[0], atom.xyz[1], atom.xyz[2]))

        except:
            try:
                os.unlink(fname_atomic_contacts)
                os.unlink(fname_interface_residue_atoms)
            except:
                pass
            raise

        print "Processed {}, {} chains, {} contacts".format(code, processed_chains, n_contacts)
        return n_contacts
