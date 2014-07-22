#!/usr/bin/env python
"""
1. Load all interface definitions from *.int files
2. Extract residue lists for each side of the interface
3. Process all pair of residues for each interface and save the following data:
    - calculate Ca-Ca distance for each type of residue-residue contacts
    - calculate orientation of each Ca frame
    - calculate relative orientation of Ca frame 2 relative to Ca frame 1. Save as quaternion (w, (x,y,z))
    - Calculate the location points of two virtual Cbeta atoms (2.4A from Ca). For Glycine take two points.
    - Calculate the distance between virtual Cbeta. For glycine take the shortest distance out of 4 possible distances
"""

import random
import numpy as np
import fnmatch
import os

from itertools import islice
from collections import defaultdict, namedtuple

from interactome.structures.nr_filter import NRFilter
from interactome.quaternions.algebra import *


aa = ["ALA", "LEU", "PRO", "GLY", "ASP", "ASN", "TYR", "HIS", "GLU", "CYS", "PHE", "VAL", "ILE", "ARG", "THR", "LYS", "SER", "GLN", "MET", "TRP"]
aa_short = "ALPGDNYHECFVIRTKSQMW"
backbone = ("N", "CA", "C", "O", "OXT")
Residue = namedtuple("Residue", "pdb chain resi resn Ca N C O")


# def write_interfaces():
#     NR = NRFilter()

#     with open("interfaces.tab", 'w') as o:

#         for root, dirnames, filenames in os.walk(get_root() + "/Interfaces/"):
#             for filename in fnmatch.filter(filenames, '*_atomic_contacts_5.0A.tab'):
#                 pdb = ""
#                 pdb, _ = os.path.basename(filename).lower().split(".", 1)
#                 fname_int = root + "/" + filename

#                 if not NR.isNR(pdb):
#                     # print "NR skip", pdb
#                     continue

#                 # print pdb

#                 interfaces = defaultdict(set)

#                 # print fname_int
#                 with open(fname_int, 'r') as f:

#                     res = ""
#                     dCA12 = 0.0
#                     for i, line in enumerate(f):
#                         # print "int"
#                         if i == 0:
#                             continue
#                         try:
#                             chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12 = line.strip().split()
#                         except:
#                             print "Error while parsing the line: ", line.strip()
#                             continue

#                         d12 = float(d12)

#                         # if d12 > 4.0: continue   # FILTER FOR 4A heavy-heavy threshold

#                         if not NR.isNR(pdb, chain1):
#                             continue
#                         if resn1 not in aa:
#                             print "skip unknown ", resn1
#                             continue
#                         if resn2 not in aa:
#                             print "skip unknown ", resn2
#                             continue

#                         interfaces[chain1].add((resi1, resn1))
#                         interfaces[chain2].add((resi2, resn2))

#                 for chain, residues in interfaces.iteritems():
#                     for (resi, resn) in sorted(residues):
#                         o.write("{}\t{}\t{}\t{}\n".format(pdb, chain, resn, resi))

# def read_interfaces():
#     interfaces = defaultdict(lambda: defaultdict(set))
#     with open("interfaces.tab", 'r') as f:
#         for line in f:
#             pdb, chain, resn, resi = line.strip().split("\t")
#             # print pdb,chain, resn, resi
#             interfaces[pdb][chain].add((resi, resn))
#     return interfaces

###############################
def residue_iterator():
    NR = NRFilter()

    for root, dirnames, filenames in os.walk(get_root() + "/Interfaces/"):
        for filename in fnmatch.filter(filenames, '*_interface_residue_atoms_5.0A.tab'):
            pdb = ""
            pdb = os.path.basename(filename).lower().split("_", 1)[0]
            # print os.path.basename(filename).lower()
            #.split(".", 1)[0].split("_",1)[0]
            fname_atoms = root + "/" + filename

            if not NR.isNR(pdb):
                # print "NR skip", pdb
                continue
            print "processing PDB", pdb

            residue = None
            with open(fname_atoms, 'r') as f:
                for line in islice(f, 1, None):
                    chain, resn, resi, atomn, x, y, z = line.strip().split()
                    coord = np.array(map(float, (x, y, z)))

                    if residue is not None:
                        if residue["resi"] != resi or residue["chain"] != chain:
                            yield(residue)
                            residue = None

                    if residue is None:
                        residue = {"pdb": pdb, "chain": chain, "resi": resi, "resn": resn}
                    if atomn == "CA":
                        residue["Ca"] = coord
                    elif atomn == "O" or atomn == "OXT":
                        residue["O"] = coord
                    elif atomn == "C":
                        residue["C"] = coord
                    elif atomn == "N":
                        residue["N"] = coord
                    elif atomn == "CB":
                        residue["Cb"] = coord
                    # print residue
                if residue is not None:
                    yield residue


def process_interfaces(fname):
    def out(m):
        return ",".join(["{:.3f}".format(x) for x in m])

    # interfaces = read_interfaces()
    with open(fname, 'w') as o:
        for r in residue_iterator():
            Ca = np.array([.0, .0, .0])
            Cb = np.array([.0, .0, .0])
            vCb = np.array([.0, .0, .0])
            C = np.array([.0, .0, .0])
            # O = np.array([.0, .0, .0])
            quaternion = np.array([.0, .0, .0, .0])
            try:
                Ca = r["Ca"]
                N = r["N"]
                C = r["C"]
                # O = r["O"]
            except:
                print "missing atoms, skipping residue", r["pdb"], r["chain"], r["resn"], r["resi"]
                continue
            # From R. Hansen et al 2011 Proteins, p. 2173 and A. Hanson JMGM 2012, p. 274
            # Calpha frame:
            X = (C-Ca) / np.linalg.norm(C-Ca)
            # R. Hansen:
            # U = np.cross(X, (N-Ca))
            # Z = U / np.linalg.norm(U)
            # A. Hansen:
            U = (N-Ca) / np.linalg.norm(N-Ca)
            Z = np.cross(X, U) / np.linalg.norm(np.cross(X, U))

            Y = np.cross(Z, X)
            F = np.matrix([X, Y, Z])
            quaternion = mat_to_quat(F)
            # print F.trace()
            # print np.linalg.det(F)
            if "Cb" in r:
                Cb = r["Cb"]
                vCb = Ca + 2.4 * ((Cb - Ca)/np.linalg.norm(Cb - Ca))
            # print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
            #     r["pdb"], r["chain"], r["resn"], r["resi"], out(Ca), out(Cb), out(vCb), out(quaternion))
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                r["pdb"], r["chain"], r["resn"], r["resi"], out(Ca), out(Cb), out(vCb), out(quaternion)))


def distributions_iterator(fname_in):
    with open(fname_in, 'r') as f:
        structure_pdb = None
        for line in  f:
            pdb, chain, resn, resi, sCa, sCb, svCb, squaternion = line.strip().split("\t")
            if structure_pdb is not None and structure_pdb != pdb:
                    yield structure_pdb, structure
                    structure_pdb = None
                    structure = None

            if structure_pdb is None:
                structure_pdb = pdb
                structure = defaultdict(list)

            Ca = np.array(map(float, sCa.split(",")))
            Cb = np.array(map(float, sCb.split(",")))
            vCb = np.array(map(float, svCb.split(",")))
            quaternion = np.array(map(float, squaternion.split(",")))

            structure[chain].append((resn, resi, Ca, Cb, vCb, quaternion))
        if structure_pdb is not None:
            yield structure_pdb, structure


def calc_distributions(fname_in, fname_out, shuffled=False):
    distance_threshold_min = 2.0
    distance_threshold_max = 15.0
    empty = np.array([.0, .0, .0])
    # set visibility scope for these vars
    k_resn = None
    l_resn = None

    c1 = 0
    with open(fname_out, 'w') as o:
        for pdb, structure in distributions_iterator(fname_in):            
            c1 += 1
            print c1, pdb,
            # if c % 100:
            #     print ".",
            c2 = 0
            for i, chain1 in enumerate(structure.iterkeys()):

                # TEMP! Normalizing for huge structures with many chains that contribute to bias in statistics
                if i > 0: continue

                for j, chain2 in enumerate(structure.iterkeys()):
                    if j > i:


                        #  # ##### Prepare lists of residues for random choice of residue #####                        
                        # if shuffled:
                        #     k_set = set()
                        #     l_set = set()

                        #     for k, residue1 in enumerate(structure[chain1]):
                        #         for l, residue2 in enumerate(structure[chain2]):
                        #             resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1
                        #             resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                        #             dCa = np.linalg.norm(Ca2 - Ca1)
                        #             if dCa < distance_threshold_min: continue
                        #             if dCa > distance_threshold_max: continue
                        #             k_set.add(k)
                        #             l_set.add(l)

                        #     k_resn = [structure[chain1][k][0] for k in k_set]
                        #     l_resn = [structure[chain2][l][0] for l in l_set]
                        #     if len(k_resn) < 5 or len(l_resn) < 5:
                        #         continue
                        # ##########################

                        residues_chain1 = structure[chain1]
                        residues_chain2 = structure[chain2]
                        if len(residues_chain1) < 5 or len(residues_chain2) < 5:
                            continue

                        residue_names1 = [r[0] for r in residues_chain1]
                        residue_names2 = [r[0] for r in residues_chain2]

                        # if len(residue_names1) < 15:
                        #     print "R1", residue_names1

                        if shuffled:
                            random.shuffle(residue_names1)
                            random.shuffle(residue_names2)

                        # if len(residue_names1) < 15:
                        #     print "R1S", residue_names1

                        for k, residue1 in enumerate(residues_chain1):
                            resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1

                            resn1_shuffl = residue_names1[k]

                            for l, residue2 in enumerate(residues_chain2):
                                resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                                # print Ca2, Ca1
                                resn2_shuffl = residue_names2[l]

                                # check if the residue name has not been changed by shuffling
                                if shuffled and resn1 == resn1_shuffl and resn2 == resn2_shuffl:
                                    continue

                                if shuffled:
                                    resn1 = residue_names1[k]
                                    resn2 = residue_names2[l]
                                res = (resn1, resn2) if aa_short.index(resn1) <= aa_short.index(resn2) else (resn2, resn1)
                                res = "{}_{}".format(*res)

                                dCa = np.linalg.norm(Ca2 - Ca1)
                                if dCa < distance_threshold_min:
                                    print "Too small distance between Ca", pdb, chain1, resi1, "and", chain2, resi2, "d=",dCa
                                    continue
                                if dCa > distance_threshold_max:
                                    # print "Too long distance between Ca", pdb, chain1, chain2, dCa
                                    continue

                                dCb = .0
                                dvCb = .0
                                cos_theta = .0
                                theta = .0
                                if not ((Cb1 == empty).all() or (Cb2 == empty).all()): 
                                    dCb = np.linalg.norm(Cb2 - Cb1)
                                    dvCb = np.linalg.norm(vCb2 - vCb1)

                                    # theta - angle between Ca - Cb vectors of two amino acids
                                    u = (Cb1 - Ca1) / np.linalg.norm(Cb1 - Ca1)
                                    v = (Cb2 - Ca2) / np.linalg.norm(Cb2 - Ca2)
                                    cos_theta = np.dot(u, v)
                                    # print cos_theta, u, v
                                    theta = np.arccos(cos_theta)
                                    if np.isnan(theta):
                                        if (u == v).all(): theta = 0
                                        else: theta = np.pi

                                    dtheta = np.linalg.norm(v-u)

                                # omega - angle between quaternions: orientations of the two amino acids
                                # quaternions are normalized already
                                cos_omega = np.dot(quaternion1, quaternion2)
                                # print cos_omega, quaternion1, quaternion2
                                omega = 2.0 * np.arccos(cos_omega)
                                if np.isnan(omega):
                                    omega = 0 
                                    # if (quaternion1 == quaternion2).all(): omega = 0
                                    # else: omega = np.pi * 2 # could be wrong. need to check what would be the angle between opposite direction quats

                                domega = np.linalg.norm(quaternion2 - quaternion1)

                                # add quaternion calculations here
                                # q_diff = a(-1) * b
                                # rotation from a to b

                                # if res == "ASP_ASP":
                                #     print "======= ASP-ASP ========"
                                #     print pdb, chain 
                                #     print residue1
                                #     print residue2
                                #     print res, round(dCa,2) #, dCb, dvCb, cos_theta, theta
                                # continue

                                # print res, dCa, dCb, dvCb, cos_theta, theta
                                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    res, dCa, dCb, dvCb, cos_theta, theta, dtheta, cos_omega, omega, domega))
                                
                                c2 += 1
            print c2


"""
def contacts_iterator():
    NR = NRFilter()

    for root, dirnames, filenames in os.walk(get_root() + "/Interfaces/"):
        for filename in fnmatch.filter(filenames, '*_atomic_contacts_5.0A.tab'):
            pdb = os.path.basename(filename).lower().split("_", 1)[0]
            fname_atoms = root + "/" + filename

            if not NR.isNR(pdb):
                continue
            print "processing PDB", pdb

            structure = defaultdict(dict())
            residue = None
            with open(fname_atoms, 'r') as f:
                for line in islice(f, 1, None):
                    chain, resn1, resi1, atomn1, chain2, resn2, resi2, atomn2, d12 = line.strip().split()
                    try:
                        d12 = float(d12)
                    except:
                        continue
                    residue = {"resn1": resn1, "resn2": resn2, "d": d12}
            yield structure
            # [chain1|chain2] = [((), ()), ...]


def calc_distributions2(fname_out, fname_shuffl_out):
    empty = np.array([.0, .0, .0])
    # set visibility scope for these vars
    k_resn = None
    l_resn = None

    with open(fname_out, 'w') as o1, open(fname_shuffl_out, 'w') as o2:
        for res in contacts_iterator():
            for i, chain1 in enumerate(structure.iterkeys()):
                for j, chain2 in enumerate(structure.iterkeys()):
                    if j > i:

                        # ##### Prepare lists of residues for random choice of residue #####                        
                        if shuffled:
                            k_set = set()
                            l_set = set()

                            for k, residue1 in enumerate(structure[chain1]):
                                for l, residue2 in enumerate(structure[chain2]):
                                    resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1
                                    resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                                    dCa = np.linalg.norm(Ca2 - Ca1)
                                    if dCa < 2.0: continue
                                    if dCa > 30.0: continue
                                    k_set.add(k)
                                    l_set.add(l)

                            k_resn = [structure[chain1][k][0] for k in k_set]
                            l_resn = [structure[chain2][l][0] for l in l_set]
                        ##########################

                        for k, residue1 in enumerate(structure[chain1]):
                            for l, residue2 in enumerate(structure[chain2]):
                                resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1
                                resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                                # print Ca2, Ca1

                                dCa = np.linalg.norm(Ca2 - Ca1)

                                if dCa < 2:
                                    print "Too small distance between Ca", pdb, chain1, resi1, "and", chain2, resi2, "d=",dCa
                                    continue

                                if dCa > 30.0:
                                    # print "Too long distance between Ca", pdb, chain1, chain2, dCa
                                    continue

                                ##### Randomize amino acid type #####
                                if shuffled:
                                    resn1 = random.choice(k_resn)
                                    resn2 = random.choice(l_resn)
                                #####################################

                                dCb = .0
                                dvCb = .0
                                cos_theta = .0
                                theta = .0
                                if not ((Cb1 == empty).all() or (Cb2 == empty).all()): 
                                    dCb = np.linalg.norm(Cb2 - Cb1)
                                    dvCb = np.linalg.norm(vCb2 - vCb1)

                                    # theta - angle between Ca - Cb vectors of two amino acids
                                    u = (Cb1 - Ca1) / np.linalg.norm(Cb1 - Ca1)
                                    v = (Cb2 - Ca2) / np.linalg.norm(Cb2 - Ca2)
                                    cos_theta = np.dot(u, v)
                                    # print cos_theta, u, v
                                    theta = np.arccos(cos_theta)
                                    if np.isnan(theta):
                                        if (u == v).all(): theta = 0
                                        else: theta = np.pi

                                    dtheta = np.linalg.norm(v-u)

                                # omega - angle between quaternions: orientations of the two amino acids
                                # quaternions are normalized already
                                cos_omega = np.dot(quaternion1, quaternion2)
                                # print cos_omega, quaternion1, quaternion2
                                omega = 2.0 * np.arccos(cos_omega)
                                if np.isnan(omega):
                                    omega = 0 
                                    # if (quaternion1 == quaternion2).all(): omega = 0
                                    # else: omega = np.pi * 2 # could be wrong. need to check what would be the angle between opposite direction quats

                                domega = np.linalg.norm(quaternion2 - quaternion1)

                                # add quaternion calculations here
                                # q_diff = a(-1) * b
                                # rotation from a to b

                                res1, res2 = (resn1, resn2) if aa.index(resn1) <= aa.index(resn2) else (resn2, resn1)
                                res = "{}_{}".format(res1, res2)

                                # if res == "ASP_ASP":
                                #     print "======= ASP-ASP ========"
                                #     print pdb, chain 
                                #     print residue1
                                #     print residue2
                                #     print res, round(dCa,2) #, dCb, dvCb, cos_theta, theta
                                # continue

                                # print res, dCa, dCb, dvCb, cos_theta, theta
                                o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    res, dCa, dCb, dvCb, cos_theta, theta, dtheta, cos_omega, omega, domega))
"""



def get_root():
    return "/Users/agoncear/projects/Interactome/Workflow"


if __name__ == '__main__':
    # write_interfaces()
    # write_coordinates()
    interface_residues = get_root() + "/Potential/nr_interface_residues.tab"
    # process_interfaces(interface_residues)

    # calc_distributions(interface_residues, get_root() + "/Potential/distance_stats.tab")
    calc_distributions(interface_residues, get_root() + "/Potential/distance_stats_shuffled.tab", shuffled=True)
