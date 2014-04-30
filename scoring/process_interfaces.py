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

import math
import random
import numpy as np
import fnmatch
import itertools
import os
import string
from collections import defaultdict, namedtuple
from non_redundant_filter import NRFilter


aa = ["ALA", "LEU", "PRO", "GLY", "ASP", "ASN", "TYR", "HIS", "GLU", "CYS", "PHE", "VAL", "ILE", "ARG", "THR", "LYS", "SER", "GLN", "MET", "TRP"]
backbone = ("N","CA","C","O","OXT")
Residue = namedtuple("Residue", "pdb chain resi resn Ca N C O")


def write_interfaces():
    NR = NRFilter()

    with open("interfaces.tab", 'w') as o:

        for root, dirnames, filenames in os.walk("results/"):
            for filename in fnmatch.filter(filenames, '*.int'):
                pdb = ""
                pdb, _ = os.path.basename(filename).lower().split(".", 1)
                fname_int = root + "/" + filename

                if not NR.isNR(pdb):
                    # print "NR skip", pdb
                    continue

                # print pdb

                interfaces = defaultdict(set)

                # print fname_int
                with open(fname_int, 'r') as f:

                    res = ""
                    dCA12 = 0.0
                    for i, line in enumerate(f):
                        # print "int"
                        if i == 0:
                            continue
                        try:
                            pdb, chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12, dCA12 = line.strip().split()
                        except:
                            print "Error while parsing the line: ", line.strip()
                            continue

                        d12 = float(d12)
                        dCA12 = float (dCA12)

                        # if d12 > 4.0: continue   # FILTER FOR 4A heavy-heavy threshold

                        if not NR.isNR(pdb, chain1):
                            continue
                        if resn1 not in aa:
                            print "skip unknown ", resn1
                            continue
                        if resn2 not in aa:
                            print "skip unknown ", resn2
                            continue

                        interfaces[chain1].add((resi1, resn1))
                        interfaces[chain2].add((resi2, resn2))

                for chain, residues in interfaces.iteritems():
                    for (resi, resn) in sorted(residues):
                        o.write("{}\t{}\t{}\t{}\n".format(pdb, chain, resn, resi))


def read_interfaces():
    interfaces = defaultdict(lambda: defaultdict(set))
    with open("interfaces.tab", 'r') as f:
        for line in f:
            pdb, chain, resn, resi = line.strip().split("\t")
            # print pdb,chain, resn, resi
            interfaces[pdb][chain].add((resi, resn))
    return interfaces


def write_coordinates():

    interfaces = read_interfaces()

    with open("interface_residues.tab", 'w') as o:
        for pdb, chains in interfaces.iteritems():
            with open("results/{}/{}.pdb".format(pdb[1:3], pdb), 'r') as f:
                chain_suffix = 0
                
                for line in f:
                    if line.startswith("MODEL"):
                        model = int(line.strip().split()[1])
                        # print "Model {}".format(model)
                        chain_suffix = model - 1

                    if line.startswith("ATOM"):
                        pdb_atomn = line[13:16].strip() # 'CA' # questions about strip!!! Calcium short vs CA
                        pdb_atomi = int(line[4:12].strip())
                        pdb_resn = line[17:20].strip()
                        pdb_chain = line[21:22] 
                        pdb_resi = line[22:26].strip()
                        pdb_element = line[76:78].strip()

                        # if pdb_element == "H": continue # skip all hydrogen atoms
                        # if pdb_resn == "UNK": continue  # skip unknown residues
                        if line[16] in string.letters: continue # skip atoms with letters in atom numbers
                        if line[26] in string.letters: continue # skip atoms with letters in residue numbers

                        chain = pdb_chain
                        if chain_suffix > 0:
                            chain = "{}_{}".format(pdb_chain, chain_suffix)

                        # print pdb, chain, pdb_resi, pdb_resn, chains[chain]

                        if chain not in chains:
                            # print "sk", pdb, chain
                            continue

                        if (pdb_resi, pdb_resn) in chains[chain]:
                            # print "MATCH", pdb, chain, pdb_resi, pdb_resn
                            # print "sk", pdb, pdb_resi, pdb_resn
                            x, y, z = [float(line[30+8*i:38+8*i]) for i in range(3)]
                            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pdb, chain, pdb_resi, pdb_resn, pdb_atomn, x, y, z))


                        # if pdb_atomn == "CA":
                        #     residues_with_CA.add((chain, pdb_resi))

                        # # atom = (pdb, chain, pdb_resi, pdb_resn, pdb_atomi, pdb_atomn)
                        # atom = (pdb, chain, pdb_resi, pdb_resn, pdb_atomn) # Do not change the order of attributes in the tuple
                        # xyz = [float(line[30+8*i:38+8*i]) for i in range(3)]
                        # # print atom
                        # # print xyz
                        # # [[l[12:26],[float(l[30+8*i:38+8*i]) for i in range(3)]]
                        # chains[chain]["atom"].append(atom)
                        # chains[chain]["xyz"].append(xyz)


def residue_iterator():
    residue = None
    with open("interface_residues.tab", 'r') as f:
        for line in f:
            pdb, chain, resi, resn, atomn, x, y, z = line.strip().split()
            x, y, z = map(float, (x, y, z))

            if residue is not None:
                if residue["resi"] != resi or residue["chain"] != chain:
                    yield(residue)
                    residue = None
            
            if residue is None:
                residue = {"pdb":pdb, "chain": chain, "resi": resi, "resn": resn}
                
            if atomn == "CA":
                residue["Ca"] = np.array([x, y, z])
            elif atomn == "O" or atomn == "OXT":
                residue["O"] = np.array([x, y, z])
            elif atomn == "C":
                residue["C"] = np.array([x, y, z])
            elif atomn == "N":
                residue["N"] = np.array([x, y, z])
            elif atomn == "CB":
                residue["Cb"] =np.array([x, y, z])
            # print residue

        if residue is not None:
            yield residue


def mat_to_quat(m):
    """
    Based on http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche52.html
    """
    def sign(x):
        if x >= 0.0: return +1.0
        else: return -1.0

    r11 = m[0,0]
    r22 = m[1,1]
    r33 = m[2,2]
    r12 = m[0,1]
    r21 = m[1,0]
    r13 = m[0,2]
    r31 = m[2,0]
    r32 = m[2,1]
    r23 = m[1,2]

    q0 = ( r11 + r22 + r33 + 1.0) / 4.0
    q1 = ( r11 - r22 - r33 + 1.0) / 4.0
    q2 = (-r11 + r22 - r33 + 1.0) / 4.0
    q3 = (-r11 - r22 + r33 + 1.0) / 4.0
    if q0 < 0.0: q0 = 0.0
    if q1 < 0.0: q1 = 0.0
    if q2 < 0.0: q2 = 0.0
    if q3 < 0.0: q3 = 0.0
    q0 = math.sqrt(q0)
    q1 = math.sqrt(q1)
    q2 = math.sqrt(q2)
    q3 = math.sqrt(q3)
    if q0 >= q1 and q0 >= q2 and q0 >= q3:
        q0 *= +1.0
        q1 *= sign(r32 - r23)
        q2 *= sign(r13 - r31)
        q3 *= sign(r21 - r12)
    elif q1 >= q0 and q1 >= q2 and q1 >= q3:
        q0 *= sign(r32 - r23)
        q1 *= +1.0
        q2 *= sign(r21 + r12)
        q3 *= sign(r13 + r31)
    elif q2 >= q0 and q2 >= q1 and q2 >= q3:
        q0 *= sign(r13 - r31)
        q1 *= sign(r21 + r12)
        q2 *= +1.0
        q3 *= sign(r32 + r23)
    elif q3 >= q0 and q3 >= q1 and q3 >= q2:
        q0 *= sign(r21 - r12)
        q1 *= sign(r31 + r13)
        q2 *= sign(r32 + r23)
        q3 *= +1.0
    else:
        print "Coding error"

    q = np.array([q0, q1, q2, q3])
    q /= np.linalg.norm(q)
    return q


def process_interfaces():
    def out(m):
        return ",".join(["{:.3f}".format(x) for x in m])

    interfaces = read_interfaces()
    
    with open("interface_residues_and_more.tab", 'w') as o:
        for r in residue_iterator():
            Ca = np.array([.0, .0, .0])
            Cb = np.array([.0, .0, .0])
            vCb = np.array([.0, .0, .0])
            quaternion = np.array([.0, .0, .0, .0])

            try:
                Ca = r["Ca"]
                N = r["N"]
                C = r["C"]
                O = r["O"]
            except:
                print "missing atoms, skipping residue", r["pdb"], r["chain"], r["resn"], r["resi"]
                continue

            X = (C-Ca) / np.linalg.norm(C-Ca)
            U = (N-Ca) / np.linalg.norm(N-Ca)
            Z = np.cross(X, U) / np.linalg.norm(np.cross(X, U))
            Y = np.cross(Z, X)

            F = np.matrix([X, Y, Z])
            # print F

            # Q = F.trace()
            # print Q
            # print np.linalg.det(F)
            # print F.trace()
            # if F.trace()
            quaternion = mat_to_quat(F)

            if "Cb" in r:
                Cb = r["Cb"]
                vCb = Ca + 2.4 * ((Cb - Ca)/np.linalg.norm(Cb - Ca))

            # print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
            #     r["pdb"], r["chain"], r["resn"], r["resi"], out(Ca), out(Cb), out(vCb), out(quaternion))
            
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                r["pdb"], r["chain"], r["resn"], r["resi"], out(Ca), out(Cb), out(vCb), out(quaternion)))


    # for pdb in
    # Cb vCb quaternion 

def distributions_iterator():
    with open("interface_residues_and_more.tab", 'r') as f:
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

            

def calc_distributions(shuffled = False):
    
    empty = np.array([.0, .0, .0])

    fname = "distance_stats_2.tab"
    if shuffled:
        fname = "distance_stats_2_shuffled.tab"

    # set visibility scope for these vars
    k_resn = None
    l_resn = None

    with open(fname, 'w') as o:
        for pdb, structure in distributions_iterator():
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


if __name__ == '__main__':
    pass
    # write_interfaces()
    # write_coordinates()
    # process_interfaces()
    # calc_distributions()

    # calc_distributions()
    # calc_distributions(shuffled = True)
