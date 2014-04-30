"""
Process all interactions
Write a dataset for residue-residue contacts, where A-B and B-A represent the same thing
Calculate the number of contacts beteen heavy atoms given Calpha threshold.

1. number of contacts vs. Calpha distance
2. min, max, avg distance between heavy atoms vs Calpha distance
3. distance between heavy atoms vs Calpha distance

"""
import fnmatch
import itertools
import os
from collections import defaultdict
from non_redundant_filter import NRFilter

aa = ["ALA", "LEU", "PRO", "GLY", "ASP", "ASN", "TYR", "HIS", "GLU", "CYS", "PHE", "VAL", "ILE", "ARG", "THR", "LYS", "SER", "GLN", "MET", "TRP"]
backbone = ("N","CA","C","O","OXT")

def process_interactions():
    NR = NRFilter()

    mode = 'w' # a
    with open("aa-contacts.tab", mode) as o1, open("aa-distances.tab", mode) as o2, open("binary_interactions.tab", mode) as o3:

        for root, dirnames, filenames in os.walk("results/"):
            for filename in fnmatch.filter(filenames, '*.int'):
                pdb = ""
                pdb, _ = os.path.basename(filename).lower().split(".", 1)
                fname_int = root + "/" + filename

                if not NR.isNR(pdb):
                    # print "NR skip", pdb
                    continue

                # print pdb

                # print fname_int
                with open(fname_int, 'r') as f:
                    # demo calc:
                    # A - X  1 .
                    # A - X  2 .
                    # A - Z  1 (AX = 2) write =>
                    # B - X  1 (AZ = 1) =>
                    # B - Y  1 (BX = 1) =>
                    # C - X  1 (BY = 1) =>
                    #          (CX = 1) =>
                    prev_residue = None
                    contact_counter = defaultdict(lambda: [0, 0])
                    interact_sidechain_list = []
                    interact_backbone_list = []
                    interact_mixed_list = []
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

                        if d12 > 4.0: continue

                        if not NR.isNR(pdb, chain1):
                            continue
                        if resn1 not in aa:
                            print "skip unknown ", resn1
                            continue
                        if resn2 not in aa:
                            print "skip unknown ", resn2
                            continue

                        idx1 = aa.index(resn1)
                        idx2 = aa.index(resn2)
                        res = ""
                        if idx1 <= idx2:
                            res = "{}_{}".format(resn1, resn2)
                        else:
                            res = "{}_{}".format(resn2, resn1)

                        interacting_residues = (chain1, resn1, resi1, chain2, resn2, resi2)

                        if atm1 in backbone or atm2 in backbone:
                            # could be a mixed case or pure backbone case

                            if atm1 in backbone and atm2 in backbone:
                                # backbone-backbone interactions only
                                if interacting_residues not in interact_backbone_list:
                                    o3.write("backbone\t{}\t{}\n".format(res, dCA12))
                                    interact_backbone_list.append(interacting_residues)
                            else:
                                # mixed backbone and sidechain interactions only
                                if interacting_residues not in interact_mixed_list:
                                    o3.write("mixed\t{}\t{}\n".format(res, dCA12))
                                    interact_mixed_list.append(interacting_residues)
                        else:
                            # sidechain-sidechain only
                            if interacting_residues not in interact_sidechain_list:
                                o3.write("sidechain\t{}\t{}\n".format(res, dCA12))
                                interact_sidechain_list.append(interacting_residues)

                            # if res == "ASP_ASP":
                            #     # print "ASP", dCA12, line,
                            #     if 10.0 < dCA12 < 11.0:
                            #         print line,


                        # if res == "ASP_LYS":
                        #     atm = None
                        #     if (atm1 == "NZ" and atm2 == "OD2") or (atm1 == "OD2" and atm2 == "NZ"):
                        #         atm = "O_N"
                        #     if (atm1 == "CB" and atm2 == "CB"): #or (atm1 == "OD2" and atm2 == "NZ"):
                        #         atm = "CB_CB"
                        #     if atm is not None:
                        #         o3.write("{}\t{}\t{}\t{}\n".format(res, atm, d12, dCA12))

                        residue = (pdb, chain1, resn1, resi1)
                        # check if the residue is the same and count contacts
                        if prev_residue is not None and prev_residue != residue:
                            # o1.write("{}\t{}\t{}\n".format(res, contact_counter, dCA12))
                            for residue_interaction, (counter, d) in contact_counter.iteritems():
                                if counter > 0:
                                    pair = residue_interaction[0]
                                    o1.write("{}\t{}\t{}\n".format(pair, counter, d))

                            contact_counter = defaultdict(lambda: [0, 0])
                        else:
                            residue_interaction = (res, chain1, resn1, resi1, chain2, resn2, resi2)
                            contact_counter[residue_interaction][0] += 1
                            contact_counter[residue_interaction][1] = dCA12
                        prev_residue = residue[:]

                        # print res, d12, dCA12
                        o2.write("{}\t{}\t{}\n".format(res, d12, dCA12))
                    
                    # repeat writing counter for the last residue, if there is anything there
                    for residue_interaction, (counter, d) in contact_counter.iteritems():
                        if counter > 0:
                            pair = residue_interaction[0]
                            o1.write("{}\t{}\t{}\n".format(pair, counter, d))
                        # o1.write("{}\t{}\t{}\n".format(res, contact_counter, dCA12)) 
                        # pass


if __name__ == '__main__':
    process_interactions()
