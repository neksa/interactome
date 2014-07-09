"""
Scoring the protein structures (complexes) with a PPI potential

0. Load the structure
1. Load the potential
2. Identify the interface in the structure: 5A heavy atom distance
3. Record all the amino acids (Calphas) in the interface, their positions and residue types
4. Summ up over all amino acid pairs and compute the LOG-Likelihood score

"""

import math
import random
import signal
import numpy as np
from collections import defaultdict
from process_interfaces import distributions_iterator, aa

def bin47(ca):
    """
    Starts with bin 0, which corresponds to 2.0-2.5
    1: 2.5 - 3.0
    2: 3.0 - 3.5
    ... so on ...
    46: 25.0+
    """
    if ca < 2.0: ca = 2.0
    if ca > 25.0: ca = 25.0
    ca -= 2.0
    bin = int(math.floor(ca/0.5))
    return bin

def bin20(ca):
    if ca < 2.0: ca = 2.0
    if ca > 22.0: ca = 22.0
    ca -= 2.0
    bin = int(math.floor(ca/1.0))
    return bin

def get_bin_table(binfunc):
    table = []
    for i in range(31*2):
        ca = i / 2.0 # Half of Angstrem - to bin - table
        bin = binfunc(ca)
        # print i, ca, bin
        table.append(bin)
    return table


def get_hist(fname, bintable):
    # nbins = binfunc(100.0) + 1
    nbins = bintable[int(math.floor(30.0/0.5))] + 1
    hist = defaultdict(lambda: [1]*nbins) # Laplace smoothing: + 1
    with open(fname) as f:
        for line in f:
            aa, ca, _ = line.split("\t", 2)
            ca = float(ca)
            # bin = binfunc(ca)
            # each bin in the table is half an angstrem - mapped to different bins in the potential
            bin = bintable[int(math.floor(ca/0.5))]
            hist[aa][bin] += 1
    # normalize
    sums = {aa: float(sum(h)) for aa, h in hist.iteritems()}
    # print "SUMS", sums["AA"]
    # print "HIST", hist["AA"]
    # print "MIN", min(sums.values())

    for aa, h in hist.iteritems():
        for bin, v in enumerate(h):
            hist[aa][bin] = v / sums[aa]

    epsilons = {aa: math.pow(10, math.ceil(math.log10(1.0 / sums[aa]))) for aa in sums.iterkeys()}
    return hist, epsilons


def get_potential(bintable):
    nbins = bintable[int(math.floor(30.0/0.5))] + 1
    h_nat, e_nat = get_hist("distance_stats_2.tab", bintable)
    h_back, e_back = get_hist("distance_stats_2_shuffled.tab", bintable)
    epsilons = {aa: max(e_nat[aa], e_back[aa]) for aa in e_nat.iterkeys()}
    # print "NAT:", h_nat["AA"]
    # print "BACK:", h_back["AA"]
    # print "NAT_EPS:", e_nat["AA"]
    # print "BACK_EPS:", e_back["AA"]
    p = defaultdict(lambda: [0.0]*nbins)
    for aa in h_nat.iterkeys():
        for bin in range(len(h_nat[aa])):
            v_nat = h_nat[aa][bin]
            v_back = h_back[aa][bin]
            if v_nat > epsilons[aa] and v_back > epsilons[aa]:
                # p[aa][bin] = - kB*300.0*round(math.log(v_nat / v_back), 5)
                p[aa][bin] = round(math.log(v_nat / v_back), 5)
    return p


def score_interfaces(fname, bintable, potential, shuffled = False):
    empty = np.array([.0, .0, .0])

    # set visibility scope for these vars
    k_resn = None
    l_resn = None

    with open(fname, 'w') as o:
        for pdb, structure in distributions_iterator():
            for i, chain1 in enumerate(structure.iterkeys()):
                for j, chain2 in enumerate(structure.iterkeys()):
                    if j > i:
                        # CALC SCORE FOR EACH INTERFACE (A PAIR OF INTERACTING CHAINS)

                        # ##### Prepare lists of residues for random choice of residue #####                        
                        if shuffled:
                            k_set = set() 
                            l_set = set()
    
                            for k, residue1 in enumerate(structure[chain1]):
                                for l, residue2 in enumerate(structure[chain2]):
                                    resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1
                                    resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                                    dCa = np.linalg.norm(Ca2 - Ca1)
                                    if not 2.0 < dCa < 30.0: continue
                                    k_set.add(k)
                                    l_set.add(l)

                            k_resn = [structure[chain1][k][0] for k in k_set]
                            l_resn = [structure[chain2][l][0] for l in l_set]
                        ##########################

                        nres = len(structure[chain1]) + len(structure[chain2])
                        score = 0.0
                        pairwise_distances = []

                        if len(structure[chain1]) < 5 or len(structure[chain2]) < 5: continue

                        for k, residue1 in enumerate(structure[chain1]):
                            for l, residue2 in enumerate(structure[chain2]):
                                resn1, resi1, Ca1, Cb1, vCb1, quaternion1 = residue1
                                resn2, resi2, Ca2, Cb2, vCb2, quaternion2 = residue2
                                # print Ca2, Ca1
                                dCa = np.linalg.norm(Ca2 - Ca1)

                                pairwise_distances.append(dCa)
                                # S = 0.0 # Ca2 - Ca1 score

                                #     print "Too small distance between Ca", pdb, chain1, resi1, "and", chain2, resi2, "d=",dCa
                                #     # print "Too long distance between Ca", pdb, chain1, chain2, dCa
                                if not 2.0 < dCa < 30.0: continue

                                ##### Randomize amino acid type #####
                                if shuffled:
                                    resn1 = random.choice(k_resn)
                                    resn2 = random.choice(l_resn)
                                #####################################
                                # res1, res2 = (resn1, resn2) if aa.index(resn1) <= aa.index(resn2) else (resn2, resn1)
                                # res = "{}_{}".format(res1, res2)

                                # bin = binfunc(dCa)
                                bin = bintable[int(math.floor(dCa/0.5))]
                                r1, r2 = (resn1, resn2) if aa.index(resn1) <= aa.index(resn2) else (resn2, resn1)
                                pair = "{}_{}".format(r1, r2)
                                S = potential[pair][bin]

                                # print "AA_pair", pdb, chain1, resn1, resi1, chain2, resn2, resi2, "CA=", dCa, "S=", S, "score=",score
                                score += S
                                # if res == "ASP_ASP":
                                #     print "======= ASP-ASP ========"
                                #     print pdb, chain 
                                #     print residue1
                                #     print residue2
                                #     print res, round(dCa,2) #, dCb, dvCb, cos_theta, theta
                                # continue

                        # FOR EACH PPI INTERFACE WRITE DATA
                        # print res, dCa, dCb, dvCb, cos_theta, theta
                        if nres > 0:
                            o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, nres, np.mean(dCa), score))


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


if __name__ == '__main__':

    # import multiprocessing as mp
    # import time
    # pool = mp.Pool(8, init_worker)

    bintable = get_bin_table(bin47)
    print "Loading potential..."
    potential = get_potential(bintable)
    print "Applying F1"
    score_interfaces("real_scores_bin47.tab", bintable, potential)
    # pool.apply_async(score_interfaces, args=("real_scores_bin47.tab", bintable, potential))
    print "Applying F2"
    score_interfaces("shuffl_scores_bin47.tab", bintable, potential, True)
    # pool.apply_async(score_interfaces, args=("shuffl_scores_bin47.tab", bintable, potential, True))


    bintable = get_bin_table(bin20)
    print "Loading potential..."
    potential = get_potential(bintable)
    print "Applying F3"
    score_interfaces("real_scores_bin20.tab", bintable, potential)
    # pool.apply_async(score_interfaces, args=("real_scores_bin20.tab", bintable, potential))
    print "Applying F4"
    score_interfaces("shuffl_scores_bin20.tab", bintable, potential, True)
    # pool.apply_async(score_interfaces, args=("shuffl_scores_bin20.tab", bintable, potential, True))

    # try:
    #     while(True):
    #         print "Watchdog... every 60 seconds (Ctrl-C to interrupt)"
    #         time.sleep(60)
    # except KeyboardInterrupt:
    #     print "Caught KeyboardInterrupt, terminating workers"
    #     pool.terminate()
    # else:
    #     pool.close()
    # finally:
    #     pool.join()

