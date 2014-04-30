"""
Read distance stats for pairs of amino acids in natural and shuffled PPI interfaces
calculate Log odds (LOD)
Save LOD to AAIndex matrix format:
ARNDCQEGHILKMFPSTWYV
"""
from collections import defaultdict
import math
from scipy.constants import k as kB

AAindex = {}
aa_list = "ARNDCQEGHILKMFPSTWYV"
AA_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
for i, a in enumerate(aa_list):
    for j, b in enumerate(aa_list):
        # the order of residues and lower diagonal of the matrix according to AAIndex DB
        if j > i: continue
        aa = a+b
        AA = AA_list[i]+"_"+AA_list[j]
        AAindex[AA] = aa
        AA = AA_list[j]+"_"+AA_list[i]
        AAindex[AA] = aa

# def AAIndex(AA):
#     return AA

def get_hist(fname):
    hist = defaultdict(lambda: [1]*20)
    with open(fname) as f:
        for line in f:
            AA, ca, _ = line.split("\t", 2)
            # aa = AAIndex(AA)
            aa = AAindex[AA]
            ca = float(ca)
            # ca_round = int(round(ca, 0)) 
            
            # if ca_round > 20: continue
            # if ca_round < 2: continue
            
            # it should be 18 bins  0:2ca-3ca, 1:3ca-4ca, 2:4ca-5ca, ... 19-20

            # if 0.0 < ca <= 2.0: bin = 0
            # if 2.0 < ca <= 4.0: bin = 1
            # if 4.0 < ca <= 6.0: bin = 2
            # if 6.0 < ca <= 8.0: bin = 3
            # if 8.0 < ca <= 10.0: bin = 4
            # if 10.0 < ca <= 12.0: bin = 5
            # if 12.0 < ca <= 14.0: bin = 6
            # if 14.0 < ca <= 16.0: bin = 7
            # if 16.0 < ca <= 18.0: bin = 8
            # if 18.0 < ca <= 20.0: bin = 9
            # if 20.0 < ca       : bin = 10

            # 0-7, 7-14, 14-20, 21-
            bin = int(math.floor(ca/3.0))
            if ca > 20: bin = 10

            # bin = int(math.floor(ca/10.0))
            # if ca > 20: bin = 2

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

def get_potential():
    h_nat, e_nat = get_hist("distance_stats_2.tab")
    h_back, e_back = get_hist("distance_stats_2_shuffled.tab")
    epsilons = {aa: max(e_nat[aa], e_back[aa]) for aa in e_nat.iterkeys()}
    # print "NAT:", h_nat["AA"]
    # print "BACK:", h_back["AA"]
    # print "NAT_EPS:", e_nat["AA"]
    # print "BACK_EPS:", e_back["AA"]
    p = defaultdict(lambda: [0.0]*20)
    for aa in h_nat.iterkeys():
        for bin in range(len(h_nat[aa])):
            v_nat = h_nat[aa][bin]
            v_back = h_back[aa][bin]
            if v_nat < epsilons[aa] or v_back < epsilons[aa]:
                p[aa][bin] = 0.0
            else:
                # p[aa][bin] = - kB*300.0*round(math.log(v_nat / v_back), 5)
                p[aa][bin] = round(math.log(v_nat / v_back), 5)
    return p

def write_potential_as_matrix(p, fname=""):
    with open(fname, 'w') as o:
        for bin in range(len(p["AA"])):
            o.write("H GONA0000{}\n".format(bin))
            o.write("D Ca-Ca distance potential for PPI interface residues\n")
            o.write("A Goncearenco, A.\n")
            o.write("M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV\n")
            for i, a in enumerate(aa_list):
                l = []
                for j, b in enumerate(aa_list):
                    if j > i: continue
                    l.append(str(p[a+b][bin]))
                o.write("  " + "  ".join(l) + "\n")
            o.write("//\n")

if __name__ == '__main__':

    p = get_potential()
    print "AA:", p["AA"]
    write_potential_as_matrix(p, "potential_3.index")




