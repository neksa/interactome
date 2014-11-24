"""
Read distance stats for pairs of amino acids in natural and shuffled PPI interfaces
calculate Log odds (LOD)
Save LOD to AAIndex matrix format:
ARNDCQEGHILKMFPSTWYV
"""
import math
# import sys
from collections import defaultdict
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

# print AAindex
# def AAIndex(AA):
#     return AA

def bin_1(ca):
    b = int(round(ca)) - 2
    if b < 0:
        b = 0
    if b > 15-2:
        b = 15-2
    return b


def get_hist(fname, binfunc):
    hist = defaultdict(lambda: [1]*20)
    with open(fname) as f:
        for line in f:
            fields = line.split("\t", 2)
            if len(fields) != 3:
                continue
            AA, ca, _ = fields
            aa = AA[0] + AA[2]  # A_A -> AA
            ca = float(ca)

            # aa = AAindex[AA]
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
            # bin = int(math.floor(ca/3.0))
            # if ca > 20: bin = 10

            # 2-3, 3-4, 4-5
            b = binfunc(b)
            # bin = int(math.floor(ca/10.0))
            # if ca > 20: bin = 2
            hist[aa][b] += 1

    # print hist
    # normalize
    sums = {aa: float(sum(h)) for aa, h in hist.iteritems()}
    # print "SUMS", sums["AA"]
    # print "HIST", hist["AA"]
    # print "MIN", min(sums.values())
    for aa, h in hist.iteritems():
        for b, v in enumerate(h):
            hist[aa][b] = v / sums[aa]
    epsilons = {aa: math.pow(10, math.ceil(math.log10(5.0 / sums[aa]))) for aa in sums.iterkeys()}
    return hist, epsilons


def get_potential(fname1, fname2):
    h_nat, e_nat = get_hist(fname1)
    h_back, e_back = get_hist(fname2)
    epsilons = {aa: max(e_nat[aa], e_back[aa]) for aa in e_nat.iterkeys()}
    # print "NAT:", h_nat["AA"]
    # print "BACK:", h_back["AA"]
    # print "NAT_EPS:", e_nat["AA"]
    # print "BACK_EPS:", e_back["AA"]
    p = defaultdict(lambda: [0.0]*20)
    for aa in h_nat.iterkeys():
        for b in range(len(h_nat[aa])):
            v_nat = h_nat[aa][b]
            v_back = h_back[aa][b]
            # print v_nat, v_back
            if v_nat < epsilons[aa] or v_back < epsilons[aa]:
                p[aa][b] = 0.0
            else:
                p[aa][b] = round(math.log(v_nat / v_back), 6)
            # else:
            #     # p[aa][bin] = - kB*300.0*round(math.log(v_nat / v_back), 5)
            #     # p[aa][b] = round(math.log(v_nat / v_back), 10)
            # print aa, b, v_nat, v_back, p[aa][b]
            # print p[aa][b]
    return p


def write_potential_as_matrix(p, fname=""):
    with open(fname, 'w') as o:
        for bin_num in range(len(p["AA"])):
            o.write("H GONA0000{}\n".format(bin_num))
            o.write("D Ca-Ca distance potential for PPI interface residues ({}-{}A)\n".format(bin_num + 1.5, bin_num + 2.5))
            o.write("A Goncearenco, A.\n")
            o.write("M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV\n")
            for i, a in enumerate(aa_list):
                l = []
                for j, b in enumerate(aa_list):
                    if j > i:
                        continue
                    v = p.get(a+b)
                    if v is None:
                        v = p.get(b+a)
                    # v[bin_num]
                    l.append(str(v[bin_num]))
                o.write("  " + "  ".join(l) + "\n")
            o.write("//\n")

    # symmetric, redundant, easy to load
    with open(fname + '-warp', 'w') as o:
        for i, a in enumerate(aa_list):
            for j, b in enumerate(aa_list):
                v = p.get(a+b)
                if v is None:
                    v = p.get(b+a)
                l = []
                l.append(a+b)
                for bin_num in range(len(p["AA"])):
                    l.append('{0: >9}'.format(v[bin_num]))
                o.write("\t".join(l) + "\n")


if __name__ == '__main__':
    def get_root():
        return "/Users/agoncear/projects/Interactome/Workflow"

    f1 = get_root() + "/Potential/distance_stats.tab"
    f2 = get_root() + "/Potential/distance_stats_shuffled.tab"
    f_pot = get_root() + "/Potential/potential_1.index"
    p = get_potential(f1, f2)
    # print "A-C:", p["AC"]
    write_potential_as_matrix(p, f_pot)
