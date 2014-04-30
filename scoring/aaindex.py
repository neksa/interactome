"""
AAIndex database parser

Data Format of AAindex1
************************************************************************
*                                                                      *
* Each entry has the following format.                                 *
*                                                                      *
* H Accession number                                                   *
* D Data description                                                   *
* R LITDB entry number                                                 *
* A Author(s)                                                          *
* T Title of the article                                               *
* J Journal reference                                                  *
* * Comment or missing                                                 *
* C Accession numbers of similar entries with the correlation          *
*   coefficients of 0.8 (-0.8) or more (less).                         *
*   Notice: The correlation coefficient is calculated with zeros       *
*   filled for missing values.                                         *
* I Amino acid index data in the following order                       *
*   Ala    Arg    Asn    Asp    Cys    Gln    Glu    Gly    His    Ile *
*   Leu    Lys    Met    Phe    Pro    Ser    Thr    Trp    Tyr    Val *
* //                                                                   *
************************************************************************
(Data Format of AAindex2 and AAindex3)
************************************************************************
*                                                                      *
* Each entry has the following format.                                 *
*                                                                      *
* H Accession number                                                   *
* D Data description                                                   *
* R LITDB entry number                                                 *
* A Author(s)                                                          *
* T Title of the article                                               *
* J Journal reference                                                  *
* * Comment or missing                                                 *
* M rowS[0] = ARNDCQEGHILKMFPSTWYV, colS[0] = ARNDCQEGHILKMFPSTWYV     *
*   AA                                                                 *
*   AR RR                                                              *
*   AN RN NN                                                           *
*   AD RD ND DD                                                        *
*   AC RC NC DC CC                                                     *
*   AQ RQ NQ DQ CQ QQ                                                  *
*   AE RE NE DE CE QE EE                                               *
*   AG RG NG DG CG QG EG GG                                            *
*   AH RH NH DH CH QH EH GH HH                                         *
*   AI RI NI DI CI QI EI GI HI II                                      *
*   AL RL NL DL CL QL EL GL HL IL LL                                   *
*   AK RK NK DK CK QK EK GK HK IK LK KK                                *
*   AM RM NM DM CM QM EM GM HM IM LM KM MM                             *
*   AF RF NF DF CF QF EF GF HF IF LF KF MF FF                          *
*   AP RP NP DP CP QP EP GP HP IP LP KP MP FP PP                       *
*   AS RS NS DS CS QS ES GS HS IS LS KS MS FS PS SS                    *
*   AT RT NT DT CT QT ET GT HT IT LT KT MT FT PT ST TT                 *
*   AW RW NW DW CW QW EW GW HW IW LW KW MW FW PW SW TW WW              *
*   AY RY NY DY CY QY EY GY HY IY LY KY MY FY PY SY TY WY YY           *
*   AV RV NV DV CV QV EV GV HV IV LV KV MV FV PV SV TV WV YV VV        *
* //                                                                   *
************************************************************************
"""

from scipy import stats
from scipy.stats import spearmanr, pearsonr
import numpy as np
import operator
import sys

AAINDEX_PATH = "/Users/agoncear/data/AAIndex/"


def aaindex1_parser():
    """
    State machine parser. Vector of values
    """

    def next_state(S):
        states = "HDAI01/"
        i = states.index(S[0])
        j = (i + 1) % len(states)
        # print "Next", S[0], i, j
        S[0] = states[j]
        return i > 0 and j == 0
    
    def state(S, item, l):
        aa = "ARNDCQEGHILKMFPSTWYV"
        header = l[0]
        data = l[2:].strip()

        if header == S[0]:
            if S[0] == "I":
                item[S[0]] = {}
            elif S[0] != "/":
                item[S[0]] = data
            return next_state(S)
        elif header == " ":
            if S[0] in ("0", "1"):
                for i, v in enumerate(data.split(), start=int(S[0]) * 10):
                    # print i, v,
                    # print float(v)
                    item["I"][aa[i]] = float(v) if v != "NA" else None
                return next_state(S)
        else:
            return False # ignore other lines

    with open(AAINDEX_PATH + "aaindex1") as f:
        S = ["H",]
        item = {}
        for line in f:
            if state(S, item, line):
                yield item
                item = {}


def aaindex23_parser(fname):
    """
    State machine parser. Matrix of values
    """

    def next_state(S):
        states = "HDAM0/"
        i = states.index(S[0])
        j = (i + 1) % len(states)
        # print "Next", S[0], i, j, states[j]
        S[0] = states[j]
        return i > 0 and j == 0
    
    def state(S, item, l):
        header = l[0]
        data = l[2:].strip()

        # print "State", S[0], header, data
        if header == S[0]:
            if S[0] == "M":
                # rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
                rows = data.split(",")[0].split()[-1]
                cols = data.split(",")[1].split()[-1]
                item[S[0]] = (rows, cols)
                item["0"] = []
            elif S[0] != "/":
                item[S[0]] = data
            return next_state(S)
        elif S[0] == "0" and header == " ":
            d = [float(v) if v not in ("-", "NA") else None for v in data.split()]
            item["0"].append(d)
            if len(item["0"]) < 20:
                return False
            return next_state(S)
        else:
            return False # ignore other lines

    # aa = "ARNDCQEGHILKMFPSTWYV"
    with open(fname) as f:
        S = ["H",]
        item = {}
        for line in f:
            if state(S, item, line):
                item["DATA"] = {}
                rows, cols = item["M"]
                for i, values in enumerate(item["0"]):
                    for j, v in enumerate(values):
                        # print i, j, v
                        item["DATA"][rows[i] + cols[j]] = v
                del item["0"]
                yield item
                item = {}


def aaindex1_get_list():
    return [(index["H"], index["D"], index["A"]) for index in aaindex1_parser()]

def aaindex23_get_list(fname):
    return [(index["H"], index["D"], index["A"]) for index in aaindex23_parser(fname)]

def aaindex1_get_data(accession):
    for index in aaindex1_parser():
        if index["H"] == accession:
            return index["I"]

def aaindex23_get_data(accession, fname):
    for index in aaindex23_parser(fname):
        if index["H"] == accession:
            return index["DATA"]


def _get_two_arrays(data1, data2):
    keys1 = set(data1.keys())
    keys2 = set(data2.keys())
    if len(keys1) != len(keys2):
        # print "Error: the number of keys is different in data1 and data2"
        return None
    if len(keys1 - keys2) > 0:
        # print "Error: keys are different in data1 and data2"
        return None
    # assume keys are equivalent
    d1 = []
    d2 = []
    for k in keys1:
        d1.append(data1[k])
        d2.append(data2[k])
    return np.array(d1), np.array(d2)


def pearson_corr(data1, data2):
    d1, d2 = _get_two_arrays(data1, data2)
    r = pearsonr(d1, d2)
    return r

# Rank correlation
def spearman_corr(data1, data2):
    d1, d2 = _get_two_arrays(data1, data2)
    rho = spearmanr(d1, d2)
    return rho


if __name__ == '__main__':
    # for i in aaindex1_get_list():
    #     print "\t".join(i)

    # print aaindex1_get_data("CIDH920104")
    # print aaindex1_get_data("MIYS990101")
    # print aaindex1_get_data("MIYS990102")
    # print aaindex1_get_data("MIYS990103")
    # print aaindex1_get_data("MIYS990104")
    # print aaindex1_get_data("MIYS990105")

    # for i in aaindex23_get_list(2):
    #     print "\t".join(i)
    # print aaindex23_get_data("PRLA000101", 2)
    # print len(aaindex23_get_data("PRLA000101", 2))

    # check matrices type 2 and 3
    # for matrix in range(2, 4):
    #     for i in aaindex23_get_list(matrix):
    #         print matrix, i[0], len(aaindex23_get_data(i[0], matrix))

    # for i in aaindex23_get_list(3):
    #     print "\t".join(i)    
    # print aaindex23_get_data("TANS760102", 3)
    # print len(aaindex23_get_data("TANS760102", 3))

    # compare BRYS930101 MIYS850103 MIYS960102 MIYS990107 KESO980102 MOOG990101 

    def fix_missing(d):
        return {a: d[a] if d[a] is not None else 0.0 for a in d.iterkeys()}

    index2 = AAINDEX_PATH + "aaindex2"
    index3 = AAINDEX_PATH + "aaindex3"
    index_gon = "potential_3.index"

    # compare = ["BRYS930101", "MIYS850103", "MIYS960102", "MIYS990107", "KESO980102", "MOOG990101"]
    index3 = index_gon
    compare = [x[0] for x in aaindex23_get_list(index3)]
    unsorted_list = []
    for i, c1 in enumerate(compare):
        for j, c2 in enumerate(compare):
            if i >= j: continue
            # print c1, c2
            data1 = fix_missing(aaindex23_get_data(c1, index3))
            data2 = fix_missing(aaindex23_get_data(c2, index3))
            try:
                r = pearson_corr(data1, data2)[0]
                rho = spearman_corr(data1, data2)[0]
                # print c1, c2, r
                if np.isnan(r) or np.isnan(rho):
                    # print "NAN"
                    raise Exception("Not a number")
                l1 = c1
                l2 = c2
                l1 = i
                text = "{}\t{}\tR={}\trho={}".format(l1, l2, round(r, 2), round(rho, 2))
                # print text
                unsorted_list.append((abs(r), text))
            except:
                pass
    # print unsorted_list
    sorted_list = sorted(unsorted_list, key=operator.itemgetter(0), reverse=True)
    # print "Top 50 pairs of correlating potentials ordered by Pearson R:"
    for i in range(len(sorted_list)):
        print sorted_list[i][1]

    sys.exit(0)
    print "####################################"

    # accession_gon = "GONA00000"

    print "Analyzing potentials in ", index_gon
    for idx in range(7):
        accession_gon = "GONA0000{}".format(idx)
        print "ACCESSION:", accession_gon
        data_gon = aaindex23_get_data(accession_gon, index_gon)

        # accession_energies = "MIYS850103"
        # accession_energies = "BRYS930101"
        # data_energies = aaindex23_get_data(accession_energies, index3)
        # my = {aa: data_energies[aa] * data_gon[aa] for aa in data_gon.iterkeys()}
        my = data_gon

        print "compare {} from {} with all the matrices".format(accession_gon, index_gon)
        unsorted_list = []
        for i in aaindex23_get_list(index3):
            # print i[0], len(aaindex23_get_data(i[0], matrix))
            database = aaindex23_get_data(i[0], index3)
            database = fix_missing(database)
            try:
                r = pearson_corr(my, database)
                rho = spearman_corr(my, database)
                text = "{}\tR={}\trho={}\t{:80s} ({})".format(i[0], round(r[0], 2), round(rho[0], 2), i[1], i[2])
                unsorted_list.append((abs(r[0]), text))
            except Exception, e:
                pass
                # print "skipping matrix", i[0]
        sorted_list = sorted(unsorted_list, key=operator.itemgetter(0), reverse=True)
        print "Top 10 matching matrices ordered by Pearson R:"
        for i in range(10):
            print sorted_list[i][1]

