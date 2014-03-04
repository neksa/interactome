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
* M rowS[0] = ARNDCQEGHILKMFPSTWYV, colS[0] = ARNDCQEGHILKMFPSTWYV           *
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

AAINDEX_PATH = "/Users/agoncear/data/AAIndex/"


def aaindex1_parser():
    """
    State machine parser.
    States: H -> D -> A -> I -> I1 -> I2 -> H
    """
    def state(S, item, l):
        header = l[0]
        data = l[2:].strip()

        # print "State", S[0]
        # print "State", S[0], header, data
        aa = "ARNDCQEGHILKMFPSTWYV"

        def next_state(S):
            states = "HDAI01/"
            i = states.index(S[0])
            j = (i + 1) % len(states)
            # print "Next", S[0], i, j
            S[0] = states[j]
            return i > 0 and j == 0

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

        # if S[0] == "H":
        #     if header == S[0]:
        #         item[S[0]] = data
        #         S[0] = "D"
        # elif S[0] == "D":
        #     if header == S[0]:
        #         item[S[0]] = data
        #         S[0] = "A"
        # elif S[0] == "A":
        #     if header == S[0]:
        #         item[S[0]] = data
        #         S[0] = "I"
        # elif S[0] == "I":
        #     if header == S[0]:
        #         item[S[0]] = {}
        #         S[0] = "I1"
        # elif S[0] == "I1":
        #     if header == " ":
        #         # print item["H"], data
        #         for i, v in enumerate(data.split(), start=0):
        #             # print i, v,
        #             # print float(v)
        #             item["I"][aa[i]] = float(v) if v != "NA" else None
        #         S[0] = "I2"
        # elif S[0] == "I2":
        #     if header == " ":
        #         # print item["H"], data
        #         for i, v in enumerate(data.split(), start=10):
        #             # print i, v,
        #             # print float(v)
        #             item["I"][aa[i]] = float(v) if v != "NA" else None
        #         S[0] = "/"
        # elif S[0] == "/":
        #     if header == S[0]:
        #         S[0] = "H"
        #         return True
        # else:
        #     raise Exception("Undefined state in parser: %s" % S[0])
        # return False


    with open(AAINDEX_PATH + "aaindex1") as f:
        S = ["H",]
        item = {}
        for line in f:
            if state(S, item, line):
                yield item
                item = {}


def aaindex1_get_list():
    return [(index["H"], index["D"], index["A"]) for index in aaindex1_parser()]


def aaindex1_get_data(accession):
    for index in aaindex1_parser():
        if index["H"] == accession:
            return index["I"]

if __name__ == '__main__':
    for i in aaindex1_get_list():
        print "\t".join(i)

    print aaindex1_get_data("CIDH920104")
    print aaindex1_get_data("MIYS990101")
    print aaindex1_get_data("MIYS990102")
    print aaindex1_get_data("MIYS990103")
    print aaindex1_get_data("MIYS990104")
    print aaindex1_get_data("MIYS990105")



