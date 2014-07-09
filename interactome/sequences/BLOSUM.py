"""
"""

"""
Read BLOSUM matrix
"""
def read_BLOSUM_matrix(matrix):
    m = {}
    amino_acids = "ARNDCQEGHILKMFPSTWYV"
    with open(matrix) as f:
        for line in f:
            P = line[0]
            if P == "#": continue
            if P == " ": continue
            if P == "*": continue
            if P not in amino_acids: continue
            m[P] = {}
            for i in range(20):
                Q = amino_acids[i]
                # if Q == "*": continue
                value = int(line[i*3 + 1: i*3 + 4].strip())
                m[P][Q] = value
    return(m)

# def export_blosum():
#     m = read_BLOSUM_matrix("BLOSUM62")
#     with open("BLOSUM62.csv", 'w') as o:
#         for P in m.iterkeys():
#             for Q in m[P].iterkeys():
#                 o.write("{}>{}\t{}\n".format(P, Q, m[P][Q]))


# if __name__ == '__main__':
#     export_blosum()
