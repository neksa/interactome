

def templates_for_coexpression(fname):
    unique_set = set()
    out_fname = "protein_pairs_" + fname.split('.')[0] + '_unique.tab'
    with open(fname, 'r') as f, open(out_fname, 'w') as o:
        for line in f:
            fields = line.split("\t")
            a = fields[4]
            b = fields[5]

            if len(a) == 0 or len(b) == 0: continue
            # expect a valid uniprot/swissprot id
            if len(a) != 6 or len(b) != 6: continue
            if a[0] < 'A' or b[0] < 'A': continue

            if a == b: continue

            print a,b 

            pair = (a,b)
            if a > b:
                pair = (b, a)
            unique_set.add(pair)

        for pair in unique_set:
            a, b = pair
            o.write("{}\t{}\n".format(a, b))



def for_coexpression(fname):
    unique_set = set()
    out_fname = "protein_pairs_" + fname.split('.')[0] + '_unique.tab'
    with open(fname, 'r') as f, open(out_fname, 'w') as o:
        for line in f:
            a, b, _ = line.split("\t", 2)
            print a,b 

            if a == b: continue

            pair = (a,b)
            if a > b:
                pair = (b, a)
            unique_set.add(pair)

        for pair in unique_set:
            a, b = pair
            o.write("{}\t{}\n".format(a, b))



def for_interactome():
    unique_set = set()
    with open("matches_human_25.tab", 'r') as f, open("draft_interactome.tab", 'w') as o:
        for line in f:
            a = ""
            b = ""
            c = ""
            try:
                a,b,c = line.split("\t", 3)
            except:
                pass

            pair = (a,b)
            if a > b:
                pair = (b, a)
            unique_set.add(pair)

        for pair in unique_set:
            a, b = pair
            o.write("{}\t{}\n".format(a, b))


if __name__ == '__main__':
    # for_interactome()
    templates_for_coexpression("template_analysis.tab")
    # for_coexpression("matches_human_25.tab")
    # for_coexpression("matches_human_20-25.tab")
    # for_coexpression("matches_human_15-20.tab")
