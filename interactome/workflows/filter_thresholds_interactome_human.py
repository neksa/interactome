from itertools import islice
from collections import defaultdict, namedtuple
import gzip


def get_scored_pairs(matches_fname):
    # pairs = defaultdict(lambda: defaultdict(list))
    max_pairs = defaultdict(dict)
    with gzip.open(matches_fname) as f:
        for line in islice(f, 1, None):
            # print "L", len(line.strip().split("\t"))
            # for i, a in enumerate(line.strip().split("\t")):
            #     print i, "==", a, "=="
            # 11 + 2 + 11 + 11 = 35
            try:
                queryA, queryB, tpl, query_type, template_type, template_nsubunits, \
                    SC1, SC2, SC3, SC4, SC5, SC6, \
                    identicalA, poitiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                    identicalB, poitiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                    siteA, siteB = line.strip().split("\t")
                    # _, _, _, _, \

                SC1, SC2, SC3, SC4, SC5, SC6 = map(float, (SC1, SC2, SC3, SC4, SC5, SC6))
                template_nsubunits = int(template_nsubunits)

                identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A = map(
                    float, (identicalA, poitiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A))

                identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B = map(
                    float, (identicalB, poitiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B))

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                pair = (queryA, queryB)
                # if queryA > queryB:
                #     pair = (queryB, queryA)

                full_identityA = identicalA / aln_lenA
                full_identityB = identicalB / aln_lenB
                full_identity = min(full_identityA, full_identityB)  # FULL SEQUENCE ALN

                bs_similarity = min(bs_positiveA / bs_lenA, bs_positiveB / bs_lenB)

                # FILTERS:
                if template_type != query_type:
                    continue

                # if bs_similarity < 0.25:
                #     continue

                # if bs_similarity < 0.5:
                #     continue

                # Combined score
                # score = 0.067 * SC1 * 5.07 * bs_similarity
                score = (1-full_identity) * 0.067 * SC1 + (full_identity) * 5.07 * bs_similarity

                # Choose the best scoring model for each protein pair
                s = max_pairs[pair].get("score")
                if s is None or score > s:
                    max_pairs[pair]["score"] = score
                    max_pairs[pair]["score1"] = SC1
                    max_pairs[pair]["bs_similarity"] = bs_similarity
                    max_pairs[pair]["sites"] = (siteA, siteB)
                    max_pairs[pair]["tpl_nsubunits"] = template_nsubunits
                    max_pairs[pair]["tpl"] = tpl
                    max_pairs[pair]["siteA"] = siteA
                    max_pairs[pair]["siteB"] = siteB

                # print pair, pairs["model_score"][pair]
                # pairs["bs_positive_homo"][pair].append()
                #     if query_type == "Unknown" or template_type == "Unknown" or query_type == template_type:
                #         score = min(
                #             bs_positiveA / bs_lenA,
                #             bs_positiveB / bs_lenB)
                #     else:
                #         score = 0.0
            except:
                print "Error reading the whole matches file. Skipping the line with an error"

    # for scoring, pair_scores in pairs.iteritems():
    #     for pair, scores in pair_scores.iteritems():
    #         max_pairs[pair][scoring] = max(scores)
    return max_pairs


# def get_uniprot_mapping():
#     from benchmark_residues_yeast import load_uniprot
#     locus_uniprot, uniprot_nresidues = load_uniprot()
#     return locus_uniprot


# def get_uniprot_mapping():
#     mapping_fname = "/Users/agoncear/data/Ecoli/gi_uniprot_mapping.tab"
#     # swissprot_fname = "/Users/agoncear/data/Ecoli/ecoli.txt"
#     mapping = {}
#     inv_mapping = {}
#     with open(mapping_fname) as f:
#         for line in f:
#             gis, uniprot = line.strip().split()
#             for gi in gis.split(","):
#                 mapping[gi] = uniprot
#                 inv_mapping[uniprot] = gi
#     return mapping


UniprotRec = namedtuple('UniprotRec', 'id accessions gene descr')
def iter_uniprot():
    fname = "/Users/agoncear/data/Uniprot/Uniprot_Release_2014_11/uniprot_sprot.dat.gz"
    i = ""
    rec = defaultdict(str)
    with gzip.open(fname) as f:
        for l in f:
            if l.startswith("//"):
                # if rec['OX'].strip() == "NCBI_TaxID=9606;": pass
                id_str = rec['ID'].split()[0]
                # AC   P31946; A8K9K2; E1P616;
                accessions = [ac[:-1] for ac in rec['AC'].strip().split()]
                # GN   Name=YWHAB; Synonyms=YWHA1;
                gene = rec['GN'].split(';')[0][5:]
                # DE   RecName: Full=HLA class I histocompatibility antigen, A-1 alpha chain;
                descr = ''
                if 'DE' in rec:
                    for d in rec['DE'].split("\n"):
                        if d.startswith("RecName"):
                            # print "SSS", d.split("=")
                            descr = d.split("=")[1][:-1]
                uniprot_rec = UniprotRec._make((id_str, accessions, gene, descr))
                # print uniprot_rec
                yield uniprot_rec

                i = ""
                rec = defaultdict(str)
                continue

            if len(l) > 5:
                nid = l[0:2]
                data = l[5:]
                if nid != "  ":
                    i = nid
                if i not in ('DE', 'GN', 'ID', 'OX', 'AC'):
                    continue
                rec[i] += data


def iter_uniprot_mapping():
    import gzip
    fname = "/Users/agoncear/data/Uniprot/Uniprot_Release_2014_11/idmapping.kb.dat.gz"
    with gzip.open(fname) as f:
        for c, l in enumerate(f):
            if c % 1000000 == 0:
                print c
            uid, utype, uacc = l.strip().split("\t", 2)
            if utype == 'UniProtKB-ID':
                yield uid, uacc

# def iter_custom_uniprot_mapping():
#     fname = "/Users/agoncear/projects/Interactome/Workflow/Interactomes/Mappings/human.txt"
#     with open(fname) as f:
#         for l in f:
#             rec = l.strip().split("\t")
#             if len(rec) == 2:
#                 uid, uacc = rec
#                 yield uid, uacc


def uniprot_accessions():
    import itertools
    accessions = {}
    # for uid, uacc in itertools.chain(iter_custom_uniprot_mapping(), [], [], []):
    for uid, uacc in itertools.chain(iter_uniprot_mapping(), [], [], []):
        accessions[uid] = uacc
    return accessions


# def uniprot_accessions():
#     # from filter_interactome import UniprotRec, iter_uniprot
#     accessions = {}
#     for r in iter_uniprot():
#         # id accessions gene descr
#         for a in r.accessions:
#             accessions[a] = r.id
#     return accessions


def pdb_chain_protein():
    fname_pdb_proteins = "/Users/agoncear/projects/Interactome/Workflow/Structures/pdb_proteins.tab"
    pdb_protein = defaultdict(list)
    with open(fname_pdb_proteins) as f_pdb:
        for line in islice(f_pdb, 1, None):
            (pdb, chain, chain_author, uniprot,
                seq_aln_begin, seq_aln_begin_ins, seq_aln_end, seq_aln_end_ins,
                db_aln_begin, db_aln_end,
                auth_aln_begin, auth_aln_end) = line.strip().split("\t")
            if uniprot == '':
                continue
            pdb_protein[(pdb, chain)].append((uniprot, int(seq_aln_begin), int(db_aln_begin)))
    return pdb_protein


def load_interactions(fname):
    interactions = []
    with open(fname) as f:
        for line in islice(f, 1, None):
            gA, nameA, gB, nameB = line.strip().split("\t")
            # print gA, nameA, gB, nameB
            gA = int(gA)
            gB = int(gB)
            pair = (gA, gB)
            if gA > gB:
                pair = (gB, gA)
            interactions.append(pair)
    return interactions


def validation_vidal(clones_by_id, vidal_fname):
    print "Loading experimental interactome for HUMAN"
    pairs_vidal = load_interactions(vidal_fname)
    pairs = set()
    for (a, b) in pairs_vidal:
        if a not in clones_by_id:
            print "Interaction for gene", a, "not mapped to uniprot"
        if b not in clones_by_id:
            print "Interaction for gene", b, "not mapped to uniprot"

        ac_a = clones_by_id.get(a, None)
        ac_b = clones_by_id.get(b, None)
        if ac_a is not None and ac_b is not None:
            print ac_a.ac, ac_b.ac
            pairs.add((ac_a.ac, ac_b.ac))
    return pairs


CLONE = namedtuple('CLONE', "uniprot genes ac names length intact")
def get_vidal_clones(fname):
    clones = {}
    clones_by_id = {}
    with open(fname) as f:
        for line in islice(f, 1, None):
            uniprot, geneids, ac, names, length, intact = line.split("\t")
            genes = map(int, geneids.split(","))
            length = int(length)
            intact = 1 if intact != "" else 0
            clones[ac] = CLONE._make((uniprot, genes, ac, names, length, intact))
            for g in genes:
                clones_by_id[g] = clones[ac]
    return clones, clones_by_id


def savefile(infile, outfile, accessions=None, mapping=None, validation_pairs=None, clones=None):
    # m = None
    # if mapping is not None:
    #     m = mapping()
    pdb_protein = pdb_chain_protein()
    pairs = get_scored_pairs(infile)

    with open(outfile, 'w') as o:
        # o.write("A\tB\tSequence similarity\tInterface compatibility\tsiteA\tsiteB\n")
        fields = "UniprotID A\tUniprotAC A\tUniprotID B\tUniprotAC B\tScore\tBinding site seq. similarity\tInterface compatibility score\tTemplate\tTemplate complex subunits\tTemplate UniprotID A\tTemplate Uniprot AC A\tTemplate UniprotID B\tTemplate UniprotAC B\tIn PDB\tLA\tLB\tNCA\tNCB\tsiteA\tsiteB\tgenesA\tgenesB\tExperimental"
        o.write(fields+"\n")
        for pair, v in pairs.iteritems():
            A, B = pair
            # uniprotA = m.get(A, "NA")
            # uniprotB = m.get(B, "NA")
            uniprotA = A
            uniprotB = B
            uniprotAC_A = accessions.get(uniprotA, "NA")
            uniprotAC_B = accessions.get(uniprotB, "NA")
            if uniprotAC_A == "NA" or uniprotAC_B == "NA":
                continue
            print A, B
            if uniprotAC_A not in clones or uniprotAC_B not in clones:
                continue

            template = v["tpl"].replace("|", ".")
            pdb, chainA, chainB = v["tpl"].split("|")
            tpl_uniprotA = pdb_protein.get((pdb, chainA.split("_")[0]), ("NA", 0, 0))[0][0]
            tpl_uniprotB = pdb_protein.get((pdb, chainB.split("_")[0]), ("NA", 0, 0))[0][0]
            tpl_uniprotAC_A = accessions.get(tpl_uniprotA, "NA")
            tpl_uniprotAC_B = accessions.get(tpl_uniprotB, "NA")
            score = v["score"]
            similarity = v["bs_similarity"]
            compatibility = v["score1"]
            siteA, siteB = v["sites"]
            template_nsubunits = v["tpl_nsubunits"]

            cA = clones[uniprotAC_A]
            cB = clones[uniprotAC_B]
            len_prot_A = cA.length
            len_prot_B = cB.length
            genesA = cA.names
            genesB = cB.names

            print "Site"
            siteA = v["siteA"]
            resiA = filter(lambda resi: resi > 0, map(lambda s: int(s[4]), [s.split(",") for s in siteA.split(";")]))
            siteB = v["siteB"]
            resiB = filter(lambda resi: resi > 0, map(lambda s: int(s[4]), [s.split(",") for s in siteB.split(";")]))

            ss = []
            for r in [s.split(",") for s in siteA.split(";")]:
                ss.append("{}{}".format(r[3], r[4]))
            sA = ", ".join(ss)

            ss = []
            for r in [s.split(",") for s in siteB.split(";")]:
                ss.append("{}{}".format(r[3], r[4]))
            sB = ", ".join(ss)

            THRESHOLD = 20
            NCA = '-'  # Not on termini
            if min(resiA) <= THRESHOLD:
                NCA = 'N'  # On N-termini
            if len_prot_A - max(resiA) <= THRESHOLD:
                if NCA == "-":
                    NCA = ""
                NCA += 'C'  # On C-termini
            if max(resiA) > len_prot_A:
                NCA = 'X'  # Error, incorrect length
            if NCA == 'NC':
                NCA = 'B'  # On both termini

            NCB = '-'  # Not on termini
            if min(resiB) <= THRESHOLD:
                NCB = 'N'  # On N-termini
            if len_prot_B - max(resiB) <= THRESHOLD:
                if NCB == "-":
                    NCB = ""
                NCB += 'C'  # On C-termini
            if max(resiB) > len_prot_B:
                NCB = 'X'  # Error, incorrect length
            if NCB == 'NC':
                NCB = 'B'  # On both termini
            print NCA, NCB

            inPDB = 0
            if (uniprotA == tpl_uniprotA and uniprotB == tpl_uniprotB) \
               or (uniprotAC_A == tpl_uniprotAC_A and uniprotAC_B == tpl_uniprotAC_B and tpl_uniprotAC_A != 'NA' and tpl_uniprotAC_B != 'NA'):
                    inPDB = 1
                    tpl_uniprotA = ""
                    tpl_uniprotB = ""
                    tpl_uniprotAC_A = ""
                    tpl_uniprotAC_B = ""

            experimental = 0
            if validation_pairs is not None:
                if ((uniprotAC_A, uniprotAC_B) in validation_pairs or (uniprotAC_B, uniprotAC_A) in validation_pairs or
                   (uniprotA, uniprotB) in validation_pairs or (uniprotB, uniprotA) in validation_pairs):
                        experimental = 1

            strfmt = ["{}"] * len(fields.split("\t"))
            strfmt = "\t".join(strfmt) + "\n"
            o.write(strfmt.format(
                uniprotA, uniprotAC_A,
                uniprotB, uniprotAC_B,
                score, similarity, compatibility,
                template, template_nsubunits,
                tpl_uniprotA, tpl_uniprotAC_A,
                tpl_uniprotB, tpl_uniprotAC_B,
                inPDB, len_prot_A, len_prot_B, NCA, NCB, sA, sB, genesA, genesB, experimental))


if __name__ == '__main__':
    DIR = "/Users/agoncear/projects/Interactome/Workflow/"

    vidal_fname = "/Users/agoncear/data/Vidal/Human/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/Human/uniprot_mapping_unique_genes.tab"

    infile = DIR + "Alignments/matches_Hsapiens.tab.gz"
    outfile = DIR + "Interactomes/human/Hsapiens.tab"
    clones, clones_by_id = get_vidal_clones(mapping_fname)
    validation_pairs = validation_vidal(clones_by_id, vidal_fname)

    accessions = uniprot_accessions()
    print "Uniprot", len(accessions)
    savefile(infile, outfile, accessions=accessions, mapping=None, validation_pairs=validation_pairs, clones=clones)

    # infile = DIR + "Alignments/matches_Ypestis.tab"
    # outfile = DIR + "Interactomes/bacteria/ypestis.tab"
    # savefile(infile, outfile, accessions=accessions, mapping=None, validation_pairs=None)
