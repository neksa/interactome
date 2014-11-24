from itertools import islice
from collections import defaultdict, namedtuple

from benchmark import load_interactions, load_gene_protein_mapping


def get_PINA_proteins(fname):
    pairs = set()
    with open(fname) as f:
        for line in f:
            interaction = line.strip().split()
            if len(interaction) != 3:
                continue
            A, _, B = interaction
            pair = A, B
            if A > B:
                pair = B, A
            pairs.add(pair)
    return pairs


# def template_analysis():
#     complex_types = {}
#     with open("/Users/agoncear/projects/Interactome/Workflow/Structures/template_analysis.tab") as f:
#         for line in islice(f, 1, None):
#             # template\tpdb\tA\tB\tprot_A\tprot_B\tcomplex_type
#             template, pdb, A, B, prot_A, prot_B, complex_type, nA, cA, nB, cB = line.strip().split("\t")
#             complex_types[template] = complex_type
#     return complex_types


def pdb_mapping(fname):
    mmcif_pdb_chain = {}
    uniprot_pdbs = defaultdict(set)
    with open(fname, 'r') as f:
        for line in islice(f, 1, None):
            pdb, chain, chain_author, uniprot, begin, end = line.strip().split("\t")
            # pdb_chain = pdb.upper() + '|' + chain
            mmcif_pdb_chain[(pdb.lower(), chain)] = chain_author
            if len(uniprot) == 0:
                continue
            uniprot = uniprot.split("-")[0]
            uniprot = uniprot.split(".")[0]
            uniprot_pdbs[uniprot].add(pdb)
    return uniprot_pdbs, mmcif_pdb_chain


def main():
    root = "/Users/agoncear/projects/Interactome/Workflow"
    # pdb_path = root + "/Interfaces/"
    struct_path = root + "/Structures"
    # pdb_templates_fname = struct_path + "/pdb_templates_5A.tab"
    pdb_proteins_fname = struct_path + "/pdb_proteins.tab"
    # template_analysis_fname = struct_path + "/template_analysis.tab"
    fname_PINA = "/Users/agoncear/data/PINA/Homo sapiens-20140521.sif"

    uniprot_pdbs, mmcif_pdb_chain = pdb_mapping(pdb_proteins_fname)
    PINA = get_PINA_proteins(fname_PINA)

    vidal_fname = "/Users/agoncear/data/Vidal/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"
    gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    records, old_new_ids = load_uniprot_mapping()
    vidal = load_interactions(vidal_fname)
    vidal_uniprot_interactions = set()
    for gA, gB in vidal:
        pA = gene_proteins[gA]
        pB = gene_proteins[gB]
        for p in pA:
            for q in pB:
                pair = p, q
                if p > q:
                    pair = q, p
                vidal_uniprot_interactions.add(pair)

    M = namedtuple('M',
            "queryA queryB template query_type template_type " +
            "n1 m1 n2 m2 " +
            "score_template_full score_template score " +
            "scaled_score_template_full scaled_score_template scaled_score " +
            "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +
            "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +
            "siteA siteB")

    # HPA_results_fname = "/Users/agoncear/projects/Interactome/Workflow/Coexpression/human_stringent_HPA.tab"
    # aln_stringent = "/Users/agoncear/projects/Interactome/Workflow/Alignments/matches_human_05_20_stringent.tab"

    aln_fname = "/Users/agoncear/projects/Interactome/Workflow/Alignments/matches_SPN_TIGR4_PPI.tab"
    network = "/Users/agoncear/projects/Interactome/Workflow/Interactomes/SPN_TIGR4_PPI.tab"

    # aln_fname = "/Users/agoncear/projects/Interactome/Workflow/Alignments/matches_HPY_26695_PPI.tab"
    # network = "/Users/agoncear/projects/Interactome/Workflow/Interactomes/HPY_26695_PPI.tab"

    interactions = defaultdict(set)
    interaction_sites = defaultdict(lambda: [set(), set()])
    with open(aln_fname) as f, open(network, 'w') as o_net:  #, open(aln_stringent, 'w') as o_aln:
        for i, line in enumerate(f):
            # if i == 0:
            #     o_aln.write(line)
            #     continue
            fields = line.strip().split("\t")
            if len(fields) != 39:
                continue
            I = M._make(fields)

            """
            ### STRINGENT FILTER ####
            seq_id = min(float(I.identicalA) / float(I.aln_lenA), float(I.identicalB) / float(I.aln_lenB))
            bs_seq_id = min(float(I.bs_identicalA) / float(I.bs_alignedA), float(I.bs_identicalB) / float(I.bs_alignedB))
            bs_seq_aligned = min(float(I.bs_alignedA)/float(I.bs_lenA), float(I.bs_alignedB)/float(I.bs_lenB))

            # print seq_id, bs_seq_id, bs_seq_aligned, I.bs_BLOSUMA, I.bs_BLOSUMB, I.bs_score1A, I.bs_score1B, I.template_type, I.query_type

            if seq_id < 0.5:
                continue
            if bs_seq_id < 0.5:
                continue
            if bs_seq_aligned < 1.0:
                continue

            if I.template_type == 'Homo' and I.query_type != 'Homo':
                continue
            if I.template_type == 'Hetereo' and I.query_type != 'Hetero':
                continue

            if I.bs_BLOSUMA < 0 or I.bs_BLOSUMB < 0:
                continue
            if I.bs_score1A < 0 or I.bs_score1B < 0:
                continue
            """
            #########################

            siteA = set()
            for site in I.siteA.split(';'):
                s = site.split(',')
                siteA.add(s[3] + s[4])

            siteB = set()
            for site in I.siteB.split(';'):
                s = site.split(',')
                siteB.add(s[3] + s[4])

            #########################

            pair = I.queryA, I.queryB
            if I.queryA > I.queryB:
                pair = I.queryB, I.queryA

            pdb, A, B = I.template.split("|")
            A = mmcif_pdb_chain.get((pdb.lower(), A), A)
            B = mmcif_pdb_chain.get((pdb.lower(), B), B)
            tpl = pdb + A + B
            print pair, I.template, tpl
            interactions[pair].add(tpl)
            interaction_sites[pair][0] |= siteA
            interaction_sites[pair][1] |= siteB
            # o_aln.write(line)

        # HPA_results = load_HPA_filter_results(HPA_results_fname, strict=True)

        o_net.write("A\tB\tTemplates\tinPINA\tinVidal\tHPA_supported\tPVH\tinPDB\tPDB\tsiteA\tsiteB\n")
        for pair, tpl in interactions.iteritems():
            inPINA = 1 if pair in PINA else 0

            # HPA_supported = 1 if HPA_results.get(pair, False) else 0
            HPA_supported = 0

            A = uniprot_pdbs[pair[0]]
            B = uniprot_pdbs[pair[1]]
            inPDB = 1 if A and B and A & B and pair[0] != pair[1] else 0
            PDB = ''
            if inPDB == 1:
                PDB = ";".join(A & B)

            inVidal = 1 if pair in vidal_uniprot_interactions else 0

            PVH = ""
            if inPINA:
                PVH += "P"
            if inVidal:
                PVH += "V"
            if HPA_supported:
                PVH += 'H'      

            siteA, siteB = interaction_sites[pair]
            siteA = ';'.join(siteA)
            siteB = ';'.join(siteB)

            o_net.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(pair[0], pair[1], ";".join(tpl), inPINA, inVidal, HPA_supported, PVH, inPDB, PDB, siteA, siteB))


UniprotRec = namedtuple('UniprotRec', 'id accessions gene descr')
def iter_uniprot():
    fname = "/Users/agoncear/data/Uniprot/uniprot_sprot.dat"
    i = ""
    rec = defaultdict(str)
    with open(fname) as f:
        for l in f:
            if l.startswith("//"):
                if rec['OX'].strip() == "NCBI_TaxID=9606;":
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
                    print uniprot_rec
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


def load_uniprot_mapping():
    fname = "/Users/agoncear/projects/Interactome/Workflow/Mappings/uniprot_description.tab"
    old_new_ids = {}
    records = {}
    try:
        with open(fname) as f:
            pass
    except:
        with open(fname, 'w') as o:
            o.write("id\taccessions\tgene\tdescr\n")
            for rec in iter_uniprot():
                o.write("{}\t{}\t{}\t{}\n".format(rec.id, ";".join(rec.accessions), rec.gene, rec.descr))
    with open(fname) as f:
        for line in islice(f, 1, None):
            id_str, accessions, gene, descr = line.strip().split("\t")
            accessions = accessions.split(";")
            rec = UniprotRec._make((id_str, accessions, gene, descr))
            if len(accessions) > 1:
                new_id = accessions[0]
                for ac in islice(accessions, 1, None):
                    old_new_ids[ac] = new_id
            records[accessions[0]] = rec
    return records, old_new_ids


def get_clone_gene_names(fname):
    import csv
    genes = set()
    with open(fname, 'rU') as f:
        clonereader = csv.reader(f, delimiter=',', quotechar='"')
        for row in clonereader:
            symbol, aliases, llid = row[8:11]
            genes.add(symbol)
            for name in aliases.split("||"):
                if 1 < len(name) < 10:
                    genes.add(name)
    return genes


def descr():
    node_props_fname = "/Users/agoncear/projects/Interactome/Workflow/Interactomes/uniprot_nodes.tab"
    vidal_fname = "/Users/agoncear/data/Vidal/ALL.tsv"
    mapping_fname = "/Users/agoncear/data/Vidal/unique_gene_ids_mapped.tab"
    clones_fname = "/Users/agoncear/data/Uetz/human cDNA collection pDonr223 Vidal-Koegl.csv"

    gene_proteins, protein_genes = load_gene_protein_mapping(mapping_fname)
    records, old_new_ids = load_uniprot_mapping()
    vidal = load_interactions(vidal_fname)
    vidal_proteins = set()
    for gA, gB in vidal:
        pA = gene_proteins[gA]
        pB = gene_proteins[gB]
        # if len(pA) > 1:
        #     print gA, pA
        for p in pA | pB:
            new_id = old_new_ids.get(p, p)  # use id as is if no mapping
            # if new_id != p:
            #     print "ID"
            # else:
            #     print "OK"
            vidal_proteins.add(new_id)

    clones = get_clone_gene_names(clones_fname)

    with open(node_props_fname, 'w') as o:
        o.write("uniprot\tgene\tdescr\tinVidalPPI\tinClones\n")
        for ac, rec in records.iteritems():
            inVidalPPI = 1 if ac in vidal_proteins else 0
            # if rec.gene == 'RASA':
            #     print ac, rec.gene, rec.gene in clones
            inClones = 1 if rec.gene in clones else 0
            for ac in rec.accessions:
                o.write("{}\t{}\t{}\t{}\t{}\n".format(ac, rec.gene, rec.descr, inVidalPPI, inClones))


if __name__ == '__main__':
    main()
    # descr()
