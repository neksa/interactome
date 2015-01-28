"""
Protein-complex homology modeling pipeline
Using PDB, Delta-BLAST, and MODELLER

Benchmarking the quality of modeling by (re)modeling the human dimers
"""

# import signal
# import time
# import math
# import pprint
import copy
from itertools import islice
from collections import defaultdict, namedtuple

from interactome.sequences.blast import BLASTReport
from interactome.structures.complexes import Complexes  # , SiteResidue
from interactome.sequences.BLOSUM import read_BLOSUM_matrix
from interactome.workflows.temp2 import pdb_proteins


def convert_alignments(in_fname, out_fname):
    pdb_path = "/Users/agoncear/projects/Interactome/Workflow/Interfaces/"

    print "Loading PDB-Uniprot mapping..."
    pdb_proteins_fname = "/Users/agoncear/projects/Interactome/Workflow/Structures/pdb_proteins.tab"
    pdb_uniprot, chain_author_to_mmcif = pdb_proteins(pdb_path, pdb_proteins_fname)
    # print chain_author_to_mmcif
    # sys.exit(0)

    with open(in_fname, 'r') as f, open(out_fname, 'w') as o:
        # for line in islice(f, 1, None):
        #query   hit     score   bit_score       evalue  identical       positive        gaps    align_len       qfrom   qto     hfrom   hto     qseq    hseq
        for i, line in enumerate(f):
            if i == 0:
                o.write(line)
                continue
            fields = line.split("\t")
            hit = fields[1]
            pdb, chain_author = hit.split("_")
            chain_mmcif = chain_author_to_mmcif.get((pdb, chain_author), None)
            # print chain_author, chain_mmcif
            if chain_mmcif is None:
                print "PDB chain not found in Biounit PDB-MMCIF mapping. Skipping BLAST hit", pdb, chain_author
                continue
                # chain_mmcif = chain_author
            fields = list(fields)
            fields[1] = pdb + "|" + chain_mmcif
            o.write("\t".join(fields))


def align(pdbchain, qseq, hseq, qfrom, hfrom, sites):
    """
    NOTE: there is an alignment problem between SEQRES in PDB and gi SEQRES in NCBI used in BLAST/Delta-BLAST. Thus the residues may be off.
          For safety, we can completely exclude these problematic alignments. (Alignment contains @ character)
          This should be fixed by NOT using NCBI's BLAST database, but building it ourselves.

    according to alignment qseq/hseq what is the substitution of the interface residue int_resn/int_resi
    is it covered by the alignment at all? if not => '*'
    is there problem with seqres-template alignment => '@'
    if the bindig site residue is covered but there is a deletion => '-'
    if it's covered and there is substitution, what's the residue name? => 'R'

    # get the sequences and align the residues:
    # i = int(int_resi) - 1
    # q = coordinate in query, resi, starting from 1
    # h = coordinate in hit, resi, starting from 1
    # qa = amino acid in query
    # ha = amino acid in hit
    # i = is the common index in qseq and hseq, since they are aligned
    """
    seqresi = [s.seqresi for s in sites]
    aln = ['*'] * len(seqresi)  # by default - not covered
    qresi = [0] * len(seqresi)  # by default - zeros in numbers of aligned residues
    q = qfrom
    h = hfrom
    for i, qa in enumerate(qseq):
        ha = hseq[i]
        try:
            site_index = seqresi.index(h)
            if ha != '-' and ha != sites[site_index].resn:
                print "the interface residue in template hit (SEQRES)", pdbchain, h, ha, \
                      "does not match what's in the structure (ATOMRES):", sites[site_index].seqresi, sites[site_index].resi, sites[site_index].resn
                aln[site_index] = '@'
            else:
                aln[site_index] = qa
                qresi[site_index] = q
        except:
            pass
        if qa != '-':
            q += 1
        if ha != '-':
            h += 1
    return zip(aln, qresi)


def binding_site_alignments(fname_interfaces, fname_blast, fname_blast_processed):
    c = Complexes()
    tpls = defaultdict(list)
    for tpl in c.loadTemplates(fname_interfaces):
        (pdb, chain1, chain2), (site1, site2) = tpl
        site1.sort(key=lambda x: x.seqresi)
        site2.sort(key=lambda x: x.seqresi)
        chain1_real = chain1.split("_")[0]
        chain2_real = chain2.split("_")[0]
        tpls[pdb.upper() + '|' + chain1_real].append((chain1, site1))
        tpls[pdb.upper() + '|' + chain2_real].append((chain2, site2))
    print "Loaded binding sites from templates"

    c = 0
    with open(fname_blast) as f, open(fname_blast_processed, 'w') as o:
        o.write("query\ttemplate\tidentical\tpositive\taln_len\tsite\n")
        for line in islice(f, 1, None):  # None
            fields = list(line.strip().split("\t"))
            c += 1
            if c % 1000 == 0:
                print '.',
            # if c % 1000000 == 0: break
            # pdb, pdb_chain = fields[1].split("|")

            # fields[1] = fields[1].replace('|', '') # remove bar from hit name (separating PDB ID and chain)
            fields[2:5] = map(float, fields[2:5])  # score, bit_score, evalue

            # print fields[5:13]
            fields[5:13] = map(int, fields[5:13])  # identical, positive, gaps, align_len, q_from, q_to, h_from, h_to
            # print fields[5:13]

            identity = int(round(float(fields[5])*100.0 / float(fields[8]), 0))  # percent identity = identical / align_len
            fields.append(identity)
            report = BLASTReport._make(fields)
            # print report

            # THRESHOLDS:
            # if not (15 <= identity):
            #     continue
            # if fields[8] < 25:
            #     continue  # report.align_len
            # if fields[4] > 0.001:
            #     continue  # report.evalue

            sites = tpls[report.hit]
            if len(sites) == 0:
                print "bs_alignments: No templates for" + report.hit
                continue

            for chain, site in sites:
                alignment = align(report.hit, report.qseq, report.hseq, report.q_from, report.h_from, site)
                aligned_bs = []
                overlap = 0
                for i, (aln, resi) in enumerate(alignment):
                    if aln != '*' and aln != '@':
                        overlap += 1
                    s = site[i]
                    aligned_bs.append("{},{},{},{},{}".format(s.resn, s.seqresi, s.ncontacts, aln, resi))
                str_site = ";".join(aligned_bs)

                # THRESHOLD: at least 3 aligned residues on the binding site
                if overlap < 3:
                    continue

                pdb = report.hit.split('|')[0]
                o.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        report.query, pdb+'|'+chain, report.identical, report.positive, report.align_len, str_site))


BS_PROPS = namedtuple('BS_PROPS', 'bs_len bs_covered bs_aligned bs_identical bs_positive bs_contacts bs_BLOSUM bs_score1')
def bs_properties(blosum, site):
    bs_len = len(site)
    bs_err = 0  # @
    bs_covered = 0  # L,-
    bs_aligned = 0  # L
    bs_identical = 0
    bs_contacts = 0
    bs_positive = 0
    bs_score1 = 0.0
    bs_BLOSUM = 0.0

    for resn1, resi, ncont, resn2, resi2 in site:
        ncont = int(ncont)
        bs_contacts += ncont
        blos = 0.0
        if resn2 == "*":  # not aligned
            pass
            blos = -4.0
        elif resn2 == "@":  # error
            bs_err += 1
            blos = -4.0
        elif resn2 == "-":  # gap
            bs_covered += 1
            blos = -4.0
        else:
            bs_covered += 1
            bs_aligned += 1
            try:
                blos = blosum[resn1][resn2]
            except:
                blos = 0.0
            if blos > 0:
                bs_positive += 1
        if resn1 == resn2:
            bs_identical += 1
        bs_BLOSUM += blos
        bs_score1 += ncont * blos
    return BS_PROPS._make((bs_len, bs_covered, bs_aligned, bs_identical, bs_positive, bs_contacts, bs_BLOSUM, bs_score1))


def merge_bs_alignments(fname_in, fname_out):
    # complex_types = {}
    # with open("/Users/agoncear/projects/Interactome/Workflow/Structures/template_analysis.tab", 'r') as f:
    #     for line in islice(f, 1, None):
    #         # template\tpdb\tA\tB\tprot_A\tprot_B\tcomplex_type
    #         template, pdb, A, B, prot_A, prot_B, complex_type = line.strip().split("\t")
    #         complex_types[template] = complex_type

    # Important: Assume that fname_in is sorted by query protein id. We can optimize a stream processing then

    blosum = read_BLOSUM_matrix("/Users/agoncear/projects/Interactome/interactome/sequences/BLOSUM62")

    with open(fname_in, 'r') as f, open(fname_out, 'w') as o:
        o.write("query\ttemplate\tidentical\tpositive\taln_len\tbs_len\tbs_covered\tbs_aligned\tbs_identical\tbs_positive\tbs_contacts\tbs_BLOSUM\tbs_score1\tsite\n")

        def sub_protein_generator(handle):
            prevquery = None
            block = []
            for c, line in enumerate(islice(handle, 1, None)):
                fields = line.strip().split("\t")
                if c % 100000 == 0:
                    print c
                if prevquery is not None and fields[0] != prevquery:
                    yield block
                    block = []
                block.append(fields)
                prevquery = fields[0]
            if len(block) > 0:
                yield block

        for block in sub_protein_generator(f):
            merged_bs_aln = defaultdict(list)
            for fields in block:
                # 0       1         2 (0)     3 (1)   4 (2)   5 (3)
                query, template, identical, positive, aln_len, site = fields
                identical = int(identical)
                positive = int(positive)
                aln_len = int(aln_len)
                residues = []
                # print query, template,
                merged = False
                for s in site.split(";"):
                    residues.append(s.split(","))  # 0:X, 1:X', 2:N, 3:Y, 4:Y'
                # print "RRRR", residues
                for i, alignment in enumerate(merged_bs_aln[(query, template)]):
                    residue_overlaps = 0
                    different_bs = False
                    # new_aln = copy.copy(alignment)
                    new_aln = alignment

                    if len(alignment[3]) != len(residues):
                        # it's a different binding site - don't merge
                        different_bs = True
                        break

                    for j, r in enumerate(residues):
                        # compare site residues
                        # print alignment[3], j
                        # print alignment[3][j]
                        if alignment[3][j][0] != r[0] or alignment[3][j][1] != r[1]:
                            # it's a different binding site - don't merge
                            # print "diff bs",
                            different_bs = True
                            break

                        aln_char = alignment[3][j][3]  # alignment[3] is site                    
                        this_char = r[3]
                        this_resi = r[4]
                        if aln_char != '*' and this_char != '*':
                            residue_overlaps += 1
                        elif aln_char == '*' and this_char != '*' and this_char != '@':
                            new_aln[3][j][3] = this_char
                            new_aln[3][j][4] = this_resi

                    if residue_overlaps == 0 and different_bs is False:
                        merged_bs_aln[(query, template)][i] = new_aln  # copy.copy(new_aln)
                        merged_bs_aln[(query, template)][i][0] += identical
                        merged_bs_aln[(query, template)][i][1] += positive
                        merged_bs_aln[(query, template)][i][2] += aln_len
                        merged = True
                        # print "Merged"
                        break
                if not merged:
                    # print "Added"
                    merged_bs_aln[(query, template)].append([identical, positive, aln_len, residues])

            # pp = pprint.PrettyPrinter(depth=6)
            for (query, template), bs_alignments in merged_bs_aln.iteritems():
                for bs_alignment in bs_alignments:
                    # pp.pprint(tuple(bs_alignment))
                    identical, positive, aln_len, merged_site = tuple(bs_alignment)

                    bs_len, bs_covered, bs_aligned, bs_identical, bs_positive, bs_contacts, bs_BLOSUM, bs_score1 = bs_properties(blosum, merged_site)

                    m_site = ["{},{},{},{},{}".format(*s) for s in merged_site]
                    merged_site_str = ";".join(m_site)

                    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            query, template, identical, positive, aln_len,
                            bs_len, bs_covered, bs_aligned, bs_identical, bs_positive,
                            bs_contacts, bs_BLOSUM, bs_score1,
                            merged_site_str))


if __name__ == '__main__':
    pass
    # dirname = "/Users/agoncear/projects/Interactome/Workflow"
    # binding_site_alignments(dirname + "/Structures/pdb_templates_5A.tab",
    #                         dirname + "/BLAST-results/human_deltablast_report2.tab",
    #                         dirname + "/BLAST-results/human_deltablast_report_merged.tab")
    # merge_bs_alignments(dirname + "/BLAST-results/human_deltablast_report_merged.tab",
    #                     dirname + "/BLAST-results/human_deltablast_report_merged-bs.tab")
