"""
    Structural complexes and protein-protein interactions
    Interactions are calculated separately by scoring potential routines (following IBIS logic)
    Calculated interactions are spread over in *.int files
    *.int files are collected by collectTemplates() into "pdb: (A,B) (A,C)" format
    Then, from BLAST hits and pairs of interacting chains templatesComplexes() identifies templates with interactions (i.e. complexes)

    Routine alignInterfaces() is responsible for alignment of interfaces, but not for scoring
"""
import fnmatch
import os
from collections import defaultdict, namedtuple
# from sequences.blast import BLASTReport

SiteResidue = namedtuple('SiteResidue', 'resn, resi, seqresi, ncontacts')


class Complexes:

    def __init__(self):
        self.aa_dict = {"ALA": "A", "LEU": "L", "PRO": "P", "GLY": "G", "ASP": "D", "ASN": "N", "TYR": "Y",
                        "HIS": "H", "GLU": "E", "CYS": "C", "PHE": "F", "VAL": "V", "ILE": "I", "ARG": "R",
                        "THR": "T", "LYS": "K", "SER": "S", "GLN": "Q", "MET": "M", "TRP": "W"}
        self.aa = self.aa_dict.keys()
        """
        ASX (B) asparagine or aspartic acid
        GLX (Z) glutamine or glutamic acid
        CSO (C) S-hydroxycysteine
        HIP (H) ND1-phosphohistidine
        HSD (H) prototropic tautomer of histidine, H on ND1 (CHARMM)
        HSE (H) prototropic tautomer of histidine, H on NE2 (CHARMM)
        HSP (H) protonated histidine
        MSE selenomethionine
        SEC (U) selenocysteine
        SEP (S) phosphoserine
        TPO (T) phosphothreonine
        PTR (Y) O-phosphotyrosine
        XLE (J) leucine or isoleucine
        XAA (X) unspecified or unknown
        """

    def align(self, pdb, chain, qseq, hseq, qfrom, hfrom, sites):
        """
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
        q = qfrom
        h = hfrom
        for i, qa in enumerate(qseq):
            ha = hseq[i]
            try:
                site_index = seqresi.index(h)
                if ha != '-' and ha != sites[site_index].resn:
                    print("the interface residue in template hit (SEQRES)", pdb, chain, h, ha,
                          "does not match what's in the structure (ATOMRES):", sites[site_index].seqresi, sites[site_index].resi, sites[site_index].resn)
                    aln[site_index] = '@'
                else:
                    aln[site_index] = qa
            except:
                pass
            if qa != '-':
                q += 1
            if ha != '-':
                h += 1
        return aln

    # def templatesWithHits(self, template, hits, benchmark=False):
    #     """
    #     This is an inverse matching: queries are assigned to templates via BLAST hits.

    #     Take a template complex. It has two chais (interacting partners).
    #     For each partner find hits.

    #     Normal procedure:
    #         Report hits
    #     Benchmark procedure (when starting from known structures):
    #         ?? Report only hits which have the same PDB identifier

    #     """
    #     min_bs_residue_covereage = 3  # 3 residues should be covered by binding site alignment
    #     min_bs_fraction_coverage = 0.7  # At least half of the interface residues should be covered by the alignments

    #     (pdb, chain1, chain2), (site1, site2) = template
    #     site1.sort(key=lambda x: x.seqresi)
    #     site2.sort(key=lambda x: x.seqresi)
    #     pdb = pdb.upper()
    #     chain1_real = chain1.split("_")[0]
    #     chain2_real = chain2.split("_")[0]

    #     # print pdb, chain1_real, chain2_real
    #     print "T:", pdb, chain1, chain2, "real:", chain1_real, chain2_real
    #     queries = []
    #     # q1_pdb = ""
    #     # q1_chain = ""
    #     # q2_pdb = ""
    #     # q2_chain = ""

    #     # Site A
    #     # hits[hit][query] = [a,a,a]
    #     for query1, q1_records in hits.get(pdb+'|'+chain1_real, {}).iteritems():
    #         # print "A:", query1
    #         alignments = []
    #         # the records are sorted by #positive substitutions descending (see blast)
    #         for q1 in q1_records:
    #             # print "id A:", q1.identity, q1.positive
    #             alignments.append(self.align(pdb, chain1, q1.qseq, q1.hseq, q1.q_from, q1.h_from, site1))
    #         # print "---- ALN A: ----"
    #         # print "".join([s.resn for s in site1])
    #         # print "-------------"
    #         # for a in alignments: print "".join(a)

    #         # now all the query-hit matches have been merged for the same site, calculate the coverage
    #         overlap = 0
    #         tmp_sites = []
    #         # consensus = []
    #         for i, s in enumerate(site1):
    #             aln = None
    #             for j, alignment in enumerate(alignments):
    #                 if j == 0:
    #                     aln = alignment[i]
    #                 # allow chimeric allignments for now:
    #                 elif (aln == '*' or aln == '@') and alignment[i] != '*' and alignment[i] != '@':
    #                     aln = alignment[i]
    #             if aln != '*' and aln != '@':
    #                 overlap += 1
    #             tmp_sites.append("{},{},{},{}".format(s.resn, s.seqresi, s.ncontacts, aln))
    #             # consensus.append(aln)
    #         # print "-------------"
    #         # print "".join(consensus)
    #         # print "overlap", overlap, "out of", len(tmp_sites)
    #         if overlap < min_bs_residue_covereage: continue
    #         if overlap / float(len(tmp_sites)) < min_bs_fraction_coverage: continue
    #         str_site1 = ";".join(tmp_sites)

    #         # Site B ###########################################################################################
    #         # if pdb+'|'+chain2_real not in hits: continue
    #         for query2, q2_records in hits.get(pdb+'|'+chain2_real, {}).iteritems():
    #             # print "B:", query2
    #             alignments = []
    #             # the records are sorted by #positive substitutions descending (see blast)
    #             for q2 in q2_records:
    #                 # print "id B:", q2.identity
    #                 alignments.append(self.align(pdb, chain2, q2.qseq, q2.hseq, q2.q_from, q2.h_from, site2))
    #                 # print alignments
    #             # now all the query-hit matches have been merged for the same site, calculate the coverage

    #             # print "---- ALN B: ----"
    #             # print [s.resn for s in site1]
    #             # for a in alignments: print a
    #             overlap = 0
    #             tmp_sites = []
    #             # consensus = []
    #             for i, s in enumerate(site2):
    #                 aln = None
    #                 for j, alignment in enumerate(alignments):
    #                     if j == 0:
    #                         aln = alignment[i]
    #                     # allow chimeric allignments for now:
    #                     elif (aln == '*' or aln == '@') and alignment[i] != '*' and alignment[i] != '@':
    #                         aln = alignment[i]
    #                 if aln != '*' and aln != '@':
    #                     overlap += 1
    #                 tmp_sites.append("{},{},{},{}".format(s.resn, s.seqresi, s.ncontacts, aln))
    #                 # consensus.append(aln)

    #             # print "-------------"
    #             # print consensus
    #             # print "overlap", overlap, "out of", len(tmp_sites)
    #             if overlap < min_bs_residue_covereage: continue
    #             if overlap / float(len(tmp_sites)) < min_bs_fraction_coverage: continue
    #             str_site2 = ";".join(tmp_sites)

    #             params = ((q1_records, str_site1), (q2_records, str_site2))
    #             if benchmark is True:
    #                 # This only makes sense for benchmark, in real life queries do not have pdb or chain ids
    #                 # matching_templates[(pdb, chain1_real, chain2_real)].add((q1_pdb, q1_chain, q2_chain, params))
    #                 # queries.add((q1_pdb+q1_chain, q1_pdb+q2_chain, params))
    #                 pass
    #             else:
    #                 queries.append((query1, query2, params))
    #     return (pdb, chain1, chain2), queries

    def loadTemplates(self, fname_interfaces):
        """
        SEQRES - ATOM mapping is taken from SIFTS for instance
        """
        with open(fname_interfaces) as f:
            for line in f:
                pdb, chain1, chain2, str_site1, str_site2 = line.strip().split("\t")
                site1 = []
                for s in str_site1.split(';'):
                    resn, resi, v = s.split(',')
                    if len(resn) > 1:
                        resn = self.aa_dict.get(resn)
                        if resn is None:
                            continue
                    resi = int(resi)
                    v = int(v)
                    # seqresi  = resi - m1.atom_from + m1.seqres_from
                    seqresi = resi
                    site1.append(SiteResidue(resn=resn, resi=resi, seqresi=seqresi, ncontacts=v))

                site2 = []
                for s in str_site2.split(';'):
                    resn, resi, v = s.split(',')
                    if len(resn) > 1:
                        resn = self.aa_dict.get(resn)
                        if resn is None:
                            continue
                    resi = int(resi)
                    v = int(v)
                    # seqresi  = resi - m2.atom_from + m2.seqres_from
                    seqresi = resi
                    site2.append(SiteResidue(resn=resn, resi=resi, seqresi=seqresi, ncontacts=v))
                yield (pdb, chain1, chain2), (site1, site2)

    def getTemplatePDBCodes(self, templates):
        pdb_list = []
        for chains_complexes in templates.itervalues():
            for pdbs in chains_complexes.itervalues():
                for pdb, tpl_chain1, tpl_chain2 in pdbs:
                    pdb_list.append(pdb)
        return pdb_list

    def getInterface(self, fname_int, distance_threshold):
        """
        Returns interface [(chain1,chain2)]: [{(resn1, resi1): ncontacts, ...}], [[(resn2, resi2): ncontacts, ...}]
        """
        with open(fname_int, 'r') as f:
            # prev_residue = None
            contacts = defaultdict(lambda: [defaultdict(int), defaultdict(int)])  # [chain1 list, chain2 list]
            # dCA12 = 0.0
            for i, line in enumerate(f):
                if i == 0:
                    continue
                try:
                    chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12 = line.strip().split()
                except:
                    print("Error while parsing the line: ", line.strip())
                    continue
                d12 = float(d12)
                # dCA12 = float (dCA12)
                # HA contact is at most (4A, 4.5A, 5A...) between atom centers
                if d12 > distance_threshold:
                    continue
                # if resn1 not in self.aa:
                #     # print "skip unknown ", resn1
                #     continue
                # if resn2 not in self.aa:
                #     # print "skip unknown ", resn2
                #     continue
                # contacts[(chain1, chain2)][0][(self.aa_dict[resn1], resi1)] += 1
                # contacts[(chain1, chain2)][1][(self.aa_dict[resn2], resi2)] += 1

                # interacting_residues = (chain1, resn1, resi1, chain2, resn2, resi2)
                contacts[(chain1, chain2)][0][(resn1, resi1)] += 1
                contacts[(chain1, chain2)][1][(resn2, resi2)] += 1
        return contacts

    def collectTemplates(self, pdb_path, fname, min_number_of_contacts=5, distance_threshold=5.0):
        """
        1. Load interfaces
        2. Filters:
            only accept interfaces with at least 5 interacting residues on each binding site with at most 5A between heavy atoms (HA)
            only accept interfaces that feature at least one "real" chain,
                the other chain could be generated by symmetry operators (and has underscore in its name)
        3. Count the number of HA contacts for each residue
        4. Save interface in fname for each pdb and two chains. List the interface residues and the number of contacts each of them forms
        """
        cnt = 0
        try:
            with open(fname, 'w') as o:
                for root, dirnames, filenames in os.walk(pdb_path):
                    for filename in fnmatch.filter(filenames, '*_atomic_contacts_5.0A.tab'):
                        pdb, _ = os.path.basename(filename).lower().split("_", 1)
                        fname_int = root + "/" + filename
                        contacts = self.getInterface(fname_int, distance_threshold)
                        # print cnt
                        cnt += 1
                        print(cnt, pdb)
                        # if cnt > 1000: continue
                        for (chain1, chain2), (site1, site2) in contacts.iteritems():
                            # print cnt, pdb, chain1, chain2
                            # the interface should contain at least 5 interacting residues on each binding site
                            if len(site1) < min_number_of_contacts:
                                continue
                            if len(site2) < min_number_of_contacts:
                                continue
                            # at least one of the chain names should not contain underscore
                            if chain1.find('_') != -1 and chain2.find('_') != -1:
                                continue

                            # @TODO: sort by residue number here!
                            str_site1 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site1.iteritems()])
                            str_site2 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site2.iteritems()])
                            o.write("{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, str_site1, str_site2))
                            # print "{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, str_site1, str_site2)
        except Exception as e:
            os.remove(fname)
            raise e
