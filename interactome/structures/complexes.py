"""
    Structural complexes and protein-protein interactions
    Interactions are calculated separately by scoring potential routines (following IBIS logic)
    Calculated interactions are spread over in *.int files 
    *.int files are collected by collectTemplates() into "pdb: (A,B) (A,C)" format
    Then, from BLAST hits and pairs of interacting chains templatesComplexes() identifies templates with interactions (i.e. complexes)

    Routine alignInterfaces() is responsible for alignment of interfaces, but not for scoring
"""

import fnmatch
import itertools
import os
from collections import defaultdict, namedtuple
from blast import BLASTReport

SiteResidue = namedtuple('SiteResidue', 'resn, resi, seqresi, ncontacts')

class Complexes:

    def __init__(self):
        self.aa_dict = {"ALA":"A", "LEU": "L", "PRO": "P", "GLY": "G", "ASP": "D", "ASN": "N", "TYR": "Y",
                    "HIS": "H", "GLU": "E", "CYS": "C", "PHE": "F", "VAL": "V", "ILE": "I", "ARG": "R",
                    "THR": "T", "LYS": "K", "SER": "S", "GLN": "Q", "MET": "M", "TRP": "W"}
        self.aa = self.aa_dict.keys()

    # def templatesComplexes(self, templates, hits):
    #     """
    #     returns dictionary
    #     query: template/pdb
    #     """
    #     templates_complexes = defaultdict(dict)
    #     for pdb, (chain1, chain2) in pdb_chain_list.iteritems():

    #         chain1_real = chain1.replace('_', '')
    #         chain2_real = chain2.replace('_', '')

    #         # print pdb, chain1, chain2
    #         templates[pdb][(chain1, chain2)] = set()

    #         for match1 in hits[pdb+chain1]:
    #             # print match1
    #             hit_id1, hsp_identity1, hsp_gaps1, hsp_qseq1, hsp_hseq1 = match1
    #             template_pdb1 = hit_id1[0:4]
    #             template_chain1 = hit_id1[5:6]

    #             for match2 in hits[pdb+chain2]:
    #                 # print match2
    #                 hit_id2, hsp_identity2, hsp_gaps2, hsp_qseq2, hsp_hseq2 = match2
    #                 template_pdb2 = hit_id2[0:4]
    #                 template_chain2 = hit_id2[5:6]

    #                 # take only complex templates
    #                 if template_pdb1 != template_pdb2: continue
    #                 templates[pdb][(chain1, chain2)].add((template_pdb1, template_chain1, template_chain2))
    #                 # print template_pdb1
    #     return templates_complexes


    def align(self, pdb, chain, qseq, hseq, qfrom, hfrom, s):
        """
        according to alignment qseq/hseq what is the substitution of the interface residue int_resn/int_resi 
        is it covered by the alignment at all? if not => '*'
        is there problem with seqres-template alignment => '@'
        if the bindig site residue is covered but there is a deletion => '-'
        if it's covered and there is substitution, what's the residue name? => 'R'
        """
        # get the sequences and align the residues:
        # i = int(int_resi) - 1
        # q = coordinate in query, resi, starting from 1
        # h = coordinate in hit, resi, starting from 1
        # qa = amino acid in query
        # ha = amino acid in hit
        # i = is the common index in qseq and hseq, since they are aligned
        q = qfrom
        h = hfrom

        # debug = False
        # if debug:
        #     print
        #     print "DEBUG0", pdb, chain
        #     print "DEBUG1", qfrom, qseq
        #     print "DEBUG2", hfrom, hseq
        #     print "DEBUG3", s.resn, s.resi, s.seqresi

        # Some heuristics to speed the process up
        # Check if the interface residue is outside the alignment
        if not hfrom < s.seqresi < hfrom + len(hseq):
            # if debug:
            #     print "Out of scope", hfrom, s.seqresi, hfrom + len(hseq)
            #     print "DEBUG4", '*'
            return '*'

        for i, qa in enumerate(qseq):
            ha = hseq[i]
            # print "=STEP= ", i, q, h, qa, ha
            if h == s.seqresi:
                if ha != '-' and ha != s.resn:
                    # print "seqres != atomres", pdb, chain
                    print "the interface residue in template hit (SEQRES)", pdb, chain, h, ha, "does not match what's in the structure (ATOMRES):", s.resn, s.resi, s.seqresi
                    return '@'
                    # print "alignment, pos =", i
                    # print qseq
                    # print hseq
                    # print
                else:
                    # if debug: print "DEBUG4", ha, "->", qa
                    return qa

            if qa != '-': q += 1
            if ha != '-': h += 1
        # the interface residue is outside the alignment OR the query sequence SEQRES did not match the atomic sequence ATOMRES
        # if debug: print "DEBUG4", '*'
        return '*'


    def templatesWithHits(self, templates, hits, benchmark=False):
        """
        This is an inverse matching: queries are assigned to templates via BLAST hits.

        Take a template complex. It has two chais (interacting partners).
        For each partner find hits.

        Normal procedure:
            Report hits
        Benchmark procedure (when starting from known structures):
            ?? Report only hits which have the same PDB identifier

        """
        # matching_templates = defaultdict(set)

        # print hits.keys()
        # c = 0

        min_bs_alignment_covereage = 3 # 3 residues should be covered by binding site alignment

        # for (pdb, chain1, chain2), (site1, site2) in templates.iteritems():
        # print templates
        # print templates[1]
        # for (pdb, chain1, chain2), (site1, site2) in templates:
        if True:
            (pdb, chain1, chain2), (site1, site2) = templates

            pdb = pdb.upper()
            # template chains and sites
            # _real chains do not contain underscores, while generated may have A_1, A_2,... in their names
            # we remove _1, _2 from the chain names
            
            # idx = chain1.find('_')
            # chain1_real = chain1[:idx] if idx != -1 else chain1
            chain1_real = chain1.split("_")[0]

            # idx = chain2.find('_')
            # chain2_real = chain2[:idx] if idx != -1 else chain2
            chain2_real = chain2.split("_")[0]

            # print pdb, chain1_real, chain2_real
            print "T:", pdb, chain1_real, chain2_real
            queries = set()
            q1_pdb = ""
            q1_chain = ""
            q2_pdb = ""
            q2_chain = ""

            for fields1 in hits[pdb+'|'+chain1_real]:

                q1 = BLASTReport._make(fields1)

                # print q1

                if benchmark == True:
                    # ONLY in Benchmark the queries are PDB structures
                    # we do not need self-hits:
                    q1_pdb = q1.query[0:4]
                    q1_chain = q1.query[4:5]
                    # if q1_pdb == pdb: continue
                    # if the query is the same as hit, make sure the chains are the same, we don't want to predict C-D based on A-B
                    if q1_pdb == pdb and q1_chain != chain1_real: continue

                # @TODO
                # skip identical matches, where query and hit is the same struture and chain
                # DEBUG!!!! >>>>>>>>
                # if q1.query != pdb+chain1_real: continue
                # if q1.query == pdb+chain1_real: continue
                # >>>>>>>>>>>>>>>>>>

                # We are interested in at least 2 residues in each of the predicted binding sites, otherwise the interface score would be dumb
                # So we skip the cases with less than 2 non-gap interface alignments
                overlap = 0
                site = []
                for s in site1:
                    aln = self.align(pdb, chain1, q1.qseq, q1.hseq, q1.q_from, q1.h_from, s)
                    if aln != '*' and aln != '@': overlap += 1
                    site.append("{},{},{},{}".format(s.resn, s.seqresi, s.ncontacts, aln))
                if overlap < min_bs_alignment_covereage: continue
                # print "OVR1", overlap
                str_site1 = ";".join(site)

                for fields2 in hits[pdb+'|'+chain2_real]:

                    q2 = BLASTReport._make(fields2)

                    if benchmark == True:
                        q2_pdb = q2.query[0:4] #!!!!
                        q2_chain = q2.query[4:5] #!!!!
                        # if the query is the same as hit, make sure the chains are the same, we don't want to predict C-D based on A-B
                        if q2_pdb == pdb and q2_chain != chain2_real: continue

                    if benchmark == True:
                        # in Benchmark we are interested in queries that are coming from the same structure
                        if q1_pdb != q2_pdb: continue    

                    overlap = 0
                    site = []
                    for s in site2:
                        aln = self.align(pdb, chain2, q2.qseq, q2.hseq, q2.q_from, q2.h_from, s)
                        if aln != '*' and aln != '@': overlap += 1
                        site.append("{},{},{},{}".format(s.resn, s.seqresi, s.ncontacts, aln))
                    if overlap < min_bs_alignment_covereage: continue
                    # print "OVR2", overlap
                    str_site2 = ";".join(site)

                    # params = ((ident1, q1_from, q1_to, h1_from, h2_to, str_site1), (ident2, q2_from, q2_to, h2_from, h2_to, str_site2))

                    params = ((q1, str_site1), (q2, str_site2))
                    
                    if benchmark == True:
                        # This only makes sense for benchmark, in real life queries do not have pdb or chain ids
                        # matching_templates[(pdb, chain1_real, chain2_real)].add((q1_pdb, q1_chain, q2_chain, params))
                        queries.add((q1_pdb+q1_chain, q1_pdb+q2_chain, params))
                    else:
                        queries.add((q1.query, q2.query, params))
                        # print "ADD!!"
                    
                    # c += 1
                    # print c, params
                    # if c > 100: break
                    # print "Added template hit"
            # theoretically, we could yield results here for pdb, chain1_real, chain2_real
            # yield (pdb, chain1_real, chain2_real), queries
            yield (pdb, chain1, chain2), queries
        # return matching_templates


    def loadTemplates(self, fname_interfaces, mapping=None):
        """
        SEQRES - ATOM mapping is taken from SIFTS for instance
        """

        templates = dict()
        with open(fname_interfaces) as f:
            for line in f:
                pdb, chain1, chain2, str_site1, str_site2 = line.strip().split("\t")

                # m1 = mapping.get(pdb+chain1, None)
                # m2 = mapping.get(pdb+chain2, None)

                # if m1 is None or m2 is None: continue

                site1 = []
                for s in str_site1.split(';'):
                    resn, resi, v = s.split(',')
                    if len(resn) > 1:
                        resn = self.aa_dict.get(resn)
                        if resn is None: continue
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
                        if resn is None: continue
                    resi = int(resi)
                    v = int(v)
                    # seqresi  = resi - m2.atom_from + m2.seqres_from
                    seqresi  = resi
                    site2.append(SiteResidue(resn=resn, resi=resi, seqresi=seqresi, ncontacts=v))

                # templates[(pdb, chain1, chain2)] = site1, site2
                yield (pdb, chain1, chain2), (site1, site2)

        # return templates


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
            prev_residue = None
            contacts = defaultdict(lambda: [defaultdict(int), defaultdict(int)]) # [chain1 list, chain2 list] 

            # dCA12 = 0.0
            for i, line in enumerate(f):
                if i == 0:
                    continue
                try:
                    chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12 = line.strip().split()
                except:
                    print "Error while parsing the line: ", line.strip()
                    continue

                d12 = float(d12)
                # dCA12 = float (dCA12)

                # HA contact is at most (4A, 4.5A, 5A...) between atom centers
                if d12 > distance_threshold: continue

                # if resn1 not in self.aa:
                #     # print "skip unknown ", resn1
                #     continue
                # if resn2 not in self.aa:
                #     # print "skip unknown ", resn2
                #     continue

                interacting_residues = (chain1, resn1, resi1, chain2, resn2, resi2)

                # contacts[(chain1, chain2)][0][(self.aa_dict[resn1], resi1)] += 1
                # contacts[(chain1, chain2)][1][(self.aa_dict[resn2], resi2)] += 1
                contacts[(chain1, chain2)][0][(resn1, resi1)] += 1
                contacts[(chain1, chain2)][1][(resn2, resi2)] += 1
        return contacts


    def collectTemplates(self, pdb_path, fname):
        """
        1. Load interfaces
        2. Filters:
            only accept interfaces with at least 5 interacting residues on each binding site with at most 5A between heavy atoms (HA)
            only accept interfaces that feature at least one "real" chain,
                the other chain could be generated by symmetry operators (and has underscore in its name) 
        3. Count the number of HA contacts for each residue
        4. Save interface in fname for each pdb and two chains. List the interface residues and the number of contacts each of them forms
        """
        min_number_of_contacts = 5
        distance_threshold = 5.0

        cnt = 0
        try:
            with open(fname, 'w') as o:
                for root, dirnames, filenames in os.walk(pdb_path):
                    for filename in fnmatch.filter(filenames, '*_atomic_contacts_5.0A.tab'):
                        pdb, _ = os.path.basename(filename).lower().split("_", 1)
                        fname_int = root + "/" + filename
                        contacts = self.getInterface(fname_int, distance_threshold)

                        # print cnt
                        cnt +=1
                        print cnt, pdb
                        # if cnt > 1000: continue

                        for (chain1, chain2), (site1, site2) in contacts.iteritems():
                            # print cnt, pdb, chain1, chain2
                            # the interface should contain at least 5 interacting residues on each binding site
                            if len(site1) < min_number_of_contacts: continue
                            if len(site2) < min_number_of_contacts: continue
                            # at least one of the chain names should not contain underscore
                            if chain1.find('_') != -1 and chain2.find('_') != -1: continue

                            # @TODO: sort by residue number here!
                            str_site1 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site1.iteritems()])
                            str_site2 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site2.iteritems()])
                            o.write("{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, str_site1, str_site2))
                            # print "{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, str_site1, str_site2)
        except Exception, e:
            os.remove(fname)
            raise e



