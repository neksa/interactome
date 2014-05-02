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
from collections import defaultdict


class Complexes:

    def __init__(self):
        self.aa_dict = {"ALA":"A", "LEU": "L", "PRO": "P", "GLY": "G", "ASP": "D", "ASN": "N", "TYR": "Y",
                    "HIS": "H", "GLU": "E", "CYS": "C", "PHE": "P", "VAL": "V", "ILE": "I", "ARG": "R",
                    "THR": "T", "LYS": "K", "SER": "S", "GLN": "Q", "MET": "M", "TRP": "W"}
        self.aa = self.aa_dict.keys()

    def templatesComplexes(self, templates, hits):
        """
        returns dictionary
        query: template/pdb
        """
        templates_complexes = defaultdict(dict)
        for pdb, (chain1, chain2) in pdb_chain_list.iteritems():

            chain1_real = chain1.replace('_', '')
            chain2_real = chain2.replace('_', '')

            # print pdb, chain1, chain2
            templates[pdb][(chain1, chain2)] = set()

            for match1 in hits[pdb+chain1]:
                # print match1
                hit_id1, hsp_identity1, hsp_gaps1, hsp_qseq1, hsp_hseq1 = match1
                template_pdb1 = hit_id1[0:4]
                template_chain1 = hit_id1[5:6]

                for match2 in hits[pdb+chain2]:
                    # print match2
                    hit_id2, hsp_identity2, hsp_gaps2, hsp_qseq2, hsp_hseq2 = match2
                    template_pdb2 = hit_id2[0:4]
                    template_chain2 = hit_id2[5:6]

                    # take only comlex templates
                    if template_pdb1 != template_pdb2: continue
                    templates[pdb][(chain1, chain2)].add((template_pdb1, template_chain1, template_chain2))
                    # print template_pdb1
        return templates_complexes


    def templatesWithHits(self, templates, hits):
        """
        This is an inverse matching: queries are assigned to templates via BLAST hits.

        Take a template complex. It has two chais (interacting partners).
        For each partner find hits.

        Normal procedure:
            Report hits
        Benchmark procedure (when starting from known structures):
            Report only hits which have the same PDB identifier

        """
        matching_templates = defaultdict(set)

        # print hits.keys()
        for pdb, (chain1, chain2, site1, site2) in templates.iteritems():

            pdb = pdb.upper()
            # template chains and sites
            # _real chains do not contain
            idx = chain1.find('_')
            chain1_real = chain1[:idx] if idx != -1 else chain1

            idx = chain2.find('_')
            chain2_real = chain2[:idx] if idx != -1 else chain2

            # print "T:", pdb, chain1_real, chain2_real

            # matching_templates[(pdb, chain1_real, chain2_real)] = set()

            # hits[hit] = query
            for query1 in hits[pdb+chain1_real]:
                q1_id, ident1, (q1_from, q1_to), (h1_from, h1_to) = query1
                # in benchmark the queries are PDB structures
                q1_pdb = q1_id[0:4]
                q1_chain = q1_id[4:5]
                # print "Q:", q1_pdb, q1_chain

                for query2 in hits[pdb+chain2_real]:
                    q2_id, ident2, (q2_from, q2_to), (h2_from, h2_to) = query2
                    # in benchmark the hits are PDB structures
                    q2_pdb = q2_id[0:4]
                    q2_chain = q2_id[4:5]
                    # print "Q:", q2_pdb, q2_chain

                    # in Benchmark we are interested in queries that are coming from the same structure
                    if q1_pdb != q2_pdb: continue

                    str_site1 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site1.iteritems()])
                    str_site2 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site2.iteritems()])

                    params = ((ident1, q1_from, q1_to, h1_from, h2_to, str_site1), (ident2, q2_from, q2_to, h2_from, h2_to, str_site2))
                    matching_templates[(pdb, chain1_real, chain2_real)].add((q1_pdb, q1_chain, q2_chain, params))
                    # print "Added template hit"
        # print matching_templates
        return matching_templates


    def loadTemplates(self, fname):
        templates = dict()
        with open(fname) as f:
            for line in f:
                pdb, chain1, chain2, str_site1, str_site2 = line.strip().split("\t")

                site1 = {}
                for s in str_site1.split(';'):
                    resn, resi, v = s.split(',')
                    site1[(resn, resi)] = v

                site2 = {}
                for s in str_site2.split(';'):
                    resn, resi, v = s.split(',')
                    site2[(resn, resi)] = v

                templates[pdb] = [chain1, chain2, site1, site2]
        return templates


    def getTemplatePDBCodes(self, templates):
        pdb_list = []
        for chains_complexes in templates.itervalues():
            for pdbs in chains_complexes.itervalues():
                for pdb, tpl_chain1, tpl_chain2 in pdbs:
                    pdb_list.append(pdb)
        return pdb_list


    def getInterface(self, fname_int):
        """
        Returns interface [(chain1,chain2)]: [{(resn1, resi1): ncontacts, ...}], [[(resn2, resi2): ncontacts, ...}]
        """
        with open(fname_int, 'r') as f:
            prev_residue = None
            contacts = defaultdict(lambda: [defaultdict(int), defaultdict(int)]) # [chain1 list, chain2 list] 

            dCA12 = 0.0
            for i, line in enumerate(f):
                if i == 0:
                    continue
                try:
                    pdb, chain1, resn1, resi1, atm1, chain2, resn2, resi2, atm2, d12, dCA12 = line.strip().split()
                except:
                    # print "Error while parsing the line: ", line.strip()
                    continue

                d12 = float(d12)
                dCA12 = float (dCA12)

                # HA contact is at most 4A between atom centers
                if d12 > 4.0: continue

                if resn1 not in self.aa:
                    # print "skip unknown ", resn1
                    continue
                if resn2 not in self.aa:
                    # print "skip unknown ", resn2
                    continue

                interacting_residues = (chain1, resn1, resi1, chain2, resn2, resi2)

                contacts[(chain1, chain2)][0][(self.aa_dict[resn1], resi1)] += 1
                contacts[(chain1, chain2)][1][(self.aa_dict[resn2], resi2)] += 1
        return contacts


    def collectTemplates(self, pdb_path, fname):
        """
        1. Load interfaces
        2. Filters:
            only accept interfaces with at least 5 interacting residues on each binding site with at most 4A between heavy atoms (HA)
            only accept interfaces that feature at least one "real" chain,
                the other chain could be generated by symmetry operators (and has underscore in its name) 
        3. Count the number of HA contacts for each residue
        4. Save interface in fname for each pdb and two chains. List the interface residues and the number of contacts each of them forms
        """
        try:
            with open(fname, 'w') as o:
                for root, dirnames, filenames in os.walk(pdb_path):
                    for filename in fnmatch.filter(filenames, '*.int'):
                        pdb, _ = os.path.basename(filename).lower().split(".", 1)
                        fname_int = root + "/" + filename
                        contacts = self.getInterface(fname_int)

                        for (chain1, chain2), (site1, site2) in contacts.iteritems():
                            # the interface should contain at least 5 interacting residues on each binding site
                            if len(site1) < 5: continue
                            if len(site2) < 5: continue
                            # at least one of the chain names should not contain underscore
                            if chain1.find('_') != -1 and chain2.find('_') != -1: continue

                            str_site1 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site1.iteritems()])
                            str_site2 = ";".join(["{},{},{}".format(k[0], k[1], v) for k, v in site2.iteritems()])
                            o.write("{}\t{}\t{}\t{}\t{}\n".format(pdb, chain1, chain2, str_site1, str_site2))
        except Exception, e:
            os.remove(fname)
            raise e



