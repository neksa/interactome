"""
    Structural complexes and protein-protein interactions
    Interactions are calculated separately by scoring potential routines (following IBIS logic)
    Calculated interactions are spread over in *.int files 
    *.int files are collected by collectTemplates() into "pdb: (A,B) (A,C)" format
    Then, from BLAST hits and pairs of interacting chains templatesComplexes() identifies templates with interactions (i.e. complexes)

    Routine alignInterfaces() is responsible for alignment of interfaces, but not for scoring
"""


class Complexes:
    """
    returns dictionary
    query: template/pdb
    """
    def templatesComplexes(self, pdb_chain_list, hits):
        templates = defaultdict(dict)
        for pdb, (chain1, chain2) in pdb_chain_list.iteritems():
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
        return templates


    def getTemplatePDBCodes(self, templates):
        pdb_list = []
        for chains_complexes in templates.itervalues():
            for pdbs in chains_complexes.itervalues():
                for pdb, tpl_chain1, tpl_chain2 in pdbs:
                    pdb_list.append(pdb)
        return pdb_list


    def collectTemplates(self, path, fname):
        pass

