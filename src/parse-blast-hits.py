"""



"""

import xml.etree.ElementTree as ET
from collections import defaultdict
from textwrap import fill

class BLAST:
    def __init__(self, fname):

        try:
            with open(fname) as f:
                xml = f.read()
        except IOError:
            print "File not found"
            return

        self.root = ET.fromstring(xml)
        # print self.root
        self.matches = defaultdict(list)

        # ./BlastOutput/BlastOutput_iterations/
        for iteration in self.root.findall("./BlastOutput_iterations/Iteration"):
            query_id = iteration.find("Iteration_query-def").text
            query_id = query_id[0:5] # 2HLEA from 2HLEA_ATOM
            # print query_id
            # print iteration
            for hit in iteration.findall("Iteration_hits/Hit"):
                # print hit
                hit_id = hit.find("Hit_accession").text
                # print hit_id
                hit_len = int(hit.find("Hit_len").text)
                for hsp in hit.findall("Hit_hsps/Hsp"):
                    hsp_identity = int(hsp.find("Hsp_identity").text)
                    
                    if not (25 <= hsp_identity < 90): continue

                    hsp_gaps = int(hsp.find("Hsp_gaps").text)
                    hsp_align_len = int(hsp.find("Hsp_align-len").text)

                    # alignment:
                    hsp_qseq = hsp.find("Hsp_qseq").text
                    hsp_hseq = hsp.find("Hsp_hseq").text

                    self.matches[query_id].append((hit_id, hsp_identity, hsp_gaps, hsp_qseq, hsp_hseq))

        # self.tree = root.find("./Hit")


    def dimer_models(self):
        for pdb, (chain1, chain2) in self.dimers.iteritems():
            print pdb, chain1, chain2

            for match1 in self.matches[pdb+chain1]:
                # print match1
                hit_id1, hsp_identity1, hsp_gaps1, hsp_qseq1, hsp_hseq1 = match1
                template_pdb1 = hit_id1[0:4]
                template_chain1 = hit_id1[5:6]


                for match2 in self.matches[pdb+chain2]:
                    # print match2
                    hit_id2, hsp_identity2, hsp_gaps2, hsp_qseq2, hsp_hseq2 = match2
                    template_pdb2 = hit_id2[0:4]
                    template_chain2 = hit_id2[5:6]

                    
                    # take only comlex templates
                    if template_pdb1 != template_pdb2: continue

                    print "\nC; Experiment 1. Modeling based on multi-chain dimer templates"
                    # print pdb, chain1, chain2
                    # print "TEMPLATE complex", template_pdb1, template_chain1, template_chain2

                    print ">P1;"+template_pdb1+template_chain1+template_chain2
                    print "structureX:{}:{}:{}:{}:{}::::".format(template_pdb1+template_chain1+template_chain2, ".", template_chain1, ".", template_chain2)
                    # from 1:A to 100:B
                    print fill(hsp_hseq1, 60)
                    print "/"
                    print fill(hsp_hseq2, 60)
                    print "*"
                    # print

                    print ">P1;query"+pdb+chain1+chain2
                    print "sequence:query{}:{}:{}:{}:{}::::".format(pdb+chain1+chain2, ".", chain1, ".", chain2)
                    # from 1:A to 100:B
                    print fill(hsp_qseq1, 60)
                    print "/"
                    print fill(hsp_qseq2, 60)
                    print "*"
                    print

            # (hit_id, hsp_identity, hsp_gaps, hsp_qseq, hsp_hseq)

    def read_dimers(self):
        self.dimers = defaultdict(list)
        with open("../SKEMPI/dimers.csv") as f:
            for line in f:
                pdb, chain1, chain2 = line.strip().split("_", 2)
                self.dimers[pdb] = (chain1, chain2)
        # print self.dimers


if __name__ == '__main__':
    blast = BLAST("../blast-pdb/CNPENUVX01R-Alignment.xml")
    # blast = BLAST("../blast-pdb/test.xml")
    blast.read_dimers()
    blast.dimer_models()

