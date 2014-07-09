"""
Test residue mapping

"""

from pdb import SIFTS
from complexes import Complexes
from blast import BLAST

blast_report = "deltablast_report.tab"
pdb_templates = "pdb_templates.tab"
matches_fname = "matches.tab"
pdb_interfaces_path = "../scoring/results/"
mapping_fname = "../../../data/EBI/SIFTS/pdb_chain_uniprot.tsv"

sifts = SIFTS()
mapping = sifts.getPDBMapping(mapping_fname)

blast = BLAST()
by_hit, by_query = blast.readReport(blast_report)

complexes = Complexes()
templates = complexes.loadTemplates(pdb_templates, mapping)

print "Info: Number of templates (PDB+chain+chain) = ", len(templates.keys())

pdb_codes = set()
for (pdb,chain1,chain2), (site1,site2) in templates.iteritems():
    pdb_codes.add(pdb)

#     idx = chain1.find('_')
#     chain1_real = chain1[:idx] if idx != -1 else chain1

#     idx = chain2.find('_')
#     chain2_real = chain2[:idx] if idx != -1 else chain2
#     print pdb, chain1, chain2, chain1_real, chain2_real

print "PDB codes:", len(pdb_codes)
# templates[]

t = templates[('2j5q','A', 'B')]
# s = "MTKVSVVGAAGTVGAAAGYNIALRDIADEVVFVDIPDKEDDTVGQAADTNHGIAYDSNTRVRQGGYEDTAGSDVVVITAGIPRQPGQTRIDLAGDNAPIMEDIQSSLDEHNDDYISLTTSNPVDLLNRHLYEAGDRSREQVIGFGGRLDSARFRYVLSEEFDAPVQNVEGTILGEHGDAQVPVFSKVRVDGTDPEFSGDEKEQLLGDLQESAMDVIERKGATEWGPARGVAHMVEAILHDTGEVLPASVKLEGEFGHEDTAFGVPVRLGSNGVEEIVEWDLDDYEQDLMADAAEKLSDQYDKIS"
s = "MTKVSVVGAAGTVGAAAGYNIALRDIADEVVFVDIPDKEDDTVGQAADTNHGIAYDSNTRVRQGGYEDTAGSDVVVITAGIPRQPGQTRIDLAGDNAPIMEDIQSSLDEHNDDYISLTTSNPVDLLNRHLYEAGDRSREQVIGFGGRLDSARFRYVLSEEFDAPVQNVEGTILGEHGDAQVPVFSKVRVDGTDPEFSGDEKEQLLGDLQESAMDVIERKGATEWGPARGVAHMVEAILHDTGEVLPASVKLEGEFGHEDTAFGVPVRLGSNGVEEIVEWDLDDYEQDLMADAAEKLSDQYDKIS"

#    1         11        21        31        41         
for x in t[0]:
    seqresi = x.resi -1 #- 21 + 1 
    X = x.resn
    Y = s[seqresi-1]
    print '+' if X==Y else '-', X, Y, x.resi, x.seqresi, seqresi


# site2.append(Site(resn=resn, resi=resi, seqresi=seqresi, ncontacts=v))

# templates[(pdb, chain1, chain2)] = site1, site2

# matching_templates = complexes.templatesWithHits(templates, by_hit)
