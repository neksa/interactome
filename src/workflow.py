"""
Protein-complex homology modeling pipeline
Using PDB, Delta-BLAST, and MODELLER

Benchmarking the quality of modeling by (re)modeling the human dimers
"""

from itertools import islice
from collections import defaultdict

from pdb import PDB, PDBSearchResults, SKEMPI
from modellerwrapper import ModellerWrapper
from blast import BLAST
from complexes import Complexes

import math

#############################################################
#############################################################

def skempi_benchmark_workflow():
    """
    1. Take dimers from SKEMPI
    2. Download PDB files
    3. Extract sequences from PDB files 
    4. And use them as queries in BLASTp
    5. In BLAST hits find templates that are protein complexes
    6. Format hetero-dimer alignments according to MODELLER rules
    7. Run MODELLER with the alignments to produce models of complexes

    Protein complex comparative modeling workflow. From query sequences to alignments and models.

    Conventions:
    PDB code: 4 chars capitalized
        PPPP
    chains: letters capitalized, follor PDB code. Up to 2 letters
        AB

    dimer is represented by pdb code and two chain letters:
        PPPPAB

    Alignments:
        Modeller-compatible PIR format
    """

    # Step 1
    skempi = SKEMPI()
    dimers = skempi.getDimers() 
    # print dimers

    # Step 2
    pdb = PDB()
    pdb.getStructureFiles(dimers.keys())

    # Step 3
    pdb.extractSEQRES(dimers.keys())

    # Step 4 
    query_fasta = "../Workflow/BLAST-results/query_sequence.fasta"
    BLAST_xml = "../Workflow/BLAST-results/BLAST_hits.xml"
    
    blast = BLAST()
    blast.pirToQuery(dimers, query_fasta)
    xml = blast.runBLASTP(query_fasta, BLAST_xml)
    hits = blast.parseHits(xml)

    # Steps 5,6
    c = Complexes()
    tpl_complexes = c.templatesComplexes(dimers, hits)
    template_pdb_codes = c.getTemplatePDBCodes(tpl_complexes)
    print template_pdb_codes
    pdb.getStructureFiles(template_pdb_codes)
    template_seqres = pdb.extractSEQRES(template_pdb_codes)

    c.alignments(tpl_complexes, query_fasta)
    # Optional Step 7.  verify model building
    c.models(tpl_complexes)


#############################################################
def human_benchmark_workflow_step1():
    """
    Step 1:
        1.1. Preparing a benchmark with Human dimers, <2.5A resolution, with two chains in bioassembly, each chain at least 50 residues long
           Manually filter PDB entries using PDB website (Advanced search option),
           Include the following fields: PDB ID,Chain ID,Resolution,Source,Taxonomy ID,Entity ID,Sequence,Chain Length,DB ID 
           Save as ../Workflow/Searches/tabularResults.csv
        1.2. Extract sequences from the search results and format them as FASTA
        1.3. Run Delta-BLAST against PDB database with FASTA file as a query
    """
    search_results = "../Workflow/Searches/tabularResults.csv"
    query_fasta = "../Workflow/BLAST-results/query_sequence.fasta"
    BLAST_xml = "../Workflow/BLAST-results/BLAST_hits.xml"

    print "Step 1"

    try:
        with open(search_results): pass
        print "\tInfo: search results file found", search_results
    except IOError:
        print "\tError: PDB search results file", search_results, "is missing"
        return

    try:
        with open(query_fasta, 'r'): pass
        print "\tInfo: query fasta file found", query_fasta
    except IOError:
        print "\tInfo: File", query_fasta, "is missing. Processing", search_results
        sr = PDBSearchResults()
        results = sr.getTabularResults()
        results_filtered = sr.filterPDB(results)
        sr.searchToFasta(results_filtered, query_fasta)

    try:
        with open(query_fasta, 'r'): pass
    except IOError:
        print "\tError: Query fasta file", query_fasta, "is missing after search results processing"
        return

    try:
        with open(BLAST_xml): pass
        print "\tInfo: BLAST results file found", BLAST_xml
    except IOError:
        blast = BLAST()
        print "\tInfo: running Delta BLAST. It may take a while..."
        blast.runDeltaBLAST(query_fasta, BLAST_xml)

    try:
        with open(BLAST_xml): pass
    except IOError:
        print "\tError: BLAST results expected but not found", BLAST_xml
        return
    
    print "Step 1 DONE"
    # hits = blast.parseHits(xml)
    # print hits

#############################################################
def human_benchmark_workflow_step2():
    """
    Step 2.
        Convert BLAST results from XML to a tab file with required fields
    """
    BLAST_xml = "../Workflow/BLAST-results/localDeltaBLAST_hits.xml"
    blast_report = "deltablast_report.tab"
    
    print "Step 2"

    try:
        with open(BLAST_xml): pass
    except IOError:
        print "\tError: BLAST results not found", BLAST_xml
        return

    try: 
        with open(blast_report): pass
        print "\tInfo: BLAST report found", blast_report
    except IOError:
        blast = BLAST()
        print "\tInfo: BLAST report not found. Converting BLAST results to a tab report. It may take a while..."
        with open(BLAST_xml, 'r') as f, open(blast_report, 'w') as o:
            query_hits = blast.iterParseHits(f)
            o.write("query\thit\tscore\tbit_score\tevalue\tidentical\tpositive\tgaps\talign_len\tqfrom\tqto\thfrom\thto\tqseq\thseq\n")
            for query, hits in query_hits.iteritems():
                for h in hits:
                    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        query, h.hit, h.score, h.bit_score, h.evalue,
                        h.identical, h.positive, h.gaps, h.align_len,
                        h.query_from, h.query_to, h.hit_from, h.hit_to, h.qseq, h.hseq))
    print "Step 2 DONE"


#############################################################
def human_benchmark_workflow_step3():
    """
    Step 3.
        Protein complex comparative modeling workflow. From query sequences to alignments and models.

        Conventions:
        PDB code: 4 chars capitalized
            PPPP
        chains: letters capitalized, follor PDB code. Up to 2 letters
            AB

        dimer is represented by pdb code and two chain letters:
            PPPPAB

        Alignments:
            Modeller-compatible PIR format
    """
    blast_report = "deltablast_report.tab"
    pdb_templates = "pdb_templates.tab"
    matches_fname = "matches.tab"
    pdb_interfaces_path = "../scoring/results/"

    print "Step 3"

    try:
        with open(blast_report): pass
        print "\tInfo: BLAST report found", blast_report
    except IOError:
        print "\tError: BLAST report not found", blast_report
        return

    c = Complexes()
    try:
        with open(pdb_templates): pass
        print "\tInfo: PDB templates summary found", pdb_templates
    except IOError:
        c.collectTemplates(pdb_interfaces_path, pdb_templates)

    try:
        with open(pdb_templates): pass
    except IOError:
        print "\tError: No PDB templates summary found", pdb_templates
        return

    templates = c.loadTemplates(pdb_templates)

    # load processed BLAST hits 
    # Allow only hits with > 25% identity, E-value < 0.01, and at least 25 residues long alignments
    hits = defaultdict(list)
    with open(blast_report) as f:
        for line in islice(f, 1, None):

            query, hit, score, bit_score, evalue, identical, positive, gaps, align_len,\
            q_from, q_to, h_from, h_to, qseq, hseq = line.strip().split("\t")

            identity = int(round(float(identical) / float(align_len), 2) * 100.0) # percent identity
            if identity < 25: continue
            if float(align_len) < 20.0: continue
            # if math.log10(float(evalue)) > -2: continue  # E < 0.01
            if float(evalue) > 0.01: continue 
            # hits[query].append(hit)

            hit = hit.replace('_', '')
            hits[hit].append((query, identity, (q_from, q_to), (h_from, h_to)))

    matching_templates = c.templatesWithHits(templates, hits)
    with open(matches_fname, 'w') as o:
        o.write("query\ttemplate\tid_A\tq_A\tt_B\tsite_A\tid_B\tq_B\tt_B\tsite_B\n")
        for tpl, queries in matching_templates.iteritems():
            tpl_pdb, tpl_A, tpl_B = tpl
            for q_pdb, q_A, q_B, params in queries:
                identity_A, qfrom_A, qto_A, hfrom_A, hto_A, site_A = params[0]
                identity_B, qfrom_B, qto_B, hfrom_B, hto_B, site_B = params[1]
                o.write("{}\t{}\t{}\t{}-{}\t{}-{}\t{}\t{}\t{}-{}\t{}-{}\t{}\n".format(
                    q_pdb + q_A + q_B,
                    tpl_pdb + tpl_A + tpl_B,
                    identity_A, qfrom_A, qto_A, hfrom_A, hto_A, site_A,
                    identity_B, qfrom_B, qto_B, hfrom_B, hto_B, site_B))

    print "Step 3 DONE"


#############################################################
def human_benchmark_workflow_step4():
    pass
    print "Step 4"
    print "Step 4 DONE"


#############################################################
def human_benchmark():
    human_benchmark_workflow_step1()
    human_benchmark_workflow_step2()
    human_benchmark_workflow_step3()
    human_benchmark_workflow_step4()


#############################################################
if __name__ == '__main__':
    human_benchmark()


