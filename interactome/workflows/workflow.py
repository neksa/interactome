"""
Protein-complex homology modeling pipeline
Using PDB, Delta-BLAST, and MODELLER

Benchmarking the quality of modeling by (re)modeling the human dimers
"""

import signal
# import time
# import math

from multiprocessing import Pool, Manager  # , Value

# from itertools import islice
# from collections import defaultdict
# from pdb import PDB, PDBSearchResults, SKEMPI  # , SIFTS
# from modellerwrapper import ModellerWrapper
from interactome.sequences.blast import BLAST  # , BLASTReport
from interactome.structures.complexes import Complexes  # , SiteResidue


#############################################################
def get_root():
    return "/Users/agoncear/projects/Interactome/Workflow"

def step1():
    """
    Step 1:
        1.1. Preparing a benchmark with Human dimers, <2.5A resolution, with two chains in bioassembly, each chain at least 50 residues long
           Manually filter PDB entries using PDB website (Advanced search option),
           Include the following fields: PDB ID,Chain ID,Resolution,Source,Taxonomy ID,Entity ID,Sequence,Chain Length,DB ID 
           Save as ../Workflow/Searches/tabularResults.csv
        1.2. Extract sequences from the search results and format them as FASTA
        1.3. Run Delta-BLAST against PDB database with FASTA file as a query
    """
    search_results = get_root() + "/Searches/tabularResults.csv"
    query_fasta = get_root() + "/BLAST-results/query_sequence.fasta"
    # BLAST_xml = get_root() + "/BLAST-results/localDeltaBLAST_hits.xml"
    BLAST_xml = get_root() + "/BLAST-results/HumanDeltaBLAST_hits.xml"

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
def step2():
    """
    Step 2.
        Convert BLAST results from XML to a tab file with required fields
    """
    # BLAST_xml = "../Workflow/BLAST-results/localDeltaBLAST_hits.xml"
    BLAST_xml = get_root() + "/BLAST-results/HumanDeltaBLAST_hits.xml"
    blast_report = get_root() + "/BLAST-results/human_deltablast_report2.tab"

    print "Step 2"

    try:
        with open(BLAST_xml): pass
    except IOError:
        print "\tError: BLAST results not found", BLAST_xml
        return

    try:
        with open(blast_report): pass
        # print "\tInfo: BLAST report found", blast_report
    except IOError:
        blast = BLAST()
        print "\tInfo: BLAST report not found. Converting BLAST results to a tab report. It may take a while..."
        blast.makeReport(BLAST_xml, blast_report)
    print "Step 2 DONE"


#############################################################

def step3():
    global blast_hits

    """
    Step 3.
        Protein complex comparative modeling workflow. From query sequences to alignments and models.
    """
    # blast_report = "deltablast_report.tab"
    blast_report = get_root() + "/BLAST-results/human_deltablast_report2.tab"
    pdb_templates = get_root() + "/Structures/pdb_templates_5A.tab"
    matches_fname = get_root() + "/Alignments/matches_human.tab"
    pdb_interfaces_path = get_root() + "/Interfaces/"
    # mapping_fname = "../../../data/EBI/SIFTS/pdb_chain_uniprot.tsv"

    parallel = True
    # parallel = False

    print "Step 3"

    try:
        with open(blast_report): pass
        print "\tInfo: BLAST report found", blast_report
    except IOError:
        print "\tError: BLAST report not found", blast_report
        return

    complexes = Complexes()
    try:
        with open(pdb_templates): pass
        print "\tInfo: PDB templates summary found", pdb_templates
    except IOError:
        complexes.collectTemplates(pdb_interfaces_path, pdb_templates)

    try:
        with open(pdb_templates): pass
    except IOError:
        print "\tError: No PDB templates summary found", pdb_templates
        return

    try:
        with open(matches_fname): pass
    except IOError:

        # data = dict()
        # if parallel:
        #     manager = Manager()
        #     data = manager.dict()

        print "Loading templates..."
        templates = complexes.loadTemplates(pdb_templates)

        print "Loading BLAST report..."
        blast = BLAST()
        blast_hits = blast.readReport(blast_report)

        # print "DEMO"
        # print data["2E5H|A"]

        pool = None
        gen_results = None

        if parallel:
            pool = Pool(8, init_worker, 10)
            gen_results = pool.imap_unordered(func_runner, templates)  # 5000
            """
            pool.imap_unordered(extract_contacts, ifilter(isNMR, mmcif.listAll()))
            """
        else:
            init_worker()
            gen_results = (func_runner(template) for template in templates)

        try:
            # while(True):
            #     print "Watchdog... every 60 seconds (Ctrl-C to interrupt)"
            #     time.sleep(60)
            c = 0
            print "Writing results to", matches_fname
            with open(matches_fname, 'w') as o:
                o.write("query_A\tquery_B\ttemplate\tid_A\taln_len_A\tsite_A\tid_B\taln_len_B\tsite_B\n")

                for result in gen_results:
                    if c % 1000 == 0: print c
                    c += 1

                    tpl, queries = result
                    tpl_pdb, tpl_A, tpl_B = tpl
                    for query_A, query_B, params in queries:
                        (A, site_A_str), (B, site_B_str) = params
                        identity_A = int(round(sum([x.identity for x in A]) / float(len(A))))
                        identity_B = int(round(sum([x.identity for x in B]) / float(len(B))))

                        aln_len_A = max(
                            sum([x.q_to - x.q_from + 1 for x in A]),
                            sum([x.h_to - x.h_from + 1 for x in A]))
                        aln_len_B = max(
                            sum([x.q_to - x.q_from + 1 for x in B]),
                            sum([x.h_to - x.h_from + 1 for x in B]))

                        o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            query_A,
                            query_B,
                            tpl_pdb + '|' + tpl_A + '|' + tpl_B,
                            identity_A, aln_len_A, site_A_str,
                            identity_B, aln_len_B, site_B_str))

        except KeyboardInterrupt:
            print "Caught KeyboardInterrupt, terminating workers"
            if pool is not None:
                pool.terminate()
                pool.join()
                pool = None
            else:
                raise

        # pool.close()
        if pool is not None:
            pool.close()
            pool.join()
        # sifts = SIFTS()
        # mapping = sifts.getPDBMapping(mapping_fname)
        # print "Info: Mapping loaded"
    print "Step 3 DONE"


#############################################################
def step4():
    pass
    print "Step 4"
    print "Step 4 DONE"


#############################################################
def human_benchmark():
    step1()
    step2()
    step3()
    step4()


#############################################################
blast_hits = None
def func_runner(template):
    global blast_hits
    # print "DEMO runner"
    # print blast_hits["2E5H|A"]
    complexes = Complexes()
    return complexes.templatesWithHits(template, blast_hits, benchmark=False)


def init_worker():
    # print "Worker initializing..."
    # blast_hits = data[0]
    # print "DEMO init"
    # print blast_hits["2E5H|A"]
    # load processed BLAST hits
    # Allow only hits with > 25% identity, E-value < 0.01, and at least 25 residues long alignments
    # print "Info: Loaded BLAST report with", len(by_query.keys()), "queries and ", len(by_hit.keys()), "hits"
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    # print "Worker initialized."


#############################################################
if __name__ == '__main__':
    human_benchmark()


