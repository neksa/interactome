"""
Protein-complex homology modeling pipeline
Using PDB, Delta-BLAST, and MODELLER

Modeling bacterial interactomes
"""

import signal
# import time
# import math
# import cProfile

from multiprocessing import Pool  # , Manager  # , Value
# from itertools import islice
# from collections import defaultdict
# from pdb import PDB, PDBSearchResults, SKEMPI  # , SIFTS
# from modellerwrapper import ModellerWrapper
from interactome.sequences.blast import BLAST  # , BLASTReport
from interactome.structures.complexes import Complexes  # , SiteResidue

from interactome.workflows.merge_alignment import binding_site_alignments, merge_bs_alignments, convert_alignments
from predict_interactions_multiscoring import predict


#############################################################
def get_root():
    return "/Users/agoncear/projects/Interactome/Workflow"


def calc_blast(infile, outfile):
    """
    Step 1:
        1.1 Take sequences in fasta format for two genomes
        1.2. Run Delta-BLAST against PDB database with FASTA file as a query
    """
    query_fasta = get_root() + "/Sequences/" + infile
    BLAST_xml = get_root() + "/BLAST-results/" + outfile

    print "Step 1: Delta-BLAST"

    try:
        with open(query_fasta, 'r'): pass
    except IOError:
        print "\tError: Query fasta file", query_fasta, "is missing"
        return

    try:
        with open(BLAST_xml): pass
        print "\tInfo: BLAST results file found", BLAST_xml
    except IOError:
        blast = BLAST()
        print "\tInfo: running Delta BLAST. It may take a while..."
        blast.runDeltaBLAST(query_fasta, BLAST_xml, inclusion_threshold=0.01, domain_inclusion_threshold=0.01, format=5, remote=False)

    try:
        with open(BLAST_xml): pass
    except IOError:
        print "\tError: BLAST results expected but not found", BLAST_xml
        return

    print "Step 1: DONE"
    # hits = blast.parseHits(xml)
    # print hits


#############################################################
def calc_blast_report(infile, outfile):
    """
    Step 2.
        Convert BLAST results from XML to a tab file with required fields
    """
    # BLAST_xml = "../Workflow/BLAST-results/localDeltaBLAST_hits.xml"
    BLAST_xml = get_root() + "/BLAST-results/" + infile
    blast_report = get_root() + "/BLAST-results/" + outfile

    print "Step 2: Make BLAST report"

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
    print "Step 2: DONE"


#############################################################
def calc_binding_site_alignments(report, merged_report, templates_fname):
    print "Step 3: Merging binding site alignments"
    tmp_report1 = get_root() + "/BLAST-results/" + merged_report + ".tmp1"
    tmp_report2 = get_root() + "/BLAST-results/" + merged_report + ".tmp2"
    merged_report = get_root() + "/BLAST-results/" + merged_report

    try:
        with open(tmp_report1): pass
    except IOError:
        convert_alignments(get_root() + "/BLAST-results/" + report, tmp_report1)

    try:
        with open(tmp_report2): pass
    except IOError:
        binding_site_alignments(get_root() + "/Structures/" + templates_fname,
                                tmp_report1,
                                tmp_report2)

    try:
        with open(merged_report): pass
    except IOError:
        merge_bs_alignments(tmp_report2, merged_report)
    print "Step 3: DONE"


#############################################################
def calc_complex_alignments(merged_report, templates_fname, outfile, parallel=True):
    global blast_hits
    """
    Step 4.
        Protein complex comparative modeling workflow. From query sequences to alignments and models.
    """
    # blast_report = "deltablast_report.tab"
    pdb_interfaces_path = get_root() + "/Interfaces/"
    blast_report = get_root() + "/BLAST-results/" + merged_report
    pdb_templates = get_root() + "/Structures/" + templates_fname
    matches_fname = get_root() + "/Alignments/" + outfile

    print "Step 4: Mapping and scoring protein complexes and interfaces based on the hits"

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

        print "Loading templates..."
        templates = complexes.loadTemplates(pdb_templates)
        list_of_structures = []
        for tpl in templates:
            pdb, chain1, chain2 = tpl[0]
            list_of_structures.append(pdb.upper() + "|" + chain1)
            list_of_structures.append(pdb.upper() + "|" + chain2)
        set_of_structures = set(list_of_structures)
        templates = complexes.loadTemplates(pdb_templates)
        # print templates

        print "Loading BLAST report..."
        blast = BLAST()
        blast_hits = blast.readReport(blast_report, set_of_structures)
        # print blast_hits

        pool = None
        gen_results = None

        if parallel:
            # pool = Pool(8, init_worker, None, 10)
            pool = Pool(8, init_worker)
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
                    if c % 1000 == 0:
                        print c
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
    print "Step 4: DONE"


#############################################################
def calc_predict_complexes(alignment_fname, templates_fname, matches_fname):
    print "Step 4: Predict"
    fname = get_root() + "/Alignments/" + matches_fname
    # fname = "0"
    try:
        with open(fname):
            pass
    except IOError:
        predict(templates_fname, alignment_fname, matches_fname)
    print "Step 4: DONE"


# #############################################################
# def analyze_template_hits(name):
#     print "Loading hits..."

#############################################################
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


blast_hits = None
def func_runner(template):
    global blast_hits
    # print "DEMO runner"
    # print blast_hits["2E5H|A"]
    complexes = Complexes()
    return complexes.templatesWithHits(template, blast_hits, benchmark=False)


#############################################################
def workflow(name):
    calc_blast('Uniprot/' + name + ".fa", name + ".xml")
    calc_blast_report(name + ".xml", name + ".tab")
    calc_binding_site_alignments(name + ".tab", name + "_merged.tab", "pdb_templates_5A.tab")
    calc_predict_complexes(name + "_merged.tab", "pdb_templates_5A.tab", "matches2_{}.tab".format(name))
    #### calc_complex_alignments(name + "_merged.tab", "pdb_templates_5A.tab", "matches_{}.tab".format(name))


#############################################################
if __name__ == '__main__':
    # workflow("Scerevisiae")
    # workflow("Ecoli")
    # workflow("Ypestis")
    # workflow("Vcholerae")
    # workflow("Hpylori")
    # workflow("Spneumoniae")

    # analyze_template_hits("Hsapiens")
    workflow("Hsapiens")

    # workflow("ecoli")
    # cProfile.run('workflow("ecoli")')
    # workflow("yeast")
    # workflow("human")
    # workflow("HPY_26695_PPI")
    # workflow("SPN_TIGR4_PPI")
