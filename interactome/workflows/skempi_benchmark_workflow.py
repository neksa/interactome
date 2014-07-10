#!/usr/bin/env python

from pdb import PDB, SKEMPI  # , SIFTS
# from modellerwrapper import ModellerWrapper
from blast import BLAST  # , BLASTReport
from complexes import Complexes  # , SiteResidue


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
    # template_seqres = pdb.extractSEQRES(template_pdb_codes)

    c.alignments(tpl_complexes, query_fasta)
    # Optional Step 7.  verify model building
    c.models(tpl_complexes)
