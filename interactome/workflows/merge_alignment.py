"""
Protein-complex homology modeling pipeline
Using PDB, Delta-BLAST, and MODELLER

Benchmarking the quality of modeling by (re)modeling the human dimers
"""

import signal
# import time
# import math

from multiprocessing import Pool, Manager  # , Value

from itertools import islice
from collections import defaultdict
# from pdb import PDB, PDBSearchResults, SKEMPI  # , SIFTS
# from modellerwrapper import ModellerWrapper
from interactome.sequences.blast import BLAST, BLASTReport
from interactome.structures.complexes import Complexes  # , SiteResidue



def merge_alignments(fname):
    """
        Parse format saved by makeReport
        convert values to float and int from string, calculate sequence identity in %

        Returns two dictionaries indexed by different keys: by hit and by query

        Reports for the same query are joined in one list. They have to be overlapped in a common alignment.
    """

    by_hit = defaultdict(lambda: defaultdict(list))
    # by_query = defaultdict(list)

    c = 0
    with open(fname) as f:
        for line in islice(f, 1, None):  # None
            fields = list(line.strip().split("\t"))
            c += 1
            # if c % 100000 == 0:
            #     print '.',
            # if c % 1000000 == 0: break
            # pdb, pdb_chain = fields[1].split("|")

            # fields[1] = fields[1].replace('|', '') # remove bar from hit name (separating PDB ID and chain)
            fields[2:5] = map(float, fields[2:5])  # score, bit_score, evalue

            # print fields[5:13]
            fields[5:13] = map(int, fields[5:13])  # identical, positive, gaps, align_len, q_from, q_to, h_from, h_to
            # print fields[5:13]

            identity = int(round(float(fields[5])*100.0 / float(fields[8]), 0))  # percent identity = identical / align_len
            fields.append(identity)
            report = BLASTReport._make(fields)

            if not (15 <= identity):
                continue
            if fields[8] < 25:
                continue  # report.align_len
            if fields[4] > 0.001:
                continue  # report.evalue

            # print report.identity
            # by_hit[report.hit].append(report)
            # by_query[report.query].append(report)
            # by_hit[fields[1]][fields[0]].append(tuple(fields)) # report.hit
            if len(by_hit[report.hit][report.query]) > 0:
                # print by_hit[report.hit][report.query]
                # print "NEW:"
                # print report
                print report.query, report.q_from
                print report.qseq
                # ALIGNMENT

            by_hit[report.hit][report.query].append(report)  # report.hit
            # by_query[report.query].append(fields)


if __name__ == '__main__':
    merge_alignments("/Users/agoncear/projects/Interactome/Workflow/BLAST-results/human_deltablast_report2.tab")
