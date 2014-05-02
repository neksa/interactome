"""
"""
import xml.etree.ElementTree as ET
from collections import defaultdict, namedtuple
from textwrap import TextWrapper

import httplib
import csv
from itertools import islice
from cStringIO import StringIO

import glob
import subprocess

class BLAST:

    def pirToQuery(self, pdb_chain_list, query_filename):
        fnames = []
        for pdb, (chain1, chain2) in pdb_chain_list.iteritems():
            fname_format = "../Workflow/Sequences/{}{}.pir"
            fname1 = fname_format.format(pdb, chain1)
            fname2 = fname_format.format(pdb, chain2)
            try:
                f = open(fname1)
                f.close()
                f = open(fname2)
                f.close()
            except IOError:
                print "Could not read sequence for one of the chains of structure PDB {}".format(pdb)
                continue
            fnames.append(fname1)
            fnames.append(fname2)

        sequence = ""
        for fname in fnames:
            with open(fname, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("structure"): continue
                    if line == "": continue
                    if line.endswith("*"):
                        line = line[:-1]
                    if line.startswith(">"):
                        line = ">" + line[4:9]
                    sequence += line + "\n"

        with open(query_filename, 'w') as o:
            o.write(sequence)


    def runBLASTP(self, sequence_file, results_file, format=5):
        """
        Returns XML
        """

        try:
            with open(results_file, 'r') as f:
                output = f.read()
                return output
        except IOError:
            pass

        proc = subprocess.Popen([
            "/usr/local/ncbi/blast/bin/blastp",
            "-out", results_file,
            "-db", "pdb",
            "-query", sequence_file, 
            "-outfmt", str(format),
            "-evalue", "1",
            "-remote"], stdout = subprocess.PIPE, stderr=subprocess.PIPE)

        output, err = proc.communicate()[0]
        print err
        # with open(results_file, 'w') as o:
        #     o.write(output)
        # return output

    def runDeltaBLAST(self, sequence_file, results_file, format=5):
        """
        Returns XML
        """

        try:
            with open(results_file, 'r') as f:
                output = f.read()
                return output
        except IOError:
            pass

        proc = subprocess.Popen([
            "/usr/local/bin/deltablast",
            "-out", results_file,
            "-db", "pdb",
            "-query", sequence_file, 
            "-outfmt", str(format),
            "-remote",
            "-inclusion_ethresh", "0.05",
            "-domain_inclusion_ethresh", "0.05"], stdout = subprocess.PIPE, stderr=subprocess.PIPE)

        output, err = proc.communicate()[0]
        print err
        # with open(results_file, 'w') as o:
        #     o.write(output)
        # return output



    def parseHits(self, xml, identity_filter = lambda x: 25 <= x < 100):
        """
        Returns Q-templates dictionary
        """
        root = ET.fromstring(xml)
        # print self.root
        matches = defaultdict(list)

        for iteration in root.findall("./BlastOutput_iterations/Iteration"):
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
                    
                    if not identity_filter(hsp_identity): continue

                    hsp_gaps = int(hsp.find("Hsp_gaps").text)
                    hsp_align_len = int(hsp.find("Hsp_align-len").text)

                    # alignment:
                    hsp_qseq = hsp.find("Hsp_qseq").text
                    hsp_hseq = hsp.find("Hsp_hseq").text

                    matches[query_id].append((hit_id, hsp_identity, hsp_gaps, hsp_qseq, hsp_hseq))
        return matches

    def iterParseHits(self, xml):
        """
        Returns Q-templates dictionary
        """
        # root = ET.fromstring(xml)
        # print self.root
        matches = defaultdict(list)

        context = ET.iterparse(xml, events=("start", "end"))
        context = iter(context)
        event, root = context.next()

        BLASTHits = namedtuple('BLASTHits', 'hit, score, bit_score, evalue, identical, positive, gaps, align_len, qseq, hseq, query_from, query_to, hit_from, hit_to')

        cnt = 0
        for event, elem in context:
            if event == "end" and elem.tag == "Iteration":
                root.clear()
                ##### BLAST iteration parsing #####
                # for iteration in root.findall("./BlastOutput_iterations/Iteration"):
                iteration = elem
                query_id = iteration.find("Iteration_query-def").text
                query_id = query_id[0:5] # 2HLEA from 2HLEA_ATOM
                cnt += 1
                print cnt, query_id
                # if cnt > 15: break
                # print iteration
                for hit in iteration.findall("Iteration_hits/Hit"):
                    # print hit
                    hit_id = hit.find("Hit_accession").text
                    # print hit_id
                    hit_len = int(hit.find("Hit_len").text)
                    for hsp in hit.findall("Hit_hsps/Hsp"):
                        hsp_identical = int(hsp.find("Hsp_identity").text)
                        hsp_positive = int(hsp.find("Hsp_positive").text)
                        hsp_score = int(hsp.find("Hsp_score").text)
                        hsp_bit_score = float(hsp.find("Hsp_bit-score").text)
                        hsp_evalue = float(hsp.find("Hsp_evalue").text)
                        
                        # if not identity_filter(hsp_identity): continue
                        hsp_gaps = int(hsp.find("Hsp_gaps").text)
                        hsp_align_len = int(hsp.find("Hsp_align-len").text)

                        # alignment:
                        hsp_qseq = hsp.find("Hsp_qseq").text
                        hsp_hseq = hsp.find("Hsp_hseq").text

                        # coordinates in sequence:
                        hsp_query_from = hsp.find("Hsp_query-from").text
                        hsp_query_to = hsp.find("Hsp_query-to").text

                        hsp_hit_from = hsp.find("Hsp_hit-from").text
                        hsp_hit_to = hsp.find("Hsp_hit-to").text

                        match = BLASTHits._make((
                            hit_id, hsp_score, hsp_bit_score, hsp_evalue,
                            hsp_identical, hsp_positive, hsp_gaps, hsp_align_len, hsp_qseq, hsp_hseq,
                            hsp_query_from, hsp_query_to, hsp_hit_from, hsp_hit_to))
                        matches[query_id].append(match)

                elem.clear()

        # for event, elem in iterparse(xml):
        #     if elem.tag == "Iteration":
        #         # .... process
        #         elem.clear()

        return matches

