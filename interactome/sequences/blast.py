"""
"""
import xml.etree.ElementTree as ET
from collections import defaultdict, namedtuple
from itertools import islice

import subprocess

BLASTReport = namedtuple('BLASTReport', 'query, hit, score, bit_score, evalue, identical, positive, gaps, align_len, q_from, q_to, h_from, h_to, qseq, hseq, identity')
BLASTHits = namedtuple('BLASTHits', 'hit, score, bit_score, evalue, identical, positive, gaps, align_len, qseq, hseq, query_from, query_to, hit_from, hit_to')


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
                    if line.startswith("structure"):
                        continue
                    if line == "":
                        continue
                    if line.endswith("*"):
                        line = line[:-1]
                    if line.startswith(">"):
                        line = ">" + line[4:9]
                    sequence += line + "\n"

        with open(query_filename, 'w') as o:
            o.write(sequence)

    def runBLASTP(self, sequence_file, results_file, evalue=0.01, format=5, remote=True):
        """
        Returns XML
        """

        try:
            with open(results_file, 'r') as f:
                output = f.read()
                return output
        except IOError:
            pass

        params = ["/usr/local/ncbi/blast/bin/blastp",
                    "-out", results_file,
                    "-query", sequence_file,
                    "-outfmt", str(format),
                    "-evalue", "{0:.5f}".format(evalue)]
        if remote:
            params.append("-db")
            params.append("pdb")
            params.append("-remote")
        else:
            params.append("-num_threads")
            params.append("8")
            params.append("-db")
            params.append("pdbaa")

        proc = subprocess.Popen(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = proc.communicate()
        print err

    def runDeltaBLAST(self, sequence_file, results_file, inclusion_threshold=0.01, domain_inclusion_threshold=0.01, format=5, remote=True):
        """
        Returns XML
        """
        try:
            with open(results_file, 'r') as f:
                output = f.read()
                return output
        except IOError:
            pass

        params = ["/usr/local/bin/deltablast",
                    "-out", results_file,
                    "-query", sequence_file,
                    "-outfmt", str(format),
                    "-inclusion_ethresh", "{0:.5f}".format(inclusion_threshold),
                    "-domain_inclusion_ethresh", "{0:.5f}".format(domain_inclusion_threshold)]
        if remote:
            params.append("-db")
            params.append("pdb")
            params.append("-remote")
        else:
            params.append("-num_threads")
            params.append("8")
            params.append("-db")
            params.append("pdbaa")

        proc = subprocess.Popen(params, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        output, err = proc.communicate()
        print err
        # with open(results_file, 'w') as o:
        #     o.write(output)
        # return output

    def parseHits(self, xml):  # identity_filter = lambda x: 25 <= x < 100):
        """
        Returns Q-templates dictionary
        """
        root = ET.fromstring(xml)
        # print self.root
        matches = defaultdict(list)

        for iteration in root.findall("./BlastOutput_iterations/Iteration"):
            query_id = iteration.find("Iteration_query-def").text

            # ATTENTION: change this in case of new BLAST query types other than swissprot or PDB
            if query_id.startswith("sp|"):
                query_id = query_id.split("|", 2)[1]
            # query_id = query_id[0:5] # 2HLEA from 2HLEA_ATOM

            # print query_id
            # print iteration
            for hit in iteration.findall("Iteration_hits/Hit"):
                # print hit
                hit_id = hit.find("Hit_accession").text
                # print hit_id
                # hit_len = int(hit.find("Hit_len").text)
                for hsp in hit.findall("Hit_hsps/Hsp"):
                    hsp_identity = int(hsp.find("Hsp_identity").text)

                    # if not identity_filter(hsp_identity): continue

                    hsp_gaps = int(hsp.find("Hsp_gaps").text)
                    # hsp_align_len = int(hsp.find("Hsp_align-len").text)

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

        cnt = 0
        for event, elem in context:
            if event == "end" and elem.tag == "Iteration":
                root.clear()
                ##### BLAST iteration parsing #####
                # for iteration in root.findall("./BlastOutput_iterations/Iteration"):
                iteration = elem

                query_id = iteration.find("Iteration_query-def").text
                if query_id.startswith("sp|"):
                    query_id = query_id.split("|", 2)[1]
                # query_id = query_id[0:5] # 2HLEA from 2HLEA_ATOM

                cnt += 1
                print cnt, query_id
                # if cnt > 15: break
                # print iteration
                for hit in iteration.findall("Iteration_hits/Hit"):
                    # print hit
                    hit_id = hit.find("Hit_accession").text
                    # hit_id = hit_id.replace("_", "|") # this is actually done separately in conversion
                    # print hit_id
                    # hit_len = int(hit.find("Hit_len").text)
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

    def makeReport(self, BLAST_xml, blast_report):
        with open(BLAST_xml, 'r') as f, open(blast_report, 'w') as o:
            query_hits = self.iterParseHits(f)
            o.write("query\thit\tscore\tbit_score\tevalue\tidentical\tpositive\tgaps\talign_len\tqfrom\tqto\thfrom\thto\tqseq\thseq\n")
            for query, hits in query_hits.iteritems():
                for h in hits:
                    o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        query, h.hit, h.score, h.bit_score, h.evalue,
                        h.identical, h.positive, h.gaps, h.align_len,
                        h.query_from, h.query_to, h.hit_from, h.hit_to, h.qseq, h.hseq))

    def readReport(self, blast_report, set_of_structures=None):
        """
            Parse format saved by makeReport
            convert values to float and int from string, calculate sequence identity in %

            Returns two dictionaries indexed by different keys: by hit and by query

            Reports for the same query are joined in one list. They have to be overlapped in a common alignment.
        """

        by_hit = defaultdict(lambda: defaultdict(list))
        # by_query = defaultdict(list)

        c = 0
        with open(blast_report) as f:
            for line in islice(f, 1, None):  # None
                fields = list(line.strip().split("\t"))
                c += 1
                if c % 100000 == 0:
                    print '.',
                # if c % 1000000 == 0: break
                # pdb, pdb_chain = fields[1].split("|")
                structure = fields[1]
                if set_of_structures is not None and structure not in set_of_structures:
                    continue

                # fields[1] = fields[1].replace('|', '') # remove bar from hit name (separating PDB ID and chain)
                fields[2:5] = map(float, fields[2:5])  # score, bit_score, evalue

                # print fields[5:13]
                fields[5:13] = map(int, fields[5:13])  # identical, positive, gaps, align_len, q_from, q_to, h_from, h_to
                # print fields[5:13]

                identity = int(round(float(fields[5])*100.0 / float(fields[8]), 0))  # percent identity = identical / align_len
                fields.append(identity)
                report = BLASTReport._make(fields)

                # if not (15 <= identity):
                #     continue

                # if fields[8] < 25:
                #     continue  # report.align_len

                # if fields[4] > 0.01:
                #     continue  # report.evalue

                # print report.identity
                # by_hit[report.hit].append(report)
                # by_query[report.query].append(report)
                # by_hit[fields[1]][fields[0]].append(tuple(fields)) # report.hit
                by_hit[report.hit][report.query].append(report)  # report.hit
                # by_query[report.query].append(fields)

        summary = 0
        sort_by_positive = lambda x: x.positive
        for k, v in by_hit.iteritems():
            for q in v.iterkeys():
                by_hit[k][q].sort(reverse=True, key=sort_by_positive)
                summary += len(by_hit[k][q])

        print "BLAST results: loaded", summary, "records"

        return dict(by_hit)  # , by_query
