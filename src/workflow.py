"""
Protein complex comparative modeling workflow.

From query sequences to alignments and models

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

import xml.etree.ElementTree as ET
from collections import defaultdict
from textwrap import fill

import httplib
import xml.etree.ElementTree as ET
import pickle
import zlib

from modeller import * #model, environ, alignment
from modeller.automodel import automodel    # Load the automodel class

import glob
import subprocess

class SKEMPI:
    def getDimers(self):
        dimers = defaultdict(list)
        with open("../SKEMPI/dimers.csv") as f:
            for line in f:
                pdb, chain1, chain2 = line.strip().split("_", 2)
                dimers[pdb] = (chain1, chain2)
        return dimers


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


    """
    Returns XML
    """
    def runBLASTP(self, sequence_file, results_file, format=5):

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
            "-remote"], stdout = subprocess.PIPE, stderr=subprocess.PIPE)

        output, err = proc.communicate()[0]
        print err
        with open(results_file, 'w') as o:
            o.write(output)
        return output



    """
    Returns Q-templates dictionary
    """
    def parseHits(self, xml, identity_filter = lambda x: 25 <= x < 100):

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


class PDB:
    def fetchRemoteStructure(self, pdb):
        conn = httplib.HTTPConnection("www.rcsb.org", timeout = 10)
        print "Downloading PDB {}...".format(pdb)
        conn.request("GET", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="+pdb.lower())
        r1 = conn.getresponse()
        
        if r1.status == httplib.NOT_FOUND: # 404
            raise ValueError("Structure not found for PDB %s (www.rcsb.org)" % (pdb,))
            
        # print r1.status, r1.reason
        pdb_structure = r1.read()
        return pdb_structure


    def getStructureFiles(self, pdb_list):
        for pdb in pdb_list:
            fname = "../Workflow/Structures/{}.pdb".format(pdb)
            try:
                f = open(fname, 'r')
                f.close()
                continue
            except IOError:
                pass
            
            try:
                structure = self.fetchRemoteStructure(pdb)
            except ValueError:
                print "Attention: structure not found in PDB", pdb
                continue

            try:
                with open(fname, 'w') as o:
                    o.write(structure)
            except IOError:
                print "Could not save PDB structure downloded from the PDB"


    """
    Returns dictionary
    chain: (header, sequence)
    """
    def extractSEQRES(self, pdb_list):
        m = None
        for pdb in pdb_list:
            fnames_pir = "../Workflow/Sequences/{}*.pir".format(pdb)
            fname_struct = "../Workflow/Structures/{}.pdb".format(pdb)

            if len(glob.glob(fnames_pir)) > 0:
                # print "File exists. skipping", glob.glob(fnames_pir)
                continue

            if m == None:
                m = Modeller()

            try:
                f = open(fname_struct, 'r')
                f.close()
                m.writeChainSeq(pdb)

            except IOError:
                print "No structure for PDB {}".format(pdb)
                continue


    """
    Returns dictionary
    chain: (header, sequence)
    """
    def extractATOMSEQ(self):
        pass


class Modeller:

    def __init__(self):
        # log.verbose()

        self.env = environ()
        self.env.io.atom_files_directory = ['.', '../Workflow/Structures/']


    def writeChainSeq(self, pdb):
        mdl = model(self.env, file="{}.pdb".format(pdb))
        for c in mdl.chains:
            if c.filter(structure_types='structureN structureX'):
                filename = "../Workflow/Sequences/{}{}.pir".format(pdb, c.name)
                print("Wrote out " + filename)
                atom_file, align_code = c.atom_file_and_code(filename)
                c.write(filename, atom_file, align_code,
                            format='PIR',
                            chop_nonstd_termini=True)

    def align(self, query_fasta_file, query_code, tpl_code):
        aln = alignment(self.env)
        aln.append(file=query_fasta_file, alignment_format='FASTA', align_codes=(query_code))
        tpl_chain_file = "../Workflow/Sequences/{}.pir".format(tpl_code)
        aln.append(file=tpl_chain_file, alignment_format='PIR')
        aln.salign(overhang=10, gap_penalties_1d=(-450, -50), output='ALIGNMENT')
        aln.write(file='modeller_tmp_aln.fasta', alignment_format='FASTA')

        
        aln_text = ''
        with open('modeller_tmp_aln.fasta', 'r') as f:
            aln_text = f.read()
        return aln_text
        # aln_B = alignment(self.env)
        # aln_B.append(file=query_fasta_file, alignment_format='FASTA', alignment_codes=query_B_code)
        # aln_B.append(file=tpl_chain2_file, alignment_format='PIR')
        

    def model(self, ali, template, query):

        a = automodel(self.env,
                    alnfile  = ali ,
                    knowns   = template,
                    sequence = 'query'+query)

        a.starting_model = 1
        a.ending_model   = 1
        a.make()


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


    def fastaToSequences(self, fasta):
        seq_id = ""
        seq = ""
        sequences = {}
        for line in fasta.split("\n"):
            # line = line.strip()
            if line.startswith(">"):
                if seq_id != "": sequences[seq_id] = seq
                seq_id = line[1:6]
                seq = ""
            else:
                seq += line
            print line, seq_id, seq

        if seq_id != "": sequences[seq_id] = seq
        return sequences

    """
    Returns alignments:
    dictionary, chain sequences
    """
    def alignments(self, templates, query_fasta_file):
        model_num = 0
        m = Modeller()
        for query_pdb, chains_complexes in templates.iteritems():
            for (query_chain1, query_chain2), template_pdbs in chains_complexes.iteritems():
                for tpl_pdb, tpl_chain1, tpl_chain2 in template_pdbs:
                    model_num += 1
                    print "* Model alignment {}: ".format(model_num)
                    print "  - Query complex [{} {} {}]".format(query_pdb, query_chain1, query_chain2)
                    print "  - Template complex : {} {} {}".format(tpl_pdb, tpl_chain1, tpl_chain2)

                    sequences = {}
                    query_codes = (query_pdb + query_chain1, query_pdb + query_chain2)
                    tpl_codes = (tpl_pdb + tpl_chain1, tpl_pdb + tpl_chain2)
                    for i, query_code in enumerate(query_codes):
                        tpl_code = tpl_codes[i]
                        fasta = m.align(query_fasta_file, query_code, tpl_code)
                        sequences.update(self.fastaToSequences(fasta))

                    # print sequences
                    
                    ali_filename = "../Workflow/Alignments/query{}{}{}_{}{}{}.ali".format(
                        query_pdb, query_chain1, query_chain2, tpl_pdb, tpl_chain1, tpl_chain2)

                    with open(ali_filename, 'w') as o:

                            o.write(">P1;"+tpl_pdb +"\n") #+tpl_chain1+tpl_chain2+"\n")
                            # o.write("structureX:{}:{}:{}:{}:{}::::\n".format(
                            #     tpl_pdb+tpl_chain1+tpl_chain2, ".", tpl_chain1, ".", tpl_chain2))

                            o.write("structureX:{}:{}:{}:{}:{}::::\n".format(
                                tpl_pdb, ".", tpl_chain1, ".", tpl_chain2))

                            o.write(fill(sequences[tpl_pdb+tpl_chain1], 60))
                            o.write("\n/\n")
                            o.write(fill(sequences[tpl_pdb+tpl_chain2], 60))
                            o.write("\n*\n")

                            seq_id = "{}{}{}_{}{}{}".format(query_pdb, query_chain1, query_chain2, tpl_pdb, tpl_chain1, tpl_chain2)
                            o.write(">P1;query{}\n".format(seq_id))
                            o.write("sequence:query{}:{}:{}:{}:{}::::\n".format(seq_id, ".", query_chain1, ".", query_chain2))
                            # from 1:A to 100:B
                            o.write(fill(sequences[query_pdb+query_chain1], 60))
                            o.write("\n/\n")
                            o.write(fill(sequences[query_pdb+query_chain1], 60))
                            o.write("\n*\n")

    def models(self, templates):
        model_num = 0
        m = Modeller()
        for query_pdb, chains_complexes in templates.iteritems():
            for (query_chain1, query_chain2), template_pdbs in chains_complexes.iteritems():
                for tpl_pdb, tpl_chain1, tpl_chain2 in template_pdbs:
                    model_num += 1
                    print "* Model building {}: ".format(model_num)
                    print "  - Query complex [{} {} {}]".format(query_pdb, query_chain1, query_chain2)
                    print "  - Template complex : {} {} {}".format(tpl_pdb, tpl_chain1, tpl_chain2)

                    query_code = "{}{}{}_{}{}{}".format(
                        query_pdb, query_chain1, query_chain2, tpl_pdb, tpl_chain1, tpl_chain2)
                    ali_filename = "../Workflow/Alignments/query{}.ali".format(query_code)
                    
                    tpl_code = tpl_pdb # + tpl_chain1 + tpl_chain2
                    # query_code = query_pdb + query_chain1 + query_chain2

                    m.model(ali_filename, tpl_code, query_code)

#############################################################

def workflow():

    skempi = SKEMPI()
    dimers = skempi.getDimers()
    # print dimers

    pdb = PDB()
    pdb.getStructureFiles(dimers.keys())
    pdb.extractSEQRES(dimers.keys())

    query_fasta = "../Workflow/BLAST-results/query_sequence.fasta"
    BLAST_xml = "../Workflow/BLAST-results/BLAST_hits.xml"

    blast = BLAST()
    blast.pirToQuery(dimers, query_fasta)
    xml = blast.runBLASTP(query_fasta, BLAST_xml)
    hits = blast.parseHits(xml)

    c = Complexes()
    tpl_complexes = c.templatesComplexes(dimers, hits)
    template_pdb_codes = c.getTemplatePDBCodes(tpl_complexes)
    print template_pdb_codes
    pdb.getStructureFiles(template_pdb_codes)
    template_seqres = pdb.extractSEQRES(template_pdb_codes)

    c.alignments(tpl_complexes, query_fasta)
    # verify model building
    # c.models(tpl_complexes)

if __name__ == '__main__':
    workflow()
