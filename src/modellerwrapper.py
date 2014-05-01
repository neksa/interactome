"""
    Interface to MODELLER 9.13 software
"""

import xml.etree.ElementTree as ET
from collections import defaultdict, namedtuple
from textwrap import TextWrapper

import httplib
import xml.etree.ElementTree as ET
import pickle
import zlib
import csv
from itertools import islice

from modeller import * #model, environ, alignment
from modeller.automodel import automodel    # Load the automodel class

from cStringIO import StringIO

import glob
import subprocess


class ModellerWrapper:

    def __init__(self, workflow="../Workflow/"):
        # log.verbose()
        self.workflow = workflow
        self.env = environ()
        self.env.io.atom_files_directory = ['.', self.workflow + '/Structures/']


    def writeChainSeq(self, pdb):
        mdl = model(self.env, file="{}.pdb".format(pdb))
        for c in mdl.chains:
            if c.filter(structure_types='structureN structureX'):
                filename = self.workflow + "/Sequences/{}{}.pir".format(pdb, c.name)
                print("Wrote out " + filename)
                atom_file, align_code = c.atom_file_and_code(filename)
                c.write(filename, atom_file, align_code,
                            format='PIR',
                            chop_nonstd_termini=True)

    def alignSequences(self, query_fasta_file, query_code, tpl_code):
        aln = alignment(self.env)
        aln.append(file=query_fasta_file, alignment_format='FASTA', align_codes=(query_code))
        tpl_chain_file = self.workflow + "/Sequences/{}.pir".format(tpl_code)
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
        

    def makeModel(self, ali, template, query):

        a = automodel(self.env,
                    alnfile  = ali ,
                    knowns   = template,
                    sequence = 'query'+query)

        a.starting_model = 1
        a.ending_model   = 1
        a.make()


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


    def dimerAlignments(self, templates, query_fasta_file):
        """
        Returns alignments:
        dictionary, chain sequences
        """
        model_num = 0
        m = Modeller()
        
        wrapper = TextWrapper(break_on_hyphens = False, width = 60)

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
                    ali_filename = self.workflow + "/Alignments/query{}{}{}_{}{}{}.ali".format(
                        query_pdb, query_chain1, query_chain2, tpl_pdb, tpl_chain1, tpl_chain2)

                    with open(ali_filename, 'w') as o:

                            o.write(">P1;"+tpl_pdb +"\n") #+tpl_chain1+tpl_chain2+"\n")
                            # o.write("structureX:{}:{}:{}:{}:{}::::\n".format(
                            #     tpl_pdb+tpl_chain1+tpl_chain2, ".", tpl_chain1, ".", tpl_chain2))

                            o.write("structureX:{}:{}:{}:{}:{}::::\n".format(
                                tpl_pdb, ".", tpl_chain1, ".", tpl_chain2))

                            o.write(wrapper.fill(sequences[tpl_pdb+tpl_chain1]))
                            o.write("\n/\n")
                            o.write(wrapper.fill(sequences[tpl_pdb+tpl_chain2]))
                            o.write("\n*\n")

                            seq_id = "{}{}{}_{}{}{}".format(query_pdb, query_chain1, query_chain2, tpl_pdb, tpl_chain1, tpl_chain2)
                            o.write(">P1;query{}\n".format(seq_id))
                            o.write("sequence:query{}:{}:{}:{}:{}::::\n".format(seq_id, ".", query_chain1, ".", query_chain2))
                            # from 1:A to 100:B
                            o.write(wrapper.fill(sequences[query_pdb+query_chain1]))
                            o.write("\n/\n")
                            o.write(wrapper.fill(sequences[query_pdb+query_chain2]))
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
                    ali_filename = self.workflow + "/Alignments/query{}.ali".format(query_code)
                    
                    tpl_code = tpl_pdb # + tpl_chain1 + tpl_chain2
                    # query_code = query_pdb + query_chain1 + query_chain2

                    m.model(ali_filename, tpl_code, query_code)


