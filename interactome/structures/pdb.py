"""
PDB module
PDBSearchResults processing module
"""

import httplib
import csv
import glob
import re
from itertools import islice
from textwrap import TextWrapper
from collections import defaultdict, namedtuple
import os
import fnmatch
import gzip

from struct import Struct, Atom


class PDBFile(Struct):

    def __init__(self, PATH, code=None, asm=None):
        Struct.__init__(self, PATH)
        self.pdb_path = PATH + "/biounit"
        self.lines = None
        self.asm = None
        self.code = None
        if code is not None:
            self.load(code, asm)

    def listAll(self):
        for root, dirnames, filenames in os.walk(self.pdb_path):
            for filename in fnmatch.filter(filenames, '*.pdb1.gz'):
                pdb_code = ""
                asm = 0
                # print filename
                if filename.endswith("pdb1.gz"):
                    pdb_code, pdb_asm, _ = os.path.basename(filename).lower().split(".", 2)
                    asm = int(pdb_asm[3:])
                else:
                    continue
                gz_fname = root + "/" + filename
                yield pdb_code, gz_fname

    def getChains(self):
        # if asm is None: asm = 1
        chains = set([a.chain for a in self.iterAtoms()])
        print chains
        return chains

    def load(self, code, asm=None):
        if asm is None:
            asm = 1
        fname = "{}/{}/{}.pdb{}.gz".format(self.pdb_path, code[1:3], code, asm)
        with gzip.open(fname, 'r') as f:
            self.lines = f.readlines()
            self.asm = asm
            self.code = None
            return
        raise Exception("Structure {} could not be loaded from file {}".format(code, fname))

    def iterAtoms(self):
        if self.lines is None:
            raise Exception('Structure not loaded')

        asm = 1
        if self.asm is not None:
            asm = self.asm
        model = 0

        for line in self.lines:
            # print line
            if line.startswith("MODEL"):
                model = int(line.strip().split()[1])

            if line.startswith("ATOM"):
                # print atom
                pdb_atomn = line[13:16].strip()  # 'CA' # questions about strip!!! Calcium short vs CA
                pdb_atomi = int(line[4:12].strip())
                pdb_resn = line[17:20].strip()
                pdb_chain = line[21:22] 
                pdb_resi = line[22:26].strip()
                try:
                    pdb_resi = int(pdb_resi)
                except: pass
                pdb_element = line[76:78].strip()

                pdb_resn_short = self.aa_long_short.get(pdb_resn, None)
                if pdb_resn_short is None: continue

                chain = "{}_{}".format(pdb_chain, model-1) if model > 1 else pdb_chain

                xyz = [float(line[30+8*i:38+8*i]) for i in range(3)]

                atom = Atom(chain=chain, chain_real=pdb_chain, chain_author=pdb_chain, resn=pdb_resn, resi=pdb_resi, resn_short=pdb_resn_short, atomn=pdb_atomn, atomi=pdb_atomi, element=pdb_element, xyz=xyz)
                yield atom

    def fetchRemoteStructure(self, pdb):
        conn = httplib.HTTPConnection("www.rcsb.org", timeout=10)
        print "Downloading PDB {}...".format(pdb)
        conn.request("GET", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="+pdb.lower())
        r1 = conn.getresponse()

        if r1.status == httplib.NOT_FOUND:  # 404
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

            if m is None:
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
