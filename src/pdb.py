"""
PDB module
PDBSearchResults processing module
"""

import httplib
import csv
import glob
from itertools import islice
from textwrap import TextWrapper
from collections import defaultdict, namedtuple

# from cStringIO import StringIO
# import subprocess
# import xml.etree.ElementTree as ET
# import pickle
# import zlib

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


############################################

class PDBSearchResults:
    def getTabularResults(self, fname):
        records = []
        PDBSearchRecord = namedtuple('PDBSearchRecord', 'pdb, chain, resolution, source, tax, entity, sequence, chain_len, db_id')
        with open(fname) as f:
            reader = csv.reader(islice(f, 1, None), delimiter=',', quotechar='"')
            for r in reader:
                if len(r) > 0:
                    rec = PDBSearchRecord._make(r)
                    # yield rec
                    records.append(rec)
            return records    
            # for rec in map(PDBSearchRecord._make, reader):
            #     print rec.pdb, rec.chain, rec.resolution
            #     print rec.sequence
            #     print
            #     records.append(rec)

    def filterPDB(self, records):
        excluded_pdb = []
        min_chain_len = 50
        tax_id = 9606
        for rec in records:
            try:
                if int(rec.tax) != tax_id:
                    excluded_pdb.append(rec.pdb)
                    continue

                if int(rec.chain_len) < min_chain_len:
                    excluded_pdb.append(rec.pdb)
                    continue
            except:
                excluded_pdb.append(rec.pdb)
                continue                

        excluded_pdb = set(excluded_pdb)
        new_records = [rec for rec in records if rec.pdb not in excluded_pdb]
        return new_records


    def searchToFasta(self, records, fasta_file):
        wrapper = TextWrapper(break_on_hyphens = False, width = 60)
        with open(fasta_file, 'w') as o:
            for r in records:
                o.write(">{}{}\n".format(r.pdb.upper(), r.chain.upper()))
                o.write("{}\n".format(wrapper.fill(r.sequence)))


#####################################################

class SKEMPI:
    def getDimers(self, fname="../SKEMPI/dimers.csv"):
        dimers = defaultdict(list)
        with open(fname) as f:
            for line in f:
                pdb, chain1, chain2 = line.strip().split("_", 2)
                dimers[pdb] = (chain1, chain2)
        return dimers

