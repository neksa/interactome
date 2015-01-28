
# import os
# import fnmatch
from collections import namedtuple


Atom = namedtuple('Atom', 'chain, chain_real, chain_author, resn, resn_short, resi, atomn, atomi, element, xyz')
ChainProteinMap = namedtuple('ChainProteinMap', 'chain_author, uniprot, seq_aln_begin, seq_aln_begin_ins, seq_aln_end, seq_aln_end_ins, db_aln_begin, db_aln_end, auth_aln_begin, auth_aln_end')


class Struct:
    def __init__(self, PATH):
        self.aa_long_short = {"ALA": 'A', "ASP": 'D', "CYS": 'C', "GLU": 'E',
                             "PHE": 'F', "GLY": 'G', "HIS": 'H', "ILE": 'I',
                             "LYS": 'K', "LEU": 'L', "MET": 'M', "ASN": 'N',
                             "PRO": 'P', "GLN": 'Q', "ARG": 'R',
                             "SER": 'S', "THR": 'T', "VAL": 'V',
                             "TRP": 'W', "TYR": 'Y'}

        # self.aa_long_short = { "ALA": 'A', "ASP": 'D', "CYS": 'C', "GLU": 'E',
        #                      "PHE": 'F', "GLY": 'G', "HIS": 'H', "ILE": 'I',
        #                      "LYS": 'K', "LEU": 'L', "MET": 'M', "ASN": 'N',
        #                      "PYL": 'O', "PRO": 'P', "GLN": 'Q', "ARG": 'R',
        #                      "SER": 'S', "THR": 'T', "SEC": 'U', "VAL": 'V',
        #                      "TRP": 'W', "TYR": 'Y' }

        self.PATH = PATH
        # root = "/Users/agoncear/projects/Interactome/scoring/"
        # self.results_dir = root+"results/"
        # self.structures_dir = "/Users/agoncear/data/PDB/biounit/" # root+"PDB/"

    def getChains(self, code):
        pass

    def iterAtoms(self, code):
        pass

    def listAll(self):
        pass

    # def runAll(self, only=None):
    #     for (pdb_code, gz_fname, fname_pdb, fname_int) in self.dirAll(only):
    #         # print pdb_code, gz_fname, fname_pdb, fname_int
    #         extract_pdb_gz_contacts(pdb_code, gz_fname, fname_pdb, fname_int)
