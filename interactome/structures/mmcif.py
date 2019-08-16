"""
mmCIF loading using parsers from pdbx
"""
from __future__ import print_function

import os
import fnmatch
import copy
import gzip
from collections import namedtuple, defaultdict

from pdbx.reader.PdbxReader import PdbxReader
from pdbx.reader.PdbxContainers import *
from .struct import Struct, Atom, ChainProteinMap

mmCIFMeta = namedtuple('Metadata', 'id, method, title, description')


class mmCifFile(Struct):

    def __init__(self, PATH, code=None, asm=None):
        Struct.__init__(self, PATH)
        self.mmcif_path = PATH + "/mmcif"
        self.asm = None
        self.code = None
        self.block = []
        if code is not None:
            self.load(code, asm)

    def getChains(self):
        chains = set([a.chain for a in self.iterAtoms()])
        # print(chains)
        return chains

    def getPDBChainMapping(self):
        """From PDB mapping to CIF mapping"""
        # if asm is None: asm = 1
        mapping = {}
        for a in self.iterAtoms():
            # mapping[a.chain_real] = a.chain_author
            # print a.chain_author,
            mapping[a.chain_author] = a.chain_real
        return mapping

    def load(self, code, asm=None):
        """ Read the main CIF data block """
        filename = "{}/{}/{}.cif.gz".format(self.mmcif_path, code[1:3], code)
        # data container
        data = []
        self.block = None
        self.code = None
        self.asm = None
        with gzip.open(filename) as cif:
            pRd = PdbxReader(cif)
            pRd.read(data)
            self.block = data[0]  # the first container
            self.code = code
            self.asm = asm
            return
        raise Exception("Could not load CIF structure {} from file {}".format(code, filename))

    def getChainProteinMapping(self):
        """
        loop_
        _struct_ref_seq.align_id
        _struct_ref_seq.ref_id
        _struct_ref_seq.pdbx_PDB_id_code
        _struct_ref_seq.pdbx_strand_id
        _struct_ref_seq.seq_align_beg
        _struct_ref_seq.pdbx_seq_align_beg_ins_code
        _struct_ref_seq.seq_align_end
        _struct_ref_seq.pdbx_seq_align_end_ins_code
        _struct_ref_seq.pdbx_db_accession
        _struct_ref_seq.db_align_beg
        _struct_ref_seq.db_align_end
        _struct_ref_seq.pdbx_auth_seq_align_beg
        _struct_ref_seq.pdbx_auth_seq_align_end
        1  1 3AZM A 4 ? 139 ? P68431 1   136 0   135
        2  2 3AZM B 4 ? 106 ? P62805 1   103 0   102
        3  3 3AZM C 4 ? 133 ? P04908 1   130 0   129
        4  4 3AZM D 4 ? 129 ? P06899 1   126 0   125
        5  1 3AZM E 4 ? 139 ? P68431 1   136 0   135
        6  2 3AZM F 4 ? 106 ? P62805 1   103 0   102
        7  3 3AZM G 4 ? 133 ? P04908 1   130 0   129
        8  4 3AZM H 4 ? 129 ? P06899 1   126 0   125
        9  5 3AZM I 1 ? 146 ? 3AZM   1   146 1   146
        10 5 3AZM J 1 ? 146 ? 3AZM   147 292 147 292
        """

        """
        _struct_ref_seq.align_id
        _struct_ref_seq.ref_id
        _struct_ref_seq.pdbx_PDB_id_code
        _struct_ref_seq.pdbx_strand_id
        _struct_ref_seq.seq_align_beg
        _struct_ref_seq.pdbx_seq_align_beg_ins_code
        _struct_ref_seq.seq_align_end
        _struct_ref_seq.pdbx_seq_align_end_ins_code
        _struct_ref_seq.pdbx_db_accession
        _struct_ref_seq.db_align_beg
        _struct_ref_seq.db_align_end
        _struct_ref_seq.pdbx_auth_seq_align_beg
        _struct_ref_seq.pdbx_auth_seq_align_end
        1  1 3B5N A 2 ? 61 ? P31109 27  86  27  86
        2  2 3B5N B 1 ? 69 ? P32867 189 257 189 257
        3  3 3B5N C 3 ? 69 ? P40357 433 499 433 499
        4  4 3B5N D 3 ? 64 ? P40357 589 650 589 650
        5  1 3B5N E 2 ? 61 ? P31109 27  86  27  86
        6  2 3B5N F 1 ? 69 ? P32867 189 257 189 257
        7  3 3B5N G 3 ? 69 ? P40357 433 499 433 499
        8  4 3B5N H 3 ? 64 ? P40357 589 650 589 650
        9  1 3B5N I 2 ? 61 ? P31109 27  86  27  86
        10 2 3B5N J 1 ? 69 ? P32867 189 257 189 257
        11 3 3B5N K 3 ? 69 ? P40357 433 499 433 499
        12 4 3B5N L 3 ? 64 ? P40357 589 650 589 650
        """

        # This one is a very complex case with insertions/deletions:
        # http://www.pdb.org/pdb/explore/remediatedSequence.do?params.desiredSequenceRnsStr=SEQRES&params.showDbRefRulerIfPossible=false&structureId=1UX4
        # Resiues 1412 ..5 residues missing.. 1417

        """
        loop_
        _struct_ref_seq.align_id
        _struct_ref_seq.ref_id
        _struct_ref_seq.pdbx_PDB_id_code
        _struct_ref_seq.pdbx_strand_id
        _struct_ref_seq.seq_align_beg
        _struct_ref_seq.pdbx_seq_align_beg_ins_code
        _struct_ref_seq.seq_align_end
        _struct_ref_seq.pdbx_seq_align_end_ins_code
        _struct_ref_seq.pdbx_db_accession
        _struct_ref_seq.db_align_beg
        _struct_ref_seq.pdbx_db_align_beg_ins_code
        _struct_ref_seq.db_align_end
        _struct_ref_seq.pdbx_db_align_end_ins_code
        _struct_ref_seq.pdbx_auth_seq_align_beg
        _struct_ref_seq.pdbx_auth_seq_align_end
        1 1 1UX4 A 1  ? 61  ? P41832 1352 ? 1412 ? 1352 1412
        2 1 1UX4 A 62 ? 410 ? P41832 1417 ? 1765 ? 1417 1765
        3 1 1UX4 B 1  ? 61  ? P41832 1352 ? 1412 ? 2352 2412
        4 1 1UX4 B 62 ? 410 ? P41832 1417 ? 1765 ? 2417 2765
        """

        if self.block is None:
            raise Exception("CIF structure not loaded")
        # chain_protein_mapping = defaultdict(set)
        chain_protein_mapping = defaultdict(list)
        pdb_chain_mapping = self.getPDBChainMapping()
        # print pdb_chain_mapping

        ref = self.block.getObj("struct_ref_seq")
        for i in range(ref.getRowCount()):
            # ref_id = ref.getValue("ref_id", i)
            chain_author = ref.getValue("pdbx_strand_id", i)
            chain = pdb_chain_mapping.get(chain_author, None)
            if chain is None:
                # this automatically removes non-protein chain, because they don't have Atom records in this version of iterAtoms()
                continue
            # print "AUTH:", chain_author, chain
            uniprot = ref.getValue("pdbx_db_accession", i)
            if uniprot == self.code.upper():  # when does this happen really? need examples
                uniprot = ''

            # _struct_ref_seq.align_id
            # _struct_ref_seq.ref_id
            # _struct_ref_seq.pdbx_PDB_id_code
            # _struct_ref_seq.pdbx_strand_id
            # _struct_ref_seq.seq_align_beg
            # _struct_ref_seq.pdbx_seq_align_beg_ins_code
            # _struct_ref_seq.seq_align_end
            # _struct_ref_seq.pdbx_seq_align_end_ins_code
            # _struct_ref_seq.pdbx_db_accession
            # _struct_ref_seq.db_align_beg
            # _struct_ref_seq.db_align_end
            # _struct_ref_seq.pdbx_auth_seq_align_beg
            # _struct_ref_seq.pdbx_auth_seq_align_end

            seq_aln_begin = int(ref.getValue("seq_align_beg", i))
            seq_aln_begin_ins = ref.getValue("pdbx_seq_align_beg_ins_code", i)
            seq_aln_end = int(ref.getValue("seq_align_end", i))
            seq_aln_end_ins = ref.getValue("pdbx_seq_align_end_ins_code", i)
            db_aln_begin = int(ref.getValue("db_align_beg", i))
            db_aln_end = int(ref.getValue("db_align_end", i))
            auth_aln_begin = int(ref.getValue("pdbx_auth_seq_align_beg", i))
            auth_aln_end = int(ref.getValue("pdbx_auth_seq_align_end", i))
            chain_protein_mapping[chain].append(ChainProteinMap(
                chain_author=chain_author, uniprot=uniprot,
                seq_aln_begin=seq_aln_begin, seq_aln_begin_ins=seq_aln_begin_ins, seq_aln_end=seq_aln_end, seq_aln_end_ins=seq_aln_end_ins,
                db_aln_begin=db_aln_begin, db_aln_end=db_aln_end,
                auth_aln_begin=auth_aln_begin, auth_aln_end=auth_aln_end))
        return chain_protein_mapping

    def getMetaData(self):
        """
        _exptl.entry_id   2AFF
        _exptl.method     'SOLUTION NMR'
        #
        _struct.entry_id                  2AFF
        _struct.title                     'The solution structure of the Ki67FHA/hNIFK(226-269)3P complex'
        _struct.pdbx_descriptor           'Antigen KI-67/MKI67 FHA domain interacting nucleolar phosphoprotein'
        _struct.pdbx_model_details        ?
        _struct.pdbx_CASP_flag            ?
        _struct.pdbx_model_type_details   ?
        #

        ELECTRON CRYSTALLOGRAPHY    .
        ELECTRON MICROSCOPY .
        EPR EPR only as a supporting method
        FIBER DIFFRACTION   .
        FLUORESCENCE TRANSFER   FLUORESCENCE TRANSFER only as a supporting method
        INFRARED SPECTROSCOPY   IR and FTIR only as supporting methods
        NEUTRON DIFFRACTION .
        POWDER DIFFRACTION  .
        SOLID-STATE NMR .
        SOLUTION NMR    .
        SOLUTION SCATTERING .
        THEORETICAL MODEL   THEORETICAL MODEL only as a supporting method
        X-RAY DIFFRACTION   .
        """
        if self.block is None:
            raise Exception("CIF structure not loaded")

        ref = self.block.getObj("exptl")
        entry_id = ref.getValue("entry_id", 0)
        method = ref.getValue("method", 0)

        ref = self.block.getObj("struct")
        title = ref.getValue("title", 0)
        description = ref.getValue("pdbx_descriptor", 0)
        return mmCIFMeta(id=entry_id, method=method, title=title, description=description)

    def mapChainProtein(self, chain_protein_mapping, chain):
        chain = chain.split("_")[0]
        protein = chain_protein_mapping.get(chain)
        return protein

    def listAll(self):
        for root, dirnames, filenames in os.walk(self.mmcif_path):
            for filename in fnmatch.filter(filenames, '*.cif.gz'):
                pdb_code = ""
                # asm = 0
                # print filename
                if filename.endswith("cif.gz"):
                    pdb_code, pdb_asm, _ = os.path.basename(filename).lower().split(".", 2)
                else:
                    continue
                gz_fname = root + "/" + filename
                yield pdb_code, gz_fname

    def iterAtoms(self):
        """
        Applies crystallographic symmetry for X-RAY structures - saves the generated chains in models
        For NMR structures returns only atoms in the first model
        """
        if self.block is None:
            raise Exception("CIF structure not loaded")

        debug = False

        asm_id = self.asm
        if asm_id is None:
            asm_id = 1

        # meta = self.getMetaData()
        only_models = [1]
        # if meta.method.endswith("NMR"):
        #     only_models = [1]

        # with gzip.open("/Users/agoncear/data/PDB/mmCIF/17/117e.cif.gz") as cif:
        # with gzip.open("/Users/agoncear/data/PDB/mmCIF/j5/2j5q.cif.gz") as cif:
        # with gzip.open("/Users/agoncear/data/PDB/mmCIF/fy/1fyh.cif.gz") as cif:
        # with gzip.open("/Users/agoncear/data/PDB/mmCIF/f9/4f9c.cif.gz") as cif:
        # Amino acid FASTA codification dictionary

        # NMR:
        # _atom_site.pdbx_PDB_model_num

        # asm_obj = self.block.getObj("pdbx_struct_assembly")
        # print "Assemblies", asm.getRowCount()

        # for i in range(asm_obj.getRowCount()):
        #     id = asm_obj.getValue("id", i)
        #     details = asm_obj.getValue("details", i)
        #     method = asm_obj.getValue("method_details", i)
        #     oligo_count = asm_obj.getValue("oligomeric_count", i)
        #     oligo_details = asm_obj.getValue("oligomeric_details", i)
        #     # print id, details, method, oligo_count, oligo_details

        atom_site_obj = self.block.getObj("atom_site")   # atoms
        atom_site_obj_ref = copy.copy(atom_site_obj)  # copy
        attributes = atom_site_obj_ref.getAttributeList()  # copy attributes
        atom_site_obj = DataCategory("atom_site", attributes)  # new atoms with old attributes

        oper_list = self.block.getObj("pdbx_struct_oper_list")
        atomNum = 1
        modelNum = 1

        # Generate chains in the assembly
        assembly_gen = self.block.getObj("pdbx_struct_assembly_gen")
        n_gen = assembly_gen.getRowCount()
        if debug:
            print("NGEN", n_gen)

        for gen in range(n_gen):
            # asm_id is the main parameter! copied from self.asm
            assemblyId = assembly_gen.getValue("assembly_id", gen)
            if debug:
                print("ASM???", assemblyId, asm_id)
            if int(assemblyId) != asm_id:
                continue

            # Operator expression. We need to parse it
            oper_expression = assembly_gen.getValue("oper_expression", gen)
            if debug:
                print("Oper expression", oper_expression)

            # how many parentheses
            paren_count = oper_expression.count("(")

            oper1 = []
            oper2 = []
            # Handles one operation assemblies (e.g., "1")
            # if paren_count == 0: oper1.append(oper_expression)
            # Handles multiple operation assemblies, no Cartesian products (e.g., "(1-5)")
            if paren_count <= 1:
                oper1.extend(self.parseOperationExpression(oper_expression))
            # Handles Cartesian product expressions (e.g., "(X0)(1-60)")
            if paren_count == 2:
                # Break the expression into two parenthesized expressions and parse them
                first_paren_end = oper_expression.find(")")
                oper1.extend(self.parseOperationExpression(oper_expression[0:first_paren_end+1]))
                oper2.extend(self.parseOperationExpression(oper_expression[first_paren_end+1:]))

            asym_id_list = assembly_gen.getValue("asym_id_list", gen)

            if debug:
                print("Asym id list", asym_id_list)
                print("Oper1", oper1)
                print("Oper2", oper2)

            # For every operation in the first parenthesized list
            for op1 in oper1:
                # Find the index of the current operation in the oper_list category table
                op1_index = 0
                for row in range(oper_list.getRowCount()):
                    if oper_list.getValue("id", row) == op1:
                        op1_index = row
                        break

                # For every operation in the second parenthesized list. There must be one cycle even if the list is empty
                for i in range(1 if (len(oper2) < 1) else len(oper2)):
                    # Find the index of the second operation in the oper_list category table
                    op2_index = -1
                    if oper2:
                        for row in range(oper_list.getRowCount()):
                            if oper_list.getValue("id", row) == oper2[i]:
                                op2_index = row
                                break

                    # Prepare the operation matrix
                    operation = self.prepareOperation(oper_list, op1_index, op2_index)
                    if debug:
                        print("OPERATION", operation)

                    # Iterate over every atom in the atom_site reference table
                    for r in range(atom_site_obj_ref.getRowCount()):
                        # If the asym_id of the current atom is not in the asym_id list, skip to the next atom
                        if asym_id_list.find(atom_site_obj_ref.getValue("label_asym_id", r)) == -1:
                            if debug:
                                print("-" + atom_site_obj_ref.getValue("label_asym_id", r), end='')
                            continue

                        # Only the first model of the original structure.
                        # We do not need all atoms in all models, only the first one
                        model = int(atom_site_obj_ref.getValue("pdbx_PDB_model_num", r))
                        if only_models is not None:
                            if model not in only_models:
                                continue

                        # print "Processing asym_id", atom_site_ref.getValue("label_asym_id", r)
                        # Retrieve the atom's row from the atom_site reference table
                        atom = atom_site_obj_ref.getFullRow(r)

                        # print r, atomNum -1, atom
                        # Add this row to the atom_site table for this assembly
                        for s in range(len(attributes)):
                            atom_site_obj.setValue(atom[s], attributes[s], atomNum - 1)

                        # Update the atom number and model number for this row
                        atom_site_obj.setValue(str(atomNum), "id", atomNum - 1)
                        atom_site_obj.setValue(str(modelNum), "pdbx_PDB_model_num", atomNum - 1)

                        # Determine and set the new coordinates
                        coords = [float(atom[10]), float(atom[11]), float(atom[12]), 1.0]
                        for a, c in enumerate("xyz"):
                            total = sum([operation[a][b] * coords[b] for b in range(4)])
                            atom_site_obj.setValue("%.3f" % total, "Cartn_" + c, atomNum - 1)
                        atomNum += 1
                    modelNum += 1

        # Retrieve the entity category table, which contains information that will be used in the FASTA header1
        entity = self.block.getObj("entity")
        # Holds non-mandatory entity attributes that could serve as FASTA header lines, ordered preferentially
        candidates = ["pdbx_description", "details", "type"]
        headerDescriptor = ""
        # Set the header descriptor
        for c in candidates :
            if entity.hasAttribute(c) :
                headerDescriptor = c
                break
        # If none of the optional descriptors are present, just use the entity id
        if not headerDescriptor:
            headerDescriptor = "id"

        """
        entity_poly = self.block.getObj("entity_poly")
        print entity_poly.getValue("pdbx_seq_one_letter_code_can")
        # for i in range(entity_poly.getRowCount()) :
        #     code = entity_poly.getValue("pdbx_seq_one_letter_code_can", i)
        #     print code

        # Retrieve the entity_poly_seq category table, which containers the monomer sequences for entities2
        entity_poly_seq = self.block.getObj("entity_poly_seq")
        # Iterate over every row in entity_poly_seq, each containing an entity monomer
        for i in range(entity_poly_seq.getRowCount()):

            # Retrieve the monomer stored in this row
            ent_id = entity_poly_seq.getValue("entity_id", i)
            num = entity_poly_seq.getValue("num", i)
            monomer = entity_poly_seq.getValue("mon_id", i)
            # print ent_id, num, monomer
        """

        """
        _atom_site.group_PDB
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_seq_id
        _atom_site.label_alt_id
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        _atom_site.occupancy
        _atom_site.B_iso_or_equiv
        _atom_site.footnote_id
        _atom_site.auth_seq_id
        _atom_site.id
        _atom_site.pdbx_PDB_model_num
        """
        for i in range(atom_site_obj.getRowCount()):
            # row = atom_site_obj.getFullRow(i)
            # print "\t".join(row)
            if atom_site_obj.getValue("group_PDB", i) != "ATOM":
                continue

            atomn = atom_site_obj.getValue("label_atom_id", i)
            element = atom_site_obj.getValue("type_symbol", i)
            atomi = int(atom_site_obj.getValue("id", i))
            resn = atom_site_obj.getValue("label_comp_id", i)  # auth_comp_id
            resi = int(atom_site_obj.getValue("label_seq_id", i))
            chain_real = atom_site_obj.getValue("label_asym_id", i)  # auth_asym_id, label_entity_id
            chain_author = atom_site_obj.getValue("auth_asym_id", i)

            resn_short = self.aa_long_short.get(resn, None)

            # do not report atoms in unknown amino acids
            # i.e. nucleic acids, ions, and ligands will not be reported
            if resn_short is None:
                continue

            model = int(atom_site_obj.getValue("pdbx_PDB_model_num", i))
            chain = "{}_{}".format(chain_real, model - 1) if model > 1 else chain_real
            xyz = map(float, [atom_site_obj.getValue("Cartn_" + c, i) for c in "xyz"])

            atom = Atom(chain=chain, chain_real=chain_real, chain_author=chain_author, resn=resn, resn_short=resn_short, resi=resi, atomn=atomn, atomi=atomi, element=element, xyz=xyz)
            yield atom

    def parseOperationExpression(self, expression):
        operations = []
        if expression != "":
            for subexpression in expression.split(")")[0].split(","):
                a = map(int, subexpression.split("-"))
                if len(a) == 2:
                    operations.extend(map(str, range(a[0], a[1] + 1)))
                else:
                    if a[0] != "":
                        operations.append(str(a[0]))
        return operations

    def prepareOperation(self, oper_list, op1index, op2index):
        """
        Code from:
        http://mmcif.wwpdb.org/docs/sw-examples/python/html/assemblies.html
        """
        # Prepare matrices for operations 1 & 2
        op1 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
        op2 = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]

        # Fill the operation matrices for operations 1 & 2
        for i in range(3):
            op1[i][3] = float(oper_list.getValue("vector[" + str(i+1) + "]", op1index))
            if op2index != -1:
                op2[i][3] = float(oper_list.getValue("vector[" + str(i+1) + "]", op2index))
            for j in range(3):
                op1[i][j] = float(oper_list.getValue("matrix[" + str(i+1) + "][" + str(j+1) + "]", op1index))
                if op2index != -1:
                    op2[i][j] = float(oper_list.getValue("matrix[" + str(i+1) + "][" + str(j+1) + "]", op2index))

        # Handles non-Cartesian product expressions
        if op2index == -1:
            return op1

        # Handles Cartesian product expressions (4x4 matrix multiplication)
        operation = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]
        for row in range(4):
            for col in range(4):
                operation[row][col] = sum([op1[row][r] * op2[r][col] for r in range(4)])
        return operation
