"""
mmCIF loading using parsers from pdbx
"""

import copy
import gzip
from pdbx.reader.PdbxReader import PdbxReader
from pdbx.reader.PdbxContainers import *
from sys import argv, exit
from os.path import splitext

from mmcifOperation import *

# Amino acid FASTA codification dictionary
codification = { "ALA" : 'A',
                 "ASP" : 'D',
                 "CYS" : 'C',
                 "GLU" : 'E',
                 "PHE" : 'F',
                 "GLY" : 'G',
                 "HIS" : 'H',
                 "ILE" : 'I',
                 "LYS" : 'K',
                 "LEU" : 'L',
                 "MET" : 'M',
                 "ASN" : 'N',
                 "PYL" : 'O',
                 "PRO" : 'P',
                 "GLN" : 'Q',
                 "ARG" : 'R',
                 "SER" : 'S',
                 "THR" : 'T',
                 "SEC" : 'U',
                 "VAL" : 'V',
                 "TRP" : 'W',
                 "TYR" : 'Y' }


data = []
# with gzip.open("/Users/agoncear/data/PDB/mmCIF/17/117e.cif.gz") as cif:
# with gzip.open("/Users/agoncear/data/PDB/mmCIF/j5/2j5q.cif.gz") as cif:
# with gzip.open("/Users/agoncear/data/PDB/mmCIF/fy/1fyh.cif.gz") as cif:
# with gzip.open("/Users/agoncear/data/PDB/mmCIF/f9/4f9c.cif.gz") as cif:
with gzip.open("/Users/agoncear/data/PDB/mmCIF/ci/4cib.cif.gz") as cif:
    pRd = PdbxReader(cif)
    # Call the read(self, containerList) method with your list
    pRd.read(data)

# Your list is now propagated with one or more DataContainer objects, which represent data blocks. To get the first data block, just use list notation: 
block = data[0]
# print block


print "Number of blocks", len(data)
asm = block.getObj("pdbx_struct_assembly")
print "Assemblies", asm.getRowCount()
for i in range(asm.getRowCount()):
    id = asm.getValue("id", i)
    details = asm.getValue("details", i)
    method = asm.getValue("method_details", i)
    oligo_count = asm.getValue("oligomeric_count", i)
    oligo_details = asm.getValue("oligomeric_details", i)
    print id, details, method, oligo_count, oligo_details


atom_site = block.getObj("atom_site")
atom_site_ref = copy.copy(atom_site)
attributes = atom_site_ref.getAttributeList()
atom_site = DataCategory("atom_site", attributes)

oper_list = block.getObj("pdbx_struct_oper_list")
oper = []
oper2 = []
atomNum = 1
modelNum = 0

assembly_gen = block.getObj("pdbx_struct_assembly_gen")
print "GEN", assembly_gen.getRowCount()
asm_id = 0
assemblyId = assembly_gen.getValue("assembly_id", asm_id)
print "Assembly generator:", assemblyId
oper_expression = assembly_gen.getValue("oper_expression", asm_id)

print "Oper expression", oper_expression
parenCount = oper_expression.count("(")
# Handles one operation assemblies (e.g., "1")
if parenCount == 0 : oper.append(oper_expression)
# Handles multiple operation assemblies, no Cartesian products (e.g., "(1-5)")
if parenCount == 1 : oper.extend(parseOperationExpression(oper_expression))
# Handles Cartesian product expressions (e.g., "(X0)(1-60)")
if parenCount == 2 :
    # Break the expression into two parenthesized expressions and parse them
    temp = oper_expression.find(")")
    oper.extend(parseOperationExpression(oper_expression[0:temp+1]))
    oper2.extend(parseOperationExpression(oper_expression[temp+1:]))
asym_id_list = assembly_gen.getValue("asym_id_list", asm_id)
temp = (1 > len(oper2)) and 1 or len(oper2)

#FIXME:
oper = ["1", "2"]

print oper
print oper2
print oper_list
print "Asym id list", asym_id_list

# For every operation in the first parenthesized list
for op1 in oper :
    print "OPER", op1
    # Find the index of the current operation in the oper_list category table
    op1index = 0
    for row in range(oper_list.getRowCount()) : 
        if oper_list.getValue("id", row) == op1 : 
            op1index = row
            break

    # For every operation in the second parenthesized list (if there is one)
    print "TEMP", temp
    for i in range(temp) :      
        # Find the index of the second operation in the oper_list category table
        op2index = -1
        if (oper2) :
            for row in range(oper_list.getRowCount()) :
                if oper_list.getValue("id", row) == oper2[i] :
                    op2index = row
                    break

        # Prepare the operation matrix
        operation = prepareOperation(oper_list, op1index, op2index)

        # Iterate over every atom in the atom_site reference table
        for r in range(atom_site_ref.getRowCount()) :
            
            # If the asym_id of the current atom is not in the asym_id list, skip to the next atom
            if (asym_id_list.find(atom_site_ref.getValue("label_asym_id", r)) == -1) :
                continue
            
            # print "Processing asym_id", atom_site_ref.getValue("label_asym_id", r)
            # Retrieve the atom's row from the atom_site reference table
            atom = atom_site_ref.getFullRow(r)

            # print r, atomNum -1, atom
            # Add this row to the atom_site table for this assembly
            for s in range(len(attributes)) :
                atom_site.setValue(atom[s], attributes[s], atomNum - 1)

            # Update the atom number and model number for this row
            atom_site.setValue(str(atomNum), "id", atomNum - 1)
            atom_site.setValue(str(modelNum), "pdbx_PDB_model_num", atomNum - 1) 

            # Determine and set the new coordinates
            coords = [float(atom[10]), float(atom[11]), float(atom[12]), 1.0]
            sum = 0.0
            xyz = ['x', 'y', 'z']
            for a in range(3) :
                sum = 0.0
                for b in range(4) :
                    sum += (operation[a][b] * coords[b])
                atom_site.setValue("%.3f" % sum, "Cartn_" + xyz[a], atomNum - 1)
            atomNum += 1
        modelNum += 1







# Retrieve the entity category table, which contains information that will be used in the FASTA header1
entity = block.getObj("entity")
# Holds non-mandatory entity attributes that could serve as FASTA header lines, ordered preferentially
candidates = ["pdbx_description", "details", "type"]
headerDescriptor = ""
# Set the header descriptor
for c in candidates :
    if entity.hasAttribute(c) :
        headerDescriptor = c
        break
# If none of the optional descriptors are present, just use the entity id
if not headerDescriptor: headerDescriptor = "id"

"""
entity_poly = block.getObj("entity_poly")
print entity_poly.getValue("pdbx_seq_one_letter_code_can")
# for i in range(entity_poly.getRowCount()) :
#     code = entity_poly.getValue("pdbx_seq_one_letter_code_can", i)
#     print code

# Retrieve the entity_poly_seq category table, which containers the monomer sequences for entities2
entity_poly_seq = block.getObj("entity_poly_seq")
# Iterate over every row in entity_poly_seq, each containing an entity monomer
for i in range(entity_poly_seq.getRowCount()) :    

    # Retrieve the monomer stored in this row
    ent_id = entity_poly_seq.getValue("entity_id", i)
    num = entity_poly_seq.getValue("num", i)
    monomer = entity_poly_seq.getValue("mon_id", i)
    # print ent_id, num, monomer
"""


"""
# Retrieve the struct_site_gen category table, which delineates members of structurally relevant sites1
print "Relevant sites:"
struct_site_gen = block.getObj("struct_site_gen")
# Iterate over the rows in struct_site_gen, where each row delineates a member of a structurally relevant site
for i in range(struct_site_gen.getRowCount()):
    # full = struct_site_gen.getFullRow(i)
    # print full 

    site = struct_site_gen.getValue("site_id", i)
    seq_id = struct_site_gen.getValue("auth_seq_id", i)
    comp_id = struct_site_gen.getValue("auth_comp_id", i)
    asym_id = struct_site_gen.getValue("auth_asym_id", i)
    print "Site", site, seq_id, comp_id, asym_id
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
"""
for i in range(atom_site.getRowCount()):
    row = atom_site.getFullRow(i)
    if atom_site.getValue("group_PDB", i) != "ATOM": continue
    print "\t".join(row)


