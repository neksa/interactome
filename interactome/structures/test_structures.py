

import unittest
import pytest


# def pytest_generate_tests(metafunc):
#     # called once per each test function
#     funcarglist = metafunc.cls.params[metafunc.function.__name__]
#     argnames = list(funcarglist[0])
#     metafunc.parametrize(argnames, [[funcargs[name] for name in argnames]
#             for funcargs in funcarglist])

@pytest.fixture
def PATH():
    """ Path to PDB data with subdirectories biounits/ and mmcif/ """
    return "/Users/agoncear/data/pdb"

@pytest.fixture
def MMCIF(PATH):
    """ mmcif instance """
    from mmcif import mmCifFile
    return mmCifFile(PATH)

@pytest.fixture
def PDB(PATH):
    """ pdb instance """
    from pdb import PDBFile
    return PDBFile(PATH)

############################################################################
def test_mmcif_operations(MMCIF):
    assert len(MMCIF.parseOperationExpression("")) == 0
    assert len(MMCIF.parseOperationExpression("1")) == 1
    assert len(MMCIF.parseOperationExpression("1,2,3")) == 3
    assert len(MMCIF.parseOperationExpression("1-4")) == 4
    assert len(MMCIF.parseOperationExpression("1-3,8,10-12,15")) == 8



def test_protein_mapping_cif_2xb2(MMCIF):
    """
    # 
    loop_
    _pdbx_struct_assembly.id 
    _pdbx_struct_assembly.details 
    _pdbx_struct_assembly.method_details 
    _pdbx_struct_assembly.oligomeric_details 
    _pdbx_struct_assembly.oligomeric_count 
    1 author_and_software_defined_assembly PISA hexameric 6 
    2 author_and_software_defined_assembly PISA hexameric 6 
    3 author_and_software_defined_assembly PISA monomeric 1 
    # 
    loop_
    _pdbx_struct_assembly_gen.assembly_id 
    _pdbx_struct_assembly_gen.oper_expression 
    _pdbx_struct_assembly_gen.asym_id_list 
    1 1 K,P,Q,G,L,M,H,J 
    2 1 A,N,O,D,B,C,I,F 
    3 1 E               
    # 

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
    1  1 2XB2 A 1 ? 411 ? P38919 1   ? 411 ? 1   411 
    2  1 2XB2 X 1 ? 411 ? P38919 1   ? 411 ? 1   411 
    3  2 2XB2 T 1 ? 150 ? O15234 137 ? 286 ? 137 286 
    4  2 2XB2 S 1 ? 150 ? O15234 137 ? 286 ? 137 286 
    5  3 2XB2 C 1 ? 146 ? P61326 1   ? 146 ? 1   146 
    6  3 2XB2 Y 1 ? 146 ? P61326 1   ? 146 ? 1   146 
    7  4 2XB2 D 1 ? 90  ? Q9Y5S9 66  ? 155 ? 66  155 
    8  4 2XB2 Z 1 ? 90  ? Q9Y5S9 66  ? 155 ? 66  155 
    9  5 2XB2 G 1 ? 60  ? Q9BZI7 411 ? 470 ? 411 470 
    10 5 2XB2 U 1 ? 60  ? Q9BZI7 411 ? 470 ? 411 470 
    11 6 2XB2 F 1 ? 3   ? 2XB2   1   ? 3   ? 1   3   
    12 7 2XB2 E 1 ? 15  ? 2XB2   1   ? 15  ? 1   15  
    13 7 2XB2 R 1 ? 15  ? 2XB2   1   ? 15  ? 1   15 
    """
    MMCIF.load("2xb2")
    mapping = MMCIF.getPDBChainMapping()
    # print mapping

    # K,P,Q,G,L,M,H,J    
    # X ? ? ? Y Z S U 
    assert len(mapping.keys()) == 5
    assert mapping['Y'] == 'L'
    assert mapping['X'] == 'K'
    assert mapping['Z'] == 'M'
    assert mapping['S'] == 'H'
    assert mapping['U'] == 'J'
    # mapping = MMCIF.getChainProteinMapping()
    # assert len(mapping.keys()) == 


def test_protein_mapping_cif_2qqk(MMCIF):
    """
    Assembly 1 - software - hexamer
    A + A_1
    H_0 + H_1
    L_0 + L_1

    A - A
    B - H
    C - L

    loop_
    _pdbx_struct_assembly_gen.assembly_id 
    _pdbx_struct_assembly_gen.oper_expression 
    _pdbx_struct_assembly_gen.asym_id_list 
    1 1,2 A,D,E     
    1 3,4 B,C       
    2 1   A,D,E,B,C 
    #
    """
    MMCIF.load("2qqk")
    mapping = MMCIF.getPDBChainMapping()
    assert len(mapping.keys()) == 3

    chains = MMCIF.getChains()
    assert len(chains) == 6
    
    # mapping = MMCIF.getChainProteinMapping()
    # assert len(mapping.keys()) == 


def test_atom_pdb_4tim(PDB):
    ca = False
    natoms = 0
    PDB.load("4tim")
    for i, atom in enumerate(PDB.iterAtoms()):
        if i == 0:
            assert atom.resn == "SER"
            assert atom.resn_short == "S"
            assert atom.resi == 2
            assert atom.chain == "A"
        if atom.chain == "A" and atom.resi == 2:
            natoms += 1
            if atom.atomn == "CA" and atom.element == "C":
                ca = True    
    assert natoms == 6
    assert ca == True


def test_atom_mmcif_4tim(MMCIF):
    ca = False
    natoms = 0
    MMCIF.load("4tim")
    for i, atom in enumerate(MMCIF.iterAtoms()):

        if i == 0:
            assert atom.resn == "SER"
            assert atom.resn_short == "S"
            assert atom.resi == 2
            assert atom.chain == "A"
        if atom.chain == "A" and atom.resi == 2:
            natoms += 1
            if atom.atomn == "CA" and atom.element == "C":
                ca = True    
    assert natoms == 6
    assert ca == True


def test_protein_mapping_cif_3azm(MMCIF):
    """
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
    MMCIF.load("3azm")
    mapping = MMCIF.getChainProteinMapping()
    assert len(mapping.keys()) == 8
    
    assert mapping["A"].uniprot == mapping["E"].uniprot
    
    m = mapping["A"]
    assert m.uniprot == "P68431"
    assert m.begin == 0
    assert m.end == 135

    m = mapping["B"]
    assert m.uniprot == "P62805"
    assert m.begin == 0
    assert m.end == 102


def test_protein_mapping_cif_3nvk(MMCIF):
    MMCIF.load("3nvk")
    mapping = MMCIF.getChainProteinMapping()
    assert len(mapping.keys()) == 6
    # CIF chains:
    assert set("ABCDEF") == set([chain for chain, m in mapping.iteritems()])
    # Original PDB chains:
    assert set("AEFHIJ") == set([m.chain_author for chain, m in mapping.iteritems()])
    # CIF B is the same as PDB F:
    assert mapping["B"].uniprot == "Q8U4M1"


def test_NMR(MMCIF):
    code = "2aff"
    MMCIF.load(code)
    chains = MMCIF.getChains()
    assert set("AB") == set(chains)
    meta = MMCIF.getMetaData()
    assert meta.method.endswith("NMR")
    assert meta.title == "The solution structure of the Ki67FHA/hNIFK(226-269)3P complex"
    assert meta.description == "Antigen KI-67/MKI67 FHA domain interacting nucleolar phosphoprotein"
    assert len(list(MMCIF.iterAtoms())) == 2155 # taking into account atoms of 3 phosphorylated residues in chain B which are not counted


@pytest.mark.parametrize("code,nchains", [
    ("2a79", 8),
    ("4tim", 2),
    ("1pfk", 4),
])
def test_chains_pdb(PDB, code, nchains):
    PDB.load(code)
    chains = PDB.getChains()
    assert len(chains) == nchains


@pytest.mark.parametrize("code,nchains", [
    ("2a79", 8),
    ("4tim", 2),
    ("1pfk", 4),
])
def test_chains_cif(MMCIF, code, nchains):
    MMCIF.load(code)
    chains = MMCIF.getChains()
    assert len(chains) == nchains


def test_list_of_chains_cif(MMCIF):
    MMCIF.load("4pfk")
    chains = MMCIF.getChains()
    assert set(chains) == set(["A", "A_1", "A_2", "A_3"])


def test_list_all_cif(MMCIF):
    lcif = MMCIF.listAll()
    assert len([x for x in lcif]) == 99775


def test_list_all_pdb(PDB):
    lpdb = PDB.listAll()
    assert len([x for x in lpdb]) == 89414


