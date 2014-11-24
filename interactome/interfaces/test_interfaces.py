import pytest
# import unittest
# import tempfile
# import os
#tmpdir

# @pytest.fixture()
# def cleandir():
#     newpath = tempfile.mkdtemp()
#     return newpath
#     # os.chdir(newpath)


@pytest.fixture
def PDB_PATH():
    """ Path to PDB data with subdirectories biounits/ and mmcif/ """
    return "/Users/agoncear/data/pdb"


@pytest.fixture
def RESULTS_PATH():
    """ PATH to Interface results """
    return "/Users/agoncear/projects/Interactome/Workflow/Interfaces"


@pytest.fixture
def MMCIF(PDB_PATH):
    """ mmcif instance """
    from interactome.structures.mmcif import mmCifFile
    return mmCifFile(PDB_PATH)


@pytest.fixture
def IFACE(tmpdir):
    """ Interface instance """
    from interactome.interfaces.interface import Interface
    return Interface(tmpdir)

# @pytest.fixture
# def PDB(PATH):
#     """ pdb instance """
#     from pdb import PDBFile
#     return PDBFile(PATH)


############################################################################
@pytest.mark.parametrize('threshold, code, ncontacts', [
    (5.0, '4pfk', 3616),
    (5.0, '1hh1', 439),

    (4.0, '4pfk', 874),
    (4.0, '1hh1', 101)
    ])
def test_interfaces(MMCIF, tmpdir, threshold, code, ncontacts):
    from interactome.interfaces.interface import Interface
    iface = Interface(tmpdir, threshold)

    iteratoms = MMCIF.iterAtoms(code)
    n = iface.findContacts(code, iteratoms)
    assert ncontacts == n


# def test_(MMCIF):
#     lcif = MMCIF.listAll()
#     assert len([x for x in lcif]) == 99775


# def test_list_of_chains_cif(MMCIF):
#     chains = MMCIF.getChains("4pfk")
#     assert set(chains) == set(["A", "A_1", "A_2", "A_3"])
