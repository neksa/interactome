
import os
import signal
import time
import sys

import multiprocessing as mp

from interactome.structures.mmcif import mmCifFile
from interactome.interfaces.interface import Interface


def get_pdb_path():
    return "/Users/agoncear/data/pdb"


def get_results_path():
    return "/Users/agoncear/projects/Interactome/Workflow/Interfaces"


def get_threshold():
    # return 4.0
    return 5.0


def extract_contacts(params):
    code, gz_filename = params
    # print "Find contacts:", code
    mmcif = mmCifFile(get_pdb_path(), code)

    meta = mmcif.getMetaData()
    if meta.method.strip() in ('X-RAY DIFFRACTION', 'SOLID-STATE NMR', 'SOLUTION NMR'):
        iface = Interface(get_results_path(), get_threshold())
        iteratoms = mmcif.iterAtoms()
        n = iface.findContacts(code, iteratoms)

    # Extract mapping:
    chain_protein_mapping = mmcif.getChainProteinMapping()
    with open(get_mapping_fname(code), 'w') as o:
        o.write("chain\tchain_author\tuniprot\tseq_aln_begin\tseq_aln_end\tdb_aln_begin\tdb_aln_end\tauth_aln_begin\tauth_aln_end\n")
        for chain, m in chain_protein_mapping.iteritems():
            o.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chain, m.chain_author, m.uniprot, m.seq_aln_begin, m.seq_aln_end, m.db_aln_begin, m.db_aln_end, m.auth_aln_begin, m.auth_aln_end))
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chain, m.chain_author, m.uniprot, m.seq_aln_begin, m.seq_aln_end, m.db_aln_begin, m.db_aln_end, m.auth_aln_begin, m.auth_aln_end))
    print code, "OK"


def get_mapping_fname(code):
    results_dir = "{}/{}".format(get_results_path(), code[1:3])
    if not os.path.isdir(results_dir):
        try:
            os.makedirs(results_dir)
        except:
            pass
    mapping = "{}/{}_chain_protein_mapping.tab".format(results_dir, code)
    return mapping

# def extract_protein_mapping(params):
#     code, gz_filename = params
#     print "PDB: ", code
#     # print code, chain_protein_mapping
#     fname = get_mapping_fname(code)

#     # try:
#     #     with open(fname, 'r') as f: pass
#     #     print "skip"
#     #     return
#     # except:
#     #     pass

#     mmcif = mmCifFile(get_pdb_path())
#     mmcif.load(code)
#     chain_protein_mapping = mmcif.getChainProteinMapping()

#     with open(fname, 'w') as o:
#         o.write("chain\tchain_author\tuniprot\tbegin\tend\n")
#         for chain, m in chain_protein_mapping.iteritems():
#             o.write("{}\t{}\t{}\t{}\t{}\n".format(chain, m.chain_author, m.uniprot, m.begin, m.end))
#             print "{}\t{}\t{}\t{}\t{}\t{}".format(code, chain, m.chain_author, m.uniprot, m.begin, m.end)


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def isNMR(params):
    # if params[0] == "2aff": return True
    # return False
    m = isNMR.method.get(params[0])
    if m is None:
        return False
    if m == "NMR":
        return True
    return False


# def accepted_method(params):
#     """
#     ELECTRON CRYSTALLOGRAPHY    .
#     ELECTRON MICROSCOPY .
#     EPR EPR only as a supporting method
#     FIBER DIFFRACTION   .
#     FLUORESCENCE TRANSFER   FLUORESCENCE TRANSFER only as a supporting method
#     INFRARED SPECTROSCOPY   IR and FTIR only as supporting methods
#     NEUTRON DIFFRACTION .
#     POWDER DIFFRACTION  .
#     SOLID-STATE NMR .
#     SOLUTION NMR    .
#     SOLUTION SCATTERING .
#     THEORETICAL MODEL   THEORETICAL MODEL only as a supporting method
#     X-RAY DIFFRACTION
#     """
#     pass


if __name__ == '__main__':
    mmcif = mmCifFile(get_pdb_path())
    pool = mp.Pool(8, init_worker)
    # isNMR.method = {}
    # with open(get_pdb_path()+"/derived_data/pdb_entry_type.txt", 'r') as f:
    #     for line in f:
    #         code, entities, method = line.strip().split()
    #         isNMR.method[code] = method
    # pool.imap_unordered(extract_contacts, ifilter(isNMR, mmcif.listAll()))
    pool.imap_unordered(extract_contacts, mmcif.listAll())
    # pool.imap_unordered(extract_protein_mapping, mmcif.listAll())

    try:
        while(True):
            print "Watchdog... every 60 seconds (Ctrl-C to interrupt)"
            time.sleep(60)

    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        pool.terminate()
        pool.join()
        sys.exit(-1)
    pool.close()
    pool.join()
